#!/usr/bin/env Rscript
# Null model for the pan-cancer feasibility screen.
# Question: does a random 21-gene set, put through the identical screen, also select KIRC?
# Data: PanCanAtlas EB++ matrix from UCSC Xena (log2(x+1)); back-transformed to the RSEM scale
# used by code/02_pancancer_screen.R, after checking agreement with the cached cBioPortal values.
# Cohorts: the same samples and OS data as the original screen (data/<CANCER>_expr.json / _clin.json),
# restricted to the 20 hard-pass cancers. Screen rules copied from 02_pancancer_screen.R.

suppressPackageStartupMessages({
  library(data.table); library(jsonlite); library(survival); library(parallel)
})
ROOT <- "/Users/xiaolin/WorkBuddy/2026-09-04-20-00-08/paper2_triaptosis"
DATA <- file.path(ROOT, "data")
OUT  <- file.path(ROOT, "results/01_screen")
N_SETS <- as.integer(Sys.getenv("N_SETS", "200"))
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a)) b else a

TRIA21 <- c("PIK3C3","PIK3R4","MTM1","KEAP1","NFE2L2","SH3GL3","ELMO2","NCKAP1","ACTN2",
            "KXD1","WDR91","FYCO1","RAB9A","SCARB2","ANXA8","CLCN4","ATP13A2",
            "UVRAG","ATG14","NRBF2","BECN1")
screen0 <- read.csv(file.path(OUT, "pancancer_screen_matrix.csv"))
hard <- screen0$cancer[screen0$pass_hard]
cat("hard-pass cancers:", length(hard), "\n")

# ---- cohorts as in the original screen ---------------------------------------
cohorts <- lapply(hard, function(cc) {
  expr <- fromJSON(file.path(DATA, paste0(cc, "_expr.json")))
  clin <- fromJSON(file.path(DATA, paste0(cc, "_clin.json")))
  os_status <- sapply(clin, function(x) as.character(x$OS_STATUS %||% ""))
  os_months <- suppressWarnings(as.numeric(sapply(clin, function(x) x$OS_MONTHS %||% NA)))
  dead <- grepl("DECEASED|1", toupper(os_status))
  keep <- !is.na(os_months) & os_months > 0 & !is.na(dead) & (os_status != "")
  sid <- names(expr)
  pid <- substr(sid, 1, 12)
  # original: rownames -> patient ID, then intersect with clinical keep (first sample per patient wins in matching)
  common <- intersect(pid, names(clin)[keep])
  s_for_p <- sid[match(common, pid)]
  idx <- match(common, names(clin)[keep])
  list(cancer = cc, samples = s_for_p,
       time = os_months[keep][idx], status = as.numeric(dead[keep][idx]),
       cached = expr[s_for_p])
})
names(cohorts) <- hard
all_samples <- unique(unlist(lapply(cohorts, `[[`, "samples")))
cat("samples needed:", length(all_samples), "\n")

# ---- Xena matrix --------------------------------------------------------------
xf <- file.path(DATA, "pancan/EBpp_geneExp.xena.gz")
hdr <- names(fread(cmd = paste("gzip -dc", shQuote(xf), "| head -1"), header = TRUE))
keep_cols <- c(1, which(hdr %in% all_samples))
cat("Xena columns matched:", length(keep_cols) - 1, "of", length(all_samples), "\n")
X <- fread(cmd = paste("gzip -dc", shQuote(xf)), select = keep_cols, header = TRUE, showProgress = FALSE)
genes <- X[[1]]; X <- as.matrix(X[, -1]); rownames(X) <- genes
X <- X[, !duplicated(colnames(X)), drop = FALSE]
alias <- c(C19orf50 = "KXD1", KIAA0831 = "ATG14")   # older symbols in the Xena matrix
rownames(X)[rownames(X) %in% names(alias)] <- alias[rownames(X)[rownames(X) %in% names(alias)]]
stopifnot(all(TRIA21 %in% rownames(X)))
X <- 2^X - 1                                   # back to RSEM scale
X[X < 0 & X > -1e-6] <- 0
cat("Xena matrix:", nrow(X), "genes x", ncol(X), "samples\n")

# ---- agreement with cached cBioPortal values (21 genes, all cohorts) --------------
agree <- rbindlist(lapply(cohorts, function(co) {
  s <- intersect(co$samples, colnames(X))
  rbindlist(lapply(TRIA21, function(g) {
    a <- vapply(co$cached[s], function(x) as.numeric(x[[g]] %||% NA), numeric(1))
    b <- X[g, s]
    ok <- is.finite(a) & is.finite(b)
    data.table(cancer = co$cancer, gene = g, n = sum(ok),
               spearman = suppressWarnings(cor(a[ok], b[ok], method = "spearman")),
               max_rel_diff = max(abs(a[ok] - b[ok]) / pmax(abs(a[ok]), 1)))
  }))
}))
cat(sprintf("agreement: median Spearman %.5f, min %.5f; median max rel diff %.2e, worst %.2e\n",
            median(agree$spearman, na.rm = TRUE), min(agree$spearman, na.rm = TRUE),
            median(agree$max_rel_diff, na.rm = TRUE), max(agree$max_rel_diff, na.rm = TRUE)))
fwrite(agree, file.path(OUT, "screen_null_source_agreement.csv"))

# ---- screen function (rules of 02_pancancer_screen.R) ------------------------------
screen_one <- function(co, gset) {
  s <- co$samples
  inx <- s %in% colnames(X)
  mat2 <- t(X[gset, s[inx], drop = FALSE])
  sdf <- data.frame(time = co$time[inx], status = co$status[inx])
  n_sig <- 0L
  for (g in colnames(mat2)) {
    v <- mat2[, g]
    if (sum(!is.na(v)) < 2 || isTRUE(sd(v, na.rm = TRUE) == 0)) next
    fit <- tryCatch(coxph(Surv(time, status) ~ v, data = sdf), error = function(e) NULL)
    if (!is.null(fit)) {
      pv <- summary(fit)$coefficients[5]
      if (!is.na(pv) && pv < 0.05) n_sig <- n_sig + 1L
    }
  }
  set.seed(42)
  z <- scale(mat2); z[!is.finite(z)] <- 0
  km <- tryCatch(kmeans(z, centers = 2, nstart = 25), error = function(e) NULL)
  if (is.null(km)) return(data.frame(cancer = co$cancer, n_sig = n_sig, logrank_p = NA, balance = NA, pass_soft = FALSE))
  lr <- tryCatch(survdiff(Surv(time, status) ~ km$cluster, data = sdf), error = function(e) NULL)
  lp <- if (!is.null(lr)) 1 - pchisq(lr$chisq, df = 1) else NA
  bal <- min(table(km$cluster)) / length(km$cluster)
  data.frame(cancer = co$cancer, n_sig = n_sig, logrank_p = lp, balance = bal,
             pass_soft = !is.na(lp) && lp < 0.05 && n_sig >= 2 && bal > 0.15)
}
run_screen <- function(gset) {
  r <- do.call(rbind, lapply(cohorts, screen_one, gset = gset))
  r <- r[order(-r$pass_soft, -r$n_sig, r$logrank_p), ]
  r$rank <- seq_len(nrow(r))
  r
}

# ---- reproduce the real screen from Xena -------------------------------------------
real <- run_screen(TRIA21)
cmp <- merge(real, screen0[, c("cancer", "n_sig_cox", "logrank_p", "pass_soft")],
             by = "cancer", suffixes = c("_xena", "_orig"))
print(cmp[order(cmp$rank), ], digits = 3)
fwrite(cmp, file.path(OUT, "screen_null_reproduction.csv"))
cat("top cancer (Xena):", real$cancer[1], "\n")

# ---- random gene sets ---------------------------------------------------------------
Xs <- X[, intersect(all_samples, colnames(X))]
ok_gene <- !grepl("^\\?", rownames(Xs)) & rowSums(!is.finite(Xs)) == 0
ok_gene <- ok_gene & vapply(seq_len(nrow(Xs)), function(i) {
  all(vapply(cohorts, function(co) {
    v <- Xs[i, intersect(co$samples, colnames(Xs))]; sd(v) > 0
  }, logical(1)))
}, logical(1))
pool <- setdiff(rownames(Xs)[ok_gene], TRIA21)
cat("eligible gene pool:", length(pool), "\n")
rm(Xs); invisible(gc())

set.seed(2026)
sets <- replicate(N_SETS, sample(pool, 21), simplify = FALSE)
res <- mclapply(seq_along(sets), function(i) {
  r <- run_screen(sets[[i]]); r$set <- i; r
}, mc.cores = 6)
bad <- vapply(res, inherits, logical(1), "try-error")
if (any(bad)) stop(sum(bad), " sets failed: ", as.character(res[[which(bad)[1]]]))
res <- do.call(rbind, res)
fwrite(res, file.path(OUT, "screen_null_sets.csv"))
fwrite(data.table(set = seq_along(sets), genes = vapply(sets, paste, "", collapse = ";")),
       file.path(OUT, "screen_null_genesets.csv"))

top <- res[res$rank == 1, ]
kirc <- res[res$cancer == "KIRC", ]
obs <- real[real$cancer == "KIRC", ]
summ <- data.frame(
  n_sets = N_SETS,
  kirc_top_frac = mean(top$cancer == "KIRC"),
  kirc_pass_soft_frac = mean(kirc$pass_soft),
  kirc_nsig_median = median(kirc$n_sig),
  kirc_nsig_ge_obs_frac = mean(kirc$n_sig >= obs$n_sig),
  kirc_logrank_le_obs_frac = mean(kirc$logrank_p <= obs$logrank_p, na.rm = TRUE),
  n_soft_pass_median = median(tapply(res$pass_soft, res$set, sum)),
  obs_kirc_nsig = obs$n_sig, obs_kirc_logrank = obs$logrank_p
)
print(summ)
fwrite(summ, file.path(OUT, "screen_null_summary.csv"))
freq <- as.data.frame(sort(table(top$cancer), decreasing = TRUE))
names(freq) <- c("cancer", "times_top")
freq$frac <- freq$times_top / N_SETS
by_c <- aggregate(cbind(pass_soft, n_sig) ~ cancer, data = res, FUN = mean)
names(by_c) <- c("cancer", "pass_soft_frac", "mean_n_sig")
freq <- merge(by_c, freq, by = "cancer", all.x = TRUE)
freq$times_top[is.na(freq$times_top)] <- 0; freq$frac[is.na(freq$frac)] <- 0
freq <- freq[order(-freq$frac, -freq$mean_n_sig), ]
print(freq, digits = 3)
fwrite(freq, file.path(OUT, "screen_null_by_cancer.csv"))
