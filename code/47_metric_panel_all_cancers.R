#!/usr/bin/env Rscript
# 47_metric_panel_all_cancers.R — 20 癌种 × 3 指标校准面板（供 Table S6）
# 与 45/46 同数据、同抽样池、同 seed=2026（同 200 个随机集）。
# 指标（全部为已发表标准方法）：
#   n_sig（P<0.05 计数，原稿口径）/ mean(-log10 P)（Fisher 1932 形式）/ Stouffer Z（1949）
# 输出: results/01_screen/metric_panel_by_cancer.csv
#   每癌种：observed 三指标 + 在 200 集零分布中的百分位与 z

suppressPackageStartupMessages({ library(data.table); library(jsonlite); library(survival); library(parallel) })
ROOT <- "/Users/xiaolin/WorkBuddy/2026-09-04-20-00-08/paper2_triaptosis"
DATA <- file.path(ROOT, "data"); OUT <- file.path(ROOT, "results/01_screen")
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a)) b else a

TRIA21 <- c("PIK3C3","PIK3R4","MTM1","KEAP1","NFE2L2","SH3GL3","ELMO2","NCKAP1","ACTN2",
            "KXD1","WDR91","FYCO1","RAB9A","SCARB2","ANXA8","CLCN4","ATP13A2",
            "UVRAG","ATG14","NRBF2","BECN1")

# ---- 队列与表达矩阵（与 45 逐行一致） ----
screen0 <- read.csv(file.path(OUT, "pancancer_screen_matrix.csv"))
hard <- screen0$cancer[screen0$pass_hard]
cohorts <- lapply(hard, function(cc) {
  expr <- fromJSON(file.path(DATA, paste0(cc, "_expr.json")))
  clin <- fromJSON(file.path(DATA, paste0(cc, "_clin.json")))
  os_status <- sapply(clin, function(x) as.character(x$OS_STATUS %||% ""))
  os_months <- suppressWarnings(as.numeric(sapply(clin, function(x) x$OS_MONTHS %||% NA)))
  dead <- grepl("DECEASED|1", toupper(os_status))
  keep <- !is.na(os_months) & os_months > 0 & !is.na(dead) & (os_status != "")
  sid <- names(expr); pid <- substr(sid, 1, 12)
  common <- intersect(pid, names(clin)[keep])
  s_for_p <- sid[match(common, pid)]
  idx <- match(common, names(clin)[keep])
  list(cancer = cc, samples = s_for_p,
       time = os_months[keep][idx], status = as.numeric(dead[keep][idx]))
})
names(cohorts) <- hard

xf <- file.path(DATA, "pancan/EBpp_geneExp.xena.gz")
all_samples <- unique(unlist(lapply(cohorts, `[[`, "samples")))
hdr <- names(fread(cmd = paste("gzip -dc", shQuote(xf), "| head -1"), header = TRUE))
keep_cols <- c(1, which(hdr %in% all_samples))
X <- fread(cmd = paste("gzip -dc", shQuote(xf)), select = keep_cols, header = TRUE, showProgress = FALSE)
genes <- X[[1]]; X <- as.matrix(X[, -1]); rownames(X) <- genes
X <- X[, !duplicated(colnames(X)), drop = FALSE]
alias <- c(C19orf50 = "KXD1", KIAA0831 = "ATG14")
rownames(X)[rownames(X) %in% names(alias)] <- alias[rownames(X)[rownames(X) %in% names(alias)]]
stopifnot(all(TRIA21 %in% rownames(X)))
X <- 2^X - 1
X[X < 0 & X > -1e-6] <- 0

# ---- 抽样池与 200 集（与 45/46 一致） ----
ok_gene <- !grepl("^\\?", rownames(X)) & (rowSums(!is.finite(X)) == 0)
for (co_n in names(cohorts)) {
  m <- X[, intersect(cohorts[[co_n]]$samples, colnames(X)), drop = FALSE]
  ok_gene <- ok_gene & (apply(m, 1, sd) > 0)
  rm(m); invisible(gc())
}
pool <- setdiff(rownames(X)[ok_gene], TRIA21)
invisible(gc())
set.seed(2026)
sets <- replicate(200, sample(pool, 21), simplify = FALSE)

# ---- 每癌种每基因 Cox P ----
p_floor <- 1e-300
panel_one <- function(co, gset) {
  s <- co$samples; inx <- s %in% colnames(X)
  mat2 <- t(X[gset, s[inx], drop = FALSE])
  sdf <- data.frame(time = co$time[inx], status = co$status[inx])
  ps <- rep(NA_real_, ncol(mat2))
  for (g in colnames(mat2)) {
    v <- mat2[, g]
    if (sum(!is.na(v)) < 2 || isTRUE(sd(v, na.rm = TRUE) == 0)) next
    fit <- tryCatch(coxph(Surv(time, status) ~ v, data = sdf), error = function(e) NULL)
    if (!is.null(fit)) ps[g] <- summary(fit)$coefficients[5]
  }
  p <- pmin(pmax(ps[!is.na(ps)], p_floor), 1 - 1e-16); k <- length(p)
  c(n_sig = sum(p < 0.05),
    mean_nlp = mean(-log10(p)),
    stouffer_z = sum(qnorm(p, lower.tail = FALSE)) / sqrt(k))
}

run_all <- function(gset) do.call(rbind, lapply(cohorts, function(co)
  data.frame(cancer = co$cancer, t(panel_one(co, gset)))))
obs <- run_all(TRIA21)

res <- mclapply(seq_along(sets), function(i) {
  r <- run_all(sets[[i]]); r$set <- i; r
}, mc.cores = 6)
null <- do.call(rbind, res)

summ <- do.call(rbind, lapply(split(null, null$cancer), function(rc) {
  cc <- rc$cancer[1]; o <- obs[obs$cancer == cc, ]
  data.frame(cancer = cc,
    n_sig_obs = o$n_sig, n_sig_pct = mean(rc$n_sig <= o$n_sig) * 100,
    mean_nlp_obs = round(o$mean_nlp, 3), mean_nlp_pct = mean(rc$mean_nlp <= o$mean_nlp) * 100,
    stouffer_obs = round(o$stouffer_z, 2), stouffer_pct = mean(rc$stouffer_z <= o$stouffer_z) * 100)
}))
summ <- summ[order(-summ$n_sig_pct), ]
fwrite(summ, file.path(OUT, "metric_panel_by_cancer.csv"))
cat("=== 20 癌种 × 3 指标面板（triaptosis 在 200 集零分布中的百分位）===\n")
print(summ, row.names = FALSE)
cat("=== done ===\n")
