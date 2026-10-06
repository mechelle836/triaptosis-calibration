#!/usr/bin/env Rscript
# Other cell-death gene sets under the exact TAS null protocol.
# Same cohort, same 10 folds (seed 2026), same elastic-net (alpha = 0.5) nested CV as
# code/07_triaptosis_null.R; each set is ranked against the 195 random 21-gene sets there.
# The three comparison sets have 19 genes, the null sets 21.

suppressPackageStartupMessages({
  library(survival); library(jsonlite); library(glmnet)
})
P2   <- "/Users/xiaolin/WorkBuddy/2026-09-04-20-00-08/paper2_triaptosis"
KIRC <- "/Users/xiaolin/CodeBuddy/Claw/pan-cancer_disulfidptosis_project"
OUT  <- file.path(P2, "results/05_null")

TRIA21 <- c("ACTN2","ANXA8","ATG14","ATP13A2","BECN1","CLCN4","ELMO2","FYCO1",
            "KEAP1","KXD1","MTM1","NCKAP1","NFE2L2","NRBF2","PIK3C3","PIK3R4",
            "RAB9A","SCARB2","SH3GL3","UVRAG","WDR91")
SETS <- list(
  Ferroptosis = c("GPX4","SLC7A11","SLC3A2","GCLC","GCLM","GSR","GSS","TFRC","FTH1",
                  "FTL","NCOA4","HMOX1","IREB2","ACSL4","LPCAT3","ALOX15","NFE2L2","AIFM2","SAT1"),
  Apoptosis = c("BAX","BAK1","BID","BCL2","BCL2L11","BCL2L1","BCL2L2","MCL1","BAD",
                "BBC3","PMAIP1","APAF1","CASP9","CASP3","CASP7","CYCS","DIABLO","XIAP","TP53"),
  Disulfidptosis = c("SLC7A11","NCKAP1","NCF1","RAC1","ACTB","MYH9","ACTN4","FLNA","FLNB",
                     "IQGAP1","CD2AP","DSTN","LIMK1","PAK1","WASF2","NCKAP1L","CYBA","NCF2","RAC2")
)

locked <- readRDS(file.path(KIRC, "results/04_mitoxy/KIRC_model_results.rds"))
d0 <- locked$d[!duplicated(locked$d$patient), ]
mat_list <- fromJSON(file.path(KIRC, "data/KIRC_full_transcriptome.json"))
common <- intersect(names(mat_list), d0$sample)
d <- d0[match(common, d0$sample), ]
allg <- sort(unique(unlist(lapply(mat_list[common], names))))
mat <- matrix(NA_real_, length(allg), length(common), dimnames = list(allg, common))
for (i in seq_along(common)) {
  v <- mat_list[[common[i]]]; mat[names(v), i] <- unlist(v)
}
X <- t(mat)
keep <- colnames(X)[colSums(is.na(X)) == 0 & apply(X, 2, sd, na.rm = TRUE) > 0]
X <- X[, keep, drop = FALSE]
dat <- data.frame(os_t = d$os_t, os_s = d$os_s); rownames(dat) <- d$sample
dat <- dat[!is.na(dat$os_t) & !is.na(dat$os_s) & dat$os_t > 0, ]
X <- X[rownames(dat), , drop = FALSE]
cat("cohort:", nrow(X), "patients,", sum(dat$os_s), "deaths\n")

K <- 10
set.seed(2026)
FOLDS <- sample(rep(1:K, length.out = nrow(dat)))

fit_enet <- function(Z, time, status) {
  cvf <- cv.glmnet(Z, Surv(time, status), family = "cox", alpha = 0.5,
                   nfolds = 10, cox.ties = "efron")
  b <- as.numeric(coef(cvf, s = "lambda.min")); names(b) <- colnames(Z); b
}
cv_nested_enet <- function(genes) {
  pred <- rep(NA_real_, nrow(dat))
  for (k in 1:K) {
    tr <- FOLDS != k; te <- FOLDS == k
    Ztr <- scale(X[tr, genes, drop = FALSE]); Ztr[!is.finite(Ztr)] <- 0
    beta <- try(fit_enet(Ztr, dat$os_t[tr], dat$os_s[tr]), silent = TRUE)
    if (inherits(beta, "try-error")) next
    mu <- attr(Ztr, "scaled:center"); sdv <- attr(Ztr, "scaled:scale"); sdv[sdv == 0] <- 1
    Zte <- scale(X[te, genes, drop = FALSE], center = mu, scale = sdv); Zte[!is.finite(Zte)] <- 0
    pred[te] <- as.numeric(Zte %*% beta)
  }
  ok <- !is.na(pred)
  as.numeric(concordance(Surv(os_t, os_s) ~ pred, data = dat[ok, ], reverse = TRUE)$concordance)
}

null <- read.csv(file.path(OUT, "triaptosis_random_null.csv"))$cv_nested
null <- null[is.finite(null)]
cat("null sets:", length(null), "\n")

# Outer folds are fixed; the inner cv.glmnet folds are random, so each set is repeated
# over 10 seeds and summarised by the mean.
SEEDS <- 2026 + 0:9
res <- do.call(rbind, lapply(c(list(Triaptosis = TRIA21), SETS), function(g) {
  g2 <- intersect(g, colnames(X))
  cs <- vapply(SEEDS, function(s) { set.seed(s); cv_nested_enet(g2) }, numeric(1))
  c_n <- mean(cs)
  data.frame(n_genes = length(g), n_found = length(g2), cv_nested = c_n,
             cv_min = min(cs), cv_max = max(cs),
             percentile = 100 * mean(null < c_n),
             z = (c_n - mean(null)) / sd(null),
             overlap_with_triaptosis = paste(intersect(g, TRIA21), collapse = ";"))
}))
res <- cbind(gene_set = rownames(res), res, row.names = NULL)
print(res, digits = 4)
write.csv(res, file.path(OUT, "pcd_benchmark_same_protocol.csv"), row.names = FALSE)
