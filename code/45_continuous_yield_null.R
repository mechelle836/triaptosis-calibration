#!/usr/bin/env Rscript
# 45_continuous_yield_null.R
# 方法学审查 #1 的回应：用连续 per-gene yield 替代"P<0.05 二值化计数"，检验核心结论
# "KIRC 是唯一 triaptosis 基因超越随机分布的癌种"是否对二值化阈值敏感。
#
# 方法依据（非自造）：
#   - mean(-log10 P) 是 Fisher 合并 P 值法 (Fisher, 1932, Statistical Methods for
#     Research Workers) 的常用变体，广泛用作基因集层面的连续汇总统计。
#   - 避免连续变量二值化的信息损失：Royston, Altman & Sauerbrei (2006, Stat Med
#     25:127-141)。
# 设计：与 38_screen_null.R 完全相同的数据、队列、随机集（seed 2026，同 200 个集），
#   仅把筛选指标从"P<0.05 基因计数"换成"21 基因 -log10(Cox P) 的均值/中位数"。
#   triaptosis 集的统计量在 200 集零分布中的位置（百分位、z）即为连续口径下的结论。
#
# 输出: results/01_screen/screen_null_continuous_yield.csv
#       results/01_screen/screen_null_continuous_summary.csv

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

# ---- 数据加载（与 38_screen_null.R 完全一致） ----
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
X <- 2^X - 1                                   # RSEM 线性尺度，与 38 一致
X[X < 0 & X > -1e-6] <- 0
cat("Xena matrix:", nrow(X), "x", ncol(X), "\n")

# ---- 连续 yield：每基因 Cox P -> mean/median(-log10 P) ----
p_floor <- 1e-300  # 防止 -log10(0)
yield_one <- function(co, gset) {
  s <- co$samples
  inx <- s %in% colnames(X)
  mat2 <- t(X[gset, s[inx], drop = FALSE])
  sdf <- data.frame(time = co$time[inx], status = co$status[inx])
  ps <- rep(NA_real_, ncol(mat2)); names(ps) <- colnames(mat2)
  for (g in colnames(mat2)) {
    v <- mat2[, g]
    if (sum(!is.na(v)) < 2 || isTRUE(sd(v, na.rm = TRUE) == 0)) next
    fit <- tryCatch(coxph(Surv(time, status) ~ v, data = sdf), error = function(e) NULL)
    if (!is.null(fit)) ps[g] <- summary(fit)$coefficients[5]
  }
  ps <- pmax(ps, p_floor, na.rm = TRUE)
  nlp <- -log10(ps)
  data.frame(cancer = co$cancer,
             mean_nlp = mean(nlp, na.rm = TRUE),
             median_nlp = median(nlp, na.rm = TRUE),
             n_sig = sum(ps < 0.05, na.rm = TRUE))
}
run_yield <- function(gset) do.call(rbind, lapply(cohorts, yield_one, gset = gset))

# ---- triaptosis 集 ----
obs <- run_yield(TRIA21)
obs_k <- obs[obs$cancer == "KIRC", ]
cat(sprintf("KIRC (triaptosis): mean(-log10P)=%.3f  median=%.3f  n_sig=%d\n",
            obs_k$mean_nlp, obs_k$median_nlp, obs_k$n_sig))

# ---- 200 随机集（与 38 同一抽样：seed 2026, pool 排除 TRIA21） ----
Xs <- X[, intersect(all_samples, colnames(X))]
ok_gene <- !grepl("^\\?", rownames(Xs)) & rowSums(!is.finite(Xs)) == 0
ok_gene <- ok_gene & vapply(seq_len(nrow(Xs)), function(i) {
  all(vapply(cohorts, function(co) {
    v <- Xs[i, intersect(co$samples, colnames(Xs))]; sd(v) > 0
  }, logical(1)))
}, logical(1))
pool <- setdiff(rownames(Xs)[ok_gene], TRIA21)
rm(Xs); invisible(gc())
set.seed(2026)
sets <- replicate(N_SETS, sample(pool, 21), simplify = FALSE)

res <- mclapply(seq_along(sets), function(i) {
  r <- run_yield(sets[[i]]); r$set <- i; r
}, mc.cores = 6)
res <- do.call(rbind, res)
fwrite(res, file.path(OUT, "screen_null_continuous_yield.csv"))

# ---- 对比：triaptosis 在连续零分布中的位置（逐癌种 + KIRC 聚焦） ----
kirc_null <- res[res$cancer == "KIRC", ]
summ <- data.frame(
  statistic = c("mean(-log10P)", "median(-log10P)", "n_sig (binary, for reference)"),
  obs_KIRC = c(obs_k$mean_nlp, obs_k$median_nlp, obs_k$n_sig),
  null_median = c(median(kirc_null$mean_nlp), median(kirc_null$median_nlp), median(kirc_null$n_sig)),
  null_q95 = c(quantile(kirc_null$mean_nlp, 0.95), quantile(kirc_null$median_nlp, 0.95), quantile(kirc_null$n_sig, 0.95)),
  null_max = c(max(kirc_null$mean_nlp), max(kirc_null$median_nlp), max(kirc_null$n_sig)),
  percentile_of_obs = c(mean(kirc_null$mean_nlp <= obs_k$mean_nlp),
                        mean(kirc_null$median_nlp <= obs_k$median_nlp),
                        mean(kirc_null$n_sig <= obs_k$n_sig)),
  z_of_obs = c((obs_k$mean_nlp - mean(kirc_null$mean_nlp)) / sd(kirc_null$mean_nlp),
               (obs_k$median_nlp - mean(kirc_null$median_nlp)) / sd(kirc_null$median_nlp),
               (obs_k$n_sig - mean(kirc_null$n_sig)) / sd(kirc_null$n_sig))
)
print(summ, digits = 3)
fwrite(summ, file.path(OUT, "screen_null_continuous_summary.csv"))

# 逐癌种：triaptosis 的 mean_nlp 在该癌种零分布中的百分位（20 个癌种）
per_cancer <- do.call(rbind, lapply(split(res, res$cancer), function(rc) {
  cc <- rc$cancer[1]; o <- obs[obs$cancer == cc, ]
  data.frame(cancer = cc,
             obs_mean_nlp = o$mean_nlp,
             null_median = median(rc$mean_nlp),
             percentile = mean(rc$mean_nlp <= o$mean_nlp),
             z = (o$mean_nlp - mean(rc$mean_nlp)) / sd(rc$mean_nlp))
}))
per_cancer <- per_cancer[order(-per_cancer$percentile), ]
print(per_cancer, digits = 3)
fwrite(per_cancer, file.path(OUT, "screen_null_continuous_by_cancer.csv"))
cat("\n判读: 若 KIRC 的连续 yield 百分位仍显著高于其他癌种且超过零分布 95 百分位,\n")
cat("则核心结论不依赖 P=0.05 二值化阈值, 对 KEAP1 边界不敏感.\n")
cat("=== done ===\n")
