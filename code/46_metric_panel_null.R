#!/usr/bin/env Rscript
# 46_metric_panel_null.R — 多指标零分布对比（KIRC 聚焦）
# 与 45_continuous_yield_null.R 完全同数据、同抽样池、同 seed=2026（同 200 个随机集），
# 仅在 KIRC 上扩展统计量面板：
#   1) n_sig（P<0.05 计数）      — 原稿口径（二值化）
#   2) mean(-log10 P)            — Fisher 合并 P 值法的常用形式 (Fisher, 1932)
#   3) Fisher 合并 P（-log10）   — Fisher (1932), df=2k 卡方
#   4) Stouffer Z                — Stouffer et al. (1949)
#   5) ACAT 统计量               — Liu Y et al. (2019) Am J Hum Genet 104:410
#   6) median(-log10 P)          — 稳健变体
# 输出: results/01_screen/metric_panel_kirc.csv

suppressPackageStartupMessages({ library(data.table); library(jsonlite); library(survival) })
ROOT <- "/Users/xiaolin/WorkBuddy/2026-09-04-20-00-08/paper2_triaptosis"
DATA <- file.path(ROOT, "data"); OUT <- file.path(ROOT, "results/01_screen")
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a)) b else a

TRIA21 <- c("PIK3C3","PIK3R4","MTM1","KEAP1","NFE2L2","SH3GL3","ELMO2","NCKAP1","ACTN2",
            "KXD1","WDR91","FYCO1","RAB9A","SCARB2","ANXA8","CLCN4","ATP13A2",
            "UVRAG","ATG14","NRBF2","BECN1")

# ---- 队列加载（与 45 逐行一致，用于构建相同抽样池） ----
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

# ---- 抽样池与 200 集（与 45 同规则：有限 + 各癌种 sd>0 + 排除 TRIA21；低内存实现） ----
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

# ---- KIRC 每基因 Cox P ----
p_floor <- 1e-300
co <- cohorts[["KIRC"]]
gene_p <- function(gset) {
  s <- co$samples; inx <- s %in% colnames(X)
  mat2 <- t(X[gset, s[inx], drop = FALSE])
  sdf <- data.frame(time = co$time[inx], status = co$status[inx])
  ps <- rep(NA_real_, ncol(mat2)); names(ps) <- colnames(mat2)
  for (g in colnames(mat2)) {
    v <- mat2[, g]
    if (sum(!is.na(v)) < 2 || isTRUE(sd(v, na.rm = TRUE) == 0)) next
    fit <- tryCatch(coxph(Surv(time, status) ~ v, data = sdf), error = function(e) NULL)
    if (!is.null(fit)) ps[g] <- summary(fit)$coefficients[5]
  }
  pmax(ps, p_floor, na.rm = TRUE)
}

metrics <- function(p) {
  p <- pmin(pmax(p[!is.na(p)], p_floor), 1 - 1e-16); k <- length(p)
  c(n_sig      = sum(p < 0.05),
    mean_nlp   = mean(-log10(p)),
    median_nlp = median(-log10(p)),
    fisher_nlp = -log10(pchisq(-2 * sum(log(p)), df = 2 * k, lower.tail = FALSE)),
    stouffer_z = sum(qnorm(p, lower.tail = FALSE)) / sqrt(k),
    acat_stat  = mean(tan((0.5 - p) * pi)))
}

obs_p <- gene_p(TRIA21); obs <- metrics(obs_p)
null <- t(sapply(sets, function(st) metrics(gene_p(st))))

df <- data.frame(metric = names(obs), observed = round(as.numeric(obs), 3),
                 null_mean = round(colMeans(null, na.rm = TRUE), 3),
                 null_q95 = round(apply(null, 2, quantile, 0.95, na.rm = TRUE), 3),
                 null_max = round(apply(null, 2, max, na.rm = TRUE), 3),
                 pct = round(sapply(names(obs), function(m) mean(null[, m] <= obs[m]) * 100), 1),
                 z = round(sapply(names(obs), function(m) (obs[m] - mean(null[, m], na.rm = TRUE)) / sd(null[, m], na.rm = TRUE)), 2))
write.csv(df, file.path(OUT, "metric_panel_kirc.csv"), row.names = FALSE)

cat("\n=== 6 指标零分布对比（KIRC，与 45 同池同 200 集）===\n")
print(df, row.names = FALSE)
pp <- sort(obs_p)
cat(sprintf("\n21 基因 P 值分布: P<0.001: %d | 0.001-0.01: %d | 0.01-0.05: %d | >=0.05: %d\n",
    sum(pp<0.001), sum(pp>=0.001 & pp<0.01), sum(pp>=0.01 & pp<0.05), sum(pp>=0.05)))
cat("最弱 6 个显著基因 P:", paste(round(tail(pp[pp<0.05], 6), 4), collapse=", "), "\n")
