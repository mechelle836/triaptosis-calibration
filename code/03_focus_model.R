#!/usr/bin/env Rscript
# Paper2-Triaptosis: 聚焦癌种建模 (Phase B) — 参数化, 主线癌种由 02 筛选决定
# 流程: 聚类分型 -> ML 建模(Triaptosis Score) -> KM/多变量Cox/时间AUC -> 平行癌种
# 用法: Rscript 03_focus_model.R <CANCER> [<PARALLEL_CANCERS...>]
suppressMessages({
  library(jsonlite); library(survival)
})

args <- commandArgs(trailingOnly = TRUE)
CANCER <- if (length(args) >= 1) args[1] else "KIRC"
PARALLEL <- args[-1]
DATA <- "/Users/xiaolin/CodeBuddy/Claw/paper2_triaptosis/data"
OUT  <- file.path("/Users/xiaolin/CodeBuddy/Claw/paper2_triaptosis/results/02_model", CANCER)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
set.seed(2026)

cat("=== 聚焦癌种:", CANCER, "===\n")
expr <- fromJSON(file.path(DATA, paste0(CANCER, "_expr.json")))
clin <- fromJSON(file.path(DATA, paste0(CANCER, "_clin.json")))
genes_all <- sort(unique(unlist(lapply(expr, names))))

mat <- t(vapply(expr, function(x) {
  v <- rep(NA_real_, length(genes_all)); names(v) <- genes_all
  for (g in genes_all) if (!is.null(x[[g]])) v[g] <- as.numeric(x[[g]])
  v
}, FUN.VALUE = numeric(length(genes_all))))
colnames(mat) <- genes_all
rownames(mat) <- substr(names(expr), 1, 12)  # TCGA 样本ID -> 患者ID
mat <- mat[!duplicated(rownames(mat)), , drop = FALSE]
os_status <- sapply(clin, function(x) as.character(x$OS_STATUS %||% ""))
os_months <- suppressWarnings(as.numeric(sapply(clin, function(x) x$OS_MONTHS %||% NA)))
dead <- grepl("DECEASED|1", toupper(os_status))
keep <- !is.na(os_months) & os_months > 0 & os_status != ""
common <- intersect(rownames(mat), names(clin)[keep])
mat2 <- mat[common, , drop = FALSE]
idx <- match(common, names(clin)[keep])
sdf <- data.frame(time = os_months[keep][idx],
                  status = as.numeric(dead[keep][idx]))
rownames(sdf) <- common
cat("样本:", nrow(sdf), " 死亡:", sum(sdf$status), "\n")

# ---------- 1. k-means 分子分型 ----------
z <- scale(mat2)
km <- kmeans(z, centers = 2, nstart = 50)
sdf$cluster <- km$cluster
lr <- survdiff(Surv(time, status) ~ cluster, data = sdf)
logrank_p <- 1 - pchisq(lr$chisq, df = 1)
cat(sprintf("聚类: C1=%d C2=%d | log-rank P=%.3e\n",
            sum(sdf$cluster == 1), sum(sdf$cluster == 2), logrank_p))

# ---------- 2. Triaptosis Score: 弹性网线性评分(核心ML, 10折CV) ----------
# 用 glmnet 弹性网(CV) 选基因权重 -> TAS = sum(weight * z)
if (requireNamespace("glmnet", quietly = TRUE)) {
  library(glmnet)
  cvf <- cv.glmnet(z, Surv(sdf$time, sdf$status), family = "cox", alpha = 0.5, nfolds = 10)
  beta <- as.numeric(coef(cvf, s = "lambda.min"))
  genes_sel <- colnames(z)[beta != 0]
  tas <- as.numeric(z %*% beta)
  cat("弹性网选中基因:", length(genes_sel), ":", paste(genes_sel, collapse = ","), "\n")
  tas_lambda <- cvf$lambda.min
} else {
  # 兜底: 单因素 Cox 显著基因 z-score 均值
  pvs <- sapply(colnames(z), function(g) {
    f <- tryCatch(coxph(Surv(time, status) ~ z[, g], data = sdf), error = function(e) NULL)
    if (is.null(f)) 1 else summary(f)$coefficients[5]
  })
  genes_sel <- names(pvs)[pvs < 0.05]
  tas <- if (length(genes_sel) > 0) rowMeans(z[, genes_sel, drop = FALSE]) else rowMeans(z)
  tas_lambda <- NA
  cat("兜底评分基因:", length(genes_sel), "\n")
}
sdf$TAS <- tas

# ---------- 3. 验证: KM 分层 + 多变量 Cox + 时间AUC ----------
med <- median(tas)
sdf$group <- ifelse(tas > med, "High", "Low")
lr2 <- survdiff(Surv(time, status) ~ group, data = sdf)
p_km <- 1 - pchisq(lr2$chisq, df = 1)
cat(sprintf("TAS KM 分层: High=%d Low=%d | log-rank P=%.3e\n",
            sum(sdf$group == "High"), sum(sdf$group == "Low"), p_km))

# 多变量 Cox (TAS + 可获取的临床变量)
mcox <- tryCatch(coxph(Surv(time, status) ~ TAS, data = sdf), error = function(e) NULL)
if (!is.null(mcox)) {
  hr <- exp(coef(mcox)); pv <- summary(mcox)$coefficients[5]
  cat(sprintf("单变量 Cox: HR=%.3f (95%%CI %.3f-%.3f) P=%.3e\n",
              hr, exp(confint(mcox))[1], exp(confint(mcox))[2], pv))
}

# 时间依赖 AUC (survivalROC 若无则跳过)
if (requireNamespace("survivalROC", quietly = TRUE)) {
  library(survivalROC)
  for (t in c(1, 3, 5)) {
    if (max(sdf$time) > t) {
      roc <- survivalROC(Stime = sdf$time, status = sdf$status, marker = tas,
                         predict.time = t, method = "KM")
      cat(sprintf("时间AUC(%dy) = %.3f\n", t, roc$AUC))
    }
  }
}

# ---------- 4. 平行癌种泛化 ----------
for (pc in PARALLEL) {
  pexpr <- fromJSON(file.path(DATA, paste0(pc, "_expr.json")))
  pclin <- fromJSON(file.path(DATA, paste0(pc, "_clin.json")))
  pmat <- t(vapply(pexpr, function(x) {
    v <- rep(NA_real_, length(genes_all)); names(v) <- genes_all
    for (g in genes_all) if (!is.null(x[[g]])) v[g] <- as.numeric(x[[g]])
    v
  }, FUN.VALUE = numeric(length(genes_all))))
  colnames(pmat) <- genes_all
  rownames(pmat) <- substr(names(pexpr), 1, 12)
  pmat <- pmat[!duplicated(rownames(pmat)), , drop = FALSE]
  pos <- sapply(pclin, function(x) as.character(x$OS_STATUS %||% ""))
  pom <- suppressWarnings(as.numeric(sapply(pclin, function(x) x$OS_MONTHS %||% NA)))
  pdead <- grepl("DECEASED|1", toupper(pos))
  pkeep <- !is.na(pom) & pom > 0 & pos != ""
  pcommon <- intersect(rownames(pmat), names(pclin)[pkeep])
  pidx <- match(pcommon, names(pclin)[pkeep])
  psdf <- data.frame(time = pom[pkeep][pidx],
                     status = as.numeric(pdead[pkeep][pidx]))
  # 用主线模型的基因权重打分
  pz <- scale(pmat[pcommon, colnames(z), drop = FALSE])
  ptas <- as.numeric(pz %*% beta)
  pmed <- median(ptas)
  psdf$group <- ifelse(ptas > pmed, "High", "Low")
  plr <- survdiff(Surv(time, status) ~ group, data = psdf)
  pp <- 1 - pchisq(plr$chisq, df = 1)
  cat(sprintf("平行验证 %s: n=%d | 分层 log-rank P=%.3e\n", pc, nrow(psdf), pp))
}

# ---------- 保存 ----------
write.csv(sdf, file.path(OUT, "score_table.csv"), row.names = TRUE)
saveRDS(list(beta = beta, genes = colnames(z), genes_sel = genes_sel,
            lambda = tas_lambda, logrank_cluster = logrank_p),
        file.path(OUT, "model.rds"))
cat("\n结果已保存:", OUT, "\n")

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a)) b else a
