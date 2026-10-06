#!/usr/bin/env Rscript
# Paper2-Triaptosis: 泛癌可跑通性筛选 (Phase A)
# 数据驱动定主线癌种: 硬筛(样本量/事件) -> 单因素Cox -> k-means聚类 log-rank
# 原则: 全部真实数据, 不靠文献拍脑袋; 与第一篇同框架, 全新实现
suppressMessages({
  library(jsonlite); library(survival)
})

ROOT <- "/Users/xiaolin/WorkBuddy/2026-09-04-20-00-08/paper2_triaptosis"
DATA <- file.path(ROOT, "data")
OUT  <- file.path(ROOT, "results/01_screen")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
gene_cox_rows <- list()

cancers <- list.files(DATA, pattern = "_expr\\.json$")
cancers <- sub("_expr\\.json$", "", cancers)
cancers <- cancers[!grepl("immune|mtor|IMmotion", cancers)]  # 排除非泛癌缓存 (IMmotion150 为 ICI 队列)
cat("癌种数:", length(cancers), "\n")

score <- function(cancer) {
  expr <- fromJSON(file.path(DATA, paste0(cancer, "_expr.json")))
  clin <- fromJSON(file.path(DATA, paste0(cancer, "_clin.json")))
  if (is.null(clin) || length(clin) == 0) return(NULL)

  # 表达矩阵: 样本 x 基因 (统一基因顺序, 强制 numeric, 缺失补 NA)
  genes_all <- sort(unique(unlist(lapply(expr, names))))
  mat <- t(vapply(expr, function(x) {
    v <- rep(NA_real_, length(genes_all)); names(v) <- genes_all
    for (g in genes_all) if (!is.null(x[[g]])) v[g] <- as.numeric(x[[g]])
    v
  }, FUN.VALUE = numeric(length(genes_all))))
  colnames(mat) <- genes_all
  rownames(mat) <- names(expr)

  # 临床: OS 状态/时间
  os_status <- sapply(clin, function(x) as.character(x$OS_STATUS %||% ""))
  os_months <- suppressWarnings(as.numeric(sapply(clin, function(x) x$OS_MONTHS %||% NA)))
  dead <- grepl("DECEASED|1", toupper(os_status))
  keep <- !is.na(os_months) & os_months > 0 & !is.na(dead) & (os_status != "")
  os_months <- os_months[keep]; dead <- dead[keep]

  n_total <- nrow(mat)
  n_os    <- sum(keep)
  n_dead  <- sum(dead)

  # 硬筛1-2: 样本量 >=150, 死亡事件 >=30
  if (n_total < 150 || n_dead < 30) {
    return(data.frame(cancer, n_total, n_os, n_dead, n_sig_cox=NA, logrank_p=NA,
                      n_c1=NA, n_c2=NA, pass_hard=FALSE, pass_soft=FALSE))
  }

  # 对齐样本: TCGA 样本ID(带 -01 后缀) -> 患者ID(前12位)
  patient_of <- function(sid) substr(sid, 1, 12)
  rownames(mat) <- patient_of(rownames(mat))
  common <- intersect(rownames(mat), names(clin)[keep])
  mat2 <- mat[common, , drop = FALSE]
  idx <- match(common, names(clin)[keep])
  osm <- os_months[idx]   # os_months / dead are already restricted to keep
  dth <- dead[idx]
  sdf <- data.frame(time = osm, status = as.numeric(dth))

  # 指标4-5: 单因素 Cox (全部基因); 同步导出逐基因 HR/P (不改筛选规则)
  n_sig <- 0; hr_all <- c()
  for (g in colnames(mat)) {
    v <- mat2[, g]
    if (sum(!is.na(v)) < 2) next
    if (isTRUE(sd(v, na.rm = TRUE) == 0)) next
    fit <- tryCatch(coxph(Surv(time, status) ~ v, data = sdf, na.action = na.omit), error = function(e) NULL)
    if (!is.null(fit)) {
      sm <- summary(fit)
      pv <- sm$coefficients[5]
      hr <- sm$coefficients[2]
      hr_all <- c(hr_all, hr)
      gene_cox_rows[[length(gene_cox_rows) + 1]] <<- data.frame(
        cancer = cancer, gene = g, HR = as.numeric(hr), p = as.numeric(pv),
        stringsAsFactors = FALSE)
      if (!is.na(pv) && pv < 0.05) n_sig <- n_sig + 1
    }
  }

  # 指标6-7: k-means 2簇 + log-rank
  # k-means cannot take missing values; cluster the samples with complete expression
  # (STAD: 38 samples lack SH3GL3 in cBioPortal)
  cc_rows <- complete.cases(mat2)
  set.seed(42)
  z <- scale(mat2[cc_rows, , drop = FALSE])
  km <- tryCatch(kmeans(z, centers = 2, nstart = 25), error = function(e) NULL)
  if (is.null(km)) {
    return(data.frame(cancer, n_total, n_os, n_dead, n_sig_cox=n_sig, logrank_p=NA,
                      n_c1=NA, n_c2=NA, pass_hard=TRUE, pass_soft=FALSE))
  }
  sdf <- sdf[cc_rows, , drop = FALSE]
  sdf$clu <- km$cluster
  lr <- tryCatch(survdiff(Surv(time, status) ~ clu, data = sdf), error = function(e) NULL)
  logrank_p <- if (!is.null(lr)) 1 - pchisq(lr$chisq, df = 1) else NA
  n_c1 <- sum(km$cluster == 1); n_c2 <- sum(km$cluster == 2)

  pass_hard <- TRUE
  pass_soft <- (!is.na(logrank_p) && logrank_p < 0.05 && n_sig >= 2 &&
                min(n_c1, n_c2) / (n_c1 + n_c2) > 0.15)  # 簇平衡

  data.frame(cancer, n_total, n_os, n_dead, n_sig_cox=n_sig, logrank_p=logrank_p,
             n_c1=n_c1, n_c2=n_c2, pass_hard=pass_hard, pass_soft=pass_soft)
}

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a)) b else a

res <- lapply(cancers, score)
res <- do.call(rbind, res[!sapply(res, is.null)])
res <- res[order(-res$pass_soft, -res$n_sig_cox, res$logrank_p), ]
rownames(res) <- NULL

write.csv(res, file.path(OUT, "pancancer_screen_matrix.csv"), row.names = FALSE)
if (length(gene_cox_rows) > 0) {
  gcox <- do.call(rbind, gene_cox_rows)
  write.csv(gcox, file.path(OUT, "gene_cancer_cox.csv"), row.names = FALSE)
  cat("gene x cancer Cox rows:", nrow(gcox), "\n")
}
print(res, row.names = FALSE)

cat("\n=== 通过软标准(数据驱动候选主线癌种) ===\n")
hits <- res[res$pass_soft, ]
if (nrow(hits) > 0) print(hits, row.names = FALSE) else cat("无\n")
cat("\n结果已保存:", file.path(OUT, "pancancer_screen_matrix.csv"), "\n")
