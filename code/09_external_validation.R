#!/usr/bin/env Rscript
# Paper2-Triaptosis: 外部验证 — GSE22541 (独立 ccRCC 队列, DFI/DFS 终点)
# + TCGA-KIRC 6:4 时间分割内部验证 (OS 终点)
suppressMessages({
  library(jsonlite); library(survival)
})

DATA <- "/Users/xiaolin/CodeBuddy/Claw/paper2_triaptosis/data"
EXT  <- file.path(DATA, "external")
MODEL_DIR <- "/Users/xiaolin/CodeBuddy/Claw/paper2_triaptosis/results/02_model/KIRC"
OUT <- "/Users/xiaolin/CodeBuddy/Claw/paper2_triaptosis/results/06_external"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
set.seed(2026)

model <- readRDS(file.path(MODEL_DIR, "model.rds"))
beta <- model$beta; genes_all <- model$genes
cat("模型基因:", length(genes_all), "\n")

# ---------- 1. 读取 GSE22541 元数据 (已由 python 干净解析) ----------
cat("=== GSE22541 ===\n")
clin_meta <- fromJSON(file.path(EXT, "GSE22541_clin.json"))
sample_ids <- clin_meta$sample
is_primary <- clin_meta$primary
dfs <- clin_meta$dfs_months
cat("样本数:", length(sample_ids), " 原发灶:", sum(is_primary), "\n")

# ---------- 2. 表达矩阵 + 探针映射 ----------
cat("读取表达矩阵...\n")
gz <- gzfile(file.path(EXT, "GSE22541_matrix.txt.gz"), "rt")
lines <- readLines(gz, warn = FALSE); close(gz)
start <- which(grepl("!series_matrix_table_begin", lines))
end   <- which(grepl("!series_matrix_table_end", lines))
tab <- read.table(text = lines[(start + 1):(end - 1)], sep = "\t",
                  header = TRUE, row.names = 1, check.names = FALSE)
rownames(tab) <- sub('^"|"$', '', rownames(tab))
cat("矩阵:", dim(tab), "\n")

# 21 基因 -> 探针 (hgu133plus2.db)
suppressMessages(library(hgu133plus2.db))
probe2gene <- AnnotationDbi::select(hgu133plus2.db, keys = rownames(tab),
                                    columns = "SYMBOL", keytype = "PROBEID")
probe2gene <- probe2gene[!is.na(probe2gene$SYMBOL), ]
mat21 <- sapply(genes_all, function(g) {
  pr <- probe2gene$PROBEID[probe2gene$SYMBOL == g]
  if (length(pr) == 0) return(rep(NA, ncol(tab)))
  if (length(pr) == 1) return(as.numeric(tab[pr, ]))
  colMeans(as.matrix(tab[pr, ]))  # 多探针取均值
})
rownames(mat21) <- colnames(tab)
cat("21 基因覆盖:", sum(colSums(!is.na(mat21)) > 0), "/", length(genes_all), "\n")

# ---------- 3. TAS 打分 (队列内 z-score) ----------
z <- scale(mat21)
tas <- as.numeric(z %*% beta); names(tas) <- rownames(mat21)
med <- median(tas, na.rm = TRUE)
grp <- ifelse(tas > med, "High", "Low")

# ---------- 4. 原发灶 DFS 验证 (生存时间 = DFS) ----------
prim_time <- dfs
prim_keep <- is_primary & !is.na(prim_time) & !is.na(tas)
cat(sprintf("\n原发灶 DFS 验证: n=%d (High=%d, Low=%d)\n",
            sum(prim_keep), sum(grp[prim_keep] == "High"), sum(grp[prim_keep] == "Low")))
if (sum(prim_keep) >= 15) {
  d <- data.frame(dfs_time = prim_time[prim_keep], grp = factor(grp[prim_keep], levels = c("Low", "High")))
  wt <- wilcox.test(dfs_time ~ grp, data = d)
  cat(sprintf("DFS 比较: High 中位=%.1f 月, Low 中位=%.1f 月 | Wilcoxon P=%.4f\n",
              median(d$dfs_time[d$grp == "High"]), median(d$dfs_time[d$grp == "Low"]), wt$p.value))
  # KM (队列设计: 全部发生转移, 视 DFS 为事件时间, 删失=0)
  d$status <- 1
  lr <- survdiff(Surv(dfs_time, status) ~ grp, data = d)
  p_lr <- 1 - pchisq(lr$chisq, df = 1)
  cat(sprintf("DFS KM log-rank P=%.4f\n", p_lr))
  mcox <- tryCatch(coxph(Surv(dfs_time, status) ~ grp, data = d), error = function(e) NULL)
  if (!is.null(mcox)) {
    cat(sprintf("Cox HR(High vs Low)=%.3f (95%%CI %.3f-%.3f) P=%.4f\n",
        exp(coef(mcox)), exp(confint(mcox))[1], exp(confint(mcox))[2],
        summary(mcox)$coefficients[5]))
  }
  saveRDS(list(d = d, p_wilcox = wt$p.value, p_logrank = p_lr), file.path(OUT, "gse22541_dfs.rds"))
}

# ---------- 5. TCGA-KIRC 6:4 时间分割内部验证 ----------
cat("\n=== TCGA-KIRC 6:4 时间分割 ===\n")
kirc_expr <- fromJSON(file.path(DATA, "KIRC_expr.json"))
kirc_clin <- fromJSON(file.path(DATA, "KIRC_clin.json"))
kmat <- t(vapply(kirc_expr, function(x) {
  v <- rep(NA_real_, length(genes_all)); names(v) <- genes_all
  for (g in genes_all) if (!is.null(x[[g]])) v[g] <- as.numeric(x[[g]])
  v
}, FUN.VALUE = numeric(length(genes_all))))
colnames(kmat) <- genes_all
rownames(kmat) <- substr(names(kirc_expr), 1, 12)
kmat <- kmat[!duplicated(rownames(kmat)), , drop = FALSE]
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a)) b else a
os_status <- sapply(kirc_clin, function(x) as.character(x$OS_STATUS %||% ""))
os_months <- suppressWarnings(as.numeric(sapply(kirc_clin, function(x) x$OS_MONTHS %||% NA)))
dead <- grepl("DECEASED|1", toupper(os_status))
keep <- !is.na(os_months) & os_months > 0 & os_status != ""
common <- intersect(rownames(kmat), names(kirc_clin)[keep])
idx <- match(common, names(kirc_clin)[keep])
ksdf <- data.frame(time = os_months[keep][idx], status = as.numeric(dead[keep][idx]))
rownames(ksdf) <- common

set.seed(42)
train_idx <- sample(seq_len(nrow(ksdf)), size = round(0.6 * nrow(ksdf)))
test_idx <- setdiff(seq_len(nrow(ksdf)), train_idx)

# 训练集: 弹性网重拟合
suppressMessages(library(glmnet))
ztr <- scale(kmat[rownames(ksdf)[train_idx], , drop = FALSE])
cvf <- cv.glmnet(ztr, Surv(ksdf$time[train_idx], ksdf$status[train_idx]),
                 family = "cox", alpha = 0.5, nfolds = 10)
zte <- scale(kmat[rownames(ksdf)[test_idx], , drop = FALSE],
             center = attr(ztr, "scaled:center"), scale = attr(ztr, "scaled:scale"))
test_risk <- as.numeric(predict(cvf, newx = zte, s = "lambda.min"))

tmed <- median(test_risk)
tgrp <- factor(ifelse(test_risk > tmed, "High", "Low"), levels = c("Low", "High"))
tsdf <- data.frame(time = ksdf$time[test_idx], status = ksdf$status[test_idx], grp = tgrp)
cat(sprintf("测试集: n=%d (High=%d Low=%d, 死亡=%d)\n", nrow(tsdf),
            sum(tgrp == "High"), sum(tgrp == "Low"), sum(tsdf$status)))
lr2 <- survdiff(Surv(time, status) ~ grp, data = tsdf)
p2 <- 1 - pchisq(lr2$chisq, df = 1)
cat(sprintf("测试集 KM log-rank P=%.4e\n", p2))
m2 <- coxph(Surv(time, status) ~ grp, data = tsdf)
cat(sprintf("测试集 Cox HR(High vs Low)=%.3f (95%%CI %.3f-%.3f) P=%.4e\n",
    exp(coef(m2)), exp(confint(m2))[1], exp(confint(m2))[2], summary(m2)$coefficients[5]))
# 测试集时间 AUC (survivalROC)
if (requireNamespace("survivalROC", quietly = TRUE)) {
  library(survivalROC)
  for (t in c(36, 60)) {
    if (max(tsdf$time) > t) {
      roc <- survivalROC(Stime = tsdf$time, status = tsdf$status, marker = test_risk,
                         predict.time = t, method = "KM")
      cat(sprintf("测试集时间AUC(%d月)=%.3f\n", t, roc$AUC))
    }
  }
}
saveRDS(list(test = tsdf, p_logrank = p2, hr = exp(coef(m2))), file.path(OUT, "tcga_split.rds"))
cat("\n结果已保存:", OUT, "\n")
