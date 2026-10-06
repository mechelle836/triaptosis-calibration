#!/usr/bin/env Rscript
# Paper2-Triaptosis: KIRC 免疫浸润分析 (ssGSEA, Charoentong 28 免疫细胞)
# TAS 高 vs 低组免疫差异 + 与 TAS 相关性
suppressMessages({
  library(GSVA); library(jsonlite); library(pheatmap)
})

DATA <- "/Users/xiaolin/CodeBuddy/Claw/paper2_triaptosis/data"
MODEL_DIR <- "/Users/xiaolin/CodeBuddy/Claw/paper2_triaptosis/results/02_model/KIRC"
OUT <- "/Users/xiaolin/CodeBuddy/Claw/paper2_triaptosis/results/04_immune"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# 1. TAS 分组
sdf <- read.csv(file.path(MODEL_DIR, "score_table.csv"), row.names = 1)
cat("TAS 样本:", nrow(sdf), "\n")

# 2. 免疫基因表达
iexpr <- fromJSON(file.path(DATA, "KIRC_immune_expr.json"))
sigs <- fromJSON(file.path(DATA, "immune_signature.json"))
genes_all <- sort(unique(unlist(lapply(iexpr, names))))
imat <- t(vapply(iexpr, function(x) {
  v <- rep(NA_real_, length(genes_all)); names(v) <- genes_all
  for (g in genes_all) if (!is.null(x[[g]])) v[g] <- as.numeric(x[[g]])
  v
}, FUN.VALUE = numeric(length(genes_all))))
colnames(imat) <- genes_all
rownames(imat) <- substr(names(iexpr), 1, 12)
cat("免疫表达矩阵:", dim(imat), "\n")

common <- intersect(rownames(imat), rownames(sdf))
cat("对齐样本:", length(common), "\n")
imat2 <- imat[common, , drop = FALSE]
tas <- sdf[common, "TAS"]
grp <- sdf[common, "group"]

# 3. ssGSEA (GSVA >=1.50 新接口)
expr_mat <- t(imat2)  # 基因 x 样本
par <- ssgseaParam(exprData = expr_mat, geneSets = sigs)
ss <- gsva(param = par, verbose = FALSE)
cat("ssGSEA 分数矩阵:", dim(ss), "\n")

# 4. 高 vs 低组差异 (Wilcoxon) + 与 TAS 相关 (Spearman)
res <- data.frame(
  cell = rownames(ss),
  mean_high = NA, mean_low = NA,
  wilcox_p = NA,
  spearman_rho = NA, spearman_p = NA
)
for (i in seq_len(nrow(ss))) {
  v <- as.numeric(ss[i, common])
  res$mean_high[i] <- mean(v[grp == "High"], na.rm = TRUE)
  res$mean_low[i]  <- mean(v[grp == "Low"], na.rm = TRUE)
  res$wilcox_p[i]  <- tryCatch(wilcox.test(v[grp == "High"], v[grp == "Low"])$p.value,
                               error = function(e) NA)
  ct <- tryCatch(cor.test(tas, v, method = "spearman"), error = function(e) NULL)
  if (!is.null(ct)) { res$spearman_rho[i] <- ct$estimate; res$spearman_p[i] <- ct$p.value }
}
res$wilcox_fdr <- p.adjust(res$wilcox_p, method = "BH")
res$spearman_fdr <- p.adjust(res$spearman_p, method = "BH")
res <- res[order(res$wilcox_p), ]

write.csv(res, file.path(OUT, "immune_by_TAS_group.csv"), row.names = FALSE)
cat("\n=== 免疫细胞 TAS 高 vs 低 (Wilcoxon, 按 P 排序前10) ===\n")
print(head(res, 10), row.names = FALSE)
n_sig <- sum(res$wilcox_fdr < 0.05, na.rm = TRUE)
cat(sprintf("\n显著差异免疫细胞 (FDR<0.05): %d/%d\n", n_sig, nrow(res)))

# 热图 (高/低组平均)
hm <- sapply(c("High", "Low"), function(g) rowMeans(ss[, common[grp == g], drop = FALSE]))
pdf(file.path(OUT, "immune_heatmap.pdf"), width = 6, height = 10)
pheatmap(hm, cluster_cols = FALSE, main = "ssGSEA immune score by TAS group (KIRC)")
dev.off()
cat("热图已保存:", file.path(OUT, "immune_heatmap.pdf"), "\n")
