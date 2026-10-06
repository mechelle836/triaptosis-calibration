# =====================================================================
# 22_immune_checkpoint.R —— TAS × 免疫检查点/抑制分子关联
# 科学问题: 高 TAS(预后差)是否伴随检查点分子与免疫抑制分子特征?
# A: TAS × 105 免疫基因 Spearman rho 排序 (高亮检查点)
# B: 关键检查点/抑制分子 High vs Low 箱线
# 输出: Fig10_immune_checkpoint.png + results/13_checkpoint/
# =====================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(jsonlite)
  library(showtext); library(patchwork); library(ggprism)
})
ROOT <- "/Users/xiaolin/WorkBuddy/2026-09-04-20-00-08/paper2_triaptosis"
OUT  <- file.path(ROOT, "figures_v2")
RES  <- file.path(ROOT, "results", "13_checkpoint")
dir.create(RES, showWarnings = FALSE, recursive = TRUE)
source(file.path(ROOT, "code", "lib", "viz_base.R"))
showtext_auto(FALSE)   # macOS 系统 Arial 直渲, showtext 宽度计算导致 facet strip 裁字

score <- read.csv(file.path(ROOT, "results/02_model/KIRC/score_table.csv"),
                  row.names = 1)
imm <- fromJSON(file.path(ROOT, "data/KIRC_immune_expr.json"))

# ---- 1. 对齐样本 (score 无 -01 后缀, imm 有 -01) ----------------------
imm_ids <- substr(names(imm), 1, 12)   # 去 -01
names(imm) <- imm_ids
common <- intersect(rownames(score), imm_ids)
cat("common samples:", length(common), "\n")
tas <- score[common, "TAS"]
names(tas) <- common
imm_mat <- do.call(rbind, lapply(imm[common], function(x) unlist(x)))
colnames(imm_mat) <- names(imm[[common[1]]])

# ---- 2. TAS × 每基因 Spearman rho + FDR ------------------------------
res <- data.frame(gene = colnames(imm_mat), rho = NA, p = NA)
for (g in colnames(imm_mat)) {
  r <- cor.test(tas, imm_mat[, g], method = "spearman")
  res$rho[res$gene == g] <- r$estimate
  res$p[res$gene == g] <- r$p.value
}
res$fdr <- p.adjust(res$p, method = "BH")
res <- res %>% arrange(rho)
write.csv(res, file.path(RES, "tas_immune_gene_corr.csv"), row.names = FALSE)
cat("sig genes (FDR<0.05):", sum(res$fdr < 0.05), "\n")

checkpoints <- c("PDCD1", "CTLA4", "ICOS", "IL10", "FOXP3", "TNFRSF18",
                 "CD274", "HAVCR2", "LAG3")
res$is_cp <- res$gene %in% checkpoints
res$label <- ifelse(res$is_cp, res$gene, "")

# ---- 3. A: 关联排序条形 (top 15 正 + top 15 负) ----------------------
topA <- bind_rows(head(res, 15), tail(res, 15)) %>%
  mutate(gene = factor(gene, gene),
         sign = ifelse(rho < 0, "negative", "positive"))
pA <- ggplot(topA, aes(rho, gene, fill = is_cp)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = label), hjust = ifelse(topA$rho < 0, 1.1, -0.1),
            size = 2.4, fontface = "italic", family = "Arial") +
  geom_vline(xintercept = 0, colour = "grey30", linewidth = 0.4) +
  scale_fill_manual(values = c(`FALSE` = "#BDD7E7", `TRUE` = "#E64B35"),
                    guide = "none") +
  labs(x = "Spearman rho (TAS vs gene)", y = NULL,
       title = "TAS \u00d7 immune gene correlations") +
  theme_prism(base_size = 8, base_family = "Arial") +
  theme(axis.text.y = element_text(size = 5.5, face = "italic"),
        axis.text.x = element_text(size = 6.5),
        plot.title = element_text(size = 9.5, face = "bold", hjust = 0.5),
        plot.margin = margin(4, 2, 2, 2, "mm"))

# ---- 4. B: 关键检查点/抑制分子 High vs Low 箱线 ----------------------
key_genes <- intersect(c("PDCD1", "CTLA4", "FOXP3", "GZMB"),
                       colnames(imm_mat))
dfb <- data.frame(TAS = tas, group = score[common, "group"],
                  imm_mat[, key_genes]) %>%
  pivot_longer(-c(TAS, group), names_to = "gene", values_to = "expr") %>%
  mutate(expr = log1p(expr),                      # TPM 跨数量级, 标准 log1p
         group = factor(group, c("Low", "High")),
         gene = factor(gene, key_genes))
# 组间差异 (Wilcoxon)
pval <- dfb %>% group_by(gene) %>%
  summarise(p = wilcox.test(expr[group == "Low"], expr[group == "High"])$p.value) %>%
  mutate(lab = sprintf("P = %.2g", p))
pB <- ggplot(dfb, aes(group, expr, fill = group)) +
  geom_boxplot(width = 0.55, outlier.size = 0.4, linewidth = 0.3,
               show.legend = FALSE) +
  scale_fill_manual(values = cols_risk, guide = "none") +
  facet_wrap(~gene, scales = "free_y", nrow = 2,
             labeller = labeller(gene = c(PDCD1 = "PD-1", CTLA4 = "CTLA-4",
                                          FOXP3 = "FOXP3", GZMB = "GZMB"))) +
  geom_text(data = pval, aes(x = 1.5, y = Inf, label = lab),
            vjust = 2.4, size = 2.0, inherit.aes = FALSE, family = "Arial") +
  labs(x = NULL, y = "Expression (log1p TPM)",
       title = "Checkpoint molecules by TAS group") +
  theme_prism(base_size = 8, base_family = "Arial") +
  theme(axis.text.x = element_text(size = 7),
        strip.text = element_text(size = 6.5, face = "bold.italic"),
        axis.title.y = element_text(size = 7.5),
        plot.title = element_text(size = 9.5, face = "bold", hjust = 0.5),
        plot.margin = margin(4, 2, 2, 2, "mm"))

fig10 <- pA | pB + plot_layout(widths = c(1, 1.9))
ggsave(file.path(OUT, "Fig10_immune_checkpoint.png"), fig10,
       width = 210, height = 155, units = "mm", dpi = 300)
ggsave(file.path(OUT, "Fig10_immune_checkpoint.pdf"), fig10,
       width = 210, height = 155, units = "mm", device = cairo_pdf,
       )
message("saved: Fig10_immune_checkpoint")

# 关键结果输出
cat("\n=== 检查点/抑制分子关联 ===\n")
print(res %>% filter(is_cp) %>% select(gene, rho, p, fdr), row.names = FALSE)