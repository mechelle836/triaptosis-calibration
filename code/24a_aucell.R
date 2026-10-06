# =====================================================================
# 24a_aucell.R —— AUCell triaptosis 通路活性 (标准 AUCell 算法)
# 与 module score 互补: AUCell 是标准通路活性评分 (基于表达排序)
# 输出: results/15_aucell/ + Fig12A (通路活性 UMAP + 细胞类型比较)
# =====================================================================
suppressPackageStartupMessages({
  library(Seurat); library(AUCell); library(Matrix)
  library(ggplot2); library(dplyr); library(showtext)
  library(patchwork)
})
ROOT <- "/Users/xiaolin/WorkBuddy/2026-09-04-20-00-08/paper2_triaptosis"
OUT  <- file.path(ROOT, "figures_v2")
RES  <- file.path(ROOT, "results", "15_aucell")
dir.create(RES, showWarnings = FALSE, recursive = TRUE)
source(file.path(ROOT, "code", "lib", "viz_base.R"))
showtext_auto(FALSE)
font_add("Arial", regular = "/System/Library/Fonts/Supplemental/Arial.ttf")
showtext_opts(dpi = 300); showtext_auto()

seu <- readRDS(file.path(ROOT, "results/11_scrna/seurat_kirc_gse159115.rds"))
mod <- readRDS(file.path(ROOT, "results/02_model/KIRC/model.rds"))
tas10 <- intersect(mod$genes_sel, rownames(seu))
cat("TAS genes in scRNA:", length(tas10), "\n")

# ---- 1. AUCell 通路活性 ----------------------------------------------
expr <- GetAssayData(seu, layer = "counts")
genesets <- list(triaptosis_TAS = tas10)   # AUCell 接受命名 list
cells_rankings <- AUCell_buildRankings(expr, plotStats = FALSE)
cells_AUC <- AUCell_calcAUC(genesets, cells_rankings, aucMaxRank = nrow(cells_rankings) * 0.05)
auc_vec <- as.numeric(getAUC(cells_AUC)["triaptosis_TAS", ])
seu$TAS_AUC <- auc_vec
saveRDS(auc_vec, file.path(RES, "tas_aucell_scores.rds"))

# ---- 2. 与 module score 一致性 (sanity check) ------------------------
cor_m <- cor(seu$TAS_AUC, seu$TAS_score1, method = "spearman")
cat("AUCell vs module score Spearman:", round(cor_m, 3), "\n")

# ---- 3. 图: A UMAP 投射 + B 细胞类型比较 -----------------------------
pal6 <- c("Tumor cells" = "#E64B35", "Macrophages" = "#F39B7F",
          "T cells" = "#3C5488", "B/Plasma/Mast" = "#8491B4",
          "Endothelial" = "#00A087", "vSMC/Pericyte" = "#B09C85")
umap_df <- data.frame(UMAP1 = Embeddings(seu, "umap")[, 1],
                      UMAP2 = Embeddings(seu, "umap")[, 2],
                      ctype = seu$ctype, AUC = seu$TAS_AUC)
pA <- ggplot(umap_df, aes(UMAP1, UMAP2, colour = AUC)) +
  geom_point(size = 0.25, alpha = 0.7) +
  scale_colour_gradient2(low = "#4575B4", mid = "#F7F7F7", high = "#D73027",
                         midpoint = median(umap_df$AUC), name = "TAS AUCell\nactivity") +
  labs(title = "Triaptosis pathway activity (AUCell)", x = "UMAP-1", y = "UMAP-2") +
  theme_classic(base_size = 8, base_family = "Arial") +
  theme(legend.position = c(0.02, 0.02), legend.justification = c(0, 0),
        legend.title = element_text(size = 6), legend.text = element_text(size = 5.5),
        legend.key.height = unit(3.2, "mm"),
        legend.background = element_rect(fill = "white", colour = NA, alpha = 0.6),
        plot.title = element_text(size = 9.5, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 6), axis.title = element_text(size = 7.5))

pB <- ggplot(umap_df, aes(ctype, AUC, fill = ctype)) +
  geom_violin(scale = "width", linewidth = 0.3, show.legend = FALSE, alpha = 0.85) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", linewidth = 0.3) +
  scale_fill_manual(values = pal6) +
  labs(x = NULL, y = "Triaptosis pathway AUCell activity",
       title = "Pathway activity by cell type") +
  theme_prism(base_size = 8, base_family = "Arial") +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 7),
        axis.title.y = element_text(size = 7.5),
        plot.title = element_text(size = 9.5, face = "bold", hjust = 0.5))

# 肿瘤 vs 微环境显著性
tumor_auc <- umap_df$AUC[umap_df$ctype == "Tumor cells"]
env_auc <- umap_df$AUC[umap_df$ctype != "Tumor cells"]
wt <- wilcox.test(tumor_auc, env_auc)
cat("Tumor vs microenvironment AUCell: P =", format(wt$p.value, scientific = TRUE),
    "| Tumor median", round(median(tumor_auc), 3),
    "| Env median", round(median(env_auc), 3), "\n")

fig12a <- pA | pB + plot_layout(widths = c(1, 1.3))
ggsave(file.path(OUT, "Fig12A_aucell.png"), fig12a,
       width = 200, height = 95, units = "mm", dpi = 300)
message("saved: Fig12A_aucell")