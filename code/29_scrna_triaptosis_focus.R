# =====================================================================
# 29_scrna_triaptosis_focus.R —— 改进单细胞分析 (学习 Dai 2026)
# 识别 triaptosis 高活性细胞类型 + 分化阶段分析
# =====================================================================
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(dplyr); library(tidyr)
  library(showtext); library(patchwork)
})
ROOT <- "/Users/xiaolin/WorkBuddy/2026-09-04-20-00-08/paper2_triaptosis"
OUT  <- file.path(ROOT, "figures_v2")
RES  <- file.path(ROOT, "results", "17_scrna_focus")
dir.create(RES, showWarnings = FALSE, recursive = TRUE)
source(file.path(ROOT, "code", "lib", "viz_base.R"))
showtext_auto(FALSE)
pal6 <- c("Tumor cells" = "#E64B35", "Macrophages" = "#F39B7F",
          "T cells" = "#3C5488", "B/Plasma/Mast" = "#8491B4",
          "Endothelial" = "#00A087", "vSMC/Pericyte" = "#B09C85")

seu <- readRDS(file.path(ROOT, "results/11_scrna/seurat_kirc_gse159115.rds"))
mod <- readRDS(file.path(ROOT, "results/02_model/KIRC/model.rds"))
tas10 <- mod$genes_sel

# ---- 1. 识别 triaptosis 高活性细胞类型 (TAS_score 最高的细胞类型) --------
df <- data.frame(ctype = seu$ctype, TAS = seu$TAS_score1)
ctype_tas <- df %>% group_by(ctype) %>%
  summarise(mean_TAS = mean(TAS), median_TAS = median(TAS),
            sd_TAS = sd(TAS), n = n(), .groups = "drop") %>%
  arrange(-mean_TAS)
write.csv(ctype_tas, file.path(RES, "triaptosis_high_activity_celltypes.csv"),
          row.names = FALSE)
cat("triaptosis 高活性细胞类型:\n"); print(ctype_tas)

# ---- 2. TAS_score 按细胞类型分布 (小提琴图) ------------------------------
pA <- ggplot(df, aes(reorder(ctype, TAS, median), TAS, fill = ctype)) +
  geom_violin(scale = "width", alpha = 0.8, show.legend = FALSE) +
  geom_boxplot(width = 0.15, outlier.size = 0.3, show.legend = FALSE) +
  scale_fill_manual(values = pal6) +
  labs(x = NULL, y = "TAS module score",
       title = "Triaptosis activity by cell type") +
  theme_cns(base_size = 7.5, base_family = "Arial") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 7),
        plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
        plot.margin = margin(2, 2, 2, 2, "mm"))

# ---- 3. TAS 高活性细胞的分化阶段 (Tumor cells 的 TAS 分布) ----------------
# 分析 Tumor cells 的 TAS 分布 (是否有亚群差异)
tumor_df <- df %>% filter(ctype == "Tumor cells")
pB <- ggplot(tumor_df, aes(TAS)) +
  geom_histogram(bins = 40, fill = "#E64B35", colour = "white", linewidth = 0.2) +
  geom_vline(xintercept = median(tumor_df$TAS), colour = "#3C5488",
             linewidth = 0.6, linetype = "dashed") +
  annotate("text", x = median(tumor_df$TAS), y = Inf, hjust = -0.1, vjust = 1.3,
           label = sprintf("median = %.2f", median(tumor_df$TAS)),
           size = 2.2, colour = "#3C5488", family = "Arial") +
  labs(x = "TAS module score", y = "Tumor cells count",
       title = "TAS distribution in tumor cells") +
  theme_cns(base_size = 7.5, base_family = "Arial") +
  theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
        plot.margin = margin(2, 2, 2, 2, "mm"))

# ---- 4. TAS 高活性 vs 低活性 Tumor cells 的差异基因 (top markers) ---------
seu$TAS_group_tumor <- ifelse(seu$ctype == "Tumor cells" &
                               seu$TAS_score1 >= quantile(seu$TAS_score1[seu$ctype == "Tumor cells"], 0.75), "TAS-high",
                        ifelse(seu$ctype == "Tumor cells" &
                               seu$TAS_score1 <= quantile(seu$TAS_score1[seu$ctype == "Tumor cells"], 0.25), "TAS-low", "mid"))
Idents(seu) <- "TAS_group_tumor"
tumor_deg <- FindMarkers(seu, ident.1 = "TAS-high", ident.2 = "TAS-low",
                         min.pct = 0.1, logfc.threshold = 0.25, verbose = FALSE)
tumor_deg$gene <- rownames(tumor_deg)
write.csv(tumor_deg, file.path(RES, "tas_high_low_tumor_deg.csv"), row.names = FALSE)
top_deg <- tumor_deg %>% arrange(p_val_adj) %>% head(10)
cat("TAS-high vs TAS-low Tumor cells top markers:\n")
print(top_deg[, c("gene", "avg_log2FC", "p_val_adj")])

# ---- 组合: A 小提琴 + B 分布 ----------------------------------------------
fig <- (pA | pB) + plot_layout(widths = c(1.4, 1))
ggsave(file.path(OUT, "Fig9_scrna_triaptosis_focus.png"), fig,
       width = 180, height = 120, units = "mm", dpi = 300)
showtext_auto()
ggsave(file.path(OUT, "Fig9_scrna_triaptosis_focus.pdf"), fig,
       width = 180, height = 120, units = "mm", device = cairo_pdf)
showtext_auto(FALSE)
message("saved: Fig9_scrna_triaptosis_focus")