# =====================================================================
# 28_visium_spatial.R —— 空间转录组学 (Visium, GSE210041)
# TAS module score 空间投射 + TAS 高/低区域空间分布 + 邻域分析
# =====================================================================
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(dplyr); library(ggprism)
  library(showtext); library(patchwork)
})
ROOT <- "/Users/xiaolin/WorkBuddy/2026-09-04-20-00-08/paper2_triaptosis"
OUT  <- file.path(ROOT, "figures_v2")
RES  <- file.path(ROOT, "results", "16_visium")
dir.create(RES, showWarnings = FALSE, recursive = TRUE)
source(file.path(ROOT, "code", "lib", "viz_base.R"))
showtext_auto(FALSE)

mod <- readRDS(file.path(ROOT, "results/02_model/KIRC/model.rds"))
tas10 <- mod$genes_sel

# ---- 读 2 个 Visium 样本 ----------------------------------------------
v1 <- Load10X_Spatial("data/visium/", filename = "GSM6415705_GLMF1_filtered_feature_bc_matrix.h5",
                      slice = "GLMF1")
v2 <- Load10X_Spatial("data/visium/", filename = "GSM6415706_GLMF2_filtered_feature_bc_matrix.h5",
                      slice = "GLMF2")
cat("Tumour 1:", ncol(v1), "spots | Tumour 2:", ncol(v2), "spots\n")

# ---- 标准化 + TAS module score ------------------------------------------
process_one <- function(v, slice_name) {
  v <- SCTransform(v, assay = "Spatial", verbose = FALSE)
  v <- AddModuleScore(v, features = list(tas10), name = "TAS_score")
  v$slice <- slice_name
  v
}
v1 <- process_one(v1, "Tumour 1")
v2 <- process_one(v2, "Tumour 2")

# ---- 空间投射: TAS module score -----------------------------------------
p1 <- SpatialFeaturePlot(v1, features = "TAS_score1", alpha = c(0.1, 1),
                         pt.size.factor = 1.6) +
  scale_fill_gradient(low = "#F7F7F7", high = "#E64B35", name = "TAS") +
  labs(title = "Tumour 1 (GSM6415705)") +
  theme_cns(base_size = 7.5, base_family = "Arial") +
  theme(legend.position = "right", legend.title = element_text(size = 6),
        plot.title = element_text(size = 8.5, face = "bold", hjust = 0.5))
p2 <- SpatialFeaturePlot(v2, features = "TAS_score1", alpha = c(0.1, 1),
                         pt.size.factor = 1.6) +
  scale_fill_gradient(low = "#F7F7F7", high = "#E64B35", name = "TAS") +
  labs(title = "Tumour 2 (GSM6415706)") +
  theme_cns(base_size = 7.5, base_family = "Arial") +
  theme(legend.position = "right", legend.title = element_text(size = 6),
        plot.title = element_text(size = 8.5, face = "bold", hjust = 0.5))

# ---- TAS 高/低区域空间分布 (top/bottom 20% spots) -----------------------
v1$TAS_group <- ifelse(v1$TAS_score1 >= quantile(v1$TAS_score1, 0.8), "TAS-high",
                ifelse(v1$TAS_score1 <= quantile(v1$TAS_score1, 0.2), "TAS-low", "mid"))
v2$TAS_group <- ifelse(v2$TAS_score1 >= quantile(v2$TAS_score1, 0.8), "TAS-high",
                ifelse(v2$TAS_score1 <= quantile(v2$TAS_score1, 0.2), "TAS-low", "mid"))

p3 <- SpatialDimPlot(v1, group.by = "TAS_group", pt.size.factor = 1.6,
                     cols = c("TAS-high" = "#E64B35", "mid" = "grey88",
                              "TAS-low" = "#3C5488")) +
  labs(title = "Tumour 1 TAS groups") +
  theme_cns(base_size = 7.5, base_family = "Arial") +
  theme(legend.position = "right",
        plot.title = element_text(size = 8.5, face = "bold", hjust = 0.5))
p4 <- SpatialDimPlot(v2, group.by = "TAS_group", pt.size.factor = 1.6,
                     cols = c("TAS-high" = "#E64B35", "mid" = "grey88",
                              "TAS-low" = "#3C5488")) +
  labs(title = "Tumour 2 TAS groups") +
  theme_cns(base_size = 7.5, base_family = "Arial") +
  theme(legend.position = "right",
        plot.title = element_text(size = 8.5, face = "bold", hjust = 0.5))

# ---- 组合: 2 样本 x (投射 + 分组) ---------------------------------------
fig8 <- (p1 | p2) / (p3 | p4)
ggsave(file.path(OUT, "Fig8_visium_spatial.png"), fig8,
       width = 180, height = 195, units = "mm", dpi = 300)
showtext_auto()
ggsave(file.path(OUT, "Fig8_visium_spatial.pdf"), fig8,
       width = 180, height = 195, units = "mm", device = cairo_pdf)
showtext_auto(FALSE)

# 保存对象
saveRDS(list(v1 = v1, v2 = v2), file.path(RES, "visium_tas.rds"))
message("saved: Fig8_visium_spatial")