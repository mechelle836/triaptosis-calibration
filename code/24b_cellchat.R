# =====================================================================
# 24b_cellchat.R —— CellChat 细胞通讯 (TAS 高 vs 低肿瘤细胞)
# 核心问题: 高 TAS 肿瘤细胞如何通过配体-受体与免疫细胞通讯 (免疫抑制机制)
# 输出: results/15_cellchat/ + Fig12B (通讯网络)
# =====================================================================
suppressPackageStartupMessages({
  library(Seurat); library(CellChat); library(dplyr); library(tidyr)
  library(ggplot2); library(showtext); library(patchwork)
})
ROOT <- "/Users/xiaolin/WorkBuddy/2026-09-04-20-00-08/paper2_triaptosis"
OUT  <- file.path(ROOT, "figures_v2")
RES  <- file.path(ROOT, "results", "15_cellchat")
dir.create(RES, showWarnings = FALSE, recursive = TRUE)
source(file.path(ROOT, "code", "lib", "viz_base.R"))
showtext_auto(FALSE)
font_add("Arial", regular = "/System/Library/Fonts/Supplemental/Arial.ttf")
showtext_opts(dpi = 300); showtext_auto()

seu <- readRDS(file.path(ROOT, "results/11_scrna/seurat_kirc_gse159115.rds"))

# ---- 1. 下采样 (保留比例, 控制计算量) ----------------------------------
set.seed(42)
n_keep <- 8000
keep <- sample(colnames(seu), min(n_keep, ncol(seu)))
seu_sub <- seu[, keep]

# ---- 2. 分组标签: Tumor cells 按 TAS 中位数分高/低 ---------------------
tas_med <- median(seu_sub$TAS_score1)
labels <- as.character(seu_sub$ctype)
labels[labels == "Tumor cells" & seu_sub$TAS_score1 >= tas_med] <- "Tumor_TAShigh"
labels[labels == "Tumor cells" & seu_sub$TAS_score1 < tas_med]  <- "Tumor_TASlow"
labels <- factor(labels)
seu_sub$labels <- labels
cat("cell labels:\n"); print(table(labels))

# ---- 3. CellChat 标准流程 --------------------------------------------
data.input <- GetAssayData(seu_sub, layer = "data")   # 细胞 x 基因 normalized
meta <- data.frame(labels = labels, row.names = colnames(seu_sub))

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat, do.fast = FALSE)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
saveRDS(cellchat, file.path(RES, "cellchat_TAS.rds"))
cat("CellChat done, interactions:", nrow(cellchat@net$count), "\n")

# ---- 4. 通讯强度对比: Tumor_TAShigh vs Tumor_TASlow outgoing ----------
net_count <- cellchat@net$count   # 信号数量矩阵
net_weight <- cellchat@net$weight # 信号强度矩阵
grp <- rownames(net_count)
if ("Tumor_TAShigh" %in% grp && "Tumor_TASlow" %in% grp) {
  high_out <- net_weight["Tumor_TAShigh", setdiff(grp, "Tumor_TAShigh")]
  low_out  <- net_weight["Tumor_TASlow",  setdiff(grp, "Tumor_TASlow")]
  cmp <- data.frame(target = names(high_out),
                    high = as.numeric(high_out),
                    low  = as.numeric(low_out)) %>%
    filter(target %in% c("Macrophages", "T cells", "B/Plasma/Mast",
                         "Endothelial", "vSMC/Pericyte")) %>%
    pivot_longer(-target, names_to = "grp", values_to = "strength")
  write.csv(cmp, file.path(RES, "tas_high_low_outgoing.csv"), row.names = FALSE)
  pC <- ggplot(cmp, aes(target, strength, fill = grp)) +
    geom_col(position = position_dodge(0.7), width = 0.6) +
    scale_fill_manual(values = c(high = "#E64B35", low = "#3C5488"),
                      name = NULL) +
    labs(x = NULL, y = "Outgoing signal strength",
         title = "TAS-high vs TAS-low tumour outgoing signals") +
    theme_prism(base_size = 8, base_family = "Arial") +
    theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 7),
          legend.position = "top",
          plot.title = element_text(size = 9.5, face = "bold", hjust = 0.5))
} else {
  pC <- ggplot() + theme_void()
}

# ---- 5. 通讯网络图 (netVisual circle) ---------------------------------
grDevices::png(file.path(OUT, "Fig12B_cellchat_network.png"),
               width = 1600, height = 1600, res = 300)
netVisual_circle(net_weight, weight.scale = TRUE,
                 vertex.label.cex = 0.7, title.name = "Cell-cell communication")
grDevices::dev.off()

# ---- 6. 汇总图: network + 差异对比 ------------------------------------
fig12b <- pC
ggsave(file.path(OUT, "Fig12B_cellchat_diff.png"), fig12b,
       width = 170, height = 110, units = "mm", dpi = 300)
message("saved: Fig12B_cellchat")