# =====================================================================
# 27e_fig6_composite.R —— Fig6 单细胞定位（投稿版）
# A UMAP + 细胞类型图例   B TAS 签名（只留一张）
# C 10 基因点图           D 肿瘤 vs 微环境平均表达（不用伪倍数）
# E CellChat 高低组差异最大的 8 条配体-受体
# AUCell UMAP 与按靶细胞汇总的 outgoing 条形进补充图。
# 输出: figures_pub/Fig7_single_cell.{png,pdf}（稿件 Figure 7）
# =====================================================================
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(dplyr); library(tidyr)
  library(CellChat); library(showtext); library(patchwork)
})
ROOT <- "/Users/xiaolin/WorkBuddy/2026-09-04-20-00-08/paper2_triaptosis"
OUT  <- file.path(ROOT, "figures_pub")
RES  <- file.path(ROOT, "results")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
source(file.path(ROOT, "code", "lib", "viz_base.R"))

pal6 <- c("Tumor cells" = "#E64B35", "Macrophages" = "#F39B7F",
          "T cells" = "#3C5488", "B/Plasma/Mast" = "#8491B4",
          "Endothelial" = "#00A087", "vSMC/Pericyte" = "#B09C85")

seu <- readRDS(file.path(RES, "11_scrna/seurat_kirc_gse159115.rds"))
auc_vec <- readRDS(file.path(RES, "15_aucell/tas_aucell_scores.rds"))
dp <- read.csv(file.path(RES, "11_scrna/scrna_dotplot_data.csv"))
mod <- readRDS(file.path(RES, "02_model/KIRC/model.rds"))
gsel <- mod$genes_sel

umap_df <- data.frame(UMAP1 = Embeddings(seu, "umap")[, 1],
                      UMAP2 = Embeddings(seu, "umap")[, 2],
                      ctype = seu$ctype, TAS = seu$TAS_score1,
                      AUC = auc_vec)
n_cells <- nrow(umap_df)

theme_umap <- function() {
  theme_classic(base_size = 7.5, base_family = "Arial") +
    theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 6.5),
          plot.margin = margin(2, 2, 2, 2, "mm"))
}

# ---- A: cell types, legend on ------------------------------------------------
pA <- ggplot(umap_df, aes(UMAP1, UMAP2, colour = ctype)) +
  geom_point(size = 0.15, alpha = 0.65) +
  scale_colour_manual(values = pal6, name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 2.2, alpha = 1))) +
  labs(x = "UMAP-1", y = "UMAP-2", tag = "A",
       title = sprintf("GSE159115 (n=%s cells)", format(n_cells, big.mark = ","))) +
  theme_umap() +
  theme(legend.position = "right")

# ---- B: one TAS UMAP ---------------------------------------------------
# TAS 是有符号的 module score, 0 有意义, 所以用以 0 为中点的双向色标
# (与 Fig8 的空间图同一套配色), 单向 white->red 会把"低于平均"的细胞画成白色。
tas_lim <- max(abs(quantile(umap_df$TAS, c(0.02, 0.98), na.rm = TRUE)))
pB <- ggplot(umap_df, aes(UMAP1, UMAP2, colour = TAS)) +
  geom_point(size = 0.15, alpha = 0.7) +
  scale_colour_gradient2(low = "#3C5488", mid = "#F7F7F7", high = "#E64B35",
                         midpoint = 0, name = "TAS",
                         limits = c(-tas_lim, tas_lim),
                         oob = scales::squish, guide = guide_cbar()) +
  labs(x = "UMAP-1", y = "UMAP-2", title = "TAS signature", tag = "B") +
  theme_umap() +
  theme(legend.position = "right")

# ---- C: dot plot -------------------------------------------------------
dp2 <- dp %>%
  mutate(gene = factor(gene, rev(gsel)),
         ctype = factor(ctype, names(pal6)))
# 色条与尺寸图例并排时总宽超出版面, "% expressed" 的键被裁掉,
# 只剩一个孤立小点; 改为上下两行并给出显式的尺寸刻度。
# x 轴标签与 A 的图例用同一套名称, 原先 "B/Mast"/"vSMC" 与 A 不一致。
pC <- ggplot(dp2, aes(ctype, gene)) +
  geom_point(aes(size = pct_expr, colour = scaled)) +
  scale_colour_gradient(low = "#F7F7F7", high = "#B2182B", name = "Scaled expr",
                        guide = guide_cbar(order = 1)) +
  scale_size_continuous(range = c(0.6, 4.2), name = "% expressing cells",
                        breaks = c(1, 10, 20, 30), limits = c(0, 40),
                        guide = guide_legend(order = 2, nrow = 1)) +
  scale_x_discrete(labels = c("Tumor cells" = "Tumor",
                              "Macrophages" = "Macrophage",
                              "T cells" = "T cell",
                              "B/Plasma/Mast" = "B/Plasma/Mast",
                              "Endothelial" = "Endothelial",
                              "vSMC/Pericyte" = "vSMC/Pericyte")) +
  labs(x = NULL, y = NULL, title = "10-gene expression", tag = "C") +
  theme_cns(base_size = 7.5) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 6),
        axis.text.y = element_text(face = "italic", size = 6.5),
        legend.position = "bottom", legend.box = "vertical",
        legend.box.just = "left",
        # 图例标题压到键的上方, 否则 "标题 + 键" 一整行比面板还宽, 右端的键会被裁掉
        legend.title.position = "top",
        legend.spacing.y = unit(0.6, "mm"),
        legend.margin = margin(0.5, 0, 0, 0, "mm"),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 5.6),
        legend.key.height = unit(2.6, "mm"),
        legend.key.width = unit(3.4, "mm"),
        plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
        plot.margin = margin(2, 2, 2, 2, "mm"))

# ---- D: tumor vs microenvironment mean (no fold axis) -----------------
pair <- dp %>%
  group_by(gene) %>%
  summarise(
    tumor = mean_expr[ctype == "Tumor cells"],
    tumor_pct = pct_expr[ctype == "Tumor cells"],
    env = weighted.mean(mean_expr[ctype != "Tumor cells"],
                        n_cells[ctype != "Tumor cells"]),
    env_pct = weighted.mean(pct_expr[ctype != "Tumor cells"],
                            n_cells[ctype != "Tumor cells"]),
    .groups = "drop") %>%
  mutate(gene = factor(gene, rev(gsel)),
         note = ifelse(env == 0, "not detected outside tumor", ""),
         delta = tumor - env,
         dir = ifelse(delta >= 0, "Higher in tumor", "Higher in microenvironment"))
write.csv(pair, file.path(OUT, "fig6_tumor_vs_env_means.csv"), row.names = FALSE)
pD <- ggplot(pair, aes(delta, gene, fill = dir)) +
  geom_col(width = 0.62, show.legend = TRUE) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey30") +
  geom_text(data = pair %>% filter(note != ""),
            aes(x = delta, y = gene, label = note), inherit.aes = FALSE,
            hjust = 0, nudge_x = 0.008, size = 1.9, colour = "grey35",
            family = "Arial") +
  scale_fill_manual(values = c("Higher in tumor" = "#3C5488",
                               "Higher in microenvironment" = "#8491B4"),
                    name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.08, 0.28))) +
  labs(x = "Mean expression, tumor \u2212 microenvironment", y = NULL, tag = "D",
       title = "Tumor vs microenvironment") +
  theme_cns(base_size = 7.5) +
  theme(axis.text.y = element_text(face = "italic", size = 6.5),
        legend.position = "bottom",
        plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
        plot.margin = margin(2, 2, 2, 2, "mm"))

# ---- E: top 8 ligand-receptor deltas ----------------------------------
cellchat <- readRDS(file.path(RES, "15_cellchat/cellchat_TAS.rds"))
prob <- cellchat@net$prob
targets <- setdiff(dimnames(prob)[[1]], c("Tumor_TAShigh", "Tumor_TASlow"))
lr_rows <- lapply(targets, function(tg) {
  hv <- prob["Tumor_TAShigh", tg, ]
  lv <- prob["Tumor_TASlow", tg, ]
  keep <- (hv + lv) > 0
  if (!any(keep)) return(NULL)
  data.frame(target = tg, lr = names(hv)[keep],
             high = as.numeric(hv[keep]), low = as.numeric(lv[keep]),
             stringsAsFactors = FALSE)
})
lr <- bind_rows(lr_rows) %>%
  mutate(delta = high - low) %>%
  arrange(desc(abs(delta)))
top8 <- lr %>% slice_head(n = 8) %>%
  mutate(lab = sprintf("%s \u2192 %s", lr, target),
         lab = factor(lab, rev(lab)),
         direction = ifelse(delta >= 0, "Higher in TAS-high", "Higher in TAS-low"))
write.csv(top8, file.path(OUT, "fig6_cellchat_top8.csv"), row.names = FALSE)
pE <- ggplot(top8, aes(delta, lab, fill = direction)) +
  geom_col(width = 0.62, show.legend = TRUE) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey30") +
  scale_fill_manual(values = c("Higher in TAS-high" = "#E64B35",
                               "Higher in TAS-low" = "#3C5488"),
                    name = NULL) +
  labs(x = "Communication probability (TAS-high \u2212 TAS-low)", y = NULL, tag = "E",
       title = "Largest outgoing differences") +
  theme_cns(base_size = 7.5) +
  theme(axis.text.y = element_text(size = 5.8),
        legend.position = "bottom",
        plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
        plot.margin = margin(2, 2, 2, 2, "mm"))

fig6 <- (pA | pB | pC) / (pD | pE) +
  plot_layout(heights = c(1.05, 1)) &
  theme(plot.tag = element_text(face = "bold", size = 11))
save_pub(fig6, "Fig7_single_cell", OUT, w = 180, h = 185)

# ---- supplements: AUCell UMAP; outgoing strength by target ------------
pS <- ggplot(umap_df, aes(UMAP1, UMAP2, colour = AUC)) +
  geom_point(size = 0.15, alpha = 0.7) +
  scale_colour_gradient(low = "#F7F7F7", high = "#B2182B", name = "AUCell",
                        limits = c(quantile(umap_df$AUC, 0.05, na.rm = TRUE),
                                   quantile(umap_df$AUC, 0.95, na.rm = TRUE)),
                        oob = scales::squish, guide = guide_cbar()) +
  labs(x = "UMAP-1", y = "UMAP-2",
       title = "Triaptosis gene-set activity (AUCell)") +
  theme_umap()
save_pub(pS, "FigS_aucell_umap", OUT, w = 90, h = 80)

outg <- read.csv(file.path(RES, "15_cellchat/tas_high_low_outgoing.csv"))
outg$grp <- factor(outg$grp, c("low", "high"))
pOut <- ggplot(outg, aes(target, strength, fill = grp)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.62) +
  scale_fill_manual(values = c(low = "#3C5488", high = "#E64B35"),
                    labels = c(low = "TAS-low tumor", high = "TAS-high tumor"),
                    name = NULL) +
  labs(x = NULL, y = "Outgoing signal strength",
       title = "Aggregated tumor outgoing signals") +
  theme_cns(base_size = 8) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "top",
        plot.title = element_text(size = 9, face = "bold", hjust = 0.5))
save_pub(pOut, "FigS_cellchat_outgoing", OUT, w = 120, h = 80)

message("Fig6 cells: ", n_cells)
message("SH3GL3 note: ", pair$note[pair$gene == "SH3GL3"])
message("top8 LR:\n", paste(sprintf("%s  d=%.4f", top8$lab, top8$delta), collapse = "\n"))
