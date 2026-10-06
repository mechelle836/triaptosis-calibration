# =====================================================================
# 27d_fig7_composite.R —— Fig7 边界（投稿版）
# 主图只留零模型直方图和 IMmotion150。
# OR 与 HR 分面，轴标题为 effect size。
# GO / 通路网络、药敏进补充图。网络图不画布局坐标。
# =====================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr)
  library(tidygraph); library(ggraph)
  library(showtext); library(patchwork)
})
ROOT <- "/Users/xiaolin/WorkBuddy/2026-09-04-20-00-08/paper2_triaptosis"
OUT  <- file.path(ROOT, "figures_pub")
RES  <- file.path(ROOT, "results")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
source(file.path(ROOT, "code", "lib", "viz_base.R"))

nv <- read.csv(file.path(RES, "05_null/triaptosis_random_null.csv"))
ici <- read.csv(file.path(RES, "03_ici/ici_null_forest.csv"))
pcd_all <- read.csv(file.path(RES, "05_null/pcd_benchmark_same_protocol.csv"))
# 头对头比较必须用同一协议: 10 seed 均值, 与铁死亡/凋亡/双硫死亡同一套外折。
tri <- pcd_all %>% filter(gene_set == "Triaptosis")
pcd <- pcd_all %>% filter(gene_set != "Triaptosis")
par <- read.csv(file.path(RES, "02_model/parallel_cancer_hr.csv"))
tas_c <- tri$cv_nested
pct <- tri$percentile
z <- tri$z

# 标签带放在柱子上方, 四个基因集各占一行, 谁和谁重合一目了然。
bins <- 18
brk <- seq(min(nv$cv_nested), max(nv$cv_nested), length.out = bins + 1)
ymax <- max(table(cut(nv$cv_nested, brk, include.lowest = TRUE)))
# 四个基因集的 C-index 都落在轴的右段, 标签统一朝左写才不会冲出面板;
# 百分位并进 Triaptosis 自己的标签, 省掉一个会压住柱子的独立注释块。
marks <- bind_rows(
  data.frame(set = "Triaptosis", value = tas_c, col = "#E64B35", lvl = 1.50,
             lab = sprintf("Triaptosis %.3f (%.1f percentile)", tas_c, pct)),
  pcd %>% transmute(set = gene_set, value = cv_nested,
                    col = c(Ferroptosis = "#00A087", Apoptosis = "#3C5488",
                            Disulfidptosis = "#8491B4")[gene_set],
                    lvl = c(Ferroptosis = 1.34, Apoptosis = 1.18,
                            Disulfidptosis = 1.02)[gene_set],
                    lab = sprintf("%s %.3f", gene_set, cv_nested))
) %>% mutate(y = lvl * ymax, xlab = value - 0.004)
pA <- ggplot(nv, aes(cv_nested)) +
  geom_histogram(bins = bins, fill = "#D9D9D9", colour = "white", linewidth = 0.2) +
  geom_segment(data = marks, aes(x = value, xend = value, y = 0, yend = y,
                                 colour = set, linetype = set == "Triaptosis"),
               inherit.aes = FALSE, linewidth = 0.4, show.legend = FALSE) +
  geom_point(data = marks, aes(x = value, y = y, colour = set), inherit.aes = FALSE,
             size = 1.3, show.legend = FALSE) +
  geom_text(data = marks, aes(x = xlab, y = y, label = lab, colour = set),
            inherit.aes = FALSE, size = 2.0, family = "Arial", hjust = 1, vjust = 0.4,
            show.legend = FALSE) +
  scale_colour_manual(values = setNames(marks$col, marks$set)) +
  scale_linetype_manual(values = c(`TRUE` = "solid", `FALSE` = "22")) +
  scale_y_continuous(limits = c(0, 1.6 * ymax), expand = expansion(mult = c(0, 0.02))) +
  scale_x_continuous(expand = expansion(mult = c(0.04, 0.03))) +
  labs(x = "Nested CV C-index of random 21-gene sets", y = "Count", tag = "A",
       title = sprintf("Null audit (n=%d random sets)", nrow(nv))) +
  theme_cns(base_size = 8) +
  theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5))

# HR 与 OR 都是零假设为 1 的比值, 可以共用一条对数轴; 原来分两个 facet,
# 上面那个只有 PFS 一行, 大半是空白, 而且两个 facet 的 x 刻度不同容易误读。
# 标签在点上方居中时, 靠近面板边缘的点会被裁掉; 按点在轴上的相对位置
# 自动决定 hjust (最左 = 0 往右写, 最右 = 1 往左写)。
hjust_in_panel <- function(v, lo, hi) {
  pmin(pmax((log10(v) - log10(lo)) / (log10(hi) - log10(lo)), 0), 1)
}
ici_lim <- c(0.3, 3.2)
ici2 <- ici %>%
  mutate(metric = ifelse(grepl("hazard", endpoint, ignore.case = TRUE), "HR", "OR"),
         endpoint = sub(" \\(.*\\)", "", endpoint),
         lab = sprintf("%s %.2f (%.2f\u2013%.2f)", metric, est, lo, hi),
         row = dplyr::case_when(
           grepl("PFS", endpoint, ignore.case = TRUE) ~ "PFS",
           grepl("ORR|objective", endpoint, ignore.case = TRUE) ~ "ORR",
           TRUE ~ "Clinical benefit"),
         hj = hjust_in_panel(est, ici_lim[1], ici_lim[2])) %>%
  arrange(est) %>%
  mutate(row = factor(row, row))
pB <- ggplot(ici2, aes(est, row)) +
  annotate("rect", xmin = 0.80, xmax = 1.25, ymin = -Inf, ymax = Inf,
           fill = "grey92", colour = NA) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.3, colour = "grey40") +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y",
                width = 0.16, linewidth = 0.45, colour = "#8491B4") +
  geom_point(size = 2.3, colour = "#8491B4") +
  geom_text(aes(x = est, label = lab, hjust = hj), vjust = -1.5, size = 2.0) +
  scale_x_log10(limits = ici_lim, breaks = c(0.5, 1, 2),
                labels = c("0.5", "1", "2")) +
  labs(x = "Effect size (log scale); grey = 0.80\u20131.25", y = NULL, tag = "B",
       title = sprintf("IMmotion150 (%d endpoints, all null)", nrow(ici2))) +
  theme_cns(base_size = 8) +
  theme(axis.text.y = element_text(size = 6.8),
        plot.title = element_text(size = 8.5, face = "bold", hjust = 0.5),
        plot.margin = margin(2, 3, 2, 2, "mm"))

# 纵轴标签带上 n 与事件数, 否则读者无法判断 CI 宽窄的来源
par_lim <- c(0.7, 8)
par2 <- par %>%
  mutate(hit = cancer == "KIRC",
         cancer = factor(cancer, rev(c("KIRC", "STAD", "SARC", "LUSC"))),
         row = sprintf("%s\nn=%d, %d events", cancer, n, events),
         lab = sprintf("%.2f (%.2f\u2013%.2f)", HR, lo, hi),
         hj = hjust_in_panel(HR, par_lim[1], par_lim[2]))
par2$row <- factor(par2$row, par2$row[order(par2$cancer)])
pC <- ggplot(par2, aes(y = row)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.3, colour = "grey40") +
  geom_segment(aes(x = 1, xend = HR, colour = hit),
               linewidth = 2.6, lineend = "butt", show.legend = FALSE) +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y",
                width = 0.18, linewidth = 0.4, colour = "grey25") +
  geom_point(aes(x = HR, colour = hit), size = 1.7, show.legend = FALSE) +
  scale_colour_manual(values = c(`TRUE` = "#E64B35", `FALSE` = "#8491B4")) +
  geom_text(aes(x = HR, label = lab, hjust = hj), vjust = -1.55, size = 2.0) +
  scale_x_log10(limits = par_lim, breaks = c(1, 2, 4),
                labels = c("1", "2", "4")) +
  labs(x = "HR, high vs low (log scale)", y = NULL, tag = "C",
       title = "Soft-pass cancers") +
  theme_cns(base_size = 8) +
  theme(axis.text.y = element_text(size = 6.2, lineheight = 0.95),
        plot.title = element_text(size = 8.5, face = "bold", hjust = 0.5),
        plot.margin = margin(2, 3, 2, 2, "mm"))

fig7 <- (pA | pB | pC) + plot_layout(widths = c(1.45, 1, 1)) &
  theme(plot.tag = element_text(face = "bold", size = 11))
save_pub(fig7, "Fig6_boundary", OUT, w = 180, h = 90)
message(sprintf("null: TAS %.4f, percentile %.2f, z %.2f, n=%d",
                tas_c, pct, z, nrow(nv)))

# ---- supplement: GO + network without layout axes ----------------------
bp <- read.csv(file.path(RES, "12_mechanism/go_bp_enrichment.csv"))
cc <- read.csv(file.path(RES, "12_mechanism/go_cc_enrichment.csv"))
mod <- readRDS(file.path(RES, "02_model/KIRC/model.rds"))
gsel <- mod$genes_sel
beta10 <- as.numeric(mod$beta[mod$genes %in% gsel])
names(beta10) <- gsel

bp2 <- bp %>% arrange(p.adjust) %>% slice_head(n = 8) %>%
  mutate(Description = factor(Description, rev(Description)),
         ratio = Count / 21)
pG <- ggplot(bp2, aes(ratio, Description)) +
  geom_point(aes(size = Count, colour = -log10(p.adjust))) +
  scale_colour_gradient(low = "#3C5488", high = "#E64B35",
                        name = expression(-log[10](FDR)), guide = guide_cbar()) +
  scale_size_continuous(range = c(2, 6), name = "Genes") +
  labs(x = "Gene ratio (of 21)", y = NULL, title = "GO biological process") +
  theme_cns(base_size = 8) +
  theme(axis.text.y = element_text(size = 6.5),
        legend.position = "right",
        plot.title = element_text(size = 9, face = "bold", hjust = 0.5))

cc2 <- cc %>% arrange(p.adjust) %>% slice_head(n = 6) %>%
  mutate(Description = case_when(
    grepl("phosphatidylinositol 3-kinase", Description) ~ "PI3K complex",
    grepl("transferase complex", Description) ~ "PI transferase complex",
    TRUE ~ Description
  ))
edges <- bind_rows(lapply(seq_len(nrow(cc2)), function(i) {
  genes <- intersect(gsel, strsplit(cc2$geneID[i], "/")[[1]])
  if (!length(genes)) return(NULL)
  data.frame(from = genes, to = cc2$Description[i], stringsAsFactors = FALSE)
}))
nodes <- data.frame(name = unique(c(edges$from, edges$to)), stringsAsFactors = FALSE) %>%
  mutate(kind = ifelse(name %in% gsel, "gene", "pathway"),
         beta = beta10[name])
g <- as_tbl_graph(edges) %>%
  activate(nodes) %>%
  left_join(nodes, by = "name")
pN <- ggraph(g, layout = "fr") +
  geom_edge_link(colour = "grey75", linewidth = 0.3) +
  geom_node_point(aes(colour = kind, size = kind)) +
  geom_node_text(aes(label = name), repel = TRUE, size = 2.1, family = "Arial",
                 point.padding = 0.25, box.padding = 0.35, max.overlaps = Inf) +
  scale_colour_manual(values = c(gene = "#3C5488", pathway = "#00A087"), name = NULL) +
  scale_size_manual(values = c(gene = 2.5, pathway = 3.5), guide = "none") +
  labs(title = "Cellular-component terms") +
  theme_void(base_family = "Arial") +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 9, face = "bold", hjust = 0.5, family = "Arial"),
        plot.margin = margin(2, 2, 2, 2, "mm"))
figS <- (pG | pN) &
  theme(plot.tag = element_text(face = "bold", size = 11))
save_pub(figS, "FigS_pathway", OUT, w = 180, h = 110)

# ---- supplement: named drugs with FDR < 0.05, grouped ------------------
drug <- read.csv(file.path(RES, "14_drug/drug_sensitivity_diff.csv"))
catalog <- grepl("^[A-Z][0-9]|_|[0-9]{3,}", drug$drug)
named <- drug %>%
  filter(fdr < 0.05, !catalog) %>%
  mutate(class = case_when(
    drug %in% c("Axitinib", "Sunitinib", "Pazopanib", "Sorafenib", "Cabozantinib") ~ "VEGFR TKI",
    drug %in% c("Romidepsin", "Vorinostat", "Panobinostat", "Entinostat") ~ "HDAC",
    drug %in% c("Camptothecin", "SN-38", "Topotecan", "Irinotecan") ~ "Topoisomerase",
    drug %in% c("Tanespimycin") ~ "HSP90",
    drug %in% c("AZD6738") ~ "ATR",
    drug %in% c("XAV939") ~ "Tankyrase",
    drug %in% c("Pevonedistat") ~ "NAE",
    drug %in% c("GSK2606414") ~ "PERK",
    TRUE ~ "Other named agent"
  )) %>%
  group_by(class) %>%
  slice_min(fdr, n = 2, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(class, diff) %>%
  mutate(drug = factor(drug, drug))
pDrug <- ggplot(named, aes(diff, drug, fill = class)) +
  geom_col(width = 0.65) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey30") +
  scale_fill_manual(values = c("VEGFR TKI" = "#E64B35", HDAC = "#3C5488",
                               Topoisomerase = "#00A087", HSP90 = "#4DBBD5",
                               ATR = "#F39B7F", Tankyrase = "#8491B4",
                               NAE = "#91D1C2", PERK = "#7E6148",
                               "Other named agent" = "#B0B0B0"),
                    name = NULL) +
  labs(x = "Predicted IC50, TAS-high minus TAS-low", y = NULL,
       title = "Named agents with FDR < 0.05") +
  theme_cns(base_size = 8) +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 9, face = "bold", hjust = 0.5))
n_sig <- sum(drug$fdr < 0.05); n_neg <- sum(drug$fdr < 0.05 & drug$diff < 0)
pDist <- ggplot(drug, aes(diff, fill = fdr < 0.05)) +
  geom_histogram(bins = 40, colour = "white", linewidth = 0.15) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey30") +
  scale_fill_manual(values = c(`TRUE` = "#3C5488", `FALSE` = "#C8C8C8"),
                    labels = c(`TRUE` = "FDR < 0.05", `FALSE` = "FDR \u2265 0.05"), name = NULL) +
  labs(x = "Predicted IC50, TAS-high minus TAS-low", y = "Compounds",
       title = sprintf("All %d compounds: %d of %d significant are lower in TAS-high",
                       nrow(drug), n_neg, n_sig)) +
  theme_cns(base_size = 8) +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 8.5, face = "bold", hjust = 0.5))
figDrug <- (pDist / pDrug) + plot_layout(heights = c(0.7, 1)) &
  theme(plot.tag = element_text(face = "bold", size = 11))
save_pub(figDrug, "FigS_drug", OUT, w = 150, h = 190)
message("drugs plotted: ", paste(named$drug, collapse = ", "))
