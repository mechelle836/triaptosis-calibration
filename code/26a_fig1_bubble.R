# =====================================================================
# 26a_fig1_bubble.R —— Fig1 升级版: 泛癌筛选气泡热图 (仙桃"气泡热图"类别)
# 21 基因 x 20 癌种 Cox HR, 气泡大小=-log10(P), 颜色=HR 方向
# =====================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(showtext); library(patchwork)
})
ROOT <- "/Users/xiaolin/WorkBuddy/2026-09-04-20-00-08/paper2_triaptosis"
OUT  <- file.path(ROOT, "figures_v2")
RES  <- file.path(ROOT, "results")
source(file.path(ROOT, "code", "lib", "viz_base.R"))

gc <- read.csv(file.path(RES, "01_screen/gene_cancer_cox.csv"))

# 排序: 癌种按显著基因数, 基因按 KIRC 显著性
cancer_sig <- gc %>% filter(p < 0.05) %>% count(cancer, sort = TRUE)
cancer_order <- cancer_sig$cancer
gene_order <- gc %>% filter(cancer == "KIRC") %>% arrange(p) %>% pull(gene)
gc2 <- gc %>%
  mutate(cancer = factor(cancer, levels = cancer_order),
         gene = factor(gene, levels = rev(gene_order)),
         log10p = -log10(p),
         direction = case_when(p >= 0.05 ~ "ns",
                               HR > 1 ~ "risk (HR>1)",
                               TRUE ~ "protective (HR<1)"))

p <- ggplot(gc2, aes(cancer, gene)) +
  geom_point(aes(size = log10p, colour = direction), alpha = 0.9) +
  scale_colour_manual(values = c("risk (HR>1)" = "#E64B35",
                                 "protective (HR<1)" = "#3C5488",
                                 "ns" = "grey88"), name = NULL) +
  scale_size_continuous(range = c(0.2, 6), name = "-log10(P)",
                        breaks = c(1.3, 2, 3, 5, 10)) +
  labs(x = NULL, y = NULL,
       title = "Pan-cancer prognostic screen of triaptosis genes") +
  theme_cns(base_size = 8, base_family = "Arial") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        axis.text.y = element_text(face = "italic", size = 7.5),
        legend.position = "right",
        plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        plot.margin = margin(4, 4, 4, 4, "mm")) +
  guides(colour = guide_legend(order = 1),
         size = guide_legend(order = 2))

ggsave(file.path(OUT, "Fig1_bubble_heatmap.png"), p,
       width = 190, height = 150, units = "mm", dpi = 300)
showtext_auto()
ggsave(file.path(OUT, "Fig1_bubble_heatmap.pdf"), p,
       width = 190, height = 150, units = "mm", device = cairo_pdf)
showtext_auto(FALSE)
message("saved: Fig1_bubble_heatmap")