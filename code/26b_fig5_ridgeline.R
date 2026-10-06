# =====================================================================
# 26b_fig5_ridgeline.R —— Fig5 升级版: 免疫浸润山峦图 (仙桃"山峦图"类别)
# 28 免疫细胞 x High/Low 组, ridge plot 展示 ssGSEA score 分布
# =====================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(ggridges)
  library(showtext); library(patchwork)
})
ROOT <- "/Users/xiaolin/WorkBuddy/2026-09-04-20-00-08/paper2_triaptosis"
OUT  <- file.path(ROOT, "figures_v2")
RES  <- file.path(ROOT, "results")
source(file.path(ROOT, "code", "lib", "viz_base.R"))

ss <- read.csv(file.path(RES, "04_immune/ssgsea_per_sample.csv"), row.names = 1)

# 长格式: cell x group x ssGSEA_score (ss 已含 group + patientId 列)
df <- ss %>%
  mutate(group = factor(group, c("Low", "High"))) %>%
  pivot_longer(-c(patientId, group), names_to = "cell", values_to = "score")
# 按 High/Low 中位数差排序 (差异大的在上)
ord <- df %>% group_by(cell) %>%
  summarise(d = median(score[group == "High"]) - median(score[group == "Low"])) %>%
  arrange(d) %>% pull(cell)
df <- df %>% mutate(cell = factor(cell, ord),
                    cell = gsub("\\.", " ", cell))

p <- ggplot(df, aes(score, cell, fill = group)) +
  geom_density_ridges(alpha = 0.75, scale = 0.95, rel_min_height = 0.01,
                      linewidth = 0.3, colour = "white") +
  scale_fill_manual(values = c(Low = "#3C5488", High = "#E64B35"), name = NULL) +
  labs(x = "ssGSEA enrichment score", y = NULL,
       title = "Immune cell infiltration distributions by TAS group") +
  theme_cns(base_size = 8, base_family = "Arial") +
  theme(axis.text.y = element_text(size = 6.5),
        axis.title.x = element_text(size = 7.5),
        legend.position = "top",
        plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        plot.margin = margin(4, 6, 4, 4, "mm"))

ggsave(file.path(OUT, "Fig5_ridgeline.png"), p,
       width = 150, height = 175, units = "mm", dpi = 300)
showtext_auto()
ggsave(file.path(OUT, "Fig5_ridgeline.pdf"), p,
       width = 150, height = 175, units = "mm", device = cairo_pdf)
showtext_auto(FALSE)
message("saved: Fig5_ridgeline")