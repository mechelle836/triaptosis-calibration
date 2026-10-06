# =====================================================================
# 27c_fig5_composite.R —— Fig5 免疫表型（投稿版）
# A 28 细胞 High−Low 哑铃（FDR 星号）
# B |delta| 最大的 6 种细胞分布
# C TAS 与免疫基因 Spearman：显著用强调色，不显著灰色
# 红/蓝只表示 High/Low，不表示相关方向。
# 输出: figures_pub/Fig5_immune_suppressive.{png,pdf}
# =====================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr)
  library(showtext); library(patchwork)
})
ROOT <- "/Users/xiaolin/WorkBuddy/2026-09-04-20-00-08/paper2_triaptosis"
OUT  <- file.path(ROOT, "figures_pub")
RES  <- file.path(ROOT, "results")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
source(file.path(ROOT, "code", "lib", "viz_base.R"))

imm <- read.csv(file.path(RES, "04_immune/immune_by_TAS_group.csv"))
ss  <- read.csv(file.path(RES, "04_immune/ssgsea_per_sample.csv"), check.names = FALSE)
cp  <- read.csv(file.path(RES, "13_checkpoint/tas_immune_gene_corr.csv"))

imm <- imm %>%
  mutate(delta = mean_high - mean_low,
         star = ifelse(wilcox_fdr < 0.05, "*", "")) %>%
  arrange(delta) %>%
  mutate(cell = factor(cell, cell))
n_sig <- sum(imm$wilcox_fdr < 0.05)
n_pt  <- nrow(ss)

# ---- A: dumbbell -------------------------------------------------------
pts <- imm %>%
  select(cell, Low = mean_low, High = mean_high) %>%
  pivot_longer(c(Low, High), names_to = "group", values_to = "mean") %>%
  mutate(group = factor(group, c("Low", "High")))
# 星号原先跟着每行哑铃的右端走, 在空白里散成一片, 无法对行阅读;
# 改为固定在面板右侧的一列, 并把 FDR 说明写进列头。
rng <- range(c(imm$mean_low, imm$mean_high))
star_x <- rng[2] + diff(rng) * 0.10
pA <- ggplot(imm, aes(y = cell)) +
  geom_segment(aes(x = mean_low, xend = mean_high, yend = cell),
               colour = "grey70", linewidth = 0.35) +
  geom_point(data = pts, aes(x = mean, colour = group), size = 1.7) +
  geom_text(aes(x = star_x, label = star), size = 3.0, colour = "grey20",
            hjust = 0.5, vjust = 0.72) +
  annotate("text", x = star_x, y = nrow(imm) + 1.6, label = "FDR\n< 0.05",
           size = 1.8, colour = "grey35", family = "Arial", lineheight = 1,
           hjust = 0.5, vjust = 0.5) +
  scale_colour_manual(values = c(Low = "#3C5488", High = "#E64B35"), name = NULL) +
  scale_x_continuous(breaks = seq(-0.2, 0.6, 0.2),
                     expand = expansion(mult = c(0.05, 0.04))) +
  scale_y_discrete(expand = expansion(add = c(0.6, 2.8))) +
  coord_cartesian(clip = "off") +
  labs(x = "Mean ssGSEA score", y = NULL, tag = "A",
       title = sprintf("Immune-cell enrichment (n=%d; %d of 28 FDR<0.05)", n_pt, n_sig)) +
  theme_cns(base_size = 7.5) +
  theme(axis.text.y = element_text(size = 6),
        # 28 行 + 横跨整幅的空白, 没有行导引线时无法把标签和哑铃对上
        panel.grid.major.y = element_line(colour = "grey93", linewidth = 0.25),
        legend.position = "top",
        plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
        plot.margin = margin(2, 4, 2, 2, "mm"))

# ---- B: top 6 |delta| distributions -----------------------------------
top6 <- imm %>% slice_max(abs(delta), n = 6) %>% arrange(delta) %>% pull(cell)
df6 <- ss %>%
  select(all_of(c(as.character(top6), "group"))) %>%
  pivot_longer(-group, names_to = "cell", values_to = "score") %>%
  mutate(group = factor(group, c("Low", "High")),
         cell = factor(cell, as.character(top6)))
# 小提琴图必须带检验结果, 否则读者无法判断这 6 个偏移是否显著
fdr6 <- imm %>% filter(cell %in% top6) %>%
  transmute(cell = factor(as.character(cell), as.character(top6)),
            lab = ifelse(wilcox_fdr < 0.001, "FDR < 0.001",
                         sprintf("FDR = %.3f", wilcox_fdr)))
ytop <- df6 %>% group_by(cell) %>% summarise(y = max(score, na.rm = TRUE), .groups = "drop")
fdr6 <- fdr6 %>% left_join(ytop, by = "cell")
pB <- ggplot(df6, aes(group, score, fill = group)) +
  geom_violin(width = 0.85, colour = NA, alpha = 0.35, show.legend = FALSE) +
  geom_boxplot(width = 0.18, outlier.size = 0.25, fill = "white",
               colour = "grey25", linewidth = 0.25, show.legend = FALSE) +
  geom_text(data = fdr6, aes(x = 1.5, y = y, label = lab), inherit.aes = FALSE,
            size = 1.8, colour = "grey25", family = "Arial", vjust = -0.3) +
  facet_wrap(~cell, ncol = 2, scales = "free_y") +
  scale_fill_manual(values = c(Low = "#3C5488", High = "#E64B35")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.16))) +
  labs(x = NULL, y = "ssGSEA score", tag = "B",
       title = sprintf("Largest High \u2212 Low shifts (n=%d)", n_pt)) +
  theme_cns(base_size = 7.5) +
  theme(axis.text.x = element_text(size = 6.5),
        strip.text = element_text(size = 6.5, face = "bold"),
        plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
        plot.margin = margin(2, 2, 2, 2, "mm"))

# ---- C: Spearman; one colour for FDR<0.05, grey for ns ----------------
canon <- c("PDCD1", "CTLA4", "FOXP3", "GZMB")
sig <- cp %>% filter(fdr < 0.05) %>% slice_max(abs(rho), n = 10)
named <- cp %>% filter(gene %in% canon)
ns  <- cp %>% filter(fdr >= 0.05) %>% slice_max(abs(rho), n = 3)
show_cp <- bind_rows(sig, named, ns) %>%
  distinct(gene, .keep_all = TRUE) %>%
  arrange(rho) %>%
  mutate(gene = factor(gene, gene),
         sig = fdr < 0.05)
pC <- ggplot(show_cp, aes(rho, gene, colour = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3, colour = "grey40") +
  geom_segment(aes(x = 0, xend = rho, yend = gene), linewidth = 0.3, show.legend = FALSE) +
  geom_point(size = 1.8) +
  scale_colour_manual(values = c(`TRUE` = "#00A087", `FALSE` = "#B0B0B0"),
                      labels = c(`TRUE` = "FDR < 0.05", `FALSE` = "ns"),
                      name = NULL) +
  labs(x = "Spearman rho with TAS", y = NULL, tag = "C",
       title = sprintf("Immune-gene correlation (n=%d)", n_pt)) +
  theme_cns(base_size = 7.5) +
  theme(axis.text.y = element_text(size = 6, face = "italic"),
        legend.position = "bottom",
        plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
        plot.margin = margin(2, 2, 2, 2, "mm"))

fig5 <- pA / (pB | pC) +
  plot_layout(heights = c(1.25, 1.15)) &
  theme(plot.tag = element_text(face = "bold", size = 11))

save_pub(fig5, "Fig5_immune_suppressive", OUT, w = 180, h = 210)
message(sprintf("immune stars %d/28; top6: %s", n_sig, paste(top6, collapse = ", ")))
message("checkpoint panel genes: ", paste(show_cp$gene, collapse = ", "))
