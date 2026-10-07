#!/usr/bin/env Rscript
# =====================================================================
# 27d2_fig6_v2.R —— Figure 6 美学精修版（替换 27d 的 A/B/C 三面板）
# 数据与 27d 完全相同（不改动任何数值），仅重排视觉语言：
#   A 零分布直方图：标注改为右上对齐图例块，避免与红色标线重叠
#   B IMmotion150：标准森林图语法，效应量文本固定成列；灰带减淡
#   C Soft-pass 癌种：与 B 同语法（点+须线），KIRC 红色高亮，标签固定成列
# 输出: figures_pub/Fig6_boundary.{png,pdf}（600 dpi，save_pub 默认已改）
# =====================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(showtext); library(patchwork)
})
ROOT <- "/Users/xiaolin/WorkBuddy/2026-09-04-20-00-08/paper2_triaptosis"
OUT  <- file.path(ROOT, "figures_pub")
RES  <- file.path(ROOT, "results")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
source(file.path(ROOT, "code", "lib", "viz_base.R"))

nv  <- read.csv(file.path(RES, "05_null/triaptosis_random_null.csv"))
ici <- read.csv(file.path(RES, "03_ici/ici_null_forest.csv"))
pcd_all <- read.csv(file.path(RES, "05_null/pcd_benchmark_same_protocol.csv"))
par <- read.csv(file.path(RES, "02_model/parallel_cancer_hr.csv"))
tri <- pcd_all %>% filter(gene_set == "Triaptosis")
pcd <- pcd_all %>% filter(gene_set != "Triaptosis")

RED  <- "#C0392B"   # 主强调色（KIRC / triaptosis）
NAVY <- "#3C5488"   # 次色
TEAL <- "#00A087"
GREY <- "#8491B4"

fig_theme <- function(base = 8) {
  theme_cns(base_size = base) +
    theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
          plot.margin = margin(2.5, 3, 2.5, 2.5, "mm"))
}

# ============================== A. 零分布直方图 ==============================
bins <- 18
brk <- seq(min(nv$cv_nested), max(nv$cv_nested), length.out = bins + 1)
ymax <- max(table(cut(nv$cv_nested, brk, include.lowest = TRUE)))

marks <- bind_rows(
  data.frame(set = "Triaptosis", value = tri$cv_nested, col = RED,
             lab = sprintf("Triaptosis %.3f (%.1f)", tri$cv_nested, tri$percentile)),
  pcd %>% transmute(set = gene_set, value = cv_nested,
                    col = c(Ferroptosis = TEAL, Apoptosis = NAVY,
                            Disulfidptosis = GREY)[gene_set],
                    lab = sprintf("%s %.3f", gene_set, cv_nested))
)
leg_y <- ymax * 1.30

pA <- ggplot(nv, aes(cv_nested)) +
  geom_histogram(bins = bins, fill = "#E4E4E4", colour = "white", linewidth = 0.25) +
  geom_segment(data = marks,
               aes(x = value, xend = value, y = 0, yend = ymax * 1.10, colour = set,
                   linetype = set == "Triaptosis"),
               inherit.aes = FALSE, linewidth = 0.45, show.legend = FALSE) +
  geom_point(data = marks, aes(x = value, y = ymax * 1.10, colour = set),
             inherit.aes = FALSE, size = 1.6, show.legend = FALSE) +
  geom_point(data = marks, aes(x = -Inf, y = -Inf, colour = lab),
             inherit.aes = FALSE, size = 2.2) +
  scale_colour_manual(values = setNames(marks$col, marks$lab), name = NULL) +
  scale_linetype_manual(values = c(`TRUE` = "solid", `FALSE` = "33")) +
  scale_y_continuous(limits = c(0, leg_y + ymax * 0.06),
                     expand = expansion(mult = c(0, 0.01))) +
  scale_x_continuous(expand = expansion(mult = c(0.04, 0.06))) +
  labs(x = "Nested CV C-index of random 21-gene sets", y = "Count", tag = "A",
       title = sprintf("Random-set null (n = %d)", nrow(nv))) +
  fig_theme() +
  theme(legend.position = "inside", legend.position.inside = c(0.97, 0.78),
        legend.justification = c(1, 0.5),
        legend.text = element_text(size = 6, family = "Arial"),
        legend.key.spacing.x = unit(0.8, "mm"),
        legend.background = element_blank())

# ============================== B. IMmotion150 ==============================
ici2 <- ici %>%
  mutate(metric = ifelse(grepl("hazard", endpoint, ignore.case = TRUE), "HR", "OR"),
         endpoint = sub(" \\(.*\\)", "", endpoint),
         lab = sprintf("%s %.2f (%.2f\u2013%.2f)", metric, est, lo, hi),
         row = dplyr::case_when(
           grepl("PFS", endpoint, ignore.case = TRUE) ~ "PFS",
           grepl("ORR|objective", endpoint, ignore.case = TRUE) ~ "ORR",
           TRUE ~ "Clinical benefit"),
         row = factor(row, c("PFS", "ORR", "Clinical benefit")))
lab_col <- 2.40   # 文本列固定位置

pB <- ggplot(ici2, aes(est, row)) +
  annotate("rect", xmin = 0.80, xmax = 1.25, ymin = -Inf, ymax = Inf,
           fill = "grey95", colour = NA) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.3, colour = "grey45") +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y",
                width = 0.14, linewidth = 0.5, colour = NAVY) +
  geom_point(size = 2.4, colour = NAVY) +
  geom_text(aes(x = lab_col, label = lab), hjust = 0, vjust = 0.4,
            size = 2.1, family = "Arial", colour = "grey20") +
  scale_x_log10(limits = c(0.3, 30), breaks = c(0.5, 1, 2),
                labels = c("0.5", "1", "2")) +
  labs(x = "Effect size, high vs low (log scale)", y = NULL, tag = "B",
       title = "IMmotion150 outcomes") +
  fig_theme() +
  theme(axis.text.y = element_text(size = 7))

# ============================== C. Soft-pass 癌种 ==============================
par2 <- par %>%
  mutate(hit = cancer == "KIRC",
         row = sprintf("%s (n = %d, %d events)", cancer, n, events),
         lab = sprintf("%.2f (%.2f\u2013%.2f)", HR, lo, hi),
         hj = ifelse(HR < 1.3, 0, 0.5))
par2$row <- factor(par2$row, rev(par2$row))
lab_col2 <- 5.6

pC <- ggplot(par2, aes(HR, row)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.3, colour = "grey45") +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y",
                width = 0.16, linewidth = 0.5,
                colour = ifelse(par2$hit, RED, NAVY)) +
  geom_point(size = 2.4, colour = ifelse(par2$hit, RED, NAVY)) +
  geom_text(aes(label = lab, hjust = hj), vjust = -1.1,
            size = 2.1, family = "Arial", colour = "grey20") +
  scale_x_log10(limits = c(0.7, 12), breaks = c(1, 2, 4),
                labels = c("1", "2", "4")) +
  labs(x = "HR, high vs low (log scale)", y = NULL, tag = "C",
       title = "Soft-pass cancers") +
  fig_theme() +
  theme(axis.text.y = element_text(size = 6.4))

fig6v2 <- (pA | pB | pC) + plot_layout(widths = c(1.10, 1.30, 1.10)) &
  theme(plot.tag = element_text(face = "bold", size = 11))
save_pub(fig6v2, "Fig6_boundary", OUT, w = 180, h = 95)
message("Fig6 v2 saved: TAS ", round(tri$cv_nested, 4), " (", round(tri$percentile, 1), "th pct)")
