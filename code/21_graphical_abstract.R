# =====================================================================
# 21_graphical_abstract.R —— 图文摘要（投稿版）
# 四步各一个图形和一个数字。无甜甜圈，无色条标题。
# 免疫面板同时画出升高和降低，不画成六根同向红条。
# =====================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(survival)
  library(ggsurvfit); library(showtext); library(patchwork)
})
ROOT <- "/Users/xiaolin/WorkBuddy/2026-09-04-20-00-08/paper2_triaptosis"
OUT  <- file.path(ROOT, "figures_pub")
RES  <- file.path(ROOT, "results")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
source(file.path(ROOT, "code", "lib", "viz_base.R"))

screen <- read.csv(file.path(RES, "01_screen/pancancer_screen_matrix.csv"))
gc <- read.csv(file.path(RES, "01_screen/gene_cancer_cox.csv"))
sc <- read.csv(file.path(RES, "02_model/KIRC/score_table.csv"), row.names = 1)
ext <- read.csv(file.path(RES, "07_external_os/external_OS_TAS_summary.csv"))
imm <- read.csv(file.path(RES, "04_immune/immune_by_TAS_group.csv"))

n_sig <- sum(gc$cancer == "KIRC" & gc$p < 0.05)
cnt <- data.frame(
  step = factor(c("All", "Hard", "Soft", "KIRC"),
                levels = c("All", "Hard", "Soft", "KIRC")),
  n = c(nrow(screen), sum(screen$pass_hard), sum(screen$pass_soft), 1)
)
p1 <- ggplot(cnt, aes(step, n, fill = step == "KIRC")) +
  geom_col(width = 0.7, show.legend = FALSE) +
  geom_text(aes(label = n), vjust = -0.3, size = 2.6, family = "Arial") +
  scale_fill_manual(values = c(`TRUE` = "#E64B35", `FALSE` = "#8491B4")) +
  scale_y_continuous(limits = c(0, 42), expand = c(0, 0)) +
  labs(x = NULL, y = "Cancer types",
       title = "Screen",
       # 必须写清数据来源: cBioPortal 矩阵是 17/21, Xena 矩阵是 18/21 (见 Fig 1D/1E)
       subtitle = sprintf("%d/21 genes Cox-significant\nin KIRC (cBioPortal)", n_sig)) +
  theme_cns(base_size = 8) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6.5),
        plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 6.5, hjust = 0.5, colour = "grey30"))

sc$group <- factor(sc$group, c("Low", "High"))
ci <- summary(coxph(Surv(time, status) ~ group, data = sc))$conf.int
hr <- unname(ci[1, 1]); lo <- unname(ci[1, 3]); hi <- unname(ci[1, 4])
p2 <- survfit2(Surv(time, status) ~ group, data = sc) |>
  ggsurvfit(linewidth = 0.6) +
  scale_ggsurvfit() +
  scale_colour_manual(values = c(Low = "#3C5488", High = "#E64B35"), name = NULL) +
  scale_x_continuous(breaks = c(0, 50, 100, 150)) +
  annotate("text", x = 0, y = 0.15, hjust = 0,
           label = sprintf("High vs Low\nHR %.2f (%.2f\u2013%.2f)", hr, lo, hi),
           size = 2.3, family = "Arial") +
  labs(x = "Months", y = "Overall survival", title = "10-gene TAS") +
  theme_cns(base_size = 8) +
  theme(legend.position = "top",
        plot.title = element_text(size = 9, face = "bold", hjust = 0.5))

p3 <- ggplot(ext, aes(x = 1, y = cohort)) +
  geom_text(aes(label = sprintf("%.3f", Cindex)), size = 6.2, fontface = "bold",
            colour = "#4DBBD5", family = "Arial", vjust = 0.15) +
  geom_text(aes(label = sprintf("n = %d, %d deaths", n, events)),
            size = 2.1, colour = "grey35", family = "Arial", vjust = 2.4) +
  scale_y_discrete(limits = rev(ext$cohort)) +
  labs(x = "C-index", y = NULL, title = "External OS") +
  theme_cns(base_size = 8) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        plot.title = element_text(size = 9, face = "bold", hjust = 0.5))

top <- imm %>%
  mutate(delta = mean_high - mean_low) %>%
  slice_max(abs(delta), n = 6) %>%
  arrange(delta) %>%
  mutate(cell = factor(cell, cell),
         direction = ifelse(delta >= 0, "Higher in TAS-high", "Higher in TAS-low"))
# 两行图例原先占掉整张图文摘要的下四分之一, 改成面板内直接标注方向
p4 <- ggplot(top, aes(delta, cell, fill = direction)) +
  geom_col(width = 0.65, show.legend = FALSE) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey30") +
  annotate("text", x = max(top$delta), y = nrow(top) + 0.75, hjust = 1,
           label = "Higher in TAS-high", colour = "#E64B35",
           size = 2.2, family = "Arial") +
  annotate("text", x = min(top$delta), y = 0.18, hjust = 0,
           label = "Higher in TAS-low", colour = "#3C5488",
           size = 2.2, family = "Arial") +
  scale_fill_manual(values = c("Higher in TAS-high" = "#E64B35",
                               "Higher in TAS-low" = "#3C5488"),
                    name = NULL) +
  scale_y_discrete(expand = expansion(add = c(1.05, 1.25))) +
  labs(x = "Mean ssGSEA, high \u2212 low", y = NULL, title = "Immune shift") +
  theme_cns(base_size = 8) +
  theme(axis.text.y = element_text(size = 6),
        plot.title = element_text(size = 9, face = "bold", hjust = 0.5))

ga <- (p1 | p2 | p3 | p4) +
  plot_layout(widths = c(0.8, 1.15, 0.9, 1.15)) +
  plot_annotation(
    # 单行标题在 180 mm 宽度上会被左右截断, 必须手动折成两行
    title = paste0("Random-gene-calibrated screen: the triaptosis gene signal concentrates in ccRCC,\n",
                   "but the resulting score is prognostic rather than triaptosis-specific"),
    theme = theme(plot.title = element_text(size = 8.5, face = "bold", hjust = 0.5,
                                           lineheight = 1.15, family = "Arial",
                                           margin = margin(b = 1.5, unit = "mm")))
  )
save_pub(ga, "GraphicalAbstract", OUT, w = 180, h = 88)
message(sprintf("GA: %d/21 significant; binary HR %.2f (%.2f-%.2f)", n_sig, hr, lo, hi))
