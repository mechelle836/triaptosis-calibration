# =====================================================================
# 27b_fig2_composite.R —— Fig2 TAS 模型（投稿版）
# A 10 基因系数  B KM + at risk + High vs Low HR (95% CI)
# C 表观 timeROC（survivalROC, KM 法）
# 508 列热图改为高低组均值，放补充图 FigS_risk_triptych。
# 不把连续 TAS 的 HR 3.52 标在二分类曲线上。
# =====================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(survival)
  library(jsonlite); library(ggsurvfit); library(showtext); library(patchwork)
})
ROOT <- "/Users/xiaolin/WorkBuddy/2026-09-04-20-00-08/paper2_triaptosis"
OUT  <- file.path(ROOT, "figures_pub")
RES  <- file.path(ROOT, "results")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
source(file.path(ROOT, "code", "lib", "viz_base.R"))

mod <- readRDS(file.path(RES, "02_model/KIRC/model.rds"))
sc  <- read.csv(file.path(RES, "02_model/KIRC/score_table.csv"), row.names = 1)
gsel <- mod$genes_sel
beta10 <- mod$beta[mod$genes %in% gsel]
names(beta10) <- gsel

# ---- A: coefficients ---------------------------------------------------
bd <- data.frame(gene = names(beta10), beta = as.numeric(beta10),
                 stringsAsFactors = FALSE) %>%
  arrange(beta) %>%
  mutate(gene = factor(gene, gene),
         sign = ifelse(beta > 0, "risk", "protective"))
pA <- ggplot(bd, aes(beta, gene, fill = sign)) +
  geom_col(width = 0.65, show.legend = FALSE) +
  geom_vline(xintercept = 0, colour = "grey30", linewidth = 0.3) +
  scale_fill_manual(values = c(risk = "#E64B35", protective = "#3C5488")) +
  labs(x = "Elastic-net coefficient", y = NULL, title = "10-gene TAS", tag = "A") +
  theme_cns(base_size = 8) +
  theme(axis.text.y = element_text(face = "italic", size = 7),
        plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
        plot.margin = margin(2, 3, 2, 2, "mm"))

# ---- B: KM, binary HR --------------------------------------------------
sc$group <- factor(sc$group, c("Low", "High"))
fit <- coxph(Surv(time, status) ~ group, data = sc)
ci <- summary(fit)$conf.int
hr <- unname(ci[1, 1]); lo <- unname(ci[1, 3]); hi <- unname(ci[1, 4])
lr <- survdiff(Surv(time, status) ~ group, data = sc)
lrp <- 1 - pchisq(lr$chisq, 1)
n_g <- table(sc$group)
ev_g <- tapply(sc$status, sc$group, sum)
# 事件数放注释块, 图例只留组名 + n, 否则图例在面板宽度内会被截断
hr_lab <- sprintf("%d vs %d deaths\nHR = %.2f (%.2f\u2013%.2f)\nlog-rank P = %.1e",
                  as.integer(ev_g["High"]), as.integer(ev_g["Low"]), hr, lo, hi, lrp)
pB <- survfit2(Surv(time, status) ~ group, data = sc) |>
  ggsurvfit(linewidth = 0.6) +
  add_confidence_interval(alpha = 0.12) +
  add_risktable(risktable_stats = "n.risk",
                stats_label = list(n.risk = "At risk"), size = 2.2,
                theme = km_risktable_theme()) +
  km_x_scale(breaks = seq(0, 140, 40)) +
  km_y_scale() +
  scale_colour_manual(values = c(Low = "#3C5488", High = "#E64B35"),
                      labels = sprintf("%s TAS (n=%d)", names(n_g), as.integer(n_g)),
                      name = NULL) +
  scale_fill_manual(values = c(Low = "#3C5488", High = "#E64B35"), guide = "none") +
  annotate("text", x = 3, y = 0.05, hjust = 0, vjust = 0,
           label = hr_lab, size = 2.3, family = "Arial", lineheight = 1.05) +
  labs(x = "Time (months)", y = "Overall survival",
       title = sprintf("TCGA-KIRC (n=%d)", nrow(sc)), tag = "B") +
  theme_cns(base_size = 8) +
  theme(legend.position = "top",
        legend.key.width = unit(5, "mm"),
        legend.text = element_text(size = 6.5),
        plot.title = element_text(size = 9, face = "bold", hjust = 0.5))
pB <- wrap_elements(full = ggsurvfit_build(pB))
message(sprintf("binary HR %.3f (%.3f-%.3f), logrank %.3g", hr, lo, hi, lrp))

# ---- C: apparent timeROC ----------------------------------------------
if (!requireNamespace("survivalROC", quietly = TRUE)) {
  stop("survivalROC is required so Fig2 AUC matches the KM-method estimate")
}
times <- c(12, 36, 60)
rocs <- lapply(times, function(tt) {
  survivalROC::survivalROC(Stime = sc$time, status = sc$status,
                           marker = sc$TAS, predict.time = tt, method = "KM")
})
roc_df <- bind_rows(lapply(seq_along(rocs), function(i) {
  data.frame(fpr = rocs[[i]]$FP, tpr = rocs[[i]]$TP,
             horizon = sprintf("%d mo (AUC %.3f)", times[i], rocs[[i]]$AUC))
}))
roc_df$horizon <- factor(roc_df$horizon, unique(roc_df$horizon))
aucs <- vapply(rocs, function(r) r$AUC, numeric(1))
message(sprintf("apparent AUC 12/36/60 = %.3f / %.3f / %.3f", aucs[1], aucs[2], aucs[3]))
# ROC 的两轴同为 0-1 概率, 必须等比; 原来横跨整幅宽度被压成 2.8:1
pC <- ggplot(roc_df, aes(fpr, tpr, colour = horizon)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.3) +
  geom_line(linewidth = 0.6) +
  scale_colour_manual(values = c("#3C5488", "#00A087", "#E64B35"), name = NULL) +
  scale_x_continuous(breaks = seq(0, 1, 0.25), expand = expansion(mult = 0.02)) +
  scale_y_continuous(breaks = seq(0, 1, 0.25), expand = expansion(mult = 0.02)) +
  coord_equal() +
  labs(x = "False positive rate", y = "True positive rate",
       title = "Apparent time-dependent ROC", tag = "C") +
  theme_cns(base_size = 8) +
  theme(legend.position = "bottom", legend.direction = "vertical",
        legend.key.height = unit(2.8, "mm"),
        legend.spacing.y = unit(0.2, "mm"),
        legend.margin = margin(0.5, 0, 0, 0, "mm"),
        legend.text = element_text(size = 6.2),
        plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
        plot.margin = margin(2, 2, 2, 2, "mm"))

fig2 <- (pA | pB | pC) +
  plot_layout(widths = c(0.72, 1.18, 1.0)) &
  theme(plot.tag = element_text(face = "bold", size = 11))
save_pub(fig2, "Fig2_model_construction", OUT, w = 180, h = 92)

# ---- supplement: rank curve, status, group-mean expression ------------
sc2 <- sc %>% arrange(TAS) %>% mutate(rank = row_number(),
                                     group = factor(group, c("Low", "High")))
cut_at <- sum(sc2$group == "Low")
pS1 <- ggplot(sc2, aes(rank, TAS, colour = group)) +
  geom_point(size = 0.35, show.legend = FALSE) +
  geom_vline(xintercept = cut_at + 0.5, linetype = "dashed", linewidth = 0.3, colour = "grey40") +
  scale_colour_manual(values = c(Low = "#3C5488", High = "#E64B35")) +
  labs(x = NULL, y = "TAS", title = "Patients ordered by TAS") +
  theme_cns(base_size = 8) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        plot.title = element_text(size = 9, face = "bold", hjust = 0.5))
pS2 <- ggplot(sc2, aes(rank, time, colour = factor(status, levels = c(0, 1),
                                                   labels = c("Censored", "Death")))) +
  geom_point(size = 0.35) +
  geom_vline(xintercept = cut_at + 0.5, linetype = "dashed", linewidth = 0.3, colour = "grey40") +
  scale_colour_manual(values = c(Censored = "#8491B4", Death = "#E64B35"), name = NULL) +
  labs(x = NULL, y = "Follow-up (months)", title = "Survival status") +
  theme_cns(base_size = 8) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = "top",
        plot.title = element_text(size = 9, face = "bold", hjust = 0.5))

ej <- fromJSON(file.path(ROOT, "data/KIRC_expr.json"))
expr <- do.call(cbind, lapply(ej, function(x) unlist(x)))
expr <- expr[intersect(gsel, rownames(expr)), , drop = FALSE]
colnames(expr) <- substr(colnames(expr), 1, 12)
expr <- expr[, intersect(colnames(expr), rownames(sc2)), drop = FALSE]
z <- t(scale(t(expr)))
means <- data.frame(
  gene = rownames(z),
  Low = rowMeans(z[, colnames(z) %in% rownames(sc2)[sc2$group == "Low"], drop = FALSE]),
  High = rowMeans(z[, colnames(z) %in% rownames(sc2)[sc2$group == "High"], drop = FALSE])
) %>%
  pivot_longer(-gene, names_to = "group", values_to = "z") %>%
  mutate(gene = factor(gene, rev(gsel)),
         group = factor(group, c("Low", "High")))
pS3 <- ggplot(means, aes(group, gene, fill = z)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  scale_fill_gradient2(low = "#3C5488", mid = "white", high = "#E64B35",
                       midpoint = 0, name = "Mean z", guide = guide_cbar()) +
  labs(x = NULL, y = NULL, title = "Mean expression by TAS group") +
  theme_cns(base_size = 8) +
  theme(axis.text.y = element_text(face = "italic", size = 7),
        plot.title = element_text(size = 9, face = "bold", hjust = 0.5))
figS <- pS1 / pS2 / pS3 +
  plot_layout(heights = c(0.7, 0.8, 1.1)) &
  theme(plot.tag = element_text(face = "bold", size = 11))
save_pub(figS, "FigS_risk_triptych", OUT, w = 120, h = 180)
