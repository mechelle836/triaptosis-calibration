# =====================================================================
# 27f_fig3_utility.R —— Fig3 预测效用（投稿版，2×2）
# A OOF timeROC（曲线由已保存的 OOF 风险重算，不重新拟合模型）
# B IPCW C-index 与 time-dependent AUC（表内数字）
# C IPCW Brier score（不再叫 calibration error）
# D 36 个月校准
# DCA 放补充图，并标明 treat-all / treat-none。
# =====================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(survival)
  library(timeROC); library(showtext); library(patchwork)
})
ROOT <- "/Users/xiaolin/WorkBuddy/2026-09-04-20-00-08/paper2_triaptosis"
OUT  <- file.path(ROOT, "figures_pub")
RES  <- file.path(ROOT, "results")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
source(file.path(ROOT, "code", "lib", "viz_base.R"))

oof <- read.csv(file.path(RES, "05_validation/oof_scores.csv"), row.names = 1)
ipcw_c <- read.csv(file.path(RES, "05_validation/ipcw_cindex.csv"))
ipcw_a <- read.csv(file.path(RES, "05_validation/auc_ipcw.csv"))
brier <- read.csv(file.path(RES, "05_validation/brier.csv"))
calib <- read.csv(file.path(RES, "05_validation/calibration.csv"))

roc <- timeROC(T = oof$time, delta = oof$status, marker = oof$risk,
               cause = 1, weighting = "marginal",
               times = c(12, 36, 60), ROC = TRUE, iid = FALSE)
message(sprintf("OOF timeROC AUC %.3f / %.3f / %.3f ; table %.3f / %.3f / %.3f",
                roc$AUC[1], roc$AUC[2], roc$AUC[3],
                ipcw_a$AUC[1], ipcw_a$AUC[2], ipcw_a$AUC[3]))
if (max(abs(roc$AUC - ipcw_a$AUC)) > 0.02) {
  warning("recomputed OOF AUC differs from auc_ipcw.csv by >0.02")
}
roc_df <- bind_rows(lapply(seq_len(3), function(i) {
  data.frame(fpr = roc$FP[, i], tpr = roc$TP[, i],
             horizon = sprintf("%d mo (AUC %.3f)", c(12, 36, 60)[i], roc$AUC[i]))
}))
roc_df$horizon <- factor(roc_df$horizon, unique(roc_df$horizon))
pA <- ggplot(roc_df, aes(fpr, tpr, colour = horizon)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.3) +
  geom_line(linewidth = 0.6) +
  scale_colour_manual(values = c("#3C5488", "#00A087", "#E64B35"), name = NULL) +
  labs(x = "False positive rate", y = "True positive rate",
       title = "OOF time-dependent ROC", tag = "A") +
  theme_cns(base_size = 8) +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 9, face = "bold", hjust = 0.5))

disc <- bind_rows(
  ipcw_c %>% transmute(times, value = C, metric = "IPCW C-index"),
  ipcw_a %>% transmute(times, value = AUC, metric = "Time-dependent AUC")
)
pB <- ggplot(disc, aes(times, value, colour = metric)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey60", linewidth = 0.3) +
  geom_line(linewidth = 0.55) +
  geom_point(size = 2) +
  geom_text(aes(label = sprintf("%.3f", value),
                vjust = ifelse(metric == "IPCW C-index", 1.9, -1.1)),
            size = 2.1, show.legend = FALSE) +
  scale_colour_manual(values = c("IPCW C-index" = "#3C5488",
                                 "Time-dependent AUC" = "#00A087"), name = NULL) +
  scale_x_continuous(breaks = c(12, 36, 60), expand = expansion(mult = 0.12)) +
  coord_cartesian(ylim = c(0.45, 0.95)) +
  labs(x = "Time (months)", y = "Out-of-fold estimate",
       title = "Discrimination (5-fold IPCW)", tag = "B") +
  theme_cns(base_size = 8) +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 9, face = "bold", hjust = 0.5))

# 数值标签一律朝内侧 (hjust 随时间点变化), 原先 12 月的 KM 标签被左轴切掉
brier <- brier %>%
  # 12 月两条线几乎重合, 标签放点的左侧 (x 轴左侧留了 12% 扩展量) 才不压点
  mutate(lab_hjust = ifelse(times == 12, 1.15, ifelse(times == 60, 1.2, 0.5)))
pC <- ggplot(brier, aes(times, Brier, colour = model, group = model)) +
  geom_line(linewidth = 0.55) +
  geom_point(size = 2) +
  geom_text(data = function(d) filter(d, model == "KM"),
            aes(label = sprintf("%.3f", Brier), hjust = lab_hjust),
            vjust = -1.4, size = 2.1, show.legend = FALSE) +
  geom_text(data = function(d) filter(d, model == "TAS"),
            aes(label = sprintf("%.3f", Brier), hjust = lab_hjust),
            vjust = 2.4, size = 2.1, show.legend = FALSE) +
  scale_colour_manual(values = c(TAS = "#E64B35", KM = "#8491B4"),
                      labels = c(KM = "Kaplan-Meier reference", TAS = "TAS"),
                      name = NULL) +
  scale_x_continuous(breaks = c(12, 36, 60), expand = expansion(mult = 0.12)) +
  scale_y_continuous(expand = expansion(mult = c(0.18, 0.14))) +
  labs(x = "Time (months)", y = "IPCW Brier score",
       title = "IPCW Brier score (lower is better)", tag = "C") +
  theme_cns(base_size = 8) +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 9, face = "bold", hjust = 0.5))

# 只有 4 个点又没有不确定度, 审稿人无法判断偏离对角线是否有意义。
# 按上游 07_model_validation.R 的同一分组重算 36 月 KM 的 Greenwood 95% CI。
oof_q <- quantile(oof$risk, probs = seq(0, 1, 0.25), na.rm = TRUE)
oof$grp <- cut(oof$risk, breaks = oof_q, include.lowest = TRUE, labels = paste0("Q", 1:4))
cal_ci <- do.call(rbind, lapply(levels(oof$grp), function(g) {
  s <- summary(survfit(Surv(time, status) ~ 1, data = oof[oof$grp == g, ]), times = 36)
  data.frame(group = g, obs = s$surv, obs_lo = s$lower, obs_hi = s$upper)
}))
calib <- calib %>% left_join(cal_ci, by = "group")
if (max(abs(calib$obs_36 - calib$obs)) > 1e-6) {
  warning("recomputed 36-month KM differs from calibration.csv; check the quartile split")
}
pD <- ggplot(calib, aes(pred_36, obs_36)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.3) +
  geom_line(colour = "grey75", linewidth = 0.4) +
  geom_errorbar(aes(ymin = obs_lo, ymax = obs_hi), width = 0.012,
                linewidth = 0.35, colour = "#3C5488") +
  geom_point(size = 2.2, colour = "#3C5488") +
  # Q1/Q2 预测值相近, 统一朝一侧会互相压字, 逐点给方向
  geom_text(aes(label = group,
                x = pred_36 + c(Q1 = 0.022, Q2 = 0.0, Q3 = -0.022, Q4 = -0.022)[group],
                y = obs_36 + c(Q1 = -0.055, Q2 = 0.075, Q3 = 0.0, Q4 = 0.0)[group],
                hjust = c(Q1 = 0, Q2 = 0.5, Q3 = 1, Q4 = 1)[group]),
            size = 2.2) +
  annotate("text", x = 0.46, y = 0.99, hjust = 0, vjust = 1, size = 2.0,
           colour = "grey35", family = "Arial", lineheight = 1.1,
           label = sprintf("n = %d per quartile\ndashed line = perfect calibration\nbars = 95%% CI of observed KM",
                           calib$n[1])) +
  coord_cartesian(xlim = c(0.45, 1), ylim = c(0.42, 1)) +
  labs(x = "Predicted 36-month survival", y = "Observed 36-month survival",
       title = "Calibration by predicted-risk quartile", tag = "D") +
  theme_cns(base_size = 8) +
  theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5))

fig3 <- (pA | pB) / (pC | pD) &
  theme(plot.tag = element_text(face = "bold", size = 11))
save_pub(fig3, "Fig3_utility", OUT, w = 180, h = 165)

if (requireNamespace("dcurves", quietly = TRUE)) {
  dca_dat <- data.frame(time = oof$time, status = oof$status, risk = 1 - oof$S36)
  dca_obj <- dcurves::dca(Surv(time, status) ~ risk, data = dca_dat,
                          time = 36, thresholds = seq(0.05, 0.90, 0.01))
  nb <- dca_obj$dca %>%
    mutate(label = dplyr::recode(variable,
                                 all = "Treat all", none = "Treat none", risk = "TAS")) %>%
    filter(label %in% c("TAS", "Treat all", "Treat none"))
  pF <- ggplot(nb, aes(threshold, net_benefit, colour = label, linetype = label)) +
    geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
    geom_line(linewidth = 0.6) +
    scale_colour_manual(values = c(TAS = "#E64B35", "Treat all" = "#8491B4",
                                   "Treat none" = "black"), name = NULL) +
    scale_linetype_manual(values = c(TAS = "solid", "Treat all" = "dashed",
                                     "Treat none" = "dotted"), name = NULL) +
    coord_cartesian(ylim = c(-0.05, 0.35)) +
    labs(x = "Threshold probability", y = "Net benefit",
         title = "Decision curve, 36-month OS (OOF)") +
    theme_cns(base_size = 8) +
    theme(legend.position = "bottom",
          plot.title = element_text(size = 9, face = "bold", hjust = 0.5))
  save_pub(pF, "FigS_dca", OUT, w = 120, h = 90)
} else {
  message("dcurves not installed; FigS_dca skipped")
}
