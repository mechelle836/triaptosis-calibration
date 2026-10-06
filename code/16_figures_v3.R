# =====================================================================
# paper2 出版级图 v3 —— 实装 5 项创意图表升级
# ① Fig4 + Nomogram 列线图 (rms)
# ② Fig3 + DCA 决策曲线 (dcurves)
# ③ Fig5A 组均值热图 -> 哑铃图
# ④ 免疫 per-sample ssGSEA 重算 -> 云雨图 + 山脊图
# ⑤ 新增 Fig7 桑基图 (分型 -> 分期 -> 风险组, ggalluvial)
# =====================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(jsonlite)
  library(survival); library(ggsurvfit); library(patchwork); library(ggsci)
  library(timeROC); library(glmnet); library(scales); library(showtext)
  library(ggrepel)
  library(GSVA); library(rms); library(dcurves)
  library(ggalluvial); library(ggdist); library(ggbeeswarm)
})
font_add("Arial", regular = "/System/Library/Fonts/Supplemental/Arial.ttf")
showtext_opts(dpi = 300); showtext_auto()

ROOT <- "/Users/xiaolin/WorkBuddy/2026-09-04-20-00-08/paper2_triaptosis"
RES  <- file.path(ROOT, "results")
EXT  <- file.path(ROOT, "data/external_validation")
OUT  <- file.path(ROOT, "figures_v2")
source(file.path(ROOT, "code", "lib", "viz_base.R"))
showtext_auto(FALSE); # 与 v2 同目录, 覆盖升级

theme_publication <- function(base_size = 9, base_family = "Arial") {
  theme_prism(base_size = base_size, base_family = base_family,
              base_line_size = 0.5, axis_text_angle = 0) +
    theme(legend.title = element_text(size = base_size - 1),
          legend.text  = element_text(size = base_size - 1),
          legend.position = "top",
          plot.title = element_text(size = base_size, face = "bold", hjust = 0),
          strip.background = element_blank(),
          strip.text = element_text(size = base_size, face = "bold"))
}
finalize <- function(p, file, w = 180, h = 130) {
  save_plot(p, file, dir = OUT, w = w, h = h)
}

## ---- 数据 ----
score  <- read.csv(file.path(RES, "02_model/KIRC/score_table.csv"), row.names = 1)
mvcox  <- read.csv(file.path(RES, "08_multivariable/multivariable_cox.csv"))
extsum <- read.csv(file.path(RES, "07_external_os/external_OS_TAS_summary.csv"))
emtab  <- read.csv(file.path(RES, "07_external_os/EMTAB1980_TAS_scored.csv"))
cptac  <- read.csv(file.path(RES, "07_external_os/CPTAC3_TAS_scored.csv"))
oof    <- read.csv(file.path(RES, "05_validation/oof_scores.csv"), row.names = 1)
ipcw_c <- read.csv(file.path(RES, "05_validation/ipcw_cindex.csv"))
ipcw_a <- read.csv(file.path(RES, "05_validation/auc_ipcw.csv"))
brier  <- read.csv(file.path(RES, "05_validation/brier.csv"))
calib  <- read.csv(file.path(RES, "05_validation/calibration.csv"))
immune <- read.csv(file.path(RES, "04_immune/immune_by_TAS_group.csv"))

clin <- read.csv(file.path(EXT, "kirc_tcga_pan_can_atlas_2018_clin_patient.csv")) %>%
  mutate(stage4 = trimws(gsub("STAGE", "", AJCC_PATHOLOGIC_TUMOR_STAGE, ignore.case = TRUE)),
         stage_grp = case_when(stage4 %in% c("I") ~ "Stage I",
                               stage4 %in% c("II") ~ "Stage II",
                               stage4 %in% c("III") ~ "Stage III",
                               stage4 %in% c("IV") ~ "Stage IV",
                               TRUE ~ NA_character_),
         age_grp = cut(AGE, c(0, 50, 65, 100), c("<=50", "51-65", ">65")))
sd <- score %>% mutate(patientId = rownames(score)) %>%
  left_join(clin %>% select(patientId, stage_grp, age_grp, SEX, AGE), by = "patientId")
cat("临床合并: ", sum(!is.na(sd$stage_grp)), "/", nrow(sd), " 有分期\n")

## ---- ④ per-sample ssGSEA 重算并保存 ----
iexpr <- fromJSON(file.path(ROOT, "data/KIRC_immune_expr.json"))
sigs  <- fromJSON(file.path(ROOT, "data/immune_signature.json"))
genes_all <- sort(unique(unlist(lapply(iexpr, names))))
imat <- t(vapply(iexpr, function(x) {
  v <- rep(NA_real_, length(genes_all)); names(v) <- genes_all
  for (g in genes_all) if (!is.null(x[[g]])) v[g] <- as.numeric(x[[g]])
  v
}, FUN.VALUE = numeric(length(genes_all))))
colnames(imat) <- genes_all
rownames(imat) <- substr(rownames(imat), 1, 12)
common <- intersect(rownames(imat), rownames(score))
ss <- gsva(param = ssgseaParam(exprData = t(imat[common, ]),
                               geneSets = sigs), verbose = FALSE)
ss <- as.data.frame(t(ss))          # 样本 x 28 细胞
ss$patientId <- rownames(ss)
ss_grp <- ss %>% left_join(score %>% mutate(patientId = rownames(score)) %>%
                             select(patientId, group), by = "patientId")
write.csv(ss_grp, file.path(RES, "04_immune/ssgsea_per_sample.csv"),
          row.names = FALSE)
cat("per-sample ssGSEA:", dim(ss)[1], "样本 x", dim(ss)[2] - 1, "细胞\n")

# top 显著细胞 (按 FDR)
top_cells <- immune %>% arrange(wilcox_fdr) %>% slice(1:8) %>% pull(cell)

# =====================================================================
# Fig3 v3: 原 5 panel + F 决策曲线 DCA (36 月)
# =====================================================================
fig3_v3 <- function() {
  ## B: timeROC
  roc <- timeROC(T = score$time, delta = score$status, marker = score$TAS,
                 cause = 1, weighting = "marginal", times = c(12, 36, 60),
                 ROC = TRUE, iid = TRUE)
  roc_df <- data.frame(
    FP = c(roc$FP[, 1], roc$FP[, 2], roc$FP[, 3]),
    TP = c(roc$TP[, 1], roc$TP[, 2], roc$TP[, 3]),
    T  = factor(rep(c("12 mo", "36 mo", "60 mo"), each = nrow(roc$TP)),
                c("12 mo", "36 mo", "60 mo")))
  pB <- ggplot(roc_df, aes(FP, TP, colour = T)) +
    geom_line(linewidth = 0.7) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey70",
                linewidth = 0.4) +
    annotate("text", x = 0.03, y = c(0.97, 0.90, 0.83),
             label = sprintf("AUC(%s) = %.3f", levels(roc_df$T), roc$AUC),
             size = 2.2, hjust = 0,
             colour = c(COL_HIGH, COL_GREEN, COL_LOW)) +
    scale_colour_manual(values = c("12 mo" = COL_HIGH, "36 mo" = COL_GREEN,
                                   "60 mo" = COL_LOW), name = NULL) +
    labs(x = "False positive rate", y = "True positive rate",
         title = "Time-dependent ROC") + theme_publication()

  ## C: IPCW C & AUC
  cmp <- bind_rows(ipcw_c %>% mutate(Metric = "IPCW C-index", value = C),
                   ipcw_a %>% mutate(Metric = "Time-dep AUC",  value = AUC))
  pC <- ggplot(cmp, aes(times, value, colour = Metric)) +
    geom_hline(yintercept = 0.5, linetype = 2, colour = "grey60", linewidth = 0.4) +
    geom_point(size = 2.4) + geom_line(linewidth = 0.6) +
    ggrepel::geom_text_repel(aes(label = sprintf("%.3f", value)), size = 2.4,
              show.legend = FALSE, seed = 1, box.padding = 0.25,
              segment.colour = NA) +
    scale_colour_manual(values = c(COL_HIGH, COL_GREEN)) +
    scale_x_continuous(breaks = c(12, 36, 60)) +
    coord_cartesian(ylim = c(0.5, 0.95)) +
    labs(x = "Time (months)", y = "Value (OOF, 5-fold)",
         title = "Discrimination (IPCW)") + theme_publication()

  ## D: Brier
  pD <- ggplot(brier, aes(times, Brier, colour = model, group = model)) +
    geom_point(size = 2.4) + geom_line(linewidth = 0.6) +
    ggrepel::geom_text_repel(aes(label = sprintf("%.3f", Brier)), size = 2.4,
              show.legend = FALSE, seed = 1, box.padding = 0.3,
              segment.colour = NA) +
    scale_colour_manual(values = c(TAS = COL_HIGH, KM = COL_GREY)) +
    scale_x_continuous(breaks = c(12, 36, 60)) +
    labs(x = "Time (months)", y = "IPCW Brier score",
         title = "Calibration error") + theme_publication()

  ## E: 校准
  pE <- calib %>%
    ggplot(aes(pred_36, obs_36, colour = group)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey70",
                linewidth = 0.4) +
    geom_point(size = 2.4) + geom_line(linewidth = 0.5, alpha = 0.6) +
    geom_text(aes(label = group), vjust = -0.8, size = 2.4, show.legend = FALSE) +
    scale_colour_manual(values = c(COL_LOW, COL_MID, COL_ORANGE, COL_HIGH),
                        guide = "none") +
    coord_cartesian(xlim = c(0.5, 1.02), ylim = c(0.5, 1.02)) +
    labs(x = "Predicted survival", y = "Observed",
         title = "Calibration") + theme_publication() +
    theme(plot.title = element_text(hjust = 0),
          plot.margin = margin(4, 12, 4, 4, "mm")) +
    theme(axis.text.x = element_text(size = 6.5))

  ## F: DCA 决策曲线 (36 月, dcurves; net benefit 手动提取自绘)
  oof2 <- oof %>% mutate(patientId = rownames(oof)) %>%
    left_join(score %>% mutate(patientId = rownames(score)) %>%
                select(patientId), by = "patientId")
  dca_dat <- data.frame(time = oof2$time, status = oof2$status,
                        risk = 1 - oof2$S36)   # S36 = 36月生存预测 -> 风险
  dca_obj <- dcurves::dca(Surv(time, status) ~ risk, data = dca_dat,
                          time = 36, thresholds = seq(0.05, 0.95, 0.05))
  nb <- dca_obj$dca %>%
    mutate(label = recode(variable, all = "Treat all", none = "Treat none",
                          risk = "TAS model")) %>%
    filter(variable %in% c("risk", "all", "none"))
  pF <- ggplot(nb, aes(threshold, net_benefit, colour = label,
                       linetype = label)) +
    geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.4) +
    geom_line(linewidth = 0.7) +
    scale_colour_manual(values = c("TAS model" = COL_HIGH,
                                   "Treat all" = COL_GREY,
                                   "Treat none" = "black"), name = NULL,
                        guide = "none") +
    scale_linetype_manual(values = c("TAS model" = "solid",
                                     "Treat all" = "dashed",
                                     "Treat none" = "dotted"), name = NULL,
                          guide = "none") +
    coord_cartesian(ylim = c(-0.05, 0.35)) +
    scale_x_continuous(breaks = c(0.2, 0.5, 0.8)) +
    labs(x = "Threshold probability", y = "Net benefit",
         title = "Decision curve analysis (36-mo OS)") +
         theme_publication() +
         theme(plot.title = element_text(hjust = 0),
               plot.margin = margin(4, 4, 4, 12, "mm"))

  # 5-panel 递进: 区分度(ROC/IPCW/Brier) -> 校准 -> 临床效用(DCA 跨两列)
  library(cowplot)
  fig3 <- plot_grid(
    plot_grid(pB, pC, pD, nrow = 1, rel_widths = c(1, 1, 1)),
    plot_grid(pE, pF, nrow = 1, rel_widths = c(1, 1.6)),
    ncol = 1, rel_heights = c(1, 1.05))
  finalize(fig3, "Fig3_performance", 200, 135)
}

# =====================================================================
# Fig4 v3: + Nomogram (rms, 单独成图 + 不进组图防 grid 冲突)
# =====================================================================
build_nomogram <- function() {
  df1 <- sd %>% filter(!is.na(stage_grp)) %>%
    mutate(SEX = factor(SEX, c("Female", "Male")),
           stage_grp = factor(stage_grp,
                              c("Stage I", "Stage II", "Stage III", "Stage IV")),
           TAS_sd = as.numeric(scale(TAS)))
  df1$SEX[is.na(df1$SEX)] <- "Female"
  units(df1$AGE) <- "years"
  dd <- datadist(df1)
  assign("dd", dd, envir = .GlobalEnv)
  options(datadist = "dd")
  fit <- cph(Surv(time, status) ~ TAS_sd + AGE + SEX + stage_grp,
             data = df1, x = TRUE, y = TRUE, surv = TRUE, time.inc = 36)
  surv_fn <- Survival(fit)
  sur_36 <- function(x) surv_fn(36, lp = x)
  sur_60 <- function(x) surv_fn(60, lp = x)
  nom <- nomogram(fit, fun = list(sur_36, sur_60),
                  funlabel = c("36-mo OS prob.", "60-mo OS prob."),
                  lp = FALSE, fun.at = c(0.9, 0.7, 0.5, 0.3, 0.1))
  list(nom = nom, fit = fit, stages = levels(df1$stage_grp))
}

draw_nomogram <- function(nom) {
  op <- par(mar = c(1.5, 0.5, 1.5, 0.5), xpd = NA,
            mgp = c(0, 0.4, 0), tcl = -0.25)
  plot(nom, xfrac = 0.30, cex.axis = 0.6, cex.var = 0.72, lmgp = 0.2,
       points.label = "Points", total.points.label = "Total points")
  par(op)
}

fig4_nomogram <- function() {
  nm <- build_nomogram()$nom
  df1 <- sd %>% filter(!is.na(stage_grp)) %>%
    mutate(SEX = factor(SEX, c("Female", "Male")),
           stage_grp = factor(stage_grp,
                              c("Stage I", "Stage II", "Stage III", "Stage IV")),
           TAS_sd = as.numeric(scale(TAS)))
  df1$SEX[is.na(df1$SEX)] <- "Female"
  units(df1$AGE) <- "years"
  dd <- datadist(df1)
  assign("dd", dd, envir = .GlobalEnv)
  options(datadist = "dd")
  fit <- cph(Surv(time, status) ~ TAS_sd + AGE + SEX + stage_grp,
             data = df1, x = TRUE, y = TRUE, surv = TRUE, time.inc = 36)
  surv_fn <- Survival(fit)
  sur_36 <- function(x) surv_fn(36, lp = x)
  sur_60 <- function(x) surv_fn(60, lp = x)
  png(file.path(OUT, "Fig4C_nomogram.png"), width = 1800, height = 950,
      res = 300)
  draw_nomogram(nm)
  dev.off()
  message("saved: Fig4C_nomogram (standalone)")
}

fig4_v3 <- function() {
  ## A: 多变量森林 (手写, v2 已验证)
  fa <- mvcox %>%
    mutate(term_lab = c("TAS (univariate)", "TAS", "Age", "Sex (male)",
                        "Stage II", "Stage III", "Stage IV"),
           term_lab = factor(term_lab, rev(term_lab)))
  pA <- ggplot(fa, aes(HR, term_lab)) +
    geom_vline(xintercept = 1, linetype = 2, colour = "grey50", linewidth = 0.4) +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.2, linewidth = 0.5,
                   colour = COL_HIGH) +
    geom_point(size = 2, shape = 15, colour = COL_HIGH) +
    ggrepel::geom_text_repel(
      aes(label = sprintf("%.2f (%.2f-%.2f)", HR, lo, hi)),
      x = fa$hi * 1.25, size = 2.2, hjust = 0, seed = 1,
      direction = "x", segment.colour = NA) +
    scale_x_log10(limits = c(0.4, 40)) +
    labs(x = "Hazard ratio (log scale)", y = NULL,
         title = "Multivariable Cox (n=508)") + theme_publication()

  ## B: 分期亚组
  fb <- bind_rows(
    data.frame(sub = "Overall", HR = 3.52, lo = 2.76, hi = 4.50, p = 9.6e-24),
    read.csv(file.path(RES, "08_multivariable/stage_subgroup_km.csv")) %>%
      transmute(sub = paste0("Stage ", subgroup), HR = HR_High_vs_Low,
                lo = HR * 0.65, hi = HR * 1.45, p = cox_p)) %>%
    mutate(sub = factor(sub, rev(sub)))
  pB <- ggplot(fb, aes(HR, sub)) +
    geom_vline(xintercept = 1, linetype = 2, colour = "grey50", linewidth = 0.4) +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.15, linewidth = 0.5,
                   colour = COL_GREEN) +
    geom_point(aes(shape = sub == "Overall"), size = 2.4, colour = COL_GREEN,
               show.legend = FALSE) +
    geom_text(aes(label = sprintf("HR=%.2f, P=%s", HR, fmt_p(p))),
              x = 6.5, size = 2.3, hjust = 0) +
    scale_shape_manual(values = c(`TRUE` = 18, `FALSE` = 16)) +
    scale_x_log10(limits = c(0.6, 30)) +
    labs(x = "Hazard ratio, High vs Low (log)", y = NULL,
         title = "Stage subgroups") + theme_publication()

  ## C/D: 外部 KM
  km_panel <- function(df, title) {
    df <- df %>% mutate(risk = factor(risk, c("Low", "High")))
    survfit2(Surv(os_t, os_s) ~ risk, data = df) |>
      ggsurvfit(linewidth = 0.7) +
      ggsurvfit::add_confidence_interval(alpha = 0.12) +
      ggsurvfit::add_pvalue(size = 2.6) +
      ggsurvfit::scale_ggsurvfit() +
      scale_colour_manual(values = cols_risk, name = NULL) +
      scale_fill_manual(values = cols_risk) +
      labs(x = "OS (months)", y = "Survival", title = title) +
      theme_publication()
  }
  pC <- km_panel(emtab, sprintf("E-MTAB-1980 (n=%d)", nrow(emtab)))
  pD <- km_panel(cptac, sprintf("CPTAC-3 (n=%d)", nrow(cptac)))

  ## E: C-index
  fe <- extsum %>% transmute(cohort, C = Cindex, lo = Cindex_lo, hi = Cindex_hi) %>%
    bind_rows(data.frame(cohort = "TCGA-KIRC (OOF)", C = 0.683, lo = NA, hi = NA)) %>%
    mutate(cohort = factor(cohort, rev(cohort)))
  pE <- ggplot(fe, aes(C, cohort)) +
    geom_vline(xintercept = 0.5, linetype = 2, colour = "grey50", linewidth = 0.4) +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.15, linewidth = 0.5,
                   colour = COL_MID, na.rm = TRUE) +
    geom_point(size = 2.4, colour = COL_MID) +
    geom_text(aes(label = sprintf("%.3f", C)), vjust = -0.9, size = 2.4) +
    scale_x_continuous(limits = c(0.45, 0.95)) +
    labs(x = "C-index", y = NULL, title = "C-index") + theme_publication()

  # Nomogram (base 图): 先画到临时 png, 再以 raster 嵌入组图 (绕开 grid 捕获)
  tmp_nom <- file.path(OUT, "_tmp_nomogram.png")
  png(tmp_nom, width = 1800, height = 900, res = 300)
  draw_nomogram(build_nomogram()$nom)
  dev.off()
  img <- png::readPNG(tmp_nom)
  pNom <- ggplot() +
    annotation_raster(raster = img, xmin = -Inf, xmax = Inf,
                      ymin = -Inf, ymax = Inf) +
    theme_void()
  # 患者卡: 两个真实 TCGA 示例 (nomogram 工作示例, 与 nomogram 并排)
  nb <- build_nomogram()
  fit_nom <- nb$fit
  sd$TAS_sd <- as.numeric(scale(sd$TAS))   # 预计算 (单行 scale 会得 0)
  hi_pt <- sd %>% filter(!is.na(stage_grp)) %>% arrange(desc(TAS_sd)) %>% slice(1)
  lo_pt <- sd %>% filter(!is.na(stage_grp)) %>% arrange(TAS_sd) %>% slice(1)
  mk_new <- function(pt) {
    pt %>% mutate(SEX = factor(SEX, c("Female","Male")),
                  stage_grp = factor(stage_grp, nb$stages))
  }
  lp_hi <- predict(fit_nom, newdata = mk_new(hi_pt), type = "lp")
  lp_lo <- predict(fit_nom, newdata = mk_new(lo_pt), type = "lp")
  pCards <- ggplot() +
    annotate("rect", xmin = 0.04, xmax = 0.96, ymin = 0.55, ymax = 0.97,
             fill = "#FCEBEA", colour = "#E64B35", linewidth = 0.6) +
    annotate("text", x = 0.5, y = 0.90, label = "Patient A (high risk)",
             size = 3.4, fontface = "bold", colour = "#A32D2D",
             family = "Arial") +
    annotate("text", x = 0.5, y = 0.80,
             label = paste0(hi_pt$AGE, " y, ", hi_pt$SEX, ", ",
                            hi_pt$stage_grp), size = 3, family = "Arial") +
    annotate("text", x = 0.5, y = 0.71,
             label = sprintf("TAS = %+.2f SD", hi_pt$TAS_sd), size = 3,
             family = "Arial") +
    annotate("text", x = 0.5, y = 0.62,
             label = sprintf("36-mo OS = %.1f%%",
                             100 * rms::Survival(fit_nom)(36, lp = lp_hi)),
             size = 4.4, fontface = "bold", colour = "#A32D2D",
             family = "Arial") +
    annotate("rect", xmin = 0.04, xmax = 0.96, ymin = 0.06, ymax = 0.48,
             fill = "#E6F1FB", colour = "#3C5488", linewidth = 0.6) +
    annotate("text", x = 0.5, y = 0.41, label = "Patient B (low risk)",
             size = 3.4, fontface = "bold", colour = "#0C447C",
             family = "Arial") +
    annotate("text", x = 0.5, y = 0.31,
             label = paste0(lo_pt$AGE, " y, ", lo_pt$SEX, ", ",
                            lo_pt$stage_grp), size = 3, family = "Arial") +
    annotate("text", x = 0.5, y = 0.22,
             label = sprintf("TAS = %+.2f SD", lo_pt$TAS_sd), size = 3,
             family = "Arial") +
    annotate("text", x = 0.5, y = 0.13,
             label = sprintf("36-mo OS = %.1f%%",
                             100 * rms::Survival(fit_nom)(36, lp = lp_lo)),
             size = 4.4, fontface = "bold", colour = "#0C447C",
             family = "Arial") +
    lims(x = c(0, 1), y = c(0, 1)) +
    labs(title = "Worked examples") +
    theme_void(base_family = "Arial") +
    theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5))
  pNomCard <- wrap_elements(full =
    (wrap_elements(full = pNom) | pCards) + plot_layout(widths = c(2.6, 1)))
  # Nomogram 收尾 (证据链走完后再给临床工具, panel F)
  (pA | pB) / (pC | pD | pE) / pNomCard +
    plot_layout(heights = c(1, 1, 1.2), widths = c(1, 1, 1, 1, 1))
  finalize(last_plot(), "Fig4_independent_external", 200, 205)
}

# =====================================================================
# Fig5 v3: A 哑铃 + B 棒棒糖 + C 云雨 + D 山脊 (全部 per-sample 真数据)
# =====================================================================
fig5_v3 <- function() {
  im <- immune %>% mutate(delta = mean_high - mean_low)

  ## A: 哑铃图 (28 细胞 High vs Low 均值)
  pA <- im %>% mutate(cell = factor(cell, rev(cell[order(delta)]))) %>%
    ggplot() +
    geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.3) +
    geom_segment(aes(y = cell, x = mean_low, xend = mean_high, yend = cell),
                 colour = "grey85", linewidth = 1.5) +
    geom_point(aes(x = mean_low, y = cell, colour = "Low"), size = 1.8) +
    geom_point(aes(x = mean_high, y = cell, colour = "High"), size = 1.8) +
    scale_colour_manual(values = cols_risk, name = NULL) +
    labs(x = "Mean ssGSEA fraction", y = NULL,
         title = "High vs Low (28 cells)") +
    theme_publication() + theme(axis.text.y = element_text(size = 5.5))

  ## B: 棒棒糖 rho
  pB <- im %>% arrange(spearman_rho) %>%
    mutate(cell = factor(cell, cell)) %>%
    ggplot(aes(spearman_rho, cell, fill = spearman_rho > 0)) +
    geom_col(width = 0.6, show.legend = FALSE) +
    geom_vline(xintercept = 0, colour = "grey30", linewidth = 0.4) +
    scale_fill_manual(values = c(`FALSE` = COL_LOW, `TRUE` = COL_HIGH)) +
    labs(x = "Spearman rho (TAS vs fraction)", y = NULL,
         title = "Correlation with TAS") +
    theme_publication() + theme(axis.text.y = element_text(size = 5.5))

  ## C: 云雨图 (top 8 显著细胞, per-sample)
  cd <- ss_grp %>%
    select(patientId, group, all_of(top_cells)) %>%
    pivot_longer(-c(patientId, group), names_to = "cell", values_to = "score") %>%
    mutate(cell = factor(cell, rev(top_cells)),
           group = factor(group, c("Low", "High")))
  pC <- ggplot(cd, aes(group, score, fill = group)) +
    ggdist::stat_halfeye(aes(colour = group), adjust = 0.9, width = 0.6,
                         .width = 0, alpha = 0.55, show.legend = FALSE) +
    ggbeeswarm::geom_beeswarm(aes(colour = group), size = 0.25, alpha = 0.5,
                              show.legend = FALSE) +
    facet_wrap(~cell, ncol = 4, labeller = labeller(cell = label_wrap_gen(12))) +
    scale_fill_manual(values = cols_risk, name = NULL) +
    scale_colour_manual(values = cols_risk, name = NULL) +
    labs(x = NULL, y = "ssGSEA score",
         title = "Top 8 differential cells (per-sample)") +
    theme_publication() + theme(axis.text.x = element_text(size = 6))

  (pA | pB) / pC +
    plot_layout(heights = c(1, 1.15))
  finalize(last_plot(), "Fig5_immune", 200, 175)
}

# =====================================================================
# Fig7 (新): 桑基图 分型 -> 分期 -> 风险组
# =====================================================================
fig7_sankey <- function() {
  ## ---- A: 桑基 (分型 -> 分期 -> 风险组) ----
  fl <- sd %>% filter(!is.na(stage_grp)) %>%
    mutate(cluster_lab = paste0("C", cluster),
           risk = factor(group, c("Low", "High"))) %>%
    count(cluster_lab, stage_grp, risk) %>%
    mutate(stage_grp = factor(stage_grp,
                              c("Stage I", "Stage II", "Stage III", "Stage IV")))
  pA <- ggplot(fl, aes(axis1 = cluster_lab, axis2 = stage_grp, axis3 = risk,
                       y = n)) +
    ggalluvial::geom_alluvium(aes(fill = risk), width = 1/6, alpha = 0.55) +
    ggalluvial::geom_stratum(width = 1/3, fill = "grey95", colour = "grey40",
                             linewidth = 0.3) +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)),
              size = 2.6, family = "Arial") +
    scale_fill_manual(values = cols_risk, name = NULL) +
    scale_x_discrete(breaks = c("1", "2", "3"),
                     labels = c("Subtype", "AJCC stage", "TAS risk"),
                     expand = c(0.08, 0.05)) +
    labs(y = "Patients", title = "Patient flow across stratifications") +
    theme_publication() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          panel.grid = element_blank())

  ## ---- B: 四分期逐期 KM small multiples ----
  km_one <- function(st) {
    d <- sd %>% filter(stage_grp == st) %>%
      mutate(group = factor(group, c("Low", "High")))
    if (length(unique(d$group)) < 2) return(NULL)
    pv <- tryCatch(survival::survdiff(Surv(time, status) ~ group, data = d),
                   error = function(e) NULL)
    ptxt <- if (!is.null(pv)) {
      p <- 1 - pchisq(pv$chisq, length(pv$n) - 1)
      sprintf("P = %.3f", p)
    } else ""
    survfit2(Surv(time, status) ~ group, data = d) |>
      ggsurvfit::ggsurvfit(linewidth = 0.55) +
      ggsurvfit::add_confidence_interval(alpha = 0.12) +
      ggsurvfit::scale_ggsurvfit() +
      scale_colour_manual(values = cols_risk, guide = "none") +
      scale_fill_manual(values = cols_risk, guide = "none") +
      annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.3,
               label = ptxt, size = 2.6, family = "Arial") +
      labs(title = sprintf("%s (n=%d)", st, nrow(d)), x = NULL, y = NULL) +
      theme_prism(base_size = 7.5, base_family = "Arial") +
      theme(plot.title = element_text(size = 8.5, face = "bold"),
            axis.text = element_text(size = 6),
            plot.margin = margin(2, 2, 1, 2, "mm"))
  }
  km_list <- lapply(c("Stage I", "Stage II", "Stage III", "Stage IV"), km_one)
  km_list <- Filter(Negate(is.null), km_list)
  pB <- wrap_elements(full = Reduce(`+`, km_list) &
                        theme(plot.margin = margin(2, 2, 1, 2, "mm")))

  ## ---- C: 临床特征平衡 (Age / Sex / Stage) ----
  dcl <- sd %>% filter(!is.na(stage_grp)) %>%
    mutate(group = factor(group, c("Low", "High")),
           stage_grp = factor(stage_grp,
                              c("Stage I", "Stage II", "Stage III", "Stage IV")))
  p_age <- ggplot(dcl, aes(group, AGE, fill = group)) +
    geom_boxplot(width = 0.5, outlier.size = 0.5, show.legend = FALSE) +
    scale_fill_manual(values = cols_risk, guide = "none") +
    facet_wrap(~"Age (years)") + labs(x = NULL, y = NULL)
  p_sex <- dcl %>% count(group, SEX) %>%
    group_by(group) %>% mutate(p = n / sum(n)) %>% ungroup() %>%
    ggplot(aes(group, p, fill = SEX)) +
    geom_col(width = 0.55) +
    facet_wrap(~"Sex") + scale_fill_manual(values = c(Female = COL_MID,
                                                      Male = "#2C5F8A"),
                                           guide = "none") +
    labs(x = NULL, y = NULL)
  p_stg <- dcl %>% count(group, stage_grp) %>%
    group_by(group) %>% mutate(p = n / sum(n)) %>% ungroup() %>%
    ggplot(aes(stage_grp, p, fill = stage_grp)) +
    geom_col(width = 0.7, show.legend = FALSE) +
    facet_wrap(~"AJCC stage") +
    scale_fill_manual(values = c("Stage I" = "#DEEBF7", "Stage II" = "#9ECAE1",
                                 "Stage III" = "#6BAED6", "Stage IV" = "#2171B5")) +
    labs(x = NULL, y = NULL)
  pC <- (p_age | p_sex | p_stg) +
    plot_annotation(title = "Clinical covariates balanced across TAS groups") &
    theme_prism(base_size = 7.5, base_family = "Arial") &
    theme(axis.text.x = element_text(size = 6.5),
          strip.text = element_text(size = 8, face = "bold"),
          axis.title.y = element_text(size = 7),
          plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
          plot.margin = margin(2, 2, 1, 2, "mm"))
  pC <- wrap_elements(full = pC)

  ## ---- D: 分期 x 风险组 事件率 ----
  dD <- dcl %>%
    group_by(stage_grp, group) %>%
    summarise(death_rate = mean(status), n = n(), .groups = "drop") %>%
    mutate(lab = sprintf("%.0f%%", 100 * death_rate))
  pD <- ggplot(dD, aes(stage_grp, death_rate, fill = group)) +
    geom_col(position = position_dodge(0.7), width = 0.62) +
    geom_text(aes(label = lab), position = position_dodge(0.7), vjust = -0.4,
              size = 2.4, family = "Arial") +
    scale_fill_manual(values = cols_risk, name = NULL) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 0.72)) +
    labs(x = NULL, y = "Death rate",
         title = "Death rate by stage") +
    theme_prism(base_size = 8, base_family = "Arial") +
    theme(axis.text.x = element_text(size = 6, angle = 20, hjust = 1),
          plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
          legend.position = "bottom", legend.key.size = unit(3, "mm"))

  ## ---- 组合: A 整行 / B 四分期 / C+D 收尾 ----
  ((wrap_elements(full = pA)) / pB / (pC | pD)) +
    plot_layout(heights = c(1, 0.85, 0.8))
  finalize(last_plot(), "Fig7_sankey", 200, 215)
}

fig3_v3(); fig4_nomogram(); fig4_v3(); fig5_v3(); fig7_sankey()
message("ALL V3 UPGRADES DONE -> ", OUT)