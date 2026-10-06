# =====================================================================
# 27g_fig4_transport.R —— Fig4 独立性与外部验证（投稿版）
# A 多变量 Cox（表内 HR/CI）
# B 分期亚组，使用表内 CI，不再用 HR×0.65/1.45
# C–D 外部 KM，纵轴 0–1，标题含 n 与事件数
# E C-index；TCGA 36 个月 OOF 点注明 CI not estimated
# 列线图单独成补充图，不附 0% OS 示例框。
# =====================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(survival)
  library(ggsurvfit); library(rms); library(showtext); library(patchwork)
})
ROOT <- "/Users/xiaolin/WorkBuddy/2026-09-04-20-00-08/paper2_triaptosis"
OUT  <- file.path(ROOT, "figures_pub")
RES  <- file.path(ROOT, "results")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
source(file.path(ROOT, "code", "lib", "viz_base.R"))

mvcox <- read.csv(file.path(RES, "08_multivariable/multivariable_cox.csv"))
grade_mv <- read.csv(file.path(RES, "08_multivariable/grade_multivariable_cox.csv"))
stage <- read.csv(file.path(RES, "08_multivariable/stage_subgroup_km.csv"))
extsum <- read.csv(file.path(RES, "07_external_os/external_OS_TAS_summary.csv"))
emtab <- read.csv(file.path(RES, "07_external_os/EMTAB1980_TAS_scored.csv"))
cptac <- read.csv(file.path(RES, "07_external_os/CPTAC3_TAS_scored.csv"))
ipcw_c <- read.csv(file.path(RES, "05_validation/ipcw_cindex.csv"))
sc <- read.csv(file.path(RES, "02_model/KIRC/score_table.csv"), row.names = 1)

# Panel A is the grade-adjusted model (n=500). Univariable TAS is refit on the same 500.
labs_g <- c(TAS = "TAS per SD, adjusted", AGE = "Age (per year)", SEXMALE = "Sex (male)",
            STAGEII = "Stage II", STAGEIII = "Stage III", STAGEIV = "Stage IV",
            GRADEG3 = "Grade 3", GRADEG4 = "Grade 4")
gr <- read.csv(file.path(ROOT, "data/external_validation/kirc_tcga_grade.csv"))
gr$GRADE <- ifelse(gr$GRADE %in% c("G1", "G2"), "G1-2", ifelse(gr$GRADE %in% c("G3", "G4"), gr$GRADE, NA))
clin <- read.csv(file.path(ROOT, "data/external_validation/kirc_tcga_pan_can_atlas_2018_clin_patient.csv"))
st <- gsub("^STAGE\\s*", "", toupper(trimws(clin$AJCC_PATHOLOGIC_TUMOR_STAGE)))
st[!st %in% c("I", "II", "III", "IV")] <- NA
dA <- sc %>% mutate(patientId = rownames(sc), TAS_sd = as.numeric(scale(TAS))) %>%
  inner_join(clin %>% transmute(patientId, STAGE = st, AGE = as.numeric(AGE), SEX), by = "patientId") %>%
  inner_join(gr %>% select(patientId, GRADE), by = "patientId") %>%
  filter(!is.na(STAGE), !is.na(GRADE), !is.na(AGE), time > 0)
uni <- summary(coxph(Surv(time, status) ~ TAS_sd, data = dA, ties = "efron"))$conf.int
fa <- bind_rows(
  data.frame(term_lab = "TAS per SD, univariable", HR = uni[1, 1], lo = uni[1, 3], hi = uni[1, 4]),
  grade_mv %>% transmute(term_lab = labs_g[term], HR, lo, hi)
) %>% mutate(term_lab = factor(term_lab, rev(c("TAS per SD, univariable", unname(labs_g)))))
pA <- ggplot(fa, aes(HR, term_lab)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.3, colour = "grey40") +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y",
                width = 0.18, linewidth = 0.4, colour = "grey30") +
  geom_point(size = 2, colour = "#E64B35") +
  geom_text(aes(x = hi, label = sprintf("%.2f (%.2f\u2013%.2f)", HR, lo, hi)),
            hjust = -0.08, size = 2.3) +
  scale_x_log10(limits = c(0.45, 40)) +
  labs(x = "Hazard ratio (log scale)", y = NULL,
       title = sprintf("With histologic grade (n=%d)", nrow(dA)), tag = "A") +
  theme_cns(base_size = 8) +
  theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
        plot.margin = margin(2, 8, 2, 2, "mm"))

fb <- stage %>%
  transmute(sub = sprintf("Stage %s (n=%d)", subgroup, n),
            HR = HR_High_vs_Low, lo = lo, hi = hi) %>%
  mutate(sub = factor(sub, rev(sub)))
pB <- ggplot(fb, aes(HR, sub)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.3, colour = "grey40") +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y",
                width = 0.18, linewidth = 0.45, colour = "#E64B35") +
  geom_point(size = 2.2, colour = "#E64B35") +
  geom_text(aes(x = hi, label = sprintf("%.2f (%.2f\u2013%.2f)", HR, lo, hi)),
            hjust = -0.08, size = 2.4) +
  scale_x_log10(limits = c(0.7, 14)) +
  labs(x = "HR, High vs Low (log scale)", y = NULL,
       title = "Stage subgroups", tag = "B") +
  theme_cns(base_size = 8) +
  theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
        plot.margin = margin(2, 6, 2, 2, "mm"))

# 外部 KM 原先没有 at-risk 表, 纵轴又用 0-1 小数, 与 Fig1C/Fig2B 不一致;
# 这里统一为百分数 + at-risk 表, 并用 wrap_elements 固定成块, 免除面板上方的大片空白。
km_panel <- function(df, meta, cohort, tag = NULL, step = 40) {
  df <- df %>% mutate(risk = factor(risk, c("Low", "High")))
  info <- meta %>% filter(cohort == !!cohort)
  br <- seq(0, ceiling(max(df$os_t, na.rm = TRUE) / step) * step, step)
  p <- survfit2(Surv(os_t, os_s) ~ risk, data = df) |>
    ggsurvfit(linewidth = 0.6) +
    add_confidence_interval(alpha = 0.12) +
    add_risktable(risktable_stats = "n.risk",
                  stats_label = list(n.risk = "At risk"), size = 2.1,
                  theme = km_risktable_theme(base_size = 6.5)) +
    scale_colour_manual(values = c(Low = "#3C5488", High = "#E64B35"),
                        labels = c("Low TAS", "High TAS"), name = NULL) +
    scale_fill_manual(values = c(Low = "#3C5488", High = "#E64B35"), guide = "none") +
    km_x_scale(breaks = br) +
    km_y_scale() +
    labs(x = "Time (months)", y = "Overall survival",
         title = sprintf("%s (n=%d, %d events)", cohort, info$n, info$events),
         tag = tag) +
    theme_cns(base_size = 8) +
    theme(legend.position = "top",
          legend.key.width = unit(4.5, "mm"),
          legend.text = element_text(size = 6.5),
          plot.title = element_text(size = 8, face = "bold", hjust = 0.5))
  wrap_elements(full = ggsurvfit_build(p))
}
pC <- km_panel(emtab, extsum, "E-MTAB-1980", tag = "C", step = 50)
pD <- km_panel(cptac, extsum, "CPTAC-3", tag = "D", step = 20)

# 没有 CI 的那一点改用空心符号并进图例, 原先的灰色小字 "CI not estimated"
# 会压在纵轴标签上。
c36 <- ipcw_c$C[ipcw_c$times == 36]
fe <- extsum %>%
  transmute(cohort, C = Cindex, lo = Cindex_lo, hi = Cindex_hi, has_ci = TRUE,
            lab = sprintf("%.3f", Cindex)) %>%
  bind_rows(data.frame(cohort = "TCGA OOF, 36 mo", C = c36, lo = NA_real_, hi = NA_real_,
                       has_ci = FALSE, lab = sprintf("%.3f", c36))) %>%
  mutate(cohort = factor(cohort, rev(c("E-MTAB-1980", "CPTAC-3", "TCGA OOF, 36 mo"))))
pE <- ggplot(fe, aes(y = cohort)) +
  geom_segment(aes(x = 0.5, xend = C), colour = "#4DBBD5", linewidth = 0.85,
               lineend = "butt") +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y",
                width = 0.18, linewidth = 0.4, colour = "#4DBBD5", na.rm = TRUE) +
  geom_point(aes(x = C, shape = has_ci), size = 2.6, colour = "#4DBBD5",
             fill = "white", stroke = 0.65) +
  geom_text(aes(x = C, label = lab), vjust = -1.35, hjust = 0.5, size = 2.3) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21),
                     breaks = c(TRUE, FALSE),
                     labels = c(`TRUE` = "with 95% CI", `FALSE` = "CI not estimated"),
                     name = NULL) +
  scale_x_continuous(limits = c(0.5, 1.02), breaks = seq(0.5, 1.0, 0.1),
                     expand = expansion(mult = c(0, 0.03))) +
  labs(x = "C-index (chance = 0.5)", y = NULL, title = "Discrimination", tag = "E") +
  theme_cns(base_size = 8) +
  theme(legend.position = "bottom", legend.direction = "vertical",
        legend.key.height = unit(2.8, "mm"),
        legend.spacing.y = unit(0.2, "mm"),
        legend.margin = margin(0.5, 0, 0, 0, "mm"),
        legend.text = element_text(size = 6.2),
        plot.title = element_text(size = 8.5, face = "bold", hjust = 0.5),
        plot.margin = margin(2, 4, 2, 2, "mm"))

fig4 <- (pA | pB) / (pC | pD | pE) +
  plot_layout(heights = c(0.85, 1)) &
  theme(plot.tag = element_text(face = "bold", size = 11),
        plot.tag.position = "topleft",
        plot.margin = margin(1, 3, 1, 1, "mm"))
save_pub(fig4, "Fig4_transport", OUT, w = 180, h = 165)

# ---- supplement nomogram, no worked-example cards ----------------------
clin <- read.csv(file.path(ROOT, "data/external_validation/kirc_tcga_pan_can_atlas_2018_clin_patient.csv")) %>%
  mutate(stage4 = trimws(gsub("STAGE", "", AJCC_PATHOLOGIC_TUMOR_STAGE, ignore.case = TRUE)),
         stage_grp = case_when(stage4 %in% "I" ~ "Stage I",
                               stage4 %in% "II" ~ "Stage II",
                               stage4 %in% "III" ~ "Stage III",
                               stage4 %in% "IV" ~ "Stage IV",
                               TRUE ~ NA_character_))
df1 <- sc %>%
  mutate(patientId = rownames(sc)) %>%
  left_join(clin %>% select(patientId, stage_grp, SEX, AGE), by = "patientId") %>%
  filter(!is.na(stage_grp), !is.na(AGE), !is.na(SEX)) %>%
  mutate(Sex = factor(SEX, c("Female", "Male")),
         Stage = factor(stage_grp, c("Stage I", "Stage II", "Stage III", "Stage IV")),
         TAS = as.numeric(scale(TAS)),
         Age = AGE)
units(df1$Age) <- "years"
dd <- datadist(df1)
assign("dd", dd, envir = .GlobalEnv)
options(datadist = "dd")
fit <- cph(Surv(time, status) ~ TAS + Age + Sex + Stage,
           data = df1, x = TRUE, y = TRUE, surv = TRUE, time.inc = 36)
surv_fn <- Survival(fit)
nom <- nomogram(fit,
                fun = list(function(x) surv_fn(36, lp = x),
                           function(x) surv_fn(60, lp = x)),
                funlabel = c("36-mo OS probability", "60-mo OS probability"),
                lp = FALSE, fun.at = c(0.9, 0.7, 0.5, 0.3, 0.1))
draw_nom <- function() {
  op <- par(mar = c(1.4, 0.5, 1.8, 0.5), xpd = NA, mgp = c(0, 0.4, 0), tcl = -0.2)
  on.exit(par(op), add = TRUE)
  plot(nom, xfrac = 0.28, cex.axis = 0.68, cex.var = 0.78, lmgp = 0.22,
       points.label = "Points", total.points.label = "Total points")
  title("Nomogram (TCGA-KIRC; TAS in SD units)", cex.main = 0.9, font.main = 1)
}
png(file.path(OUT, "FigS_nomogram.png"), width = 180, height = 110, units = "mm", res = 300)
draw_nom()
dev.off()
grDevices::quartz(file = file.path(OUT, "FigS_nomogram.pdf"), type = "pdf",
                  width = 180 / 25.4, height = 110 / 25.4)
draw_nom()
dev.off()
message(sprintf("nomogram n=%d; stage HR I-II %.2f (%.2f-%.2f); III-IV %.2f (%.2f-%.2f)",
                nrow(df1),
                stage$HR_High_vs_Low[1], stage$lo[1], stage$hi[1],
                stage$HR_High_vs_Low[2], stage$lo[2], stage$hi[2]))
message(sprintf("TCGA OOF C at 36 mo = %.3f (no CI in source table)", c36))
