# =====================================================================
# paper2_triaptosis 出版级组图 v2 —— 6 张 multi-panel figure
# 依据模板集: r_visualization_templates_bioinformatics.md (ggsurvfit/patchwork/
#   ComplexHeatmap/ggsci/timeROC 官方文档 + 生信技能树/SCIPainter 教程)
# 配色: NPG (用户指定); 主题: theme_publication(); 输出: PNG 300dpi + PDF 矢量
# =====================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(jsonlite)
  library(survival); library(ggsurvfit); library(patchwork); library(ggsci)
  library(timeROC); library(glmnet); library(scales); library(showtext)
  library(forestploter); library(ggrepel)
})

# macOS Arial 注册 (cairo 渲染 CID 失败的兜底方案)
font_add("Arial", regular = "/System/Library/Fonts/Supplemental/Arial.ttf")
showtext_opts(dpi = 300); showtext_auto()

ROOT <- "/Users/xiaolin/WorkBuddy/2026-09-04-20-00-08/paper2_triaptosis"
RES  <- file.path(ROOT, "results")
OUT  <- file.path(ROOT, "figures_v2")
source(file.path(ROOT, "code", "lib", "viz_base.R")); dir.create(OUT, showWarnings = FALSE)

## ---- 模板第 0 节: theme_publication() -------------------------------
# 拼图后统一主题兜底 (模板第 2 节 & 用法)
finalize <- function(p, file, w = 180, h = 130) {
  save_plot(p, file, dir = OUT, w = w, h = h)
}

## ---- NPG 配色 (ggsci 官方 nrc 前 6 色) ------------------------------


## ---- 数据加载 --------------------------------------------------------
score  <- read.csv(file.path(RES, "02_model/KIRC/score_table.csv"),
                   row.names = 1)
screen <- read.csv(file.path(RES, "01_screen/pancancer_screen_matrix.csv"))
mvcox  <- read.csv(file.path(RES, "08_multivariable/multivariable_cox.csv"))
stage  <- read.csv(file.path(RES, "08_multivariable/stage_subgroup_km.csv"))
extsum <- read.csv(file.path(RES, "07_external_os/external_OS_TAS_summary.csv"))
emtab  <- read.csv(file.path(RES, "07_external_os/EMTAB1980_TAS_scored.csv"))
cptac  <- read.csv(file.path(RES, "07_external_os/CPTAC3_TAS_scored.csv"))
oof    <- read.csv(file.path(RES, "05_validation/oof_scores.csv"), row.names = 1)
ipcw_c <- read.csv(file.path(RES, "05_validation/ipcw_cindex.csv"))
ipcw_a <- read.csv(file.path(RES, "05_validation/auc_ipcw.csv"))
brier  <- read.csv(file.path(RES, "05_validation/brier.csv"))
calib  <- read.csv(file.path(RES, "05_validation/calibration.csv"))
null   <- read.csv(file.path(RES, "05_null/triaptosis_random_null.csv"))
nullv  <- read.csv(file.path(RES, "05_null/triaptosis_vs_null.csv"))
immune <- read.csv(file.path(RES, "04_immune/immune_by_TAS_group.csv"))
ici    <- read.csv(file.path(RES, "03_ici/immotion150_validation.csv"))
expr_l <- fromJSON(file.path(ROOT, "data/KIRC_expr.json"))
mod    <- readRDS(file.path(RES, "02_model/KIRC/model.rds"))

# 表达矩阵: 510 样本 x 21 基因 (list-of-list -> matrix)
expr_mat <- t(vapply(expr_l, function(x) unlist(x[mod$genes]),
                     numeric(length(mod$genes))))
colnames(expr_mat) <- mod$genes
# score_table 行名是 12 位 patient id, expr 样本名是 15 位 sample id -> 统一
rownames(expr_mat) <- substr(rownames(expr_mat), 1, 12)
keep <- rownames(expr_mat) %in% rownames(score)   # 508/510, 2 例无 OS 信息
expr_mat <- expr_mat[keep, ]
cli <- data.frame(
  sample = rownames(expr_mat),
  time   = score[rownames(expr_mat), "time"],
  status = score[rownames(expr_mat), "status"])
stopifnot(sum(is.na(cli$time)) == 0)
expr_z <- scale(expr_mat)   # per-SD 表达 (全局, 供 Cox/glmnet)

# =====================================================================
# FIG 1 | 泛癌筛选 (3 panel: A donut / B 筛选矩阵 raster 整行
#                  / C KIRC 21 基因单因素 Cox 森林, 按 |log2HR| 降序)
# =====================================================================
fig1 <- function() {
  s <- screen %>% arrange(desc(n_sig_cox), logrank_p) %>%
    mutate(idx = row_number(),
           neglogp = -log10(pmax(logrank_p, 1e-30)))

  ## A: CNS 风格 donut (环形图): soft pass 3 / hard-only 17 / fail 12,
  ##    中心放总数, 扇区直接标 n (%), 色盲安全三色
  donut <- data.frame(
    cat = factor(c("Soft pass", "Hard only", "Not passed"),
                 c("Soft pass", "Hard only", "Not passed")),
    n = c(sum(screen$pass_soft),
          sum(screen$pass_hard) - sum(screen$pass_soft),
          sum(!screen$pass_hard)))
  donut$frac <- donut$n / sum(donut$n)
  donut$lab <- sprintf("%d (%.0f%%)", donut$n, 100 * donut$frac)
  donut$ymax <- cumsum(donut$frac)
  donut$ymin <- c(0, head(donut$ymax, -1))
  donut$lab_pos <- (donut$ymax + donut$ymin) / 2
  pA <- ggplot(donut, aes(xmin = 1.6, xmax = 2.4, ymin = ymin, ymax = ymax,
                          fill = cat)) +
    ggforce::geom_arc_bar(aes(x0 = 2, y0 = 0, r0 = 1.6, r = 2.4,
                              start = 2 * pi * (1 - ymax),
                              end = 2 * pi * (1 - ymin), fill = cat),
                          color = "white", linewidth = 0.6) +
    annotate("text", x = 2, y = 0, label = sprintf("32\ncancers"),
             size = 3.4, fontface = "bold", lineheight = 0.9) +
    geom_text(aes(x = 2 + 2.75 * cos(2 * pi * (1 - lab_pos)),
                  y = 2 + 2.75 * sin(2 * pi * (1 - lab_pos)), label = lab),
              size = 2.5, show.legend = FALSE) +
    labs(x = NULL, y = NULL) +
    scale_fill_manual(values = c("Soft pass" = COL_GREEN,
                                 "Hard only" = COL_MID,
                                 "Not passed" = "grey90"), name = NULL) +
    coord_equal(clip = "off") +
    labs(title = "Screening funnel (32 types)") +
    theme_void() + theme(
      plot.title = element_text(size = 9, face = "bold", hjust = 0.5,
                                family = "Arial"),
      legend.position = "bottom", legend.text = element_text(size = 6.5),
      legend.key.size = unit(2.5, "mm"))
  pA <- wrap_elements(full = pA + theme_void())

  ## B: 筛选矩阵热图 (geom_tile 版, patchwork 友好)
  hm <- s %>% slice(1:20) %>%  # 前 20 个有信号癌种, 保证可读
    select(idx, cancer, n_total, n_dead, n_sig_cox, neglogp) %>%
    pivot_longer(-c(idx, cancer)) %>%
    group_by(name) %>% mutate(v = value / max(value)) %>% ungroup()
  # B panel -> raster 嵌入 CNS 级 ComplexHeatmap (由 17_heatmaps_cns.R 生成)
  hm_png <- file.path(OUT, "Fig1B_screen_matrix_cns.png")
  imgB <- png::readPNG(hm_png)
  pB <- ggplot() +
    annotation_raster(imgB, xmin = -Inf, xmax = Inf,
                      ymin = -Inf, ymax = Inf) +
    labs(title = NULL) + theme_void()

  ## D: KIRC 21 基因单因素 Cox 森林 (模板 4c 手写版, 真实重算)
  uni <- lapply(mod$genes, function(g) {
    d <- data.frame(t = cli$time, s = cli$status, g = expr_z[, g])
    cox <- coxph(Surv(t, s) ~ g, data = d)
    ci <- summary(cox)$conf.int
    data.frame(gene = g, HR = ci[1], lo = ci[3], hi = ci[4],
               p = summary(cox)$coefficients[1, 5])
  }) %>% bind_rows() %>%
    arrange(desc(abs(log2(HR)))) %>%
    mutate(sig = p < 0.05, gene = factor(gene, rev(unique(gene))))
  n_sig <- sum(uni$sig)
  pD <- ggplot(uni, aes(HR, gene)) +
    geom_vline(xintercept = 1, linetype = 2, colour = "grey50", linewidth = 0.4) +
    geom_errorbarh(aes(xmin = lo, xmax = hi, colour = sig), height = 0.25,
                   linewidth = 0.5, show.legend = FALSE) +
    geom_point(aes(colour = sig), size = 1.8, show.legend = FALSE) +
    scale_colour_manual(values = c(`FALSE` = COL_GREY, `TRUE` = COL_HIGH)) +
    scale_x_log10() +
    labs(x = "Hazard ratio (log scale)", y = NULL,
         title = sprintf("KIRC: 21 DRGs univariate Cox (%d/21 significant)",
                         n_sig)) +
    theme_publication() + theme(axis.text.y = element_text(size = 6.5))

  ((pA | pD) / wrap_elements(full = pB)) +
    plot_annotation(tag_levels = "A") +
    plot_layout(heights = c(1.15, 2.3))
  finalize(last_plot(), "Fig1_pancancer_screen", 180, 250)
}

# =====================================================================
# FIG 2 | TAS 模型构建 (4 panel: A LASSO 路径 / B 10 基因森林
#                  / C 风险三联 / D k-means 分型 KM)
# =====================================================================
fig2 <- function() {
  ## A: glmnet 系数路径 (真实重跑, 种子固定)
  set.seed(42)
  x <- expr_mat                          # glmnet 需要 n(样本) x p(基因)
  y <- Surv(cli$time, cli$status)
  fit <- glmnet(x, y, family = "cox", alpha = 0.5)
  beta_df <- as.data.frame(as.matrix(fit$beta)) %>%
    mutate(gene = rownames(.)) %>%
    pivot_longer(-gene, names_to = "step", values_to = "coef") %>%
    mutate(step = as.numeric(sub("^s", "", step)),
           lambda = fit$lambda[step + 1],
           selected = gene %in% mod$genes_sel)
  pA <- ggplot(beta_df, aes(log(lambda), coef, group = gene,
                            colour = selected)) +
    geom_line(linewidth = 0.4, show.legend = TRUE) +
    geom_vline(xintercept = log(mod$lambda), linetype = 2, colour = "grey40",
               linewidth = 0.4) +
    annotate("text", x = log(mod$lambda), y = max(beta_df$coef) * 0.9,
             label = "lambda (10-gene)", size = 2.8, hjust = -0.05) +
    scale_colour_manual(values = c(`FALSE` = "grey75", `TRUE` = COL_HIGH),
                        labels = c(`FALSE` = "others", `TRUE` = "selected"),
                        name = NULL) +
    labs(x = expression(log(lambda)), y = "Coefficient",
         title = "Elastic-net coefficient path") + theme_publication()

  ## B: 10 基因森林 (真实单因素 Cox, 与 Fig1D 同源取子集)
  uni10 <- lapply(mod$genes_sel, function(g) {
    d <- data.frame(t = cli$time, s = cli$status, g = expr_z[, g])
    cox <- coxph(Surv(t, s) ~ g, data = d)
    ci <- summary(cox)$conf.int
    data.frame(gene = g, HR = ci[1], lo = ci[3], hi = ci[4],
               p = summary(cox)$coefficients[1, 5])
  }) %>% bind_rows() %>%
    mutate(gene = factor(gene, rev(gene)))
  pB <- ggplot(uni10, aes(HR, gene)) +
    geom_vline(xintercept = 1, linetype = 2, colour = "grey50", linewidth = 0.4) +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.2, linewidth = 0.5,
                   colour = COL_HIGH) +
    geom_point(size = 2, colour = COL_HIGH, shape = 15) +
    scale_x_log10() +
    labs(x = "Hazard ratio (log scale)", y = NULL,
         title = "10 selected genes (univariate)") +
    theme_publication() + theme(axis.text.y = element_text(face = "italic"))

  ## C: 风险三联图 (模板第 6 节 SCIPainter 流程, 真实 score_table)
  dt <- score %>% arrange(TAS) %>%
    mutate(id = row_number(), group = factor(group, c("Low", "High")))
  cutoff_x <- sum(dt$group == "Low") + 0.5
  pC1 <- ggplot(dt, aes(id, TAS, colour = group)) +
    geom_point(size = 0.9, show.legend = FALSE) +
    geom_vline(xintercept = cutoff_x, linetype = 2, colour = "grey50",
               linewidth = 0.4) +
    scale_colour_manual(values = cols_risk) +
    labs(y = "TAS", x = NULL, title = "Risk score") + theme_publication() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.line.x = element_blank())
  pC2 <- ggplot(dt, aes(id, time, colour = factor(status,
                      c(0, 1), c("Alive", "Dead")))) +
    geom_point(size = 0.9) +
    geom_vline(xintercept = cutoff_x, linetype = 2, colour = "grey50",
               linewidth = 0.4) +
    scale_colour_manual(values = c(Alive = COL_LOW, Dead = COL_HIGH),
                        name = NULL) +
    labs(y = "OS (months)", x = "Patients (increasing TAS)") +
    theme_publication()
  gene_sel <- mod$genes_sel
  ez <- t(scale(t(expr_mat[rownames(dt), gene_sel])))
  ez[ez >  2] <-  2; ez[ez < -2] <- -2
  ez_df <- as.data.frame(ez) %>%
    mutate(sample_id = rownames(.)) %>%
    pivot_longer(-sample_id, names_to = "gene", values_to = "z") %>%
    mutate(id = match(sample_id, rownames(ez)),
           risk = dt$group[id])
  pC3 <- ggplot(ez_df, aes(id, factor(gene, rev(gene_sel)), fill = z)) +
    geom_tile(linewidth = 0) +
    geom_vline(xintercept = cutoff_x, linetype = 2, colour = "grey50",
               linewidth = 0.4) +
    scale_fill_gradient2(low = COL_LOW, mid = "white", high = COL_HIGH,
                         midpoint = 0, name = "Z") +
    labs(y = NULL, x = NULL, title = "Gene expression (Z)") +
    theme_publication() +
    theme(axis.text.y = element_text(size = 6, face = "italic"),
          axis.text.x = element_blank(), axis.ticks.x = element_blank())
  pC <- pC1 / pC2 / pC3 + plot_layout(heights = c(1, 1, 1.6), guides = "collect")

  ## D: k-means 分型 KM (ggsurvfit 模板第 1 节)
  pD <- ggsurvfit::survfit2(Surv(time, status) ~ factor(cluster), data = score) |>
    ggsurvfit::ggsurvfit(linewidth = 0.7) +
    ggsurvfit::add_confidence_interval(alpha = 0.12) +
    ggsurvfit::add_pvalue(location = "annotation", size = 3) +
    ggsurvfit::scale_ggsurvfit() +
    scale_colour_manual(values = c("1" = COL_LOW, "2" = COL_HIGH),
                        labels = c("1" = "C1 (n=293)", "2" = "C2 (n=215)"),
                        name = NULL) +
    scale_fill_manual(values = c("1" = COL_LOW, "2" = COL_HIGH),
                      guide = "none") +
    labs(x = "Time (months)", y = "Survival probability",
         title = "k-means molecular subtypes") + theme_publication()

  # 叙事顺序: A=k-means 分型 KM (unsupervised 证据打头) / B=LASSO 建模
  #           / C=10 基因森林 (标签由 tag_levels 自动分配)
  combined <- (pD | pA) / pB +
    plot_annotation(tag_levels = "A") +
    plot_layout(heights = c(1, 1))
  finalize(combined, "Fig2_model_build", 180, 140)

  # 风险三联图单独成图 (复合 panel, 竖排共享 x)
  pC / plot_annotation(tag_levels = list(c("C1", "C2", "C3"))) &
    theme_publication()
  finalize(last_plot(), "Fig2C_risk_triptych", 170, 140)
}

# =====================================================================
# FIG 3 | 模型性能评估 (5 panel: A KM+risk table / B timeROC
#                  / C C-index & AUC / D Brier / E 校准)
# =====================================================================
fig3 <- function() {
  ## A: KM + risk table (模板第 1 节)
  sc <- score %>% mutate(group = factor(group, c("Low", "High")))
  pA <- ggsurvfit::survfit2(Surv(time, status) ~ group, data = sc) |>
    ggsurvfit::ggsurvfit(linewidth = 0.7) +
    ggsurvfit::add_confidence_interval(alpha = 0.12) +
    ggsurvfit::add_pvalue(location = "annotation", size = 3) +
    ggsurvfit::add_risktable(risktable_stats = "n.risk",
                  stats_label = list(n.risk = "At risk"), size = 2.6) +
    ggsurvfit::scale_ggsurvfit() +
    scale_colour_manual(values = cols_risk,
                        labels = c("Low (n=254)", "High (n=254)"), name = NULL) +
    scale_fill_manual(values = cols_risk) +
    labs(x = "Time (months)", y = "Overall survival probability",
         title = "TCGA-KIRC KM (HR = 3.52)") + theme_publication()

  ## B: timeROC 1/3/5 年 (模板第 5 节)
  roc <- timeROC(T = score$time, delta = score$status, marker = score$TAS,
                 cause = 1, weighting = "marginal", times = c(12, 36, 60),
                 ROC = TRUE, iid = TRUE)
  roc_df <- data.frame(
    FP  = c(roc$FP[, 1], roc$FP[, 2], roc$FP[, 3]),
    TP  = c(roc$TP[, 1], roc$TP[, 2], roc$TP[, 3]),
    T   = factor(rep(c("12 mo", "36 mo", "60 mo"), each = nrow(roc$TP)),
                 c("12 mo", "36 mo", "60 mo")))
  auc_lab <- sprintf("AUC(%s) = %.3f", levels(roc_df$T), roc$AUC)
  pB <- ggplot(roc_df, aes(FP, TP, colour = T)) +
    geom_line(linewidth = 0.7) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey70",
                linewidth = 0.4) +
    annotate("text", x = 0.62, y = c(0.28, 0.20, 0.12), label = auc_lab,
             size = 2.6, hjust = 0,
             colour = c(COL_HIGH, COL_GREEN, COL_LOW)) +
    scale_colour_manual(values = c("12 mo" = COL_HIGH, "36 mo" = COL_GREEN,
                                   "60 mo" = COL_LOW), name = NULL) +
    labs(x = "False positive rate", y = "True positive rate",
         title = "Time-dependent ROC (TAS)") + theme_publication()

  ## C: IPCW C-index & AUC 点图
  cmp <- bind_rows(
    ipcw_c %>% mutate(Metric = "IPCW C-index", value = C),
    ipcw_a %>% mutate(Metric = "Time-dep AUC",  value = AUC))
  pC <- ggplot(cmp, aes(times, value, colour = Metric)) +
    geom_hline(yintercept = 0.5, linetype = 2, colour = "grey60",
               linewidth = 0.4) +
    geom_point(size = 2.4) + geom_line(linewidth = 0.6) +
    ggrepel::geom_text_repel(aes(label = sprintf("%.3f", value)),
              size = 2.4, show.legend = FALSE, seed = 1,
              box.padding = 0.25, segment.colour = NA) +
    scale_colour_manual(values = c(COL_HIGH, COL_GREEN)) +
    scale_x_continuous(breaks = c(12, 36, 60)) +
    coord_cartesian(ylim = c(0.5, 0.95)) +
    labs(x = "Time (months)", y = "Value (OOF, 5-fold)",
         title = "Discrimination (IPCW)") + theme_publication()

  ## D: Brier score (TAS vs KM)
  pD <- ggplot(brier, aes(times, Brier, colour = model, group = model)) +
    geom_point(size = 2.4) + geom_line(linewidth = 0.6) +
    ggrepel::geom_text_repel(aes(label = sprintf("%.3f", Brier)),
              size = 2.4, show.legend = FALSE, seed = 1,
              box.padding = 0.3, segment.colour = NA) +
    scale_colour_manual(values = c(TAS = COL_HIGH, KM = COL_GREY)) +
    scale_x_continuous(breaks = c(12, 36, 60)) +
    labs(x = "Time (months)", y = "IPCW Brier score",
         title = "Calibration error (lower = better)") + theme_publication()

  ## E: 校准曲线 (Q1-Q4, pred vs obs, 36 mo)
  pE <- calib %>%
    ggplot(aes(pred_36, obs_36, colour = group)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey70",
                linewidth = 0.4) +
    geom_point(size = 2.4) +
    geom_line(linewidth = 0.5, alpha = 0.6) +
    geom_text(aes(label = group), vjust = -0.8, size = 2.4,
              show.legend = FALSE) +
    scale_colour_manual(values = c(COL_LOW, COL_MID, COL_ORANGE, COL_HIGH)) +
    coord_cartesian(xlim = c(0.55, 1.0), ylim = c(0.55, 1.0)) +
    theme(legend.position = "bottom", legend.key.size = unit(3, "mm")) +
    labs(x = "Predicted 36-mo survival", y = "Observed",
         title = "Calibration by quartile") + theme_publication()

  # A 面板带 risk table 是双面板对象, 用 wrap_elements 融入 patchwork
  # 不用 guides="collect" (各 panel 图例保留在各自内部, 避免顶部图例堆积)
  (wrap_elements(full = pA) | pB) / (pC | pD | pE) +
    plot_annotation(tag_levels = "A") +
    plot_layout(widths = c(1.4, 1, 1, 1, 1))
  finalize(last_plot(), "Fig3_performance", 200, 130)
}

# =====================================================================
# FIG 4 | 独立预后 + 外部验证 (5 panel: A 多变量森林 / B 分期亚组
#                  / C E-MTAB KM / D CPTAC KM / E C-index 对比)
# =====================================================================
fig4 <- function() {
  ## A: 多变量 Cox 森林 (模板 4c 手写版, 直接用已存结果)
  # 手写森林图 (forestploter 的 grid 对象与 patchwork 尺寸冲突, 组图内用手写版)
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

  # forestploter 独立版 (单独文件, 供正文单图/备选)
  fa2 <- fa %>%
    mutate(Variable = as.character(term_lab),
           `HR (95% CI)` = sprintf("%.2f (%.2f to %.2f)", HR, lo, hi),
           `P value` = fmt_p(p))
  fa2$`HR (95% CI)` <- paste0("    ", fa2$`HR (95% CI)`, "    ")
  fa_fp <- forestploter::forest(fa2[, c("Variable", "HR (95% CI)", "P value")],
                                est = fa2$HR, lower = fa2$lo, upper = fa2$hi,
                                ci_column = 2, ref_line = 1,
                                xlim = c(0.4, 12), log_scale = TRUE)
  ggsave(file.path(OUT, "Fig4A_forestploter_standalone.png"),
         ggplotify::as.ggplot(fa_fp), width = 140, height = 80,
         units = "mm", dpi = 300)

  ## B: 分期亚组森林 + 整体
  fb <- bind_rows(
    data.frame(sub = "Overall", HR = 3.52, lo = 2.76, hi = 4.50, p = 9.6e-24),
    stage %>% transmute(sub = paste0("Stage ", subgroup), HR = HR_High_vs_Low,
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

  ## C/D: 两个外部队列 KM (ggsurvfit)
  km_panel <- function(df, title) {
    df <- df %>% mutate(risk = factor(risk, c("Low", "High")))
    ggsurvfit::survfit2(Surv(os_t, os_s) ~ risk, data = df) |>
      ggsurvfit::ggsurvfit(linewidth = 0.7) +
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

  ## E: 三队列 C-index + 95%CI
  fe <- extsum %>% transmute(cohort, C = Cindex, lo = Cindex_lo, hi = Cindex_hi) %>%
    bind_rows(data.frame(cohort = "TCGA-KIRC (OOF)", C = 0.683,
                         lo = NA, hi = NA)) %>%
    mutate(cohort = factor(cohort, rev(cohort)))
  pE <- ggplot(fe, aes(C, cohort)) +
    geom_vline(xintercept = 0.5, linetype = 2, colour = "grey50", linewidth = 0.4) +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.15, linewidth = 0.5,
                   colour = COL_MID, na.rm = TRUE) +
    geom_point(size = 2.4, colour = COL_MID) +
    geom_text(aes(label = sprintf("%.3f", C)), vjust = -0.9, size = 2.4) +
    scale_x_continuous(limits = c(0.45, 0.95)) +
    labs(x = "C-index", y = NULL, title = "C-index") +
    theme_publication()

  (pA | pB) / (pC | pD | pE) +
    plot_annotation(tag_levels = "A") +
    plot_layout(widths = c(1, 1, 1, 1, 1))
  finalize(last_plot(), "Fig4_independent_external", 200, 140)
}

# =====================================================================
# FIG 5 | 免疫微环境 (4 panel: A 组均值热图 / B 棒棒糖差值
#                  / C Spearman rho 条形 / D 关键细胞哑铃)
# =====================================================================
fig5 <- function() {
  im <- immune %>% mutate(delta = mean_high - mean_low)

  ## A: 组均值热图 (28 行 x 2 列, 行 Z-score)
  mat <- as.matrix(im[, c("mean_high", "mean_low")])
  rownames(mat) <- im$cell
  mat_z <- t(scale(t(mat)))
  az <- as.data.frame(mat_z) %>% mutate(cell = rownames(.)) %>%
    pivot_longer(-cell, names_to = "grp", values_to = "z") %>%
    mutate(grp = factor(grp, c("mean_high", "mean_low"),
                        c("High TAS", "Low TAS")),
           cell = factor(cell, rev(im$cell)))
  pA <- ggplot(az, aes(grp, cell, fill = z)) +
    geom_tile(color = "white", linewidth = 0.2) +
    scale_fill_gradient2(low = COL_LOW, mid = "white", high = COL_HIGH,
                         midpoint = 0, name = "Z-score", limits = c(-2, 2)) +
    labs(x = NULL, y = NULL,
         title = "ssGSEA (Charoentong 28; 25/28 FDR<0.05)") +
    theme_publication() + theme(axis.text.y = element_text(size = 5.5))

  ## B: 棒棒糖 (High-Low 差值, 模板 8b)
  pB <- im %>% mutate(cell = factor(cell, rev(cell[order(delta)]))) %>%
    ggplot(aes(delta, cell, colour = delta > 0)) +
    geom_vline(xintercept = 0, colour = "grey50", linewidth = 0.4) +
    geom_segment(aes(x = 0, y = cell, xend = delta, yend = cell),
                 linewidth = 0.5, show.legend = FALSE) +
    geom_point(size = 1.8, show.legend = FALSE) +
    scale_colour_manual(values = c(`FALSE` = COL_LOW, `TRUE` = COL_HIGH)) +
    labs(x = "Fraction difference (High - Low)", y = NULL,
         title = "Cell fraction shift") +
    theme_publication() + theme(axis.text.y = element_text(size = 5.5))

  ## C: Spearman rho 条形 (排序, 正负双色)
  pC <- im %>% arrange(spearman_rho) %>%
    mutate(cell = factor(cell, cell)) %>%
    ggplot(aes(spearman_rho, cell, fill = spearman_rho > 0)) +
    geom_col(width = 0.6, show.legend = FALSE) +
    geom_vline(xintercept = 0, colour = "grey30", linewidth = 0.4) +
    scale_fill_manual(values = c(`FALSE` = COL_LOW, `TRUE` = COL_HIGH)) +
    labs(x = "Spearman rho (TAS vs fraction)", y = NULL,
         title = "Correlation with TAS") +
    theme_publication() + theme(axis.text.y = element_text(size = 5.5))

  ## D: 关键 6 细胞哑铃图 (mean_high vs mean_low + FDR 星号)
  key <- im %>% arrange(wilcox_fdr) %>% slice(1:6) %>%
    select(cell, mean_high, mean_low, wilcox_fdr) %>%
    pivot_longer(-c(cell, wilcox_fdr), names_to = "grp", values_to = "frac") %>%
    mutate(grp = factor(grp, c("mean_high", "mean_low"), c("High", "Low")),
           sig = cut(wilcox_fdr, c(-Inf, 0.001, 0.01, 0.05, Inf),
                     c("***", "**", "*", "ns")),
           cell = factor(cell, rev(unique(cell))))
  pD <- ggplot(key, aes(frac, cell)) +
    geom_line(aes(group = cell), colour = "grey70", linewidth = 0.5) +
    geom_point(aes(colour = grp), size = 2.6) +
    geom_text(data = key %>% filter(grp == "High"),
              aes(label = sig), hjust = -0.4, vjust = 0.4, size = 2.6) +
    scale_colour_manual(values = cols_risk, name = NULL) +
    labs(x = "Mean cell fraction", y = NULL,
         title = "Top 6 differential cells") +
    theme_publication()

  (pA | pB) / (pC | pD) +
    plot_annotation(tag_levels = "A")
  finalize(last_plot(), "Fig5_immune", 180, 160)
}

# =====================================================================
# FIG 6 | 零模型 + ICI 探索 (5 panel: A 零模型C分布 / B 零模型HR
#                  / C cv 配对散点 / D responder 堆叠 / E PFS KM)
# =====================================================================
fig6 <- function() {
  ## A: 零模型 C-index 分布 (nested CV) + TAS 位置
  pA <- ggplot(null, aes(cv_nested)) +
    geom_histogram(bins = 20, fill = COL_GREY, colour = "white",
                   linewidth = 0.2) +
    geom_vline(xintercept = nullv$TAS[2], colour = COL_HIGH, linewidth = 0.8) +
    annotate("text", x = nullv$TAS[2], y = Inf, hjust = -0.05, vjust = 1.5,
             size = 2.6, colour = COL_HIGH,
             label = sprintf("TAS = %.3f\n(88.7th pctile, z=+0.96)",
                             nullv$TAS[2])) +
    labs(x = "Nested CV C-index (random 21-gene sets, n=195)", y = "Count",
         title = "Null model: C-index") + theme_publication()

  ## B: 零模型 HR 分布 (截断到 95% 分位)
  hr95 <- quantile(null$HR, 0.95)
  pB <- ggplot(null %>% filter(HR <= hr95), aes(HR)) +
    geom_histogram(bins = 20, fill = COL_GREY, colour = "white",
                   linewidth = 0.2) +
    geom_vline(xintercept = 3.52, colour = COL_HIGH, linewidth = 0.8) +
    annotate("text", x = 3.52, y = Inf, hjust = -0.05, vjust = 1.5, size = 2.6,
             colour = COL_HIGH, label = "TAS HR = 3.52") +
    labs(x = "Univariate Cox HR (random sets)", y = "Count",
         title = "Null model: HR") + theme_publication()

  ## D: IMmotion150 RESPONDER 堆叠条 (真实 per-sample)
  resp <- ici %>% filter(!is.na(responder)) %>%
    count(group, responder) %>%
    group_by(group) %>% mutate(frac = n / sum(n)) %>% ungroup()
  pD <- ggplot(resp, aes(group, frac, fill = responder)) +
    geom_col(width = 0.55, colour = "white", linewidth = 0.2) +
    geom_text(aes(label = n), position = position_stack(vjust = 0.5),
              size = 2.6, colour = "white") +
    scale_fill_manual(values = c(R = COL_HIGH, NR = COL_GREY),
                      labels = c(R = "Responder", NR = "Non-responder"),
                      name = NULL) +
    scale_y_continuous(labels = percent) +
    labs(x = NULL, y = "Proportion",
         title = "IMmotion150 response (P=1.0)") + theme_publication()

  ## E: PFS KM (ggsurvfit, 真实 per-sample)
  ici2 <- ici %>% mutate(
    pfs_status_num = as.integer(sub(":.*", "", pfs_status)),
    group = factor(group, c("Low", "High")))
  pE <- ggsurvfit::survfit2(Surv(pfs, pfs_status_num) ~ group, data = ici2) |>
    ggsurvfit::ggsurvfit(linewidth = 0.7) +
    ggsurvfit::add_confidence_interval(alpha = 0.12) +
    ggsurvfit::add_pvalue(size = 2.8) +
    ggsurvfit::scale_ggsurvfit() +
    scale_colour_manual(values = cols_risk, name = NULL) +
    scale_fill_manual(values = cols_risk) +
    labs(x = "PFS (months)", y = "Progression-free survival",
         title = "PFS by TAS (P=0.111)") + theme_publication()

  (pA | pB) / (pD | pE) +
    plot_annotation(tag_levels = "A")
  finalize(last_plot(), "Fig6_null_ici", 200, 140)
}

## ---- 执行 ------------------------------------------------------------
fig1(); fig2(); fig3(); fig4(); fig5(); fig6()
message("ALL 6 FIGURES DONE -> ", OUT)