#!/usr/bin/env Rscript
# =====================================================================
# 18_figures_claims.R — claim-mapped main figures (no fabricated CIs)
# Fig1 screen | Fig2 modular TAS | Fig3 utility | Fig4 independence
# Fig5 immune | Fig6 boundary | GA | S1 extras | S2 immune-clinical corr
# =====================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(jsonlite)
  library(survival); library(ggsurvfit); library(patchwork); library(scales)
  library(showtext); library(ggrepel)
  library(timeROC); library(glmnet)
})

ROOT <- "/Users/xiaolin/WorkBuddy/2026-09-04-20-00-08/paper2_triaptosis"
RES  <- file.path(ROOT, "results")
OUT  <- file.path(ROOT, "figures_main")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
source(file.path(ROOT, "code/lib/viz_base.R"))

# Okabe-Ito High/Low already in viz_base after palette update
cols_mod <- c(
  "VPS34 complex"     = "#0072B2",
  "PI(3)P"            = "#009E73",
  "Endosome-lysosome" = "#E69F00",
  "Endocytosis"       = "#CC79A7",
  "Oxidative stress"  = "#56B4E9")
gene_mod <- c(
  PIK3C3 = "VPS34 complex", PIK3R4 = "VPS34 complex", ATG14 = "VPS34 complex",
  NRBF2 = "VPS34 complex", BECN1 = "VPS34 complex", UVRAG = "VPS34 complex",
  MTM1 = "PI(3)P",
  KXD1 = "Endosome-lysosome", WDR91 = "Endosome-lysosome",
  FYCO1 = "Endosome-lysosome", RAB9A = "Endosome-lysosome",
  SCARB2 = "Endosome-lysosome", ATP13A2 = "Endosome-lysosome",
  CLCN4 = "Endosome-lysosome",
  SH3GL3 = "Endocytosis", ELMO2 = "Endocytosis", NCKAP1 = "Endocytosis",
  ACTN2 = "Endocytosis", ANXA8 = "Endocytosis",
  KEAP1 = "Oxidative stress", NFE2L2 = "Oxidative stress")
immune_class <- c(
  "Macrophage" = "Myeloid", "M2 macrophage" = "Myeloid",
  "M1 macrophage" = "Myeloid", "Monocyte" = "Myeloid",
  "MDSC" = "Myeloid", "Neutrophil" = "Myeloid", "Eosinophil" = "Myeloid",
  "Mast cell" = "Myeloid",
  "NK cell" = "NK / innate", "CD56bright NK cell" = "NK / innate",
  "CD56dim NK cell" = "NK / innate", "NKT cell" = "NK / innate",
  "Activated dendritic cell" = "NK / innate",
  "Immature dendritic cell" = "NK / innate",
  "Plasmacytoid dendritic cell" = "NK / innate",
  "B cell" = "B lineage", "Memory B cell" = "B lineage",
  "Plasma cell" = "B lineage",
  "Type 1 T helper cell" = "T helper", "Type 2 T helper cell" = "T helper",
  "Type 17 T helper cell" = "T helper",
  "T follicular helper cell" = "T helper", "Regulatory T cell" = "T helper",
  "Activated CD4 T cell" = "T helper",
  "Activated CD8 T cell" = "CD8 / effector",
  "Effector memory CD8 T cell" = "CD8 / effector",
  "Central memory CD8 T cell" = "CD8 / effector",
  "Gamma delta T cell" = "CD8 / effector")
class_ord <- c("Myeloid", "NK / innate", "B lineage", "T helper", "CD8")
immune_class[immune_class == "CD8 / effector"] <- "CD8"

finalize <- function(p, file, w = 180, h = 130, apply_theme = TRUE) {
  save_plot(p, file, dir = OUT, w = w, h = h, apply_theme = apply_theme)
}

## ---- data ----
score  <- read.csv(file.path(RES, "02_model/KIRC/score_table.csv"), row.names = 1)
screen <- read.csv(file.path(RES, "01_screen/pancancer_screen_matrix.csv"))
gcox   <- read.csv(file.path(RES, "01_screen/gene_cancer_cox.csv"))
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
immune <- read.csv(file.path(RES, "04_immune/immune_by_TAS_group.csv"))
ssg    <- read.csv(file.path(RES, "04_immune/ssgsea_per_sample.csv"), check.names = FALSE)
ici    <- read.csv(file.path(RES, "03_ici/immotion150_validation.csv"))
mod    <- readRDS(file.path(RES, "02_model/KIRC/model.rds"))
clin   <- read.csv(file.path(ROOT,
  "data/external_validation/kirc_tcga_pan_can_atlas_2018_clin_patient.csv"))

clin <- clin %>%
  mutate(stage4 = toupper(trimws(gsub("^STAGE\\s*", "",
                                      AJCC_PATHOLOGIC_TUMOR_STAGE))),
         STAGE = factor(ifelse(stage4 %in% c("I", "II", "III", "IV"),
                               stage4, NA),
                        levels = c("I", "II", "III", "IV")))
sd <- score %>% mutate(patientId = rownames(score)) %>%
  left_join(clin %>% select(patientId, AGE, SEX, STAGE), by = "patientId")

beta_df <- data.frame(gene = mod$genes, beta = as.numeric(mod$beta),
                      stringsAsFactors = FALSE) %>%
  mutate(module = unname(gene_mod[gene]),
         selected = gene %in% mod$genes_sel,
         module = factor(module, names(cols_mod)))

# =====================================================================
# FIG 1  Screening decision map
# =====================================================================
fig1 <- function() {
  suppressPackageStartupMessages(library(ggtext))
  col_sel  <- viz_palette$risk_high   # KIRC / selected
  col_soft <- viz_palette$accent      # soft pass
  col_hard <- "#C7C7C7"               # hard-pass only
  y_cut <- -log10(0.05)

  s <- screen %>%
    mutate(neglogp = -log10(pmax(logrank_p, 1e-30)),
           bal = pmin(n_c1, n_c2) / (n_c1 + n_c2),
           pass_n = n_total >= 150,
           pass_ev = n_dead >= 30,
           pass_cox = !is.na(n_sig_cox) & n_sig_cox >= 2,
           pass_clu = !is.na(logrank_p) & logrank_p < 0.05 &
             !is.na(bal) & bal > 0.15,
           soft = as.logical(pass_soft),
           n_ok = as.integer(pass_n) + as.integer(pass_ev) +
             as.integer(pass_cox) + as.integer(pass_clu),
           band = factor(
             case_when(
               soft ~ "Soft pass",
               as.logical(pass_hard) ~ "Hard only",
               TRUE ~ "Failed hard"),
             c("Soft pass", "Hard only", "Failed hard")))

  ## A  decision space: legends outside the panel, KIRC called out in free space
  sh <- s %>%
    filter(as.logical(pass_hard), !is.na(logrank_p)) %>%
    mutate(cls = factor(
      case_when(
        cancer == "KIRC" ~ "Selected (KIRC)",
        soft ~ "Soft pass",
        TRUE ~ "Hard-pass only"),
      c("Selected (KIRC)", "Soft pass", "Hard-pass only")))
  lab_df <- sh %>%
    filter(cancer %in% c("SARC", "LUSC", "PAAD")) %>%
    mutate(cancer = factor(cancer, c("SARC", "LUSC", "PAAD"))) %>%
    arrange(cancer) %>%
    mutate(lab = ifelse(cancer == "PAAD", "PAAD (cluster fail)",
                        as.character(cancer)),
           lcol = c(col_soft, col_soft, "#6E6E6E"),
           nx = c(-1.35, -1.15, 2.35),
           ny = c(0.72, -0.55, 0.62))

  pA <- ggplot(sh, aes(n_sig_cox, neglogp)) +
    annotate("rect", xmin = 2, xmax = 18.6, ymin = y_cut, ymax = 5.25,
             fill = col_soft, alpha = 0.055) +
    geom_vline(xintercept = 2, linetype = 2, colour = "#C6C6C6",
               linewidth = 0.3) +
    geom_hline(yintercept = y_cut, linetype = 2, colour = "#C6C6C6",
               linewidth = 0.3) +
    annotate("text", x = 2.35, y = 5.05, label = "Soft-pass region",
             size = 2.15, colour = "#00816C", hjust = 0, fontface = "italic") +
    annotate("text", x = 17.6, y = 1.1, label = "italic(P)==0.05",
             parse = TRUE, size = 2.1, colour = "#8C8C8C", hjust = 1) +
    annotate("text", x = 2.28, y = 2.85, label = "2 genes",
             size = 2.1, colour = "#8C8C8C", hjust = 0) +
    geom_point(data = dplyr::filter(sh, cancer == "KIRC"),
               aes(size = n_dead), shape = 21, fill = NA,
               colour = col_sel, stroke = 1.15, alpha = 0.38,
               show.legend = FALSE) +
    geom_point(aes(size = n_dead, fill = cls, alpha = cls), shape = 21,
               colour = "white", stroke = 0.45) +
    ggrepel::geom_text_repel(
      data = lab_df, aes(label = lab, colour = lcol),
      size = 2.45, fontface = "bold",
      nudge_x = lab_df$nx, nudge_y = lab_df$ny,
      min.segment.length = 0.05, box.padding = 0.28,
      point.padding = 0.2, segment.colour = "#A0A0A0",
      segment.size = 0.26, seed = 3, show.legend = FALSE) +
    scale_colour_identity() +
    ggtext::geom_richtext(
      data = data.frame(x = 2.4, y = 4.4), aes(x, y), inherit.aes = FALSE,
      label = paste0("<b>KIRC</b> &nbsp;17/21 Cox-significant<br>",
                     "<i>P</i><sub>k-means</sub> = 1.1&times;10<sup>&minus;5</sup>"),
      size = 2.35, hjust = 0, vjust = 0.5, lineheight = 1.25,
      fill = "#FFF6EF", colour = col_sel,
      label.colour = alpha(col_sel, 0.5),
      label.padding = unit(c(3.6, 4.6, 3.6, 4.6), "pt"),
      label.r = unit(1.4, "pt")) +
    annotate("curve", x = 11.1, xend = 15.9, y = 4.4, yend = 4.88,
             curvature = -0.24, linewidth = 0.34, colour = col_sel,
             arrow = arrow(length = unit(2, "pt"), type = "closed")) +
    scale_fill_manual(values = c(`Selected (KIRC)` = col_sel,
                                 `Soft pass` = col_soft,
                                 `Hard-pass only` = col_hard),
                      name = NULL) +
    scale_alpha_manual(values = c(`Selected (KIRC)` = 1, `Soft pass` = 0.95,
                                  `Hard-pass only` = 0.6),
                       guide = "none") +
    scale_size_continuous(range = c(1.9, 7.2), breaks = c(100, 200, 400),
                          name = "Deaths") +
    scale_x_continuous(limits = c(-0.8, 18.6), breaks = seq(0, 18, 3),
                       expand = c(0, 0)) +
    scale_y_continuous(limits = c(-0.18, 5.25), expand = c(0, 0)) +
    labs(x = "Cox-significant genes (of 21)",
         y = expression(-log[10](italic(P)[k-means])),
         title = "Screening decision space") +
    guides(fill = guide_legend(order = 1, ncol = 1,
                               override.aes = list(size = 2.6, alpha = 1)),
           size = guide_legend(order = 2, ncol = 1,
                               override.aes = list(fill = "grey60",
                                                   colour = "white",
                                                   alpha = 1))) +
    theme_publication() +
    theme(legend.position = "top",
          legend.justification = "left",
          legend.box = "horizontal",
          legend.box.just = "top",
          legend.margin = margin(0, 0, 1, 0),
          legend.spacing.x = unit(10, "pt"),
          legend.key.size = unit(7.5, "pt"),
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 7),
          panel.grid = element_blank())

  ## B  Cox-gene yield for hard-pass cancers with ≥4 genes; remainder as footnotes
  soft_v <- s %>% filter(band == "Soft pass") %>%
    arrange(desc(cancer == "KIRC"), desc(n_sig_cox)) %>% pull(cancer)
  hard_all <- s %>% filter(band == "Hard only") %>%
    arrange(desc(n_sig_cox), tidyr::replace_na(logrank_p, 1))
  hard_v <- hard_all %>% filter(n_sig_cox >= 4) %>% pull(cancer)
  hard_rest <- hard_all %>% filter(n_sig_cox < 4) %>% pull(cancer)
  fail_v <- s %>% filter(band == "Failed hard") %>%
    arrange(cancer) %>% pull(cancer)

  pack_y <- function(ids, y0) y0 - seq_along(ids) + 1
  gap <- 2.15
  y_soft <- pack_y(soft_v, 0)
  y_hard <- pack_y(hard_v, min(y_soft) - gap)
  ymap <- setNames(c(y_soft, y_hard), c(soft_v, hard_v))
  col_band <- c(`Soft pass` = col_soft, `Hard only` = "#8A93A6")
  sB <- s %>%
    filter(cancer %in% names(ymap)) %>%
    mutate(y = unname(ymap[cancer]),
           is_k = cancer == "KIRC",
           stem = ifelse(is_k, col_sel, unname(col_band[as.character(band)])),
           cox_n = n_sig_cox)
  hdr <- data.frame(
    lab  = c("Soft pass", "Hard-pass only  (\u2265 4 genes)"),
    y    = c(max(y_soft), max(y_hard)) + 1.05,
    hcol = c(col_soft, "#5B667A"))
  bands <- data.frame(
    ymin = c(min(y_soft), min(y_hard)) - 0.40,
    ymax = c(max(y_soft), max(y_hard)) + 0.40,
    bg = c("#E7F5F2", "#F4F5F7"))
  foot_y <- min(y_hard) - 1.85
  ylim_b <- c(foot_y - 2.35, max(hdr$y) + 0.25)
  fail_l1 <- paste(fail_v[1:6], collapse = "    ")
  fail_l2 <- paste(fail_v[7:12], collapse = "    ")

  pB <- ggplot(sB, aes(cox_n, y)) +
    geom_rect(data = bands, aes(ymin = ymin, ymax = ymax),
              xmin = -Inf, xmax = Inf, fill = bands$bg,
              colour = NA, inherit.aes = FALSE) +
    geom_vline(xintercept = 2, linetype = 2, colour = "#C8C8C8",
               linewidth = 0.3) +
    geom_segment(aes(x = 0, xend = cox_n, yend = y, colour = stem),
                 linewidth = 0.6) +
    geom_point(aes(colour = stem,
                   fill = ifelse(pass_clu, stem, NA_character_)),
               shape = 21, size = 2.35, stroke = 0.55) +
    geom_text(aes(x = -0.45, label = cancer, colour = stem,
                  fontface = ifelse(is_k, "bold", "plain")),
              hjust = 1, size = 2.7) +
    geom_text(data = hdr, aes(x = -0.45, y = y, label = lab, colour = hcol),
              hjust = 1, size = 2.55, fontface = "bold", inherit.aes = FALSE) +
    geom_text(data = dplyr::filter(sB, is_k),
              aes(x = cox_n + 0.45, label = cox_n),
              hjust = 0, size = 2.5, fontface = "bold", colour = col_sel) +
    geom_text(data = dplyr::filter(sB, cancer == "PAAD"),
              aes(x = cox_n + 0.40, label = "cluster n.s."),
              hjust = 0, size = 2.15, colour = "#6E6E6E") +
    annotate("text", x = -4.4, y = foot_y + 0.70,
             label = paste0(length(hard_rest),
                            " further hard-pass cancers with 2\u20133 genes, all cluster n.s."),
             hjust = 0, size = 2.15, colour = "#6E6E6E") +
    annotate("text", x = -4.4, y = foot_y + 0.05,
             label = paste(hard_rest, collapse = "    "),
             hjust = 0, size = 2.15, colour = "#8A8A8A") +
    annotate("text", x = -4.4, y = foot_y - 0.45,
             label = "Below hard gates (12)",
             hjust = 0, size = 2.55, fontface = "bold", colour = "#9A9A9A") +
    annotate("text", x = -4.4, y = foot_y - 1.15, label = fail_l1,
             hjust = 0, size = 2.15, colour = "#9A9A9A") +
    annotate("text", x = -4.4, y = foot_y - 1.85, label = fail_l2,
             hjust = 0, size = 2.15, colour = "#9A9A9A") +
    scale_colour_identity() +
    scale_fill_identity(na.value = "white") +
    scale_x_continuous(limits = c(-4.8, 20.2), breaks = c(2, 9, 17),
                       expand = c(0, 0)) +
    coord_cartesian(ylim = ylim_b, clip = "off") +
    labs(title = "Cox-gene yield among hard-pass cancers",
         x = "Cox-significant genes (of 21)", y = NULL,
         caption = "Filled = cluster P < 0.05;  open = cluster not significant.") +
    theme_publication() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          axis.text.x = element_text(size = 6.5),
          panel.grid = element_blank(),
          plot.caption = element_text(size = 6.2, colour = "grey35", hjust = 0,
                                      margin = margin(t = 4)),
          plot.margin = margin(4, 8, 4, 6),
          legend.position = "none")


  ## C 21 x 20 -log10(P) heatmap
  hard_ord <- sh$cancer[order(-sh$n_sig_cox, sh$logrank_p)]
  gc <- gcox %>%
    mutate(module = unname(gene_mod[gene]),
           module = factor(module, names(cols_mod)),
           gene = factor(gene, rev(names(gene_mod))),
           cancer = factor(cancer, hard_ord),
           nlp = pmin(-log10(pmax(p, 1e-12)), 8),
           star = ifelse(p < 0.05, "*", ""))
  pC <- ggplot(gc, aes(cancer, gene, fill = nlp)) +
    geom_tile(colour = "white", linewidth = 0.18) +
    geom_text(aes(label = star), size = 2.1, colour = "grey15", vjust = 0.75) +
    scale_fill_gradientn(
      colours = c("#FFFFFF", "#FCE3CD", "#F6A868", "#E06C1B", "#A63C05"),
      values = scales::rescale(c(0, 1.3, 3, 5, 8)), limits = c(0, 8),
      name = expression(-log[10](italic(P)))) +
    scale_x_discrete(labels = function(x)
      ifelse(x == "KIRC",
             paste0("<b style='color:", col_sel, ";'>KIRC</b>"), x)) +
    facet_grid(module ~ ., scales = "free_y", space = "free_y") +
    labs(x = NULL, y = NULL,
         title = "Per-gene Cox signal across hard-pass cancers (* P < 0.05)") +
    theme_publication() +
    theme(axis.text.x = ggtext::element_markdown(angle = 45, hjust = 1,
                                                 size = 5.5),
          axis.text.y = element_text(size = 5.5, face = "italic"),
          axis.ticks = element_line(linewidth = 0.25),
          strip.text.y = element_text(size = 6, angle = 0),
          legend.position = "top",
          legend.justification = "left",
          legend.margin = margin(0, 0, 1, 0),
          legend.key.height = unit(5.5, "pt"),
          legend.key.width = unit(15, "pt"),
          plot.title = element_text(margin = margin(b = 2)))

  top <- (pA | pB) + plot_layout(widths = c(1.15, 1))
  fin <- (top / pC) + plot_layout(heights = c(1.32, 1)) +
    plot_annotation(tag_levels = "A")
  finalize(fin, "Fig1_screen_decision", 190, 218, apply_theme = FALSE)
}

# =====================================================================
# FIG 2  Modular TAS
# =====================================================================
fig2 <- function() {
  ba <- beta_df %>%
    arrange(module, beta) %>%
    mutate(gene = factor(gene, gene))
  pA <- ggplot(ba, aes(beta, gene, colour = module)) +
    geom_vline(xintercept = 0, colour = "grey50", linewidth = 0.35) +
    geom_segment(aes(x = 0, xend = beta, yend = gene,
                     alpha = selected), linewidth = 0.55,
                 show.legend = FALSE) +
    geom_point(aes(shape = selected, fill = module), size = 2.1,
               colour = "grey20") +
    scale_colour_manual(values = cols_mod, name = NULL) +
    scale_fill_manual(values = cols_mod, guide = "none") +
    scale_shape_manual(values = c(`TRUE` = 21, `FALSE` = 1),
                       labels = c(`TRUE` = "In TAS (10)", `FALSE` = "Dropped"),
                       name = NULL) +
    scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.35)) +
    labs(x = "Elastic-net coefficient (frozen)", y = NULL,
         title = "21-gene modules; 10 retained") +
    theme_publication() +
    theme(axis.text.y = element_text(size = 6.2, face = "italic"))

  sc <- score %>% mutate(group = factor(group, c("Low", "High")))
  pB <- survfit2(Surv(time, status) ~ group, data = sc) |>
    ggsurvfit(linewidth = 0.75) +
    add_confidence_interval(alpha = 0.12) +
    add_risktable(risktable_stats = "n.risk",
                  stats_label = list(n.risk = "At risk"), size = 2.5) +
    annotate("text", x = 10, y = 0.12, size = 2.8, hjust = 0,
             label = "log-rank P = 2.9e-14") +
    scale_ggsurvfit() +
    scale_colour_manual(values = cols_risk,
                        labels = c("Low (n=254)", "High (n=254)"),
                        name = NULL) +
    scale_fill_manual(values = cols_risk, guide = "none") +
    labs(x = "Time (months)", y = "Overall survival",
         title = "TCGA-KIRC TAS (HR = 3.52, 2.76-4.50)") +
    theme_publication()
  pB <- ggsurvfit::ggsurvfit_build(pB)

  if (!requireNamespace("survivalROC", quietly = TRUE))
    stop("survivalROC required for published time-dependent AUC")
  roc_km <- lapply(c(12, 36, 60), function(tt) {
    survivalROC::survivalROC(Stime = score$time, status = score$status,
                             marker = score$TAS, predict.time = tt,
                             method = "KM")
  })
  roc_df <- bind_rows(lapply(seq_along(roc_km), function(i) {
    data.frame(FP = roc_km[[i]]$FP, TP = roc_km[[i]]$TP,
               T = c("12 mo", "36 mo", "60 mo")[i])
  })) %>%
    mutate(T = factor(T, c("12 mo", "36 mo", "60 mo")))
  auc_lab <- sapply(roc_km, `[[`, "AUC")
  pC <- ggplot(roc_df, aes(FP, TP, colour = T)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey70",
                linewidth = 0.35) +
    geom_line(linewidth = 0.7) +
    annotate("text", x = 0.52, y = c(0.28, 0.20, 0.12),
             label = sprintf("AUC(%s) = %.3f", levels(roc_df$T), auc_lab),
             size = 2.5, hjust = 0,
             colour = c(viz_palette$risk_high, viz_palette$accent,
                        viz_palette$risk_low)) +
    scale_colour_manual(values = c("12 mo" = viz_palette$risk_high,
                                   "36 mo" = viz_palette$accent,
                                   "60 mo" = viz_palette$risk_low),
                        name = NULL) +
    labs(x = "False positive rate", y = "True positive rate",
         title = "Time-dependent ROC") +
    theme_publication()

  (pA | (wrap_elements(full = pB) / pC + plot_layout(heights = c(1.55, 0.85)))) +
    plot_annotation(tag_levels = "A") +
    plot_layout(widths = c(0.9, 1.15))
  finalize(last_plot(), "Fig2_modular_TAS", 190, 185, apply_theme = FALSE)
}

# =====================================================================
# FIG 3  Prediction-grade utility card
# =====================================================================
fig3 <- function() {
  cmp <- bind_rows(
    ipcw_c %>% transmute(times, value = C, Metric = "IPCW C-index"),
    ipcw_a %>% transmute(times, value = AUC, Metric = "AUC-IPCW"))
  pA <- ggplot(cmp, aes(times, value, colour = Metric)) +
    geom_hline(yintercept = 0.5, linetype = 2, colour = "grey60",
               linewidth = 0.35) +
    geom_line(linewidth = 0.6) + geom_point(size = 2.2) +
    geom_text(aes(label = sprintf("%.3f", value)), vjust = -0.9, size = 2.3,
              show.legend = FALSE) +
    scale_colour_manual(values = c("IPCW C-index" = viz_palette$risk_high,
                                   "AUC-IPCW" = viz_palette$accent),
                        name = NULL) +
    scale_x_continuous(breaks = c(12, 36, 60)) +
    coord_cartesian(ylim = c(0.50, 0.85)) +
    labs(x = "Time (months)", y = "OOF value (5-fold)",
         title = "Discrimination") +
    theme_publication()

  br <- brier %>%
    mutate(times = factor(times, c(12, 36, 60),
                          c("12 mo", "36 mo", "60 mo")),
           model = factor(model, c("KM", "TAS")))
  pB <- ggplot(br, aes(model, Brier, group = times, colour = times)) +
    geom_line(linewidth = 0.55) +
    geom_point(size = 2.3) +
    geom_text(aes(label = sprintf("%.3f", Brier)), size = 2.2,
              vjust = -0.9, show.legend = FALSE) +
    annotate("text", x = 1.5, y = max(br$Brier) * 0.98,
             label = "IBS relative gain 9.5%", size = 2.4, colour = "grey30") +
    scale_colour_manual(values = c("12 mo" = viz_palette$risk_high,
                                   "36 mo" = viz_palette$accent,
                                   "60 mo" = viz_palette$risk_low),
                        name = NULL) +
    labs(x = NULL, y = "IPCW Brier score", title = "Calibration error") +
    theme_publication()

  cl <- calib %>%
    pivot_longer(c(pred_12, obs_12, pred_36, obs_36, pred_60, obs_60),
                 names_to = "k", values_to = "v") %>%
    separate(k, c("kind", "mo"), sep = "_") %>%
    pivot_wider(names_from = kind, values_from = v) %>%
    mutate(mo = factor(paste0(mo, " mo"), c("12 mo", "36 mo", "60 mo")))
  pC <- ggplot(cl, aes(pred, obs, colour = group)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey70",
                linewidth = 0.35) +
    geom_point(size = 2.1) +
    facet_wrap(~mo, nrow = 1) +
    scale_colour_manual(values = c(Q1 = viz_palette$risk_low, Q2 = "#56B4E9",
                                   Q3 = "#E69F00", Q4 = viz_palette$risk_high),
                        name = NULL) +
    coord_equal(xlim = c(0.30, 1), ylim = c(0.30, 1)) +
    labs(x = "Predicted OS", y = "Observed OS",
         title = "Quartile calibration") +
    theme_publication() +
    theme(legend.position = "bottom", legend.key.size = unit(3, "mm"))

  if (requireNamespace("dcurves", quietly = TRUE)) {
    dca_dat <- data.frame(time = oof$time, status = oof$status,
                          risk = 1 - oof$S36)
    dca_obj <- dcurves::dca(Surv(time, status) ~ risk, data = dca_dat,
                            time = 36, thresholds = seq(0.05, 0.95, 0.05))
    nb <- dca_obj$dca %>%
      mutate(label = recode(variable,
                            all = "Treat all", none = "Treat none",
                            risk = "TAS model")) %>%
      filter(variable %in% c("risk", "all", "none"))
    pD <- ggplot(nb, aes(threshold, net_benefit, colour = label,
                         linetype = label)) +
      geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.35) +
      geom_line(linewidth = 0.7) +
      scale_colour_manual(values = c("TAS model" = viz_palette$risk_high,
                                     "Treat all" = viz_palette$neutral,
                                     "Treat none" = "black"), name = NULL) +
      scale_linetype_manual(values = c("TAS model" = "solid",
                                       "Treat all" = "dashed",
                                       "Treat none" = "dotted"), name = NULL) +
      coord_cartesian(ylim = c(-0.04, 0.32)) +
      labs(x = "Threshold probability", y = "Net benefit (36-mo OS)",
           title = "Decision curve") +
      theme_publication()
  } else {
    pD <- ggplot() + theme_void() +
      labs(title = "DCA skipped (dcurves not installed)")
  }

  (pA | pB) / (pC | pD) + plot_annotation(tag_levels = "A")
  finalize(last_plot(), "Fig3_utility_card", 190, 165)
}

# =====================================================================
# FIG 4  Independence + transport
# =====================================================================
surv_at <- function(fit, newdata, t) {
  sf <- survfit(fit, newdata = newdata)
  s <- summary(sf, times = t, extend = TRUE)
  if (is.null(s$surv)) return(rep(NA_real_, nrow(newdata)))
  if (is.null(dim(s$surv))) return(as.numeric(s$surv))
  as.numeric(s$surv[1, ])
}

fig4 <- function() {
  dd <- sd %>% filter(!is.na(STAGE))
  grid <- do.call(rbind, lapply(levels(dd$STAGE), function(st) {
    d1 <- dd %>% filter(STAGE == st)
    if (nrow(d1) < 20 || sum(d1$status) < 5) return(NULL)
    fit <- coxph(Surv(time, status) ~ TAS, data = d1)
    xs <- seq(quantile(d1$TAS, 0.05), quantile(d1$TAS, 0.95), length.out = 40)
    nd <- data.frame(TAS = xs)
    data.frame(STAGE = paste("Stage", st), TAS = xs,
               S36 = surv_at(fit, nd, 36),
               S60 = surv_at(fit, nd, 60), n = nrow(d1))
  }))
  gl <- grid %>%
    pivot_longer(c(S36, S60), names_to = "horizon", values_to = "surv") %>%
    mutate(horizon = recode(horizon, S36 = "36 mo", S60 = "60 mo"))
  pA <- ggplot(gl, aes(TAS, surv, colour = horizon)) +
    geom_line(linewidth = 0.7) +
    facet_wrap(~STAGE, nrow = 1) +
    scale_colour_manual(values = c("36 mo" = viz_palette$accent,
                                   "60 mo" = viz_palette$risk_high),
                        name = NULL) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "TAS (continuous)", y = "Predicted OS",
         title = "Stage-conditional survival") +
    theme_publication()

  fa <- mvcox %>% filter(model == "Multivariable") %>%
    mutate(term_lab = recode(term,
                             TAS = "TAS (per SD)",
                             AGE = "Age (per year)",
                             SEXMALE = "Sex (male)",
                             STAGEII = "Stage II vs I",
                             STAGEIII = "Stage III vs I",
                             STAGEIV = "Stage IV vs I"),
           term_lab = factor(term_lab, rev(term_lab)),
           is_tas = term == "TAS")
  pB <- ggplot(fa, aes(HR, term_lab)) +
    geom_vline(xintercept = 1, linetype = 2, colour = "grey50",
               linewidth = 0.35) +
    geom_rect(data = fa %>% filter(is_tas),
              aes(xmin = 0.35, xmax = 12, ymin = as.numeric(term_lab) - 0.4,
                  ymax = as.numeric(term_lab) + 0.4),
              fill = "#F4F1EA", colour = NA, inherit.aes = FALSE) +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.18, linewidth = 0.5,
                   colour = viz_palette$risk_high) +
    geom_point(size = 2.1, shape = 15, colour = viz_palette$risk_high) +
    geom_text(aes(x = hi, label = sprintf("%.2f (%.2f-%.2f)", HR, lo, hi)),
              hjust = -0.08, size = 2.2) +
    scale_x_log10(limits = c(0.35, 20)) +
    labs(x = "Hazard ratio (log)", y = NULL,
         title = "Multivariable Cox (n=508)") +
    theme_publication()

  set.seed(42)
  hold <- sample(seq_len(nrow(score)), size = round(0.4 * nrow(score)))
  # Higher TAS / OOF risk = worse OS, so concordance uses the negative marker
  c_hold <- survival::concordance(Surv(time, status) ~ I(-TAS),
                                  data = score[hold, ])
  c_oof  <- survival::concordance(Surv(time, status) ~ I(-risk), data = oof)
  se <- function(obj) sqrt(as.numeric(obj$var)[1])
  pass <- bind_rows(
    data.frame(cohort = "TCGA-KIRC (OOF)",
               C = c_oof$concordance,
               lo = c_oof$concordance - 1.96 * se(c_oof),
               hi = c_oof$concordance + 1.96 * se(c_oof),
               note = "RNA-seq"),
    data.frame(cohort = "TCGA hold-out 40%",
               C = c_hold$concordance,
               lo = c_hold$concordance - 1.96 * se(c_hold),
               hi = c_hold$concordance + 1.96 * se(c_hold),
               note = "RNA-seq"),
    extsum %>% transmute(cohort, C = Cindex, lo = Cindex_lo, hi = Cindex_hi,
                         note = ifelse(cohort == "E-MTAB-1980",
                                       "Agilent", "TPM"))) %>%
    mutate(cohort = factor(cohort, rev(unique(cohort))))
  pC <- ggplot(pass, aes(C, cohort)) +
    geom_vline(xintercept = 0.5, linetype = 2, colour = "grey50",
               linewidth = 0.35) +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.15, linewidth = 0.5,
                   colour = viz_palette$secondary, na.rm = TRUE) +
    geom_point(size = 2.3, colour = viz_palette$secondary) +
    geom_text(aes(x = hi, label = sprintf("%.3f  %s", C, note)),
              hjust = -0.1, size = 2.3, na.rm = TRUE) +
    scale_x_continuous(limits = c(0.45, 1.15), expand = expansion(mult = c(0.02, 0.02))) +
    labs(x = "C-index (95% CI)", y = NULL, title = "Cohort passport") +
    theme_publication()

  km_small <- function(df, title, foot = NULL) {
    df <- df %>% mutate(risk = factor(risk, c("Low", "High")))
    p <- survfit2(Surv(os_t, os_s) ~ risk, data = df) |>
      ggsurvfit(linewidth = 0.65) +
      add_confidence_interval(alpha = 0.12) +
      scale_ggsurvfit() +
      scale_colour_manual(values = cols_risk, name = NULL) +
      scale_fill_manual(values = cols_risk, guide = "none") +
      labs(x = "OS (months)", y = "Survival", title = title,
           subtitle = foot) +
      theme_publication() +
      theme(plot.subtitle = element_text(size = 7, colour = "grey30"))
    p
  }
  pD1 <- km_small(emtab, sprintf("E-MTAB-1980 (n=%d, 23 events)", nrow(emtab)),
                  "KM P = 1.8e-03")
  pD2 <- km_small(cptac, sprintf("CPTAC-3 (n=%d, 20 events)", nrow(cptac)),
                  "KM P = 0.071; continuous Cox P = 4.8e-06")

  (pA) / (pB | pC) / (pD1 | pD2) +
    plot_annotation(tag_levels = "A") +
    plot_layout(heights = c(0.85, 1, 0.95))
  finalize(last_plot(), "Fig4_independence_transport", 190, 210)
}

# =====================================================================
# FIG 5  Immune phenotype
# =====================================================================
fig5 <- function() {
  col_hi <- viz_palette$risk_high
  col_lo <- viz_palette$risk_low
  im <- immune %>%
    mutate(
      class = unname(immune_class[as.character(cell)]),
      class = factor(class, class_ord),
      dir = factor(ifelse(spearman_rho >= 0, "Enriched in high TAS",
                          "Depleted in high TAS"),
                   c("Depleted in high TAS", "Enriched in high TAS")),
      sig = spearman_fdr < 0.05)
  cell_ord <- im %>%
    arrange(class, spearman_rho) %>%
    pull(cell)
  im <- im %>% mutate(cell = factor(cell, rev(cell_ord)))

  n_sig <- sum(im$sig)
  pA <- ggplot(im, aes(spearman_rho, cell)) +
    geom_vline(xintercept = 0, colour = "#C8C8C8", linewidth = 0.35) +
    geom_segment(aes(x = 0, xend = spearman_rho, yend = cell, colour = dir,
                     alpha = sig), linewidth = 0.55, show.legend = FALSE) +
    geom_point(aes(fill = dir, colour = dir, alpha = sig),
               shape = 21, size = 2.05, stroke = 0.35) +
    geom_text(data = dplyr::filter(im, sig),
              aes(x = spearman_rho + 0.022 * sign(spearman_rho), label = "*"),
              size = 2.6, colour = "grey20", vjust = 0.35) +
    facet_grid(class ~ ., scales = "free_y", space = "free_y") +
    scale_fill_manual(values = c(`Depleted in high TAS` = col_lo,
                                 `Enriched in high TAS` = col_hi),
                      name = NULL) +
    scale_colour_manual(values = c(`Depleted in high TAS` = col_lo,
                                   `Enriched in high TAS` = col_hi),
                        guide = "none") +
    scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.28), guide = "none") +
    scale_x_continuous(limits = c(-0.42, 0.40), breaks = c(-0.3, 0, 0.3),
                       expand = c(0, 0)) +
    labs(x = expression(Spearman~rho~(TAS~vs~ssGSEA)),
         y = NULL,
         title = sprintf("28 cell types vs TAS  (%d/28 FDR < 0.05)", n_sig)) +
    guides(fill = guide_legend(nrow = 2, override.aes = list(size = 2.4, alpha = 1))) +
    theme_publication() +
    theme(axis.text.y = element_text(size = 6.1),
          axis.text.x = element_text(size = 6.5),
          strip.text.y = element_text(size = 6.2, angle = 0),
          legend.position = "bottom",
          legend.justification = "left",
          legend.margin = margin(2, 0, 0, 0),
          panel.grid = element_blank(),
          plot.margin = margin(4, 6, 2, 2))

  keep_lo <- c("M2 macrophage", "Mast cell", "Type 17 T helper cell")
  keep_hi <- c("CD56bright NK cell", "Memory B cell", "T follicular helper cell")
  pretty <- c(
    "M2 macrophage" = "M2",
    "Mast cell" = "Mast",
    "Type 17 T helper cell" = "Th17",
    "CD56bright NK cell" = "CD56br NK",
    "Memory B cell" = "Memory B",
    "T follicular helper cell" = "Tfh")
  make_pole <- function(cells, title) {
    qdat <- im %>%
      filter(as.character(cell) %in% cells) %>%
      transmute(cell_raw = as.character(cell),
                qlab = paste0("q = ", fmt_p(wilcox_fdr)))
    dd <- ssg %>%
      select(patientId, group, all_of(cells)) %>%
      pivot_longer(-c(patientId, group), names_to = "cell_raw",
                   values_to = "score") %>%
      mutate(group = factor(group, c("Low", "High")),
             cell = factor(unname(pretty[cell_raw]), pretty[cells]))
    qdat <- qdat %>%
      mutate(cell = factor(unname(pretty[cell_raw]), levels(dd$cell))) %>%
      left_join(dd %>% group_by(cell) %>%
                  summarise(y = max(score) + 0.06 * diff(range(score)),
                            .groups = "drop"), by = "cell")
    ggplot(dd, aes(group, score, fill = group)) +
      geom_violin(colour = NA, alpha = 0.32, width = 0.9) +
      geom_boxplot(width = 0.24, outlier.size = 0.22, alpha = 0.92,
                   linewidth = 0.3) +
      geom_text(data = qdat, aes(x = 1.5, y = y, label = qlab),
                inherit.aes = FALSE, size = 2.15, colour = "grey25") +
      facet_wrap(~cell, nrow = 1, scales = "free_y") +
      scale_fill_manual(values = cols_risk,
                        labels = c("Low TAS", "High TAS"), name = NULL) +
      scale_y_continuous(expand = expansion(mult = c(0.04, 0.12))) +
      labs(x = NULL, y = "ssGSEA score", title = title) +
      theme_publication() +
      theme(strip.text = element_text(size = 7.2),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            legend.position = "bottom",
            legend.justification = "left",
            panel.grid = element_blank(),
            plot.margin = margin(2, 4, 1, 2))
  }
  pB <- (make_pole(keep_lo, "Depleted in high TAS") /
         make_pole(keep_hi, "Enriched in high TAS")) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom", legend.justification = "left")

  right <- wrap_elements(pB) +
    labs(title = "Six cells named in the text") +
    theme(plot.title = element_text(size = 9, face = "bold", hjust = 0,
                                    family = "Arial", margin = margin(b = 2)))

  (pA | right) + plot_annotation(tag_levels = "A") +
    plot_layout(widths = c(1.05, 1.12))
  finalize(last_plot(), "Fig5_immune_phenotype", 190, 158, apply_theme = FALSE)
}

# =====================================================================
# FIG 6  Boundary: what TAS is / is not
# =====================================================================
or_ci <- function(tab) {
  ft <- fisher.test(tab)
  data.frame(est = unname(ft$estimate), lo = ft$conf.int[1],
             hi = ft$conf.int[2], p = ft$p.value)
}

parallel_hr <- function(cancer) {
  expr <- fromJSON(file.path(ROOT, "data", paste0(cancer, "_expr.json")))
  clinj <- fromJSON(file.path(ROOT, "data", paste0(cancer, "_clin.json")))
  genes <- mod$genes
  mat <- t(vapply(expr, function(x) {
    v <- rep(NA_real_, length(genes)); names(v) <- genes
    for (g in genes) if (!is.null(x[[g]])) v[g] <- as.numeric(x[[g]])
    v
  }, FUN.VALUE = numeric(length(genes))))
  colnames(mat) <- genes
  rownames(mat) <- substr(names(expr), 1, 12)
  mat <- mat[!duplicated(rownames(mat)), , drop = FALSE]
  `%||%` <- function(a, b) if (is.null(a) || length(a) == 0 ||
                               (length(a) == 1 && is.na(a))) b else a
  os_status <- sapply(clinj, function(x) as.character(x$OS_STATUS %||% ""))
  os_months <- suppressWarnings(as.numeric(sapply(clinj, function(x) x$OS_MONTHS %||% NA)))
  dead <- grepl("DECEASED|1", toupper(os_status))
  keep <- !is.na(os_months) & os_months > 0 & os_status != ""
  common <- intersect(rownames(mat), names(clinj)[keep])
  idx <- match(common, names(clinj)[keep])
  z <- scale(mat[common, , drop = FALSE])
  tas <- as.numeric(z %*% mod$beta)
  d <- data.frame(time = os_months[keep][idx],
                  status = as.numeric(dead[keep][idx]),
                  group = factor(ifelse(tas > median(tas), "High", "Low"),
                                 c("Low", "High")))
  fx <- coxph(Surv(time, status) ~ group, data = d)
  lr <- survdiff(Surv(time, status) ~ group, data = d)
  ci <- exp(confint(fx))
  data.frame(cancer = cancer, n = nrow(d), events = sum(d$status),
             HR = exp(coef(fx))[1], lo = ci[1], hi = ci[2],
             km_p = 1 - pchisq(lr$chisq, 1))
}

fig6 <- function() {
  mu <- mean(null$cv_nested, na.rm = TRUE)
  sdv <- sd(null$cv_nested, na.rm = TRUE)
  zcut <- mu + 1.96 * sdv
  marks <- data.frame(
    name = c("TAS", "Apoptosis", "Ferroptosis", "Mitoxyperilysis"),
    C = c(0.6788, 0.681, 0.690, 0.659),
    pct = c(88.7, 94.0, 97.3, 27.3),
    y = c(14.2, 11.2, 8.2, 14.2))
  pA <- ggplot(null, aes(cv_nested)) +
    annotate("rect", xmin = zcut, xmax = Inf, ymin = -Inf, ymax = Inf,
             fill = "#F4F1EA", colour = NA) +
    geom_density(fill = viz_palette$neutral, colour = NA, alpha = 0.45) +
    geom_vline(xintercept = zcut, linetype = 3, colour = "grey40",
               linewidth = 0.4) +
    geom_vline(data = marks, aes(xintercept = C, colour = name),
               linewidth = 0.55) +
    geom_text(data = marks,
              aes(x = C, y = y, label = sprintf("%s\n%.1f%%", name, pct),
                  colour = name),
              vjust = 0, size = 2.2, show.legend = FALSE) +
    scale_colour_manual(values = c(TAS = viz_palette$risk_high,
                                   Apoptosis = "#0072B2",
                                   Ferroptosis = "#009E73",
                                   Mitoxyperilysis = "grey40"),
                        name = NULL) +
    annotate("text", x = zcut, y = 0, label = "z = 1.96",
             hjust = -0.1, vjust = -0.4, size = 2.2, colour = "grey35") +
    labs(x = "Nested CV C-index", y = "Density",
         title = "Null model: TAS at 88.7th percentile",
         subtitle = "Random 21-gene Cox/KM significance rate = 100%") +
    theme_publication()

  ici2 <- ici %>%
    mutate(group = factor(group, c("Low", "High")),
           resp = ifelse(responder == "R", "R",
                         ifelse(responder == "NR", "NR", NA)))
  tab_orr <- table(ici2$group, ici2$resp)
  orr <- or_ci(tab_orr)
  clin_ici <- fromJSON(file.path(ROOT, "data/IMmotion150_clin.json"))
  cb <- sapply(ici2$sample, function(s) {
    v <- clin_ici[[s]]$CLINICAL_BENEFIT
    if (is.null(v) || v == "") return(NA_character_)
    ifelse(toupper(v) %in% c("TRUE", "1", "YES"), "Yes", "No")
  })
  tab_cb <- table(ici2$group[!is.na(cb)], cb[!is.na(cb)])
  cbr <- or_ci(tab_cb)
  ici2$pfs_event <- as.numeric(grepl("1|Progressed", ici2$pfs_status))
  fx <- coxph(Surv(pfs, pfs_event) ~ group, data = ici2)
  ci <- exp(confint(fx))
  ici_forest <- data.frame(
    endpoint = factor(c("ORR (odds ratio)", "Clinical benefit (OR)",
                        "PFS (hazard ratio)"),
                      rev(c("ORR (odds ratio)", "Clinical benefit (OR)",
                            "PFS (hazard ratio)"))),
    est = c(orr$est, cbr$est, exp(coef(fx))[1]),
    lo = c(orr$lo, cbr$lo, ci[1]),
    hi = c(orr$hi, cbr$hi, ci[2]),
    p = c(orr$p, cbr$p, summary(fx)$coefficients[1, 5]))
  write.csv(ici_forest, file.path(RES, "03_ici/ici_null_forest.csv"),
            row.names = FALSE)
  pB <- ggplot(ici_forest, aes(est, endpoint)) +
    geom_vline(xintercept = 1, linetype = 2, colour = "grey50",
               linewidth = 0.35) +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.15, linewidth = 0.5,
                   colour = viz_palette$neutral) +
    geom_point(size = 2.3, colour = viz_palette$neutral) +
    geom_text(aes(label = sprintf("%.2f (%.2f-%.2f), P=%s",
                                  est, lo, hi, fmt_p(p))),
              hjust = 0, nudge_x = 0.08, size = 2.3) +
    scale_x_log10(limits = c(0.35, 6)) +
    labs(x = "Effect (High vs Low, log)", y = NULL,
         title = "TAS does not predict atezolizumab benefit") +
    theme_publication()

  kirc_hr <- data.frame(cancer = "KIRC", n = 508, events = sum(score$status),
                        HR = 3.56, lo = 2.51, hi = 5.06, km_p = 2.90e-14)  # high vs low, same estimand as SARC/LUSC
  par <- bind_rows(kirc_hr, parallel_hr("SARC"), parallel_hr("LUSC")) %>%
    mutate(cancer = factor(cancer, rev(c("KIRC", "SARC", "LUSC"))))
  write.csv(par, file.path(RES, "02_model/parallel_cancer_hr.csv"),
            row.names = FALSE)
  pC <- ggplot(par, aes(HR, cancer)) +
    geom_vline(xintercept = 1, linetype = 2, colour = "grey50",
               linewidth = 0.35) +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.15, linewidth = 0.5,
                   colour = ifelse(par$cancer == "KIRC",
                                   viz_palette$risk_high,
                                   viz_palette$neutral)) +
    geom_point(size = 2.3,
               colour = ifelse(par$cancer == "KIRC",
                               viz_palette$risk_high, viz_palette$neutral)) +
    geom_text(aes(label = sprintf("HR=%.2f, P=%s (n=%d)", HR, fmt_p(km_p), n)),
              hjust = 0, nudge_x = 0.12, size = 2.3) +
    scale_x_log10(limits = c(0.5, 12)) +
    labs(x = "HR High vs Low (frozen TAS)", y = NULL,
         title = "Parallel cancers: signal is KIRC-specific") +
    theme_publication()

  (pA) / (pB | pC) + plot_annotation(tag_levels = "A") +
    plot_layout(heights = c(1.05, 1))
  finalize(last_plot(), "Fig6_boundary", 190, 170)
}

# =====================================================================
# Graphical abstract
# =====================================================================
fig_ga <- function() {
  fun <- data.frame(
    step = factor(c("32 types", "Hard 20", "Soft 3", "KIRC"),
                  c("32 types", "Hard 20", "Soft 3", "KIRC")),
    n = c(32, 20, 3, 1),
    lab = c("32", "20", "3", "1  (17/21 Cox)"))
  pL <- ggplot(fun, aes(step, n, fill = step)) +
    geom_col(width = 0.62, show.legend = FALSE) +
    geom_text(aes(label = lab), vjust = -0.25, size = 2.5) +
    scale_fill_manual(values = c("32 types" = "grey75",
                                 "Hard 20" = viz_palette$neutral,
                                 "Soft 3" = viz_palette$accent,
                                 "KIRC" = viz_palette$risk_high)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.32))) +
    labs(x = NULL, y = "Cancer types", title = "Screen") +
    theme_publication()

  util <- data.frame(
    item = factor(c("OOF C-index   0.68",
                    "E-MTAB C-index   0.69",
                    "CPTAC C-index   0.73",
                    "Stage-adjusted HR   2.78"),
                  rev(c("OOF C-index   0.68",
                        "E-MTAB C-index   0.69",
                        "CPTAC C-index   0.73",
                        "Stage-adjusted HR   2.78"))))
  pM <- ggplot(util, aes(0.5, item)) +
    geom_tile(width = 0.92, height = 0.78, fill = "#E8F2EF") +
    geom_text(aes(label = item), size = 2.7) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0.4, 4.6), expand = FALSE) +
    labs(x = NULL, y = NULL, title = "Utility") +
    theme_void() +
    theme(plot.title = element_text(size = 9, face = "bold", hjust = 0))

  bnd <- data.frame(
    lab = c("ICI association\nNone (ORR P = 1.00)",
            "Gene-set specificity\n88.7th percentile, z = 0.96"),
    y = c(2, 1))
  pR <- ggplot(bnd, aes(0.06, y)) +
    geom_point(shape = 4, size = 3.4, stroke = 1.2,
               colour = viz_palette$risk_high) +
    geom_text(aes(label = lab), hjust = 0, nudge_x = 0.08,
              size = 2.5, lineheight = 0.95) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0.35, 2.65), expand = FALSE) +
    labs(x = NULL, y = NULL, title = "Boundary") +
    theme_void() +
    theme(plot.title = element_text(size = 9, face = "bold", hjust = 0))

  (pL | pM | pR) +
    plot_layout(widths = c(1, 1.05, 1.15)) +
    plot_annotation(
      title = "Triaptosis pan-cancer screen to KIRC TAS (utility, not mechanism)")
  finalize(last_plot(), "Fig0_graphical_abstract", 210, 78, apply_theme = FALSE)
}

# =====================================================================
# S1  LASSO path, risk triptych, k-means KM, nomogram
# =====================================================================
fig_s1 <- function() {
  expr_l <- fromJSON(file.path(ROOT, "data/KIRC_expr.json"))
  expr_mat <- t(vapply(expr_l, function(x) unlist(x[mod$genes]),
                       numeric(length(mod$genes))))
  colnames(expr_mat) <- mod$genes
  rownames(expr_mat) <- substr(rownames(expr_mat), 1, 12)
  keep <- rownames(expr_mat) %in% rownames(score)
  expr_mat <- expr_mat[keep, ]
  set.seed(42)
  fit <- glmnet(expr_mat, Surv(score[rownames(expr_mat), "time"],
                               score[rownames(expr_mat), "status"]),
                family = "cox", alpha = 0.5)
  bpath <- as.data.frame(as.matrix(fit$beta)) %>%
    mutate(gene = rownames(.)) %>%
    pivot_longer(-gene, names_to = "step", values_to = "coef") %>%
    mutate(step = as.numeric(sub("^s", "", step)),
           lambda = fit$lambda[step + 1],
           selected = gene %in% mod$genes_sel)
  lam_lab <- fit$lambda[which.min(abs(fit$lambda - mod$lambda))]
  lab_df <- dplyr::filter(bpath, selected, lambda == lam_lab) %>%
    arrange(coef)
  top_g <- lab_df$gene[which.max(abs(lab_df$coef))]
  rest <- lab_df$gene != top_g
  lab_df$y_lab <- lab_df$coef
  lab_df$y_lab[rest] <- seq(-0.006, 0.020, length.out = sum(rest))
  x_end <- log(lam_lab) + 0.55
  pL <- ggplot() +
    geom_line(data = dplyr::filter(bpath, !selected),
              aes(log(lambda), coef, group = gene),
              colour = "grey65", linewidth = 0.32) +
    geom_line(data = dplyr::filter(bpath, selected),
              aes(log(lambda), coef, group = gene),
              colour = viz_palette$risk_high, linewidth = 0.6) +
    geom_vline(xintercept = log(mod$lambda), linetype = 2,
               colour = "#8A8A8A", linewidth = 0.35) +
    annotate("text", x = log(mod$lambda) - 0.12,
             y = max(bpath$coef, na.rm = TRUE) * 0.96,
             label = "lambda[min]~(10~genes)", parse = TRUE,
             size = 2.25, hjust = 1, colour = "#555555") +
    geom_segment(data = lab_df,
                 aes(x = log(lambda), xend = x_end, y = coef, yend = y_lab),
                 colour = "#D4A890", linewidth = 0.25) +
    geom_text(data = lab_df,
              aes(x = x_end + 0.06, y = y_lab, label = gene),
              hjust = 0, size = 2.15, fontface = "italic",
              colour = viz_palette$risk_high) +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.28))) +
    labs(x = expression(log(lambda)), y = "Cox coefficient",
         title = "Elastic-net path (TAS construction)",
         subtitle = "Grey: 11 genes shrunk out. Orange: 10 genes retained.") +
    theme_publication() +
    theme(legend.position = "none",
          panel.grid = element_blank(),
          plot.subtitle = element_text(size = 6.5, colour = "grey35",
                                       face = "plain"))

  sc_km <- score %>% mutate(clu = factor(cluster, levels = c(1, 2)))
  nclu <- as.integer(table(sc_km$clu))
  lr <- survdiff(Surv(time, status) ~ clu, data = sc_km)
  p_lr <- 1 - pchisq(lr$chisq, length(lr$n) - 1)
  cox_cl <- coxph(Surv(time, status) ~ clu, data = sc_km)
  hr_cl <- summary(cox_cl)$conf.int[1, c(1, 3, 4)]
  pKM <- survfit2(Surv(time, status) ~ clu, data = sc_km) |>
    ggsurvfit(linewidth = 0.75) +
    add_confidence_interval(alpha = 0.08) +
    add_risktable(risktable_stats = "n.risk",
                  stats_label = list(n.risk = "At risk"), size = 2.4) +
    annotate("text", x = 8, y = 0.22, size = 2.45, hjust = 0, parse = TRUE,
             label = {
               e <- floor(log10(p_lr))
               m <- p_lr / 10^e
               paste0("log-rank~italic(P)==", sprintf("%.1f", m),
                      "%*%10^{", e, "}")
             }) +
    annotate("text", x = 8, y = 0.10, size = 2.45, hjust = 0,
             label = sprintf("HR = %.2f (%.2f-%.2f)",
                             hr_cl[1], hr_cl[2], hr_cl[3])) +
    scale_ggsurvfit() +
    scale_colour_manual(values = c("1" = viz_palette$risk_low,
                                   "2" = viz_palette$risk_high),
                        labels = c("1" = sprintf("C1 (n=%d)", nclu[1]),
                                   "2" = sprintf("C2 (n=%d)", nclu[2])),
                        name = NULL) +
    scale_fill_manual(values = c("1" = viz_palette$risk_low,
                                 "2" = viz_palette$risk_high), guide = "none") +
    labs(x = "Time (months)", y = "Overall survival",
         title = "k-means on 21 genes (screening gate)",
         subtitle = "Unsupervised clusters on the full 21-gene matrix, not TAS.") +
    theme_publication() +
    theme(plot.subtitle = element_text(size = 6.5, colour = "grey35",
                                       face = "plain"))
  pKM <- ggsurvfit::ggsurvfit_build(pKM)

  dt <- score %>% arrange(TAS) %>%
    mutate(id = row_number(), group = factor(group, c("Low", "High")))
  cutx <- sum(dt$group == "Low") + 0.5
  p1 <- ggplot(dt, aes(id, TAS, colour = group)) +
    geom_point(size = 0.55, show.legend = FALSE) +
    geom_vline(xintercept = cutx, linetype = 2, colour = "grey50",
               linewidth = 0.35) +
    scale_colour_manual(values = cols_risk) +
    labs(y = "TAS", x = NULL, title = "Risk score") +
    theme_publication() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  p2 <- ggplot(dt, aes(id, time, colour = factor(status, c(0, 1),
                                                 c("Alive", "Dead")))) +
    geom_point(size = 0.55) +
    geom_vline(xintercept = cutx, linetype = 2, colour = "grey50",
               linewidth = 0.35) +
    scale_colour_manual(values = c(Alive = viz_palette$risk_low,
                                   Dead = viz_palette$risk_high), name = NULL) +
    labs(y = "OS (months)", x = "Patients (increasing TAS)") +
    theme_publication()
  ez <- scale(expr_mat[rownames(dt), mod$genes_sel])
  ez[ez > 2] <- 2; ez[ez < -2] <- -2
  ez_df <- as.data.frame(ez) %>%
    mutate(sample_id = rownames(.)) %>%
    pivot_longer(-sample_id, names_to = "gene", values_to = "z") %>%
    mutate(id = match(sample_id, rownames(ez)))
  p3 <- ggplot(ez_df, aes(id, factor(gene, rev(mod$genes_sel)), fill = z)) +
    geom_tile() +
    geom_vline(xintercept = cutx, linetype = 2, colour = "grey50",
               linewidth = 0.35) +
    scale_fill_gradient2(low = viz_palette$risk_low, mid = "white",
                         high = viz_palette$risk_high, name = "Z") +
    labs(y = NULL, x = NULL, title = "Selected-gene expression") +
    theme_publication() +
    theme(axis.text.y = element_text(size = 6, face = "italic"),
          axis.text.x = element_blank(), axis.ticks.x = element_blank())
  pTrip <- p1 / p2 / p3 + plot_layout(heights = c(1, 1, 1.4))

  (pL | wrap_elements(full = pKM)) / wrap_elements(full = pTrip) +
    plot_annotation(tag_levels = "A") +
    plot_layout(heights = c(1.05, 1.15))
  finalize(last_plot(), "FigS1_model_extras", 190, 215, apply_theme = FALSE)

  if (requireNamespace("rms", quietly = TRUE)) {
    suppressPackageStartupMessages(library(rms))
    df1 <- sd %>% filter(!is.na(STAGE)) %>%
      mutate(SEX = factor(SEX, c("Female", "Male")),
             STAGE = factor(STAGE, c("I", "II", "III", "IV")),
             TAS_sd = as.numeric(scale(TAS)))
    units(df1$AGE) <- "years"
    dd <- datadist(df1)
    assign("dd", dd, envir = .GlobalEnv)
    options(datadist = "dd")
    fit <- rms::cph(Surv(time, status) ~ TAS_sd + AGE + SEX + STAGE,
                    data = df1, x = TRUE, y = TRUE, surv = TRUE, time.inc = 36)
    surv_fn <- Survival(fit)
    nom <- rms::nomogram(fit,
                         fun = list(function(x) surv_fn(36, lp = x),
                                    function(x) surv_fn(60, lp = x)),
                         funlabel = c("36-mo OS", "60-mo OS"),
                         lp = FALSE, fun.at = c(0.9, 0.7, 0.5, 0.3, 0.1))
    png(file.path(OUT, "FigS1D_nomogram.png"), width = 1800, height = 950,
        res = 300)
    plot(nom, xfrac = 0.30, cex.axis = 0.6, cex.var = 0.72)
    dev.off()
    pdf(file.path(OUT, "FigS1D_nomogram.pdf"), width = 7.1, height = 3.7)
    plot(nom, xfrac = 0.30, cex.axis = 0.6, cex.var = 0.72)
    dev.off()
    message("saved: FigS1D_nomogram")
  }
}

# =====================================================================
# S2  TAS + clinical + six immune cells (Spearman)
# =====================================================================
fig_s2 <- function() {
  keep <- c("M2 macrophage", "Mast cell", "Type 17 T helper cell",
            "CD56bright NK cell", "Memory B cell", "T follicular helper cell")
  d <- ssg %>%
    left_join(sd %>% select(patientId, TAS, AGE, STAGE), by = "patientId") %>%
    mutate(stage_n = as.numeric(STAGE)) %>%
    select(TAS, AGE, stage_n, all_of(keep))
  vars <- names(d)
  grid <- expand.grid(a = vars, b = vars, stringsAsFactors = FALSE)
  cr <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
    x <- d[[grid$a[i]]]; y <- d[[grid$b[i]]]
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 10) return(data.frame(a = grid$a[i], b = grid$b[i],
                                        rho = NA, p = NA))
    ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman"))
    data.frame(a = grid$a[i], b = grid$b[i],
               rho = unname(ct$estimate), p = ct$p.value)
  }))
  cr <- cr %>%
    mutate(a = factor(a, vars), b = factor(b, rev(vars)),
           star = ifelse(p < 0.05 & as.character(a) != as.character(b),
                         "*", ""))
  p <- ggplot(cr, aes(a, b, fill = rho)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    geom_text(aes(label = ifelse(is.na(rho), "",
                                 sprintf("%.2f%s", rho, star))),
              size = 2.3) +
    scale_fill_gradient2(low = viz_palette$risk_low, mid = "white",
                         high = viz_palette$risk_high, midpoint = 0,
                         name = expression(rho)) +
    labs(x = NULL, y = NULL,
         title = "TAS, clinical variables, and six immune cells") +
    theme_publication() +
    theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 6),
          axis.text.y = element_text(size = 6.5))
  finalize(p, "FigS2_tas_immune_clinical", 165, 135)
}

fig1(); fig2(); fig3(); fig4(); fig5(); fig6()
fig_ga(); fig_s1(); fig_s2()
message("ALL CLAIM FIGURES DONE -> ", OUT)
