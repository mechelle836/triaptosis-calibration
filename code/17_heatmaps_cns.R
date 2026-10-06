# =====================================================================
# 17_heatmaps_cns.R —— CNS 级热图重做（按 ComplexHeatmap 官方规范）
# 场景 A: Fig1B 32 癌种筛选矩阵 -> 多单列热图拼接 + cell_fun 数值 + 注释条
# 场景 B: Fig2D 10 基因 x 508 样本 -> column_split + 4 层 top_annotation
# 规范依据: /tmp/cns_heatmap_spec.md (ComplexHeatmap book 第 2-5 章,
#           Nat Commun 2025 doi:10.1038/s41467-025-64666-7 等)
# =====================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(jsonlite)
  library(ComplexHeatmap); library(circlize); library(RColorBrewer)
  library(showtext)
})
font_add("Arial", regular = "/System/Library/Fonts/Supplemental/Arial.ttf")
showtext_opts(dpi = 300); showtext_auto()

ROOT <- "/Users/xiaolin/WorkBuddy/2026-09-04-20-00-08/paper2_triaptosis"
RES  <- file.path(ROOT, "results")
OUT  <- file.path(ROOT, "figures_v2")
source(file.path(ROOT, "code", "lib", "viz_base.R"))

COL_OKBLUE <- "#2166AC"; COL_OKRED <- "#B2182B"; COL_OKMID <- "#F7F7F7"

## ---- 数据 ----
screen <- read.csv(file.path(RES, "01_screen/pancancer_screen_matrix.csv"))
score  <- read.csv(file.path(RES, "02_model/KIRC/score_table.csv"), row.names = 1)
ext_clin <- read.csv(file.path(
  ROOT, "data/external_validation/kirc_tcga_pan_can_atlas_2018_clin_patient.csv"))
mod <- readRDS(file.path(RES, "02_model/KIRC/model.rds"))
expr_l <- fromJSON(file.path(ROOT, "data/KIRC_expr.json"))
expr_mat <- t(vapply(expr_l, function(x) unlist(x[mod$genes]),
                     numeric(length(mod$genes))))
colnames(expr_mat) <- mod$genes
rownames(expr_mat) <- substr(rownames(expr_mat), 1, 12)
keep <- rownames(expr_mat) %in% rownames(score)
expr_mat <- expr_mat[keep, ]
ord_score <- score[rownames(expr_mat), ]

# =====================================================================
# 场景 A: 32 癌种筛选矩阵 (CNS 规范: 每指标独立单色色带 + cell_fun
#         数值 + pass_soft 注释 + KIRC 高亮, 不聚类语义排序)
# =====================================================================
fig1B_cns <- function() {
  ## 只保留硬筛通过的癌种 (NA 灰块全部消除), 语义排序
  s <- screen %>%
    filter(pass_hard) %>%
    arrange(desc(n_sig_cox), logrank_p)
  stopifnot(nrow(s) == sum(screen$pass_hard))

  ## rank 百分位归一化 (解决线性归一化弱值全白的区分度问题)
  pct <- function(x) rank(x, na.last = "keep") / length(x)
  s <- s %>%
    mutate(
      `Total (n)`    = pct(n_total),
      `Deaths (n)`   = pct(n_dead),
      `Sig. Cox /21` = pct(n_sig_cox),
      `-log10(P)`    = pct(-log10(pmax(logrank_p, 1e-30))))
  m <- as.matrix(s[, c("Total (n)", "Deaths (n)", "Sig. Cox /21", "-log10(P)")])
  rownames(m) <- s$cancer

  ## 每指标独立单色色带 (rank 百分位 0-1)
  cfs <- list(
    colorRamp2(c(0, 0.5, 1), c("#F7FBFF", "#6BAED6", "#08519C")),
    colorRamp2(c(0, 0.5, 1), c("#F7F7F7", "#9E9AC8", "#3F3B7A")),
    colorRamp2(c(0, 0.5, 1), c("#FFF5F0", "#FB6A4A", "#A50F15")),
    colorRamp2(c(0, 0.5, 1), c("#F7FCF5", "#41AB5D", "#00441B")))

  ## cell_fun 闭包: 每格标注原始数值 (n = 整数, P = 2 位小数)
  raw <- list(
    function(i) sprintf("%d", s$n_total[i]),
    function(i) sprintf("%d", s$n_dead[i]),
    function(i) sprintf("%d", s$n_sig_cox[i]),
    function(i) sprintf("%.2f", -log10(pmax(s$logrank_p[i], 1e-30))))
  make_cf <- function(jcol) {
    force(jcol)
    rf <- raw[[jcol]]
    function(j, i, x, y, w, h, fill) {
      v <- m[i, jcol]
      if (is.na(v)) return(invisible(NULL))   # NA 格跳过 (显示 na_col 灰)
      grid::grid.text(rf(i), x, y,
                      gp = grid::gpar(fontsize = 7,
                                      col = if (v > 0.55) "white" else "grey15",
                                      fontfamily = "Arial"))
    }
  }

  ## Pass 注释条 (Yes 深绿 / No 浅灰)
  pass_lab <- ifelse(s$pass_soft, "Yes", "No")
  ha_row <- rowAnnotation(
    `Pass` = anno_simple(pass_lab, pt_size = unit(2, "mm"),
                         col = c("Yes" = "#00441B", "No" = "grey92")),
    annotation_name_side = "top",
    show_annotation_name = FALSE,
    show_legend = FALSE)

  ## 行名: KIRC 加粗红色
  name_col <- ifelse(rownames(m) == "KIRC", viz_palette$risk_high, "black")
  name_fw  <- ifelse(rownames(m) == "KIRC", "bold", "plain")

  ht_list <- NULL
  for (j in seq_len(ncol(m))) {
    ht <- Heatmap(
      m[, j, drop = FALSE], name = colnames(m)[j],
      col = cfs[[j]],
      cluster_rows = FALSE, cluster_columns = FALSE,
      cell_fun = make_cf(j),
      width = unit(1.9, "cm"),
      row_names_side = "left",
      row_names_gp = grid::gpar(fontsize = 8.5, col = name_col,
                                fontface = name_fw, fontfamily = "Arial"),
      column_names_gp = grid::gpar(fontsize = 9, fontfamily = "Arial"),
      column_names_rot = 0,
      column_title = colnames(m)[j],
      column_title_gp = grid::gpar(fontsize = 9.5, fontface = "bold",
                                   fontfamily = "Arial"),
      show_heatmap_legend = FALSE,
      na_col = "grey92",
      border = TRUE,
      right_annotation = if (j == ncol(m)) ha_row else NULL)
    ht_list <- if (is.null(ht_list)) ht else ht_list + ht
  }

  png(file.path(OUT, "Fig1B_screen_matrix_cns.png"),
      width = 2000, height = 1750, res = 300)
  draw(ht_list, heatmap_legend_side = "bottom",
       annotation_legend_side = "bottom",
       column_title = "Pan-cancer screen matrix (20 cancers passing hard filter)",
       column_title_gp = grid::gpar(fontsize = 12, fontface = "bold",
                                    fontfamily = "Arial"),
       padding = unit(c(8, 4, 4, 4), "mm"))
  dev.off()
  pdf(file.path(OUT, "Fig1B_screen_matrix_cns.pdf"),
      width = 6.7, height = 5.8, useDingbats = FALSE)
  draw(ht_list, heatmap_legend_side = "bottom",
       annotation_legend_side = "bottom",
       column_title = "Pan-cancer screen matrix (20 cancers passing hard filter)",
       column_title_gp = grid::gpar(fontsize = 12, fontface = "bold",
                                    fontfamily = "Arial"),
       padding = unit(c(8, 4, 4, 4), "mm"))
  dev.off()
  message("saved: Fig1B_screen_matrix_cns")
}

fig2D_cns <- function() {
  gene_sel <- mod$genes_sel
  ez <- t(scale(t(expr_mat[, gene_sel])))
  ez[ez > 2] <- 2; ez[ez < -2] <- -2
  ez <- ez[, gene_sel]          # 列序 = 基因顺序

  stg <- ext_clin %>% select(patientId, AJCC_PATHOLOGIC_TUMOR_STAGE) %>%
    mutate(stage = sub("STAGE ", "", AJCC_PATHOLOGIC_TUMOR_STAGE))
  meta <- ord_score %>% mutate(patientId = rownames(ord_score)) %>%
    left_join(stg %>% select(patientId, stage), by = "patientId") %>%
    mutate(group = as.character(group),
           status_f = ifelse(status == 1, "Dead", "Alive"))
  # 列排序: Low 块 TAS 升序 -> High 块 TAS 升序
  col_ord <- order(factor(meta$group, levels = c("Low", "High")),
                   meta$TAS)   # Low 块在前 (factor 级别序, 非字母序)
  ez <- ez[col_ord, ]
  meta <- meta[col_ord, ]

  col_fun <- colorRamp2(c(-2, 0, 2), c(COL_OKBLUE, COL_OKMID, COL_OKRED))
  cols4 <- c("I" = "#F7F7F7", "II" = "#BDD7E7",
             "III" = "#6BAED6", "IV" = "#2171B5")
  stg_vec <- ifelse(is.na(meta$stage), "Stage I", meta$stage)  # 3 例 NA 归 I

  ha_top <- HeatmapAnnotation(
    Risk = anno_simple(meta$group,
                       col = c("Low" = viz_palette$risk_low,
                               "High" = viz_palette$risk_high),
                       pt_size = unit(1.5, "mm"), height = unit(3, "mm")),
    `TAS` = anno_lines(meta$TAS, border = FALSE, gp = grid::gpar(lwd = 0.6),
                       height = unit(1.1, "cm"),
                       add_points = FALSE),
    Status = anno_simple(meta$status_f,
                         col = c("Alive" = viz_palette$neutral,
                                 "Dead" = viz_palette$risk_high),
                         height = unit(3, "mm")),
    Stage = anno_simple(stg_vec, col = cols4, height = unit(3, "mm")),
    annotation_name_side = "left",
    annotation_name_gp = grid::gpar(fontsize = 7.5, fontfamily = "Arial"),
    show_legend = FALSE)

  col_title <- c(sprintf("Low risk (n=%d)", sum(meta$group == "Low")),
                 sprintf("High risk (n=%d)", sum(meta$group == "High")))

  ht <- Heatmap(
    t(ez), name = "Z-score",   # 10 基因行 x 508 样本列
    col = col_fun,
    top_annotation = ha_top,
    column_split = factor(meta$group, levels = c("Low", "High")),
    # 必须 factor: 字符向量会按字母序分块 (High 在前) 导致标题错位
    column_title = col_title,
    column_title_gp = grid::gpar(fontsize = 9.5, fontface = "bold",
                                 fontfamily = "Arial"),
    column_gap = unit(2.5, "mm"),
    cluster_rows = TRUE, clustering_method_rows = "ward.D2",
    cluster_columns = FALSE,
    row_dend_width = unit(7, "mm"),
    row_names_gp = grid::gpar(fontsize = 8, fontface = "italic",
                              fontfamily = "Arial"),
    column_names_gp = grid::gpar(fontsize = 4.5),
    show_column_names = FALSE,
    border = TRUE,
    heatmap_legend_param = list(
      title = "Z-score", legend_height = unit(2.8, "cm"),
      title_gp = gpar(fontsize = 8, fontface = "bold"),
      labels_gp = gpar(fontsize = 6.5), title_position = "topleft"))

  png(file.path(OUT, "Fig2D_gene_heatmap_cns.png"),
      width = 2400, height = 1250, res = 300)
  draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom",
       column_title = "Triaptosis signature gene expression (10-gene TAS)",
       column_title_gp = grid::gpar(fontsize = 11, fontface = "bold",
                                    fontfamily = "Arial"),
       merge_legends = TRUE,
       padding = unit(c(6, 2, 2, 2), "mm"))
  dev.off()
  pdf(file.path(OUT, "Fig2D_gene_heatmap_cns.pdf"),
      width = 8, height = 4.2, useDingbats = FALSE)
  draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom",
       column_title = "Triaptosis signature gene expression (10-gene TAS)",
       column_title_gp = grid::gpar(fontsize = 11, fontface = "bold",
                                    fontfamily = "Arial"),
       merge_legends = TRUE,
       padding = unit(c(6, 2, 2, 2), "mm"))
  dev.off()
  message("saved: Fig2D_gene_heatmap_cns")
}

fig1B_cns(); fig2D_cns()
message("CNS HEATMAPS DONE")