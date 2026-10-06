# =====================================================================
# 30_spatial_tas.R
# GSE210041 两张 ccRCC FFPE Visium：同一套 10 基因 TAS 是否落在肿瘤区域。
# 标签来自 GSE159115 的 Seurat 锚点转移；CA9 作为不在 TAS 里的独立对照。
# 不做空间 CellChat，不做 RCTD，不做生存。
# 输出: results/16_spatial/spatial_tas_summary.csv
#       figures_pub/Fig8_spatial.{png,pdf}
# =====================================================================
# Visium objects on Seurat 5.5 fail AddModuleScore if created as Assay5.
options(Seurat.object.assay.version = "v3")
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})
ROOT <- "/Users/xiaolin/WorkBuddy/2026-09-04-20-00-08/paper2_triaptosis"
RAW  <- file.path(ROOT, "data/spatial/GSE210041")
RES  <- file.path(ROOT, "results/16_spatial")
OUT  <- file.path(ROOT, "figures_pub")
dir.create(RES, showWarnings = FALSE, recursive = TRUE)
source(file.path(ROOT, "code/lib/viz_base.R"))
showtext_auto(FALSE)

samples <- list(
  list(id = "GSM6415705", label = "Tumour 1",
       h5 = "GSM6415705_GLMF1_filtered_feature_bc_matrix.h5",
       tar = "GSM6415705_GLMF1_spatial.tar.gz"),
  list(id = "GSM6415706", label = "Tumour 7",
       h5 = "GSM6415706_GLMF2_filtered_feature_bc_matrix.h5",
       tar = "GSM6415706_GLMF2_spatial.tar.gz")
)

prepare_dir <- function(sm) {
  dest <- file.path(RAW, sm$id)
  dir.create(dest, showWarnings = FALSE, recursive = TRUE)
  h5_src <- file.path(RAW, sm$h5)
  if (!file.exists(h5_src)) stop("missing h5: ", h5_src)
  file.copy(h5_src, file.path(dest, "filtered_feature_bc_matrix.h5"), overwrite = TRUE)
  spat <- file.path(dest, "spatial")
  if (!file.exists(file.path(spat, "tissue_lowres_image.png")) &&
      !file.exists(file.path(spat, "tissue_hires_image.png"))) {
    tar <- file.path(RAW, sm$tar)
    if (!file.exists(tar)) stop("missing spatial tar: ", tar)
    tmp <- file.path(dest, "tar_tmp")
    unlink(tmp, recursive = TRUE)
    dir.create(tmp)
    untar(tar, exdir = tmp)
    hits <- list.files(tmp, pattern = "tissue_lowres_image.png|tissue_hires_image.png",
                       recursive = TRUE, full.names = TRUE)
    if (length(hits) == 0) stop("spatial images not found in ", tar)
    src_dir <- dirname(hits[1])
    unlink(spat, recursive = TRUE)
    dir.create(spat, showWarnings = FALSE)
    file.copy(list.files(src_dir, full.names = TRUE), spat, recursive = TRUE)
    unlink(tmp, recursive = TRUE)
  }
  dest
}

mod <- readRDS(file.path(ROOT, "results/02_model/KIRC/model.rds"))
tas_genes <- mod$genes_sel
ref <- readRDS(file.path(ROOT, "results/11_scrna/seurat_kirc_gse159115.rds"))
DefaultAssay(ref) <- "RNA"

transfer_labels <- function(vis, ref) {
  feats <- intersect(VariableFeatures(ref), rownames(vis))
  if (length(feats) < 300) feats <- intersect(rownames(ref), rownames(vis))
  feats <- head(feats, 2000)
  anchors <- FindTransferAnchors(
    reference = ref, query = vis, features = feats,
    normalization.method = "LogNormalize",
    reference.reduction = "pca", dims = 1:30, verbose = FALSE
  )
  pred <- TransferData(anchorset = anchors, refdata = ref$ctype,
                       dims = 1:30, verbose = FALSE)
  vis$predicted.id <- pred$predicted.id
  vis$prediction.score.max <- pred$prediction.score.max
  vis
}

score_one <- function(sm, ref) {
  dest <- prepare_dir(sm)
  message("loading ", sm$label)
  vis <- Load10X_Spatial(data.dir = dest, filename = "filtered_feature_bc_matrix.h5",
                         assay = "Spatial")
  vis <- NormalizeData(vis, assay = "Spatial", verbose = FALSE)
  present <- intersect(tas_genes, rownames(vis))
  if (length(present) < 5) stop("too few TAS genes in ", sm$id, ": ", paste(present, collapse = ","))
  vis <- AddModuleScore(vis, features = list(present), name = "TAS_score", assay = "Spatial")
  vis$section <- sm$label
  vis$gsm <- sm$id
  if ("CA9" %in% rownames(vis)) {
    vis$CA9 <- as.numeric(GetAssayData(vis, assay = "Spatial", layer = "data")["CA9", ])
  } else {
    vis$CA9 <- NA_real_
  }
  vis <- tryCatch(transfer_labels(vis, ref), error = function(e) {
    message("label transfer failed for ", sm$id, ": ", conditionMessage(e))
    vis$predicted.id <- NA_character_
    vis$prediction.score.max <- NA_real_
    vis
  })
  vis$is_tumor <- vis$predicted.id == "Tumor cells"
  list(vis = vis, n_genes = length(present), genes = present)
}

scored <- lapply(samples, score_one, ref = ref)
rm(ref)

summarise_one <- function(obj) {
  vis <- obj$vis
  md <- vis@meta.data
  tumor <- md$is_tumor %in% TRUE
  rest <- md$is_tumor %in% FALSE
  wt_p <- NA_real_
  if (sum(tumor) >= 5 && sum(rest) >= 5) {
    wt_p <- wilcox.test(md$TAS_score1[tumor], md$TAS_score1[rest], exact = FALSE)$p.value
  }
  sp_r <- NA_real_
  sp_p <- NA_real_
  if (sum(is.finite(md$CA9)) > 20) {
    ct <- cor.test(md$TAS_score1, md$CA9, method = "spearman", exact = FALSE)
    sp_r <- unname(ct$estimate)
    sp_p <- ct$p.value
  }
  data.frame(
    gsm = md$gsm[1],
    section = md$section[1],
    n_spots = nrow(md),
    n_tas_genes = obj$n_genes,
    tas_genes = paste(obj$genes, collapse = ";"),
    n_tumor_spots = sum(tumor),
    n_rest_spots = sum(rest),
    median_pred_score = median(md$prediction.score.max, na.rm = TRUE),
    median_pred_score_tumor = if (any(tumor)) median(md$prediction.score.max[tumor], na.rm = TRUE) else NA_real_,
    median_tas_tumor = if (any(tumor)) median(md$TAS_score1[tumor]) else NA_real_,
    median_tas_rest = if (any(rest)) median(md$TAS_score1[rest]) else NA_real_,
    wilcox_p = wt_p,
    spearman_ca9 = sp_r,
    spearman_p = sp_p,
    stringsAsFactors = FALSE
  )
}

summary_df <- bind_rows(lapply(scored, summarise_one))
write.csv(summary_df, file.path(RES, "spatial_tas_summary.csv"), row.names = FALSE)
message("spatial summary:")
print(summary_df)

spot_df <- bind_rows(lapply(scored, function(obj) {
  md <- obj$vis@meta.data
  data.frame(section = md$section, tas = md$TAS_score1, ca9 = md$CA9,
             predicted = md$predicted.id,
             compartment = ifelse(md$is_tumor %in% TRUE, "Tumor spots", "Other spots"),
             stringsAsFactors = FALSE)
}))
write.csv(spot_df, file.path(RES, "spatial_spot_scores.csv"), row.names = FALSE)

# ---- figure ----------------------------------------------------------
# 组织学图必须给比例尺。不去解析 scalefactors_json, 而是用 Visium 自身的几何
# 自校准: 相邻 spot 中心距固定为 100 um, 所以绘图坐标里量出最近邻距离
# 就得到 um/绘图单位, 无论 Seurat 对坐标做了什么缩放都成立。
plot_units_per_um <- function(vis) {
  co <- GetTissueCoordinates(vis)
  xy <- as.matrix(co[, intersect(c("imagecol", "imagerow", "x", "y"), colnames(co))])
  if (ncol(xy) < 2) return(NA_real_)
  idx <- seq_len(min(nrow(xy), 600))
  d <- as.matrix(dist(xy[idx, , drop = FALSE]))
  diag(d) <- Inf
  nn <- apply(d, 1, min)
  median(nn, na.rm = TRUE) / 100      # 100 um center-to-center
}

# SpatialFeaturePlot 的注释坐标系和它实际画出来的组织学图对不上 (试过 spot 层
# 范围和 panel_params, annotate 都落到组织外面), 而比例尺必须和组织同一坐标系。
# 所以直接自己画: 低分辨率组织图按其像素尺寸铺成 raster, spot 坐标乘 lowres
# scale factor 落在同一像素网格上, coord_fixed(expand = FALSE) 锁死 1:1。
lowres_scalef <- function(img_obj) {
  sf <- tryCatch(Seurat::ScaleFactors(img_obj)$lowres, error = function(e) NULL)
  if (is.null(sf)) sf <- img_obj@scale.factors$lowres
  sf
}

map_one <- function(obj, bar_um = 1000) {
  vis <- obj$vis
  img_obj <- vis@images[[1]]
  sf <- lowres_scalef(img_obj)
  img <- img_obj@image                      # H x W x 3, 0-1
  iw <- ncol(img); ih <- nrow(img)
  img <- img * 0.62 + 0.38                  # 压淡组织学, 让 spot 颜色读得出来

  co <- GetTissueCoordinates(vis)
  xcol <- intersect(c("imagecol", "x"), colnames(co))[1]
  ycol <- intersect(c("imagerow", "y"), colnames(co))[1]
  df <- data.frame(x = co[[xcol]] * sf,
                   y = ih - co[[ycol]] * sf,   # 图像行号自上而下, 翻成向上为正
                   TAS = vis$TAS_score1[match(rownames(co), colnames(vis))])

  upu <- plot_units_per_um(vis) * sf         # 低分辨率像素 / um
  bar_len <- bar_um * upu
  # 裁掉定位框留白, 只留 spot 范围外 3%
  pad <- 0.03 * max(diff(range(df$x)), diff(range(df$y)))
  xl <- c(min(df$x) - pad, max(df$x) + pad)
  yl <- c(min(df$y) - pad, max(df$y) + pad)
  # 点大小按 spot 间距 (100 um) 给, 相邻 spot 几乎相接; ggplot 的 size 约等于
  # 毫米直径, 所以要把像素换成成图上的毫米 (每个 panel 约 86 mm 宽)。
  # 成图 180 mm 宽, 上排两张切片受高度限制, 每张绘图区实测约 68 mm
  mm_per_px <- 68 / diff(xl)
  pt_size <- 100 * upu * mm_per_px * 0.85
  x1 <- xl[2] - 0.03 * diff(xl); x0 <- x1 - bar_len
  yb <- yl[1] + 0.045 * diff(yl)
  message(sprintf("%s: image %dx%d px, sf %.5f, %.4f px/um, bar %.1f px, pt %.2f mm, xlim[%.0f,%.0f] ylim[%.0f,%.0f]",
                  unique(vis$section), iw, ih, sf, upu, bar_len, pt_size,
                  xl[1], xl[2], yl[1], yl[2]))

  lim <- quantile(df$TAS, c(0.05, 0.95), na.rm = TRUE)
  ggplot(df, aes(x, y, fill = TAS)) +
    annotation_raster(img, xmin = 0, xmax = iw, ymin = 0, ymax = ih,
                      interpolate = TRUE) +
    geom_point(shape = 21, stroke = 0, size = pt_size, alpha = 0.85) +
    scale_fill_gradient2(low = "#3C5488", mid = "#F7F7F7", high = "#E64B35",
                         midpoint = 0,
                         limits = c(min(lim[1], -0.05), max(lim[2], 0.05)),
                         oob = scales::squish, name = "TAS",
                         guide = guide_cbar(barwidth = unit(16, "mm"),
                                            barheight = unit(2, "mm"),
                                            title.position = "left")) +
    annotate("segment", x = x0, xend = x1, y = yb, yend = yb,
             linewidth = 1.0, colour = "black", lineend = "butt") +
    annotate("text", x = (x0 + x1) / 2, y = yb, vjust = -0.6,
             label = sprintf("%g mm", bar_um / 1000),
             size = 1.9, family = "Arial", colour = "black") +
    coord_fixed(xlim = xl, ylim = yl, expand = FALSE) +
    labs(title = sub("Tumour", "Tumor", unique(vis$section))) +
    theme_void(base_size = 7.5, base_family = "Arial") +
    theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5,
                                    margin = margin(b = 1, unit = "mm")),
          legend.position = "bottom",
          legend.margin = margin(t = 0.5, unit = "mm"),
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 7, vjust = 1),
          plot.margin = margin(1, 1, 1, 1, unit = "mm"))
}
pA <- map_one(scored[[1]])
pB <- map_one(scored[[2]])

comp <- spot_df %>% filter(!is.na(predicted)) %>%
  mutate(section = sub("Tumour", "Tumor", section))
# 小提琴图原先既没有检验结果也没有每组 spot 数, 补上 Mann-Whitney P 与 n
n_lab <- comp %>% count(section, compartment) %>%
  mutate(lab = sprintf("n=%d", n))
p_lab <- summary_df %>%
  transmute(section = sub("Tumour", "Tumor", section),
            lab = ifelse(wilcox_p < 0.001, "Mann-Whitney P < 0.001",
                         sprintf("Mann-Whitney P = %.3f", wilcox_p)))
ytop <- max(comp$tas, na.rm = TRUE)
ybot <- min(comp$tas, na.rm = TRUE)
pC <- ggplot(comp, aes(compartment, tas, fill = compartment)) +
  geom_violin(scale = "width", linewidth = 0.25, colour = "grey30") +
  geom_boxplot(width = 0.18, outlier.shape = NA, fill = "white", linewidth = 0.25) +
  geom_text(data = n_lab, aes(x = compartment, y = ybot, label = lab),
            inherit.aes = FALSE, size = 1.8, colour = "grey35", family = "Arial",
            vjust = 1.4) +
  geom_text(data = p_lab, aes(x = 1.5, y = ytop, label = lab), inherit.aes = FALSE,
            size = 1.9, colour = "grey25", family = "Arial", vjust = -0.5) +
  facet_wrap(~ section, nrow = 1) +
  scale_fill_manual(values = c("Tumor spots" = "#E64B35", "Other spots" = "#3C5488"),
                    name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.13, 0.12))) +
  labs(x = NULL, y = "TAS module score", title = "Transferred tumor spots") +
  theme_cns(base_size = 7.5) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 6.5),
        plot.title = element_text(size = 9, face = "bold", hjust = 0.5))

# 两张切片叠在一张散点里, CA9 的零膨胀在 x=0 处堆成一条竖条, 看着像画错;
# 改为按切片分面, 并明确写出零值 spot 的比例和 Spearman 的 P。
sp_df2 <- spot_df %>% mutate(section = sub("Tumour", "Tumor", section))
rho_tab <- summary_df %>%
  transmute(section = sub("Tumour", "Tumor", section),
            lab = sprintf("rho = %.3f, %s", spearman_ca9,
                          ifelse(spearman_p < 0.001, "P < 0.001",
                                 sprintf("P = %.3f", spearman_p))))
zero_tab <- sp_df2 %>%
  group_by(section) %>%
  summarise(pct0 = 100 * mean(ca9 == 0, na.rm = TRUE), .groups = "drop") %>%
  mutate(lab = sprintf("CA9 = 0 in %.0f%% of spots", pct0))
xmax <- max(sp_df2$ca9, na.rm = TRUE)
ymax_d <- max(sp_df2$tas, na.rm = TRUE)
ymin_d <- min(sp_df2$tas, na.rm = TRUE)
pD <- ggplot(sp_df2, aes(ca9, tas, colour = section)) +
  geom_point(size = 0.3, alpha = 0.45) +
  geom_text(data = rho_tab, aes(x = xmax, y = ymax_d, label = lab),
            inherit.aes = FALSE, hjust = 1, vjust = -0.3,
            size = 1.9, family = "Arial", colour = "grey25") +
  geom_text(data = zero_tab, aes(x = xmax, y = ymin_d, label = lab),
            inherit.aes = FALSE, hjust = 1, vjust = 1.2,
            size = 1.8, family = "Arial", colour = "grey40") +
  facet_wrap(~ section, nrow = 1) +
  scale_colour_manual(values = c("Tumor 1" = "#0072B2", "Tumor 7" = "#E69F00"),
                      guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.12, 0.12))) +
  labs(x = expression(italic(CA9)~"(log-normalized)"), y = "TAS module score",
       title = "TAS versus CA9, all spots") +
  theme_cns(base_size = 7.5) +
  theme(strip.text = element_text(size = 6.5, face = "bold"),
        plot.title = element_text(size = 9, face = "bold", hjust = 0.5))

fig <- (pA | pB) / (pC | pD) +
  plot_layout(heights = c(1.15, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 11))
save_pub(fig, "Fig8_spatial", OUT, w = 180, h = 175)
message("done")
