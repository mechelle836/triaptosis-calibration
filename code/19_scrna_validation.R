# =====================================================================
# 19_scrna_validation.R —— 单细胞验证：10 个 TAS 基因在 ccRCC 细胞图谱
# 数据: GSE159115 (Zhang 2021 PNAS, Chinnaiyan lab)
#   - RAW.tar: 14 样本 10X H5 (7 ccRCC + chRCC + normal)
#   - 官方注释: GSE159115_ccRCC_anno.csv.gz (Tumor/Macro/Tcell/Endo/...)
#              GSE159115_normal_anno.csv.gz (PT-A/TAL/... 正常肾节段)
# 输出: Fig8 单细胞验证 DotPlot (气泡=表达%, 颜色=scaled mean expression)
#       + 数值表 results/11_scrna/scrna_dotplot_data.csv
# =====================================================================
suppressPackageStartupMessages({
  library(Seurat); library(jsonlite); library(ggplot2); library(dplyr)
  library(tidyr); library(showtext); library(patchwork); library(ggprism)
})
font_add("Arial", regular = "/System/Library/Fonts/Supplemental/Arial.ttf")
showtext_opts(dpi = 300); showtext_auto()

ROOT <- "/Users/xiaolin/WorkBuddy/2026-09-04-20-00-08/paper2_triaptosis"
SCR  <- file.path(ROOT, "data/scrna")
OUT  <- file.path(ROOT, "figures_v2")
RES  <- file.path(ROOT, "results"); dir.create(file.path(RES, "11_scrna"),
                                               showWarnings = FALSE)
source(file.path(ROOT, "code", "lib", "viz_base.R"))

mod <- readRDS(file.path(RES, "02_model/KIRC/model.rds"))
gene_sel <- mod$genes_sel          # 10 个 TAS 基因

# ---------------------------------------------------------------------
# 1. 解压 RAW.tar
# ---------------------------------------------------------------------
h5_files <- list.files(SCR, pattern = "\\.h5$", full.names = TRUE)
cat("H5 files:", length(h5_files), "\n")

# ---------------------------------------------------------------------
# 2. 读注释, 建 barcode -> cell type 映射 (ccRCC 样本)
# ---------------------------------------------------------------------
anno <- read.csv(gzfile(file.path(SCR, "GSE159115_ccRCC_anno.csv.gz")))
anno$cell <- gsub("^\"|\"$", "", anno$cell)
anno$anno <- gsub("^\"|\"$", "", anno$anno)
anno$sample <- gsub("^\"|\"$", "", anno$sample)
cat("ccRCC annotated cells:", nrow(anno), "\n")
cat("cell types:", paste(sort(unique(anno$anno)), collapse = ", "), "\n")

# 细胞类型归并为 6 大类 (CNS DotPlot 可读性)
anno$ctype <- case_when(
  anno$anno == "Tumor"          ~ "Tumor cells",
  anno$anno %in% c("Macro", "Macro_MKI67") ~ "Macrophages",
  anno$anno %in% c("Tcell", "Tcell_CD8")   ~ "T cells",
  anno$anno %in% c("Bcell", "Plasma", "Mast") ~ "B/Plasma/Mast",
  grepl("^Endo", anno$anno)     ~ "Endothelial",
  anno$anno %in% c("vSMC", "Peri") ~ "vSMC/Pericyte",
  TRUE ~ NA_character_)
cat("merged types:", paste(sort(unique(na.omit(anno$ctype))),
                           collapse = ", "), "\n")

# ---------------------------------------------------------------------
# 3. 逐样本读 H5, 提取 10 基因 + barcodes 匹配
# ---------------------------------------------------------------------
cc_samples <- unique(anno$sample)
expr_list <- list(); lib_list <- list()
for (sm in cc_samples) {
  f <- h5_files[grepl(sm, h5_files)]
  if (length(f) == 0) next
  counts <- Seurat::Read10X_h5(f)
  if (is.list(counts)) counts <- counts$`Gene Expression`
  lib_sm <- Matrix::colSums(counts)                 # 全转录组 lib size (关键!)
  names(lib_sm) <- paste0(sm, "_", colnames(counts))
  gi <- intersect(gene_sel, rownames(counts))
  sub <- counts[gi, , drop = FALSE]
  colnames(sub) <- paste0(sm, "_", colnames(sub))   # 对齐 anno cell 命名
  expr_list[[sm]] <- as.matrix(sub)
  lib_list[[sm]] <- lib_sm
  cat(" ", sm, dim(sub)[2], "cells,", length(gi), "/", length(gene_sel),
      "genes found\n")
}
expr_all <- do.call(cbind, expr_list)
lib_all <- unlist(lib_list, use.names = FALSE)

# 对齐注释
common_cells <- intersect(colnames(expr_all), anno$cell)
cat("matched cells:", length(common_cells), "\n")
E <- expr_all[, common_cells]
ct_vals <- anno$ctype[match(common_cells, anno$cell)]
ct <- ct_vals; names(ct) <- common_cells

# CP10K + log1p 归一化 (全转录组 lib size, 10 基因子集会全变 NaN)
lib <- lib_all[colnames(E)]
E_norm <- sweep(E, 2, lib, "/") * 1e4
E_norm <- log1p(E_norm)

# ---------------------------------------------------------------------
# 4. 聚合: per cell type -> mean expression (log) + percent expressed
# ---------------------------------------------------------------------
df_long <- E_norm %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "cell", values_to = "expr") %>%
  mutate(ctype = ct[cell]) %>%
  filter(!is.na(ctype))

agg <- df_long %>%
  group_by(gene, ctype) %>%
  summarise(mean_expr = mean(expr),
            pct_expr  = mean(expr > 0) * 100,
            n_cells   = n(), .groups = "drop")

agg_all <- agg %>%
  mutate(ctype = factor(ctype,
    c("Tumor cells", "Macrophages", "T cells",
      "B/Plasma/Mast", "Endothelial", "vSMC/Pericyte")),
    gene = factor(gene, rev(gene_sel)))

# ---------------------------------------------------------------------
# 5. CNS DotPlot (气泡大小 = %expressed, 颜色 = scaled mean)
# ---------------------------------------------------------------------
agg_all <- agg_all %>%
  group_by(gene) %>%
  mutate(scaled = (mean_expr - min(mean_expr)) /
         (max(mean_expr) - min(mean_expr) + 1e-9)) %>%
  ungroup()

p <- ggplot(agg_all, aes(ctype, gene)) +
  geom_point(aes(size = pct_expr, colour = scaled)) +
  scale_colour_gradient(low = "#BDD7E7", high = "#08519C",
                        name = "Scaled mean\nexpression") +
  scale_size(range = c(1.5, 6.5), name = "% expressed",
             breaks = c(5, 25, 50, 75)) +
  labs(x = NULL, y = NULL,
       title = "GSE159115: TAS signature genes across ccRCC cell atlas") +
  theme_prism(base_size = 9, base_family = "Arial", base_line_size = 0.5) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 8),
        axis.text.y = element_text(face = "italic", size = 8.5),
        legend.position = "right",
        plot.title = element_text(size = 10, face = "bold", hjust = 0))

ggsave(file.path(OUT, "Fig8_scrna_dotplot.png"), p,
       width = 180, height = 110, units = "mm", dpi = 300)
ggsave(file.path(OUT, "Fig8_scrna_dotplot.pdf"), p,
       width = 180, height = 110, units = "mm", device = cairo_pdf,
       )
message("saved: Fig8_scrna_dotplot")

# 汇总: Tumor vs 微环境均值倍数
env_types <- c("Macrophages", "T cells", "B/Plasma/Mast", "Endothelial",
               "vSMC/Pericyte")
sum_df <- agg_all %>%
  select(gene, ctype, mean_expr) %>%
  pivot_wider(names_from = ctype, values_from = mean_expr) %>%
  rowwise() %>%
  mutate(env_mean = mean(c_across(all_of(env_types)))) %>%
  ungroup() %>%
  mutate(`Tumor/env fold` = round(`Tumor cells` / (env_mean + 1e-9), 2)) %>%
  select(-env_mean)
print(as.data.frame(sum_df), width = 200)
write.csv(sum_df, file.path(RES, "11_scrna/tumor_vs_env_fold.csv"),
          row.names = FALSE)
