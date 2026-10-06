# =====================================================================
# 19b_scrna_validation.R —— 单细胞验证 (干净重写版)
# 数据: GSE159115 (Zhang 2021 PNAS) 7 个 ccRCC 样本 H5 + 官方注释 csv
# 输出: Fig8_scrna_dotplot (气泡=表达细胞%, 颜色=scaled mean expression)
#       results/11_scrna/scrna_dotplot_data.csv
# =====================================================================
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(dplyr); library(tidyr)
  library(ggplot2); library(showtext); library(patchwork)
})
font_add("Arial", regular = "/System/Library/Fonts/Supplemental/Arial.ttf")
showtext_opts(dpi = 300); showtext_auto()
ROOT <- "/Users/xiaolin/WorkBuddy/2026-09-04-20-00-08/paper2_triaptosis"
SCR  <- file.path(ROOT, "data/scrna")
OUT  <- file.path(ROOT, "figures_v2")
RES  <- file.path(ROOT, "results", "11_scrna")
dir.create(RES, showWarnings = FALSE, recursive = TRUE)
source(file.path(ROOT, "code", "lib", "viz_base.R"))

mod <- readRDS(file.path(ROOT, "results/02_model/KIRC/model.rds"))
gene_sel <- mod$genes_sel

# ---- 1. 读注释 + 归并细胞大类 ----------------------------------------
anno <- read.csv(gzfile(file.path(SCR, "GSE159115_ccRCC_anno.csv.gz")))
anno$ctype <- dplyr::case_when(
  anno$anno == "Tumor"                        ~ "Tumor cells",
  anno$anno %in% c("Macro", "Macro_MKI67")    ~ "Macrophages",
  anno$anno %in% c("Tcell", "Tcell_CD8")      ~ "T cells",
  anno$anno %in% c("Bcell", "Plasma", "Mast") ~ "B/Plasma/Mast",
  grepl("^Endo", anno$anno)                   ~ "Endothelial",
  anno$anno %in% c("vSMC", "Peri")            ~ "vSMC/Pericyte",
  TRUE                                        ~ NA_character_)
anno <- anno %>% filter(!is.na(ctype), pct_MT < 0.2)   # 官方 QC 口径
cat("annotated cells:", nrow(anno), "\n")

# ---- 2. 读 7 个 H5: 全转录组 lib size + 10 基因子阵 -------------------
h5_files <- list.files(SCR, pattern = "\\.h5$", full.names = TRUE)
stopifnot(length(h5_files) == 7)

expr_list <- vector("list", 7); lib_list <- vector("list", 7)
for (k in seq_along(h5_files)) {
  counts <- Seurat::Read10X_h5(h5_files[k])
  if (is.list(counts)) counts <- counts$`Gene Expression`
  sm <- regmatches(basename(h5_files[k]),
                   regexpr("SI_[0-9]+", basename(h5_files[k])))
  lib_list[[k]] <- setNames(Matrix::colSums(counts),
                            paste0(sm, "_", colnames(counts)))
  gi <- intersect(gene_sel, rownames(counts))
  stopifnot(length(gi) == length(gene_sel))
  sub <- counts[gi, , drop = FALSE]
  colnames(sub) <- paste0(sm, "_", colnames(sub))
  expr_list[[k]] <- as.matrix(sub)
  cat(sm, ncol(sub), "cells\n")
}
E  <- do.call(cbind, expr_list)
LB <- unlist(lib_list, use.names = FALSE)
# cbind 与 unlist 同序; 顺序对齐校验 (名字可能因 list 拼接有差异)
stopifnot(length(LB) == ncol(E))

# ---- 3. 匹配注释 + CP10K/log1p 归一化 --------------------------------
idx <- match(anno$cell, colnames(E))
anno_e <- anno[!is.na(idx), ]
idx <- idx[!is.na(idx)]
E2 <- E[, idx]
LB2 <- LB[idx]
En <- log1p(sweep(E2, 2, LB2, "/") * 1e4)
cat("cells used:", ncol(E2), " NA in En:", sum(is.na(En)), "\n")

# ---- 4. 聚合 ----------------------------------------------------------
df <- En %>% as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "cell", values_to = "expr") %>%
  mutate(ctype = anno_e$ctype[match(cell, anno_e$cell)])

agg <- df %>%
  group_by(gene, ctype) %>%
  summarise(mean_expr = mean(expr),
            pct_expr  = mean(expr > 0) * 100,
            n_cells   = n(), .groups = "drop") %>%
  mutate(gene = factor(gene, rev(gene_sel)),
         ctype = factor(ctype,
           c("Tumor cells", "Macrophages", "T cells", "B/Plasma/Mast",
             "Endothelial", "vSMC/Pericyte"))) %>%
  arrange(gene, ctype)

agg <- agg %>%
  group_by(gene) %>%
  mutate(scaled = (mean_expr - min(mean_expr)) /
         (max(mean_expr) - min(mean_expr) + 1e-9)) %>%
  ungroup()
stopifnot(sum(is.na(agg$scaled)) == 0)
write.csv(agg, file.path(RES, "scrna_dotplot_data.csv"), row.names = FALSE)

# ---- 5. CNS DotPlot + Tumor/env fold 条形列 ---------------------------
p_dot <- ggplot(agg, aes(ctype, gene)) +
  geom_point(aes(size = pct_expr, colour = scaled)) +
  scale_colour_gradient(low = "#BDD7E7", high = "#08519C",
                        name = "Scaled mean\nexpression") +
  scale_size(range = c(1.5, 6.5), name = "% expressed",
             breaks = c(5, 25, 50, 75)) +
  labs(x = NULL, y = NULL,
       title = "GSE159115: TAS signature genes across ccRCC cell atlas") +
  theme_prism(base_size = 9, base_family = "Arial", base_line_size = 0.5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        plot.margin = margin(5, 2, 2, 2, "mm"),
        axis.text.y = element_text(face = "italic", size = 8.5),
        legend.position = "right",
        plot.title = element_text(size = 10, face = "bold", hjust = 0))

# fold 条形 (与 DotPlot 基因行对齐; SH3GL3 微环境 0 表达 -> tumor-specific)
fold_df <- agg %>% group_by(gene) %>% mutate(g = as.character(gene)) %>%
  ungroup() %>%
  select(gene, ctype, mean_expr) %>%
  pivot_wider(names_from = ctype, values_from = mean_expr) %>%
  rowwise() %>%
  mutate(env_mean = mean(c_across(c(`Macrophages`, `T cells`,
                                    `B/Plasma/Mast`, `Endothelial`,
                                    `vSMC/Pericyte`)))) %>%
  ungroup() %>%
  mutate(fold = `Tumor cells` / (env_mean + 1e-9),
         fold_cap = pmin(fold, 8),
         label = ifelse(fold > 100, "tumor-specific",
                        sprintf("%.1fx", fold)))
p_fold <- ggplot(fold_df, aes(fold_cap, gene)) +
  geom_col(fill = COL_HIGH, width = 0.55, alpha = 0.85) +
  geom_text(aes(label = label, x = fold_cap), hjust = -0.12, size = 2.3,
            family = "Arial") +
  scale_x_continuous(limits = c(0, 10.5), breaks = c(0, 2, 4, 8)) +
  labs(x = "Tumor / microenvironment fold", y = NULL,
       title = "Tumour enrichment") +
  theme_prism(base_size = 8, base_family = "Arial", base_line_size = 0.5) +
  theme(axis.text.y = element_text(face = "italic", size = 8.5),
        axis.text.x = element_text(size = 7),
        plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
        plot.margin = margin(5, 4, 2, 0, "mm"),
        panel.grid.minor = element_blank())

fig8 <- p_dot | p_fold + plot_layout(widths = c(2.3, 1))
ggsave(file.path(OUT, "Fig8_scrna_dotplot.png"), fig8,
       width = 210, height = 110, units = "mm", dpi = 300)
ggsave(file.path(OUT, "Fig8_scrna_dotplot.pdf"), fig8,
       width = 210, height = 110, units = "mm", device = cairo_pdf,
       )
message("saved: Fig8_scrna_dotplot")

# ---- 6. 汇总: Tumor vs 微环境倍数 --------------------------------------
env <- c("Macrophages", "T cells", "B/Plasma/Mast", "Endothelial",
         "vSMC/Pericyte")
sum_df <- agg %>%
  select(gene, ctype, mean_expr) %>%
  pivot_wider(names_from = ctype, values_from = mean_expr) %>%
  rowwise() %>%
  mutate(env_mean = mean(c_across(all_of(env)))) %>%
  ungroup() %>%
  mutate(`Tumor/env fold` = round(`Tumor cells` / (env_mean + 1e-9), 2)) %>%
  select(-env_mean)
print(as.data.frame(sum_df), width = 220)
write.csv(sum_df, file.path(RES, "tumor_vs_env_fold.csv"), row.names = FALSE)