# =====================================================================
# 19c_scrna_umap.R —— 单细胞验证完整组图 (Seurat 标准流程)
# A UMAP by cell type | B UMAP by TAS module score | C DotPlot+fold
# D Violin (WDR91/SH3GL3) | E 细胞比例 per patient
# =====================================================================
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(dplyr); library(tidyr)
  library(ggplot2); library(showtext); library(patchwork)
  library(RColorBrewer)
})
ROOT <- "/Users/xiaolin/WorkBuddy/2026-09-04-20-00-08/paper2_triaptosis"
SCR  <- file.path(ROOT, "data/scrna")
OUT  <- file.path(ROOT, "figures_v2")
RES  <- file.path(ROOT, "results", "11_scrna")
source(file.path(ROOT, "code", "lib", "viz_base.R"))
showtext_auto(FALSE)

mod <- readRDS(file.path(ROOT, "results/02_model/KIRC/model.rds"))
gene_sel <- mod$genes_sel

# ---- 1. 官方注释 ------------------------------------------------------
anno <- read.csv(gzfile(file.path(SCR, "GSE159115_ccRCC_anno.csv.gz"))) %>%
  mutate(ctype = dplyr::case_when(
    anno == "Tumor"                        ~ "Tumor cells",
    anno %in% c("Macro", "Macro_MKI67")    ~ "Macrophages",
    anno %in% c("Tcell", "Tcell_CD8")      ~ "T cells",
    anno %in% c("Bcell", "Plasma", "Mast") ~ "B/Plasma/Mast",
    grepl("^Endo", anno)                   ~ "Endothelial",
    anno %in% c("vSMC", "Peri")            ~ "vSMC/Pericyte",
    TRUE                                   ~ NA_character_)) %>%
  filter(!is.na(ctype), pct_MT < 0.2)
cat("annotated:", nrow(anno), "\n")

# ---- 2. Seurat 标准流程 (缓存 rds) ------------------------------------
seur_rds <- file.path(RES, "seurat_kirc_gse159115.rds")
if (file.exists(seur_rds)) {
  seu <- readRDS(seur_rds)
  cat("loaded cached seurat:", ncol(seu), "cells\n")
} else {
  h5_files <- list.files(SCR, pattern = "\\.h5$", full.names = TRUE)
  obj_list <- lapply(h5_files, function(f) {
    counts <- Seurat::Read10X_h5(f)
    if (is.list(counts)) counts <- counts$`Gene Expression`
    sm <- regmatches(basename(f), regexpr("SI_[0-9]+", basename(f)))
    colnames(counts) <- paste0(sm, "_", colnames(counts))  # 防 merge 后缀
    obj <- CreateSeuratObject(counts, project = sm,
                              min.cells = 3, min.features = 200)
    obj$patient <- sm
    mt.genes <- grep("^MT-", rownames(obj), value = TRUE)
    obj$percent.mt <- if (length(mt.genes))
      Matrix::colSums(obj[mt.genes, ]) / Matrix::colSums(obj) * 100 else 0
    obj
  })
  seu <- merge(obj_list[[1]], obj_list[-1])
  seu <- JoinLayers(seu)
  # 合并对象上统一重算 QC 指标 (Seurat v5: 显式取 counts 层, 不用 seu[genes,])
  cnts <- GetAssayData(seu, layer = "counts")
  mt_i <- grep("^MT-", rownames(cnts))
  seu$percent.mt <- Matrix::colSums(cnts[mt_i, , drop = FALSE]) /
    Matrix::colSums(cnts) * 100
  keep <- seu$nFeature_RNA > 200 & seu$percent.mt < 20
  seu <- seu[, keep]
  # 官方注释映射
  idx <- match(colnames(seu), anno$cell)
  keep2 <- !is.na(idx)
  seu$ctype <- anno$ctype[idx]
  seu$ctype[is.na(seu$ctype)] <- NA
  seu <- seu[, which(!is.na(seu$ctype))]
  seu <- NormalizeData(seu, verbose = FALSE)
  seu <- FindVariableFeatures(seu, nfeatures = 2000, verbose = FALSE)
  seu <- ScaleData(seu, verbose = FALSE)
  seu <- RunPCA(seu, npcs = 30, verbose = FALSE)
  seu <- RunUMAP(seu, dims = 1:30, umap.method = "uwot",
                 n.neighbors = 30, min.dist = 0.3, verbose = FALSE)
  # TAS module score
  tas_in <- intersect(gene_sel, rownames(seu))
  seu <- AddModuleScore(seu, features = list(tas_in), name = "TAS_score")
  seu$ctype <- factor(seu$ctype,
    c("Tumor cells", "Macrophages", "T cells", "B/Plasma/Mast",
      "Endothelial", "vSMC/Pericyte"))
  saveRDS(seu, seur_rds)
  cat("seurat built:", ncol(seu), "cells\n")
}
cat("final cells:", ncol(seu), "\n")

pal6 <- c("Tumor cells" = "#E64B35", "Macrophages" = "#F39B7F",
          "T cells" = "#3C5488", "B/Plasma/Mast" = "#8491B4",
          "Endothelial" = "#00A087", "vSMC/Pericyte" = "#B09C85")
cat("== DIAG: seu$ctype table ==\n")
print(table(seu$ctype, useNA = "ifany"))
umap_df <- data.frame(UMAP1 = Embeddings(seu, "umap")[, 1],
                      UMAP2 = Embeddings(seu, "umap")[, 2],
                      ctype = seu$ctype,
                      patient = seu$patient,
                      TAS_score = seu$TAS_score1)
cat("== DIAG: umap_df ctype NA:", sum(is.na(umap_df$ctype)),
    " patient NA:", sum(is.na(umap_df$patient)), "==\n")

# ---- A: UMAP by cell type --------------------------------------------
pA <- ggplot(umap_df, aes(UMAP1, UMAP2, colour = ctype)) +
  geom_point(size = 0.25, alpha = 0.6, show.legend = TRUE) +
  scale_colour_manual(values = pal6, name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 2.5))) +
  labs(title = "Cell types", x = "UMAP-1", y = "UMAP-2") +
  theme_classic(base_size = 8, base_family = "Arial") +
  theme(legend.position = c(0.02, 0.02), legend.justification = c(0, 0),
        legend.text = element_text(size = 6),
        legend.key.size = unit(2.2, "mm"),
        legend.background = element_rect(fill = "white", colour = NA, alpha = 0.6),
        plot.title = element_text(size = 9.5, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 6), axis.title = element_text(size = 7.5),
        axis.ticks = element_line(linewidth = 0.3),
        plot.margin = margin(2, 2, 2, 2, "mm"))

# ---- B: UMAP by TAS module score -------------------------------------
pB <- ggplot(umap_df, aes(UMAP1, UMAP2, colour = TAS_score)) +
  geom_point(size = 0.25, alpha = 0.65) +
  scale_colour_gradient2(low = "#4575B4", mid = "#F7F7F7", high = "#D73027",
                         midpoint = median(umap_df$TAS_score),
                         name = "TAS module\nscore") +
  labs(title = "TAS signature projection", x = "UMAP-1", y = "UMAP-2") +
  theme_classic(base_size = 8, base_family = "Arial") +
  theme(legend.position = c(0.02, 0.02), legend.justification = c(0, 0),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 5.5),
        legend.key.height = unit(3.2, "mm"),
        legend.background = element_rect(fill = "white", colour = NA, alpha = 0.6),
        plot.title = element_text(size = 9.5, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 6), axis.title = element_text(size = 7.5),
        axis.ticks = element_line(linewidth = 0.3),
        plot.margin = margin(2, 2, 2, 2, "mm"))

# ---- C: DotPlot + fold (复用 19b 逻辑) --------------------------------
LB_lib <- NULL  # lib size per cell
counts_tas <- t(FetchData(seu, vars = gene_sel, layer = "counts"))  # 10 基因 x 细胞
lib_all <- Matrix::colSums(GetAssayData(seu, layer = "counts"))
En <- log1p(sweep(counts_tas, 2, lib_all[colnames(counts_tas)], "/") * 1e4)
cat("== DIAG counts_tas dim:", dim(counts_tas), " NA:", sum(is.na(counts_tas)),
    " lib NA:", sum(is.na(lib_all[colnames(counts_tas)])), "==\n")
df_long <- En %>% as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "cell", values_to = "expr") %>%
  mutate(ctype = seu$ctype[match(cell, colnames(seu))])
cat("== DIAG df_long rows:", nrow(df_long), " NA expr:", sum(is.na(df_long$expr)),
    " NA ctype:", sum(is.na(df_long$ctype)), "==\n")
agg <- df_long %>%
  group_by(gene, ctype) %>%
  summarise(mean_expr = mean(expr), pct_expr = mean(expr > 0) * 100,
            n_cells = n(), .groups = "drop") %>%
  filter(is.finite(mean_expr)) %>%
  mutate(gene = factor(gene, rev(gene_sel)),
         ctype = factor(ctype, names(pal6))) %>%
  arrange(gene, ctype) %>%
  group_by(gene) %>%
  mutate(scaled = (mean_expr - min(mean_expr)) /
         (max(mean_expr) - min(mean_expr) + 1e-9)) %>%
  ungroup()
stopifnot(sum(is.na(agg$scaled)) == 0)
cat("== DIAG agg rows:", nrow(agg), " genes:", length(unique(agg$gene)),
    " ctypes:", length(unique(agg$ctype)), "==\n")
cat("== DIAG mean_expr range:", range(agg$mean_expr), "==\n")

p_dot <- ggplot(agg, aes(ctype, gene)) +
  geom_point(aes(size = pct_expr, colour = scaled)) +
  scale_colour_gradient(low = "#BDD7E7", high = "#08519C",
                        name = "Scaled mean\nexpression") +
  scale_size(range = c(1.2, 5.5), name = "% expressed",
             breaks = c(10, 40, 70)) +
  labs(x = NULL, y = NULL) +
  theme_prism(base_size = 8.5, base_family = "Arial", base_line_size = 0.5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7.5),
        axis.text.y = element_text(face = "italic", size = 8),
        legend.position = "right", legend.box = "vertical",
        legend.title = element_text(size = 6.5),
        legend.text = element_text(size = 5.5),
        legend.key.size = unit(2.6, "mm"),
        plot.margin = margin(4, 2, 2, 2, "mm"))

env <- c("Macrophages", "T cells", "B/Plasma/Mast", "Endothelial",
         "vSMC/Pericyte")
fold_df <- agg %>%
  group_by(gene) %>%
  summarise(tumor = mean_expr[ctype == "Tumor cells"],
            env_mean = mean(mean_expr[ctype %in% env]), .groups = "drop") %>%
  mutate(fold = tumor / (env_mean + 1e-9),
         fold_cap = pmin(fold, 8),
         label = ifelse(fold > 100, "tumor-specific", sprintf("%.1fx", fold)),
         gene = factor(gene, rev(gene_sel)))

p_fold <- ggplot(fold_df, aes(fold_cap, gene)) +
  geom_col(fill = COL_HIGH, width = 0.55, alpha = 0.85) +
  geom_text(aes(label = label, x = fold_cap), hjust = -0.12, size = 2.2,
            family = "Arial") +
  scale_x_continuous(limits = c(0, 12.5), breaks = c(0, 2, 4, 8)) +
  labs(x = "Tumor / micro. fold", y = NULL, title = "Tumour enrichment") +
  theme_prism(base_size = 8, base_family = "Arial", base_line_size = 0.5) +
  theme(axis.text.y = element_text(face = "italic", size = 8),
        axis.text.x = element_text(size = 6.5),
        plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
        plot.margin = margin(4, 4, 2, 0, "mm"))
pC <- p_dot | p_fold + plot_layout(widths = c(2.1, 1))

# ---- D: Violin WDR91 + SH3GL3 ----------------------------------------
vln_genes <- intersect(c("WDR91", "SH3GL3"), rownames(seu))
vln_df <- FetchData(seu, vars = c(vln_genes, "ctype")) %>%
  pivot_longer(-ctype, names_to = "gene", values_to = "expr") %>%
  mutate(gene = factor(gene, vln_genes),
         ctype = factor(ctype, names(pal6)))
pD <- ggplot(vln_df, aes(ctype, expr, fill = ctype)) +
  geom_violin(scale = "width", linewidth = 0.3, show.legend = FALSE,
              alpha = 0.85) +
  scale_fill_manual(values = pal6) +
  facet_wrap(~gene, scales = "free_y") +
  stat_summary(fun = mean, geom = "point", size = 1.2, show.legend = FALSE) +
  labs(x = NULL, y = NULL,
       title = "Tumour-enriched genes (log1p CP10K)") +
  theme_prism(base_size = 8, base_family = "Arial") +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 6.5),
        strip.text = element_text(size = 8.5, face = "bold.italic"),
        axis.title.y = element_text(size = 7.5),
        plot.title = element_text(size = 9.5, face = "bold", hjust = 0.5),
        plot.margin = margin(4, 2, 2, 2, "mm"))

# ---- E: 细胞比例 per patient ------------------------------------------
prop_df <- umap_df %>%
  count(patient, ctype) %>%
  group_by(patient) %>% mutate(p = n / sum(n)) %>% ungroup() %>%
  mutate(patient = factor(patient, sort(unique(patient))),
         ctype = factor(ctype, names(pal6)))
pE <- ggplot(prop_df, aes(patient, p, fill = ctype)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.2) +
  scale_fill_manual(values = pal6, name = NULL) +
  guides(fill = guide_legend(nrow = 1, override.aes = list(size = 2))) +
  labs(x = NULL, y = "Fraction",
       title = "Cell composition per patient") +
  theme_classic(base_size = 8, base_family = "Arial") +
  theme(legend.position = "bottom", legend.text = element_text(size = 5.5),
        legend.key.size = unit(2, "mm"),
        axis.text.x = element_text(size = 6.5, angle = 30, hjust = 1),
        axis.title.y = element_text(size = 7.5),
        plot.title = element_text(size = 9.5, face = "bold", hjust = 0.5),
        plot.margin = margin(4, 2, 2, 2, "mm"))

# ---- 组合: 3 UMAP 横排 / C 跨两列 + D / E ----
fig8 <- pA | pB
fig8 <- fig8 / pC / (pD | pE) +
  plot_annotation(tag_levels = "A") +
  plot_layout(heights = c(1, 1.05, 0.95))
ggsave(file.path(OUT, "Fig8_scrna_dotplot.png"), fig8,
       width = 210, height = 190, units = "mm", dpi = 300)
ggsave(file.path(OUT, "Fig8_scrna_dotplot.pdf"), fig8,
       width = 210, height = 190, units = "mm", device = cairo_pdf,
       )
message("saved: Fig8_scrna_full")