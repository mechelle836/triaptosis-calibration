#!/usr/bin/env Rscript
# AUCell as a second single-cell scoring method for the 10 TAS genes.
# Reads the cached GSE159115 Seurat object and AUCell vector; repeats the tests of
# code/29_scrna_section_stats.R and reports agreement with the module score.

suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(tidyr)
})
ROOT <- "/Users/xiaolin/WorkBuddy/2026-09-04-20-00-08/paper2_triaptosis"
RES  <- file.path(ROOT, "results")
seu <- readRDS(file.path(RES, "11_scrna/seurat_kirc_gse159115.rds"))
auc <- readRDS(file.path(RES, "15_aucell/tas_aucell_scores.rds"))
stopifnot(length(auc) == ncol(seu))
if (!is.null(names(auc))) stopifnot(identical(names(auc), colnames(seu)))

df <- data.frame(ctype = as.character(seu$ctype), patient = as.character(seu$patient),
                 module = as.numeric(seu$TAS_score1), auc = as.numeric(auc))
df$is_tumor <- df$ctype == "Tumor cells"

by_type <- df %>% group_by(ctype) %>%
  summarise(n = n(), median_auc = median(auc), frac_nonzero = mean(auc > 0), .groups = "drop") %>%
  arrange(desc(median_auc))
wt <- wilcox.test(auc ~ is_tumor, data = df, exact = FALSE)
by_pt <- df %>% group_by(patient, is_tumor) %>%
  summarise(m = median(auc), mean_auc = mean(auc), .groups = "drop") %>%
  mutate(where = ifelse(is_tumor, "tumor", "rest")) %>% select(-is_tumor) %>%
  pivot_wider(names_from = where, values_from = c(m, mean_auc))
wt_pt <- wilcox.test(by_pt$mean_auc_tumor, by_pt$mean_auc_rest, paired = TRUE, exact = FALSE)
rho <- cor.test(df$auc, df$module, method = "spearman", exact = FALSE)

tests <- data.frame(
  contrast = c("tumor_vs_rest_cells_mean", "tumor_vs_rest_patient_means", "auc_vs_module_spearman"),
  n = c(nrow(df), nrow(by_pt), nrow(df)),
  detail = c(sprintf("mean tumor %.4f vs rest %.4f; median %.4f vs %.4f",
                     mean(df$auc[df$is_tumor]), mean(df$auc[!df$is_tumor]),
                     median(df$auc[df$is_tumor]), median(df$auc[!df$is_tumor])),
             sprintf("%d/%d patients tumor mean higher", sum(by_pt$mean_auc_tumor > by_pt$mean_auc_rest), nrow(by_pt)),
             sprintf("rho %.3f", unname(rho$estimate))),
  p = c(wt$p.value, wt_pt$p.value, rho$p.value)
)
print(by_type); print(by_pt); print(tests)
write.csv(by_type, file.path(RES, "15_aucell/aucell_by_celltype.csv"), row.names = FALSE)
write.csv(tests, file.path(RES, "15_aucell/aucell_tests.csv"), row.names = FALSE)
