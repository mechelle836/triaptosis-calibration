# =====================================================================
# 29_scrna_section_stats.R
# 只读已缓存的 GSE159115 Seurat 对象，补细胞类型和患者水平的 TAS 数字。
# 不重跑 CellChat，不改对象。
# 输出:
#   results/11_scrna/tas_by_celltype.csv
#   results/11_scrna/tas_by_patient.csv
# =====================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
})
ROOT <- "/Users/xiaolin/WorkBuddy/2026-09-04-20-00-08/paper2_triaptosis"
RES  <- file.path(ROOT, "results", "11_scrna")
seu <- readRDS(file.path(RES, "seurat_kirc_gse159115.rds"))

df <- data.frame(
  ctype   = as.character(seu$ctype),
  patient = as.character(seu$patient),
  tas     = as.numeric(seu$TAS_score1),
  stringsAsFactors = FALSE
)
df$is_tumor <- df$ctype == "Tumor cells"

by_type <- df %>%
  group_by(ctype) %>%
  summarise(
    n = n(),
    median_tas = median(tas),
    q1 = as.numeric(quantile(tas, 0.25)),
    q3 = as.numeric(quantile(tas, 0.75)),
    .groups = "drop"
  ) %>%
  arrange(desc(median_tas))

# 细胞水平：肿瘤细胞 vs 其余全部细胞（n 很大，p 值只作描述）
wt_cell <- wilcox.test(tas ~ is_tumor, data = df, exact = FALSE)
cell_row <- data.frame(
  contrast = "tumor_vs_rest_cells",
  n_tumor = sum(df$is_tumor),
  n_rest = sum(!df$is_tumor),
  median_tumor = median(df$tas[df$is_tumor]),
  median_rest = median(df$tas[!df$is_tumor]),
  statistic = unname(wt_cell$statistic),
  p = wt_cell$p.value,
  stringsAsFactors = FALSE
)

by_pt <- df %>%
  group_by(patient, is_tumor) %>%
  summarise(n = n(), median_tas = median(tas), .groups = "drop") %>%
  mutate(where = ifelse(is_tumor, "tumor", "rest")) %>%
  select(-is_tumor) %>%
  pivot_wider(names_from = where, values_from = c(n, median_tas))

wt_pt <- wilcox.test(by_pt$median_tas_tumor, by_pt$median_tas_rest,
                     paired = TRUE, exact = FALSE)
pt_row <- data.frame(
  contrast = "tumor_vs_rest_patient_medians",
  n_patients = nrow(by_pt),
  n_tumor_higher = sum(by_pt$median_tas_tumor > by_pt$median_tas_rest),
  median_of_patient_tumor = median(by_pt$median_tas_tumor),
  median_of_patient_rest = median(by_pt$median_tas_rest),
  statistic = unname(wt_pt$statistic),
  p = wt_pt$p.value,
  stringsAsFactors = FALSE
)

write.csv(by_type, file.path(RES, "tas_by_celltype.csv"), row.names = FALSE)
write.csv(by_pt, file.path(RES, "tas_by_patient.csv"), row.names = FALSE)
write.csv(rbind(
  data.frame(contrast = cell_row$contrast, n = cell_row$n_tumor + cell_row$n_rest,
             detail = sprintf("median tumor %.4f vs rest %.4f",
                              cell_row$median_tumor, cell_row$median_rest),
             statistic = cell_row$statistic, p = cell_row$p),
  data.frame(contrast = pt_row$contrast, n = pt_row$n_patients,
             detail = sprintf("%d/%d patients tumor median higher; patient-median %.4f vs %.4f",
                              pt_row$n_tumor_higher, pt_row$n_patients,
                              pt_row$median_of_patient_tumor, pt_row$median_of_patient_rest),
             statistic = pt_row$statistic, p = pt_row$p)
), file.path(RES, "tas_localization_tests.csv"), row.names = FALSE)

message("cell types:")
print(by_type)
message("patients:")
print(by_pt)
message("tests:")
print(cell_row)
print(pt_row)
