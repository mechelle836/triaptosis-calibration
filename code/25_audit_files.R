# =====================================================================
# 25_audit_files.R —— 生成 3 个可复核文件 (数字审计闭环)
# 1. ibs.csv: IBS 汇总 (07 公式梯形积分, 完整曲线)
# 2. cluster_assignment.csv: k-means 聚类分组明细
# 3. stage_four_km.csv: 四分期 KM P 值明细
# =====================================================================
suppressPackageStartupMessages({
  library(survival); library(dplyr)
})
ROOT <- "/Users/xiaolin/WorkBuddy/2026-09-04-20-00-08/paper2_triaptosis"
RES <- file.path(ROOT, "results")

# ---- 1. ibs.csv (按 07 公式: 梯形积分 / max(times)) -------------------
brier <- read.csv(file.path(RES, "05_validation/brier.csv"))
ibs_one <- function(model_name) {
  b <- brier$Brier[brier$model == model_name]
  t <- brier$times[brier$model == model_name]
  sum(diff(c(0, t)) * (c(0, b[-length(b)]) + b) / 2) / max(t)
}
ibs_tas <- ibs_one("TAS"); ibs_km <- ibs_one("KM")
ibs_df <- data.frame(
  metric = "integrated_Brier_score",
  TAS = round(ibs_tas, 4), KM = round(ibs_km, 4),
  relative_gain_pct = round(100 * (ibs_km - ibs_tas) / ibs_km, 1),
  formula = "trapezoid integral of IPCW Brier over 12-60 months / 60")
write.csv(ibs_df, file.path(RES, "05_validation/ibs.csv"), row.names = FALSE)
cat("ibs.csv: TAS=", round(ibs_tas, 4), " KM=", round(ibs_km, 4), "\n")

# ---- 2. cluster_assignment.csv (k-means 聚类分组) ----------------------
sc <- read.csv(file.path(RES, "02_model/KIRC/score_table.csv"), row.names = 1)
cl_assign <- data.frame(sample = rownames(sc), cluster = sc$cluster,
                        cluster_lab = paste0("C", sc$cluster))
write.csv(cl_assign, file.path(RES, "02_model/KIRC/cluster_assignment.csv"),
          row.names = FALSE)
cat("cluster_assignment.csv: C1=", sum(sc$cluster == 1),
    " C2=", sum(sc$cluster == 2), "\n")

# ---- 3. stage_four_km.csv (四分期 KM P 值明细) -------------------------
clin <- read.csv(file.path(ROOT, "data/external_validation/kirc_tcga_pan_can_atlas_2018_clin_patient.csv")) %>%
  mutate(stage4 = trimws(gsub("STAGE", "", AJCC_PATHOLOGIC_TUMOR_STAGE, ignore.case = TRUE)),
         stage_grp = case_when(stage4 %in% c("I") ~ "Stage I",
                               stage4 %in% c("II") ~ "Stage II",
                               stage4 %in% c("III") ~ "Stage III",
                               stage4 %in% c("IV") ~ "Stage IV",
                               TRUE ~ NA_character_))
sd <- sc %>% mutate(patientId = rownames(sc)) %>%
  left_join(clin %>% select(patientId, stage_grp), by = "patientId")
stage_km <- lapply(c("Stage I", "Stage II", "Stage III", "Stage IV"), function(st) {
  d <- sd %>% filter(stage_grp == st) %>% mutate(group = factor(group, c("Low", "High")))
  if (length(unique(d$group)) < 2) return(NULL)
  pv <- survdiff(Surv(time, status) ~ group, data = d)
  p <- 1 - pchisq(pv$chisq, length(pv$n) - 1)
  data.frame(stage = st, n = nrow(d),
             n_low = sum(d$group == "Low"), n_high = sum(d$group == "High"),
             km_logrank_p = signif(p, 3))
})
stage_km <- do.call(rbind, stage_km)
write.csv(stage_km, file.path(RES, "08_multivariable/stage_four_km.csv"),
          row.names = FALSE)
cat("stage_four_km.csv:\n"); print(stage_km)
message("3 audit files saved")