#!/usr/bin/env Rscript
# =============================================================================
# 13_multivariable_cox.R — Paper2-Triaptosis KIRC 多变量 Cox (2026-09-02)
#
# 背景: 03_focus_model.R 的"多变量 Cox"实际只拟合了 TAS 单变量。本脚本补:
#   1. TAS + AGE + SEX + AJCC 分期 的多变量 Cox (调整后 HR)
#   2. 分期亚组 KM (I-II vs III-IV 内 TAS High/Low)
#   3. 聚类亚型与 TAS 的独立性交叉表
# 临床来源: data/external_validation/KIRC_clin_ext.csv (10b 号脚本, kirc_tcga
#           或 PanCanAtlas 端点), 与 results/02_model/KIRC/score_table.csv 对接
# 输出: results/08_multivariable/
# =============================================================================
suppressPackageStartupMessages(library(survival))
P2  <- "/Users/xiaolin/WorkBuddy/2026-09-04-20-00-08/paper2_triaptosis"
OUT <- file.path(P2, "results", "08_multivariable")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

score <- read.csv(file.path(P2, "results/02_model/KIRC/score_table.csv"),
                  stringsAsFactors = FALSE, row.names = 1)
score$patientId <- substr(rownames(score), 1, 12)
# 连续 TAS 以 TCGA-KIRC 全队列 1 个 SD 为单位报告 HR
score$TAS <- as.numeric(scale(score$TAS))
cat("score_table:", nrow(score), "patients,", sum(score$status), "deaths\n")

# ---------- 临床扩展: 找到可用的抓取产物 ----------
cand <- c(file.path(P2, "data/external_validation/kirc_tcga_pan_can_atlas_2018_clin_patient.csv"),
          file.path(P2, "data/external_validation/kirc_tcga_clin_patient.csv"),
          file.path(P2, "data/external_validation/KIRC_clin_ext.csv"))
clin_f <- cand[file.exists(cand)]
if (length(clin_f) == 0) stop("未找到 KIRC 临床扩展文件, 先跑 10b_fetch_kirc_clinical_ext.py")
clin <- read.csv(clin_f[1], stringsAsFactors = FALSE)
cat("clin file:", clin_f[1], " rows:", nrow(clin), " cols:", paste(colnames(clin), collapse = ","), "\n")

if (!"AGE" %in% colnames(clin) && "PATIENT_AGE" %in% colnames(clin))
  clin$AGE <- clin$PATIENT_AGE

# 分期归一: "Stage I" / "Stage II" ... -> I/II/III/IV
stage_col <- intersect(c("AJCC_PATHOLOGIC_TUMOR_STAGE", "AJCC_TUMOR_STAGE", "TUMOR_STAGE"),
                       colnames(clin))[1]
if (!is.na(stage_col)) {
  st <- toupper(trimws(clin[[stage_col]]))
  st <- gsub("^STAGE\\s*", "", st)
  st[!st %in% c("I", "II", "III", "IV")] <- NA
  clin$STAGE <- factor(st, levels = c("I", "II", "III", "IV"))
}
clin$AGE <- suppressWarnings(as.numeric(clin$AGE))
if ("SEX" %in% colnames(clin)) clin$SEX <- factor(toupper(trimws(clin$SEX)))

d <- merge(score[, c("patientId", "time", "status", "TAS", "cluster", "group")],
           clin[, intersect(c("patientId", "AGE", "SEX", "STAGE", "GRADE"), colnames(clin))],
           by = "patientId")
cat("merged:", nrow(d), "patients,", sum(d$status), "deaths\n")
cat("非缺失: AGE", sum(!is.na(d$AGE)), " SEX", sum(!is.na(d$SEX)),
    " STAGE", sum(!is.na(d$STAGE)),
    if ("GRADE" %in% colnames(d)) sprintf(" GRADE %d", sum(!is.na(d$GRADE))) else "", "\n")

# ---------- 1. 多变量 Cox ----------
terms <- c("TAS", "AGE", "SEX", "STAGE")[c(TRUE, !all(is.na(d$AGE)),
                                          !all(is.na(d$SEX)), !all(is.na(d$STAGE)))]
cat("\n--- 单变量 vs 多变量 Cox ---\n")
rows <- list()
fit_u <- coxph(Surv(time, status) ~ TAS, data = d)
s_u <- summary(fit_u)
rows[[1]] <- data.frame(model = "Univariate", term = "TAS",
                        HR = s_u$conf.int["TAS", "exp(coef)"],
                        lo = s_u$conf.int["TAS", "lower .95"],
                        hi = s_u$conf.int["TAS", "upper .95"],
                        p = s_u$coefficients["TAS", "Pr(>|z|)"])
form <- as.formula(paste("Surv(time, status) ~", paste(terms, collapse = " + ")))
fit_m <- coxph(form, data = d)
s_m <- summary(fit_m)
for (t in terms) {
  rn <- rownames(s_m$coefficients)
  hits <- rn[grepl(paste0("^", t, "(II|III|IV|MALE|FEMALE|>|<)"), rn) | rn == t]
  for (r in hits) {
    rows[[length(rows) + 1]] <- data.frame(
      model = "Multivariable", term = r,
      HR = s_m$conf.int[r, "exp(coef)"], lo = s_m$conf.int[r, "lower .95"],
      hi = s_m$conf.int[r, "upper .95"],
      p = s_m$coefficients[r, "Pr(>|z|)"])
  }
}
res <- do.call(rbind, rows)
res$HR <- round(res$HR, 3); res$lo <- round(res$lo, 3); res$hi <- round(res$hi, 3)
res$p <- signif(res$p, 3)
print(res)
write.csv(res, file.path(OUT, "multivariable_cox.csv"), row.names = FALSE)
cat(sprintf("\n多变量模型 C-index = %.3f  (单变量 %.3f)\n",
            s_m$concordance["C"], s_u$concordance["C"]))

# ---------- 2. 分期亚组 ----------
if (!all(is.na(d$STAGE))) {
  cat("\n--- 分期亚组 (TAS 中位数分层 KM) ---\n")
  sub_rows <- list()
  for (sub in list(c("I", "II"), c("III", "IV"))) {
    dd <- d[d$STAGE %in% sub, ]
    if (nrow(dd) < 30) next
    dd$group <- factor(dd$group, levels = c("Low", "High"))  # 基线=Low
    lr <- survdiff(Surv(time, status) ~ group, data = dd)
    p <- 1 - pchisq(lr$chisq, df = 1)
    fx <- coxph(Surv(time, status) ~ group, data = dd)
    ci <- exp(confint(fx))
    sub_rows[[length(sub_rows) + 1]] <- data.frame(
      subgroup = paste(sub, collapse = "-"), n = nrow(dd),
      events = sum(dd$status),
      HR_High_vs_Low = round(exp(coef(fx)[1]), 2),
      lo = round(ci[1], 2), hi = round(ci[2], 2),
      cox_p = signif(summary(fx)$coefficients[1, 5], 3), km_p = signif(p, 3))
  }
  if (length(sub_rows)) {
    subs <- do.call(rbind, sub_rows)
    print(subs)
    write.csv(subs, file.path(OUT, "stage_subgroup_km.csv"), row.names = FALSE)
  }
}

# ---------- 3. TAS 分组 x 分期交叉 ----------
if (!all(is.na(d$STAGE))) {
  tab <- table(d$group, d$STAGE)
  cat("\n--- TAS 分组 x 分期 ---\n"); print(tab)
  write.csv(data.frame(stage = colnames(tab), Low = tab["Low", ], High = tab["High", ]),
            file.path(OUT, "group_by_stage.csv"), row.names = FALSE)
}
cat("\n>>> DONE", OUT, "\n")
