#!/usr/bin/env Rscript
# 44_logscale_sensitivity.R
# 方法学审查 #2 的回应：表达尺度敏感性分析。
# 手稿主分析使用非对数 RSEM 尺度 (§2.1 "original, non-log scale")，本脚本用
# log2(RSEM+1) 尺度重跑同一管线，检验 TRRS 构建与验证结论是否对尺度选择稳健。
#
# 方法依据（非自造）：RNA-seq 计数建模的方差稳定化标准实践
#   (Law et al. 2014, Genome Biol 15:R29, voom)；尺度敏感性属于 STROBE/STRATOS
#   框架下的标准敏感性分析 (Sauerbrei et al. 2014, Stat Med)。
#
# 关键事实：UCSC Xena EB++ 矩阵本身即 log2(RSEM+1) 格式——log 尺度数据 = Xena
#   原始值；linear 尺度 = 2^x - 1（38_screen_null.R 第 58 行的反变换）。
#
# 设计：linear 与 log 两端用完全相同的流程（21 基因 elastic-net alpha=0.5、
#   lambda.min、队列内 z-score、冻结系数外部验证、5 折 OOF 折内重选 lambda），
#   对比：选中基因、per-SD HR（单变量/多变量校正）、apparent C、OOF C、
#   E-MTAB-1980 / CPTAC-3 的 C-index。
#
# 输出: results/08_scale_sensitivity/scale_sensitivity_summary.csv
#       results/08_scale_sensitivity/scale_sensitivity_genes.csv

suppressPackageStartupMessages({
  library(data.table); library(jsonlite); library(survival); library(glmnet)
})
P2  <- "/Users/xiaolin/WorkBuddy/2026-09-04-20-00-08/paper2_triaptosis"
MP  <- "/Users/xiaolin/CodeBuddy/Claw/pan-cancer_disulfidptosis_project"
OUT <- file.path(P2, "results/08_scale_sensitivity")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
set.seed(2026)

TRIA21 <- c("PIK3C3","PIK3R4","MTM1","KEAP1","NFE2L2","SH3GL3","ELMO2","NCKAP1","ACTN2",
            "KXD1","WDR91","FYCO1","RAB9A","SCARB2","ANXA8","CLCN4","ATP13A2",
            "UVRAG","ATG14","NRBF2","BECN1")
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a)) b else a

# ---------- KIRC 临床（与 38/02 口径一致：每患者第一份样本） ----------
clin <- fromJSON(file.path(P2, "data/KIRC_clin.json"))
os_status <- sapply(clin, function(x) as.character(x$OS_STATUS %||% ""))
os_months <- suppressWarnings(as.numeric(sapply(clin, function(x) x$OS_MONTHS %||% NA)))
dead <- grepl("DECEASED|1", toupper(os_status))
keep <- !is.na(os_months) & os_months > 0 & !is.na(dead) & (os_status != "")
kclin <- data.frame(patient = names(clin)[keep],
                    time = os_months[keep],
                    status = as.numeric(dead[keep]),
                    stringsAsFactors = FALSE)

# ---------- KIRC 表达：Xena EBpp（log2(RSEM+1) 原始格式） ----------
kexpr <- fromJSON(file.path(P2, "data/KIRC_expr.json"))
sid <- names(kexpr); pid <- substr(sid, 1, 12)
common_p <- intersect(pid, kclin$patient)
s_for_p <- sid[match(common_p, pid)]
kclin <- kclin[match(common_p, kclin$patient), ]

xf <- file.path(P2, "data/pancan/EBpp_geneExp.xena.gz")
hdr <- names(fread(cmd = paste("gzip -dc", shQuote(xf), "| head -1"), header = TRUE))
keep_cols <- c(1, which(hdr %in% s_for_p))
X <- fread(cmd = paste("gzip -dc", shQuote(xf)), select = keep_cols, header = TRUE, showProgress = FALSE)
gn <- X[[1]]; X <- as.matrix(X[, -1]); rownames(X) <- gn
alias <- c(C19orf50 = "KXD1", KIAA0831 = "ATG14")
rownames(X)[rownames(X) %in% names(alias)] <- alias[rownames(X)[rownames(X) %in% names(alias)]]
stopifnot(all(TRIA21 %in% rownames(X)))
X <- X[, s_for_p, drop = FALSE]
cat("KIRC:", ncol(X), " samples, events =", sum(kclin$status), "\n")

X_log <- t(X[TRIA21, , drop = FALSE])            # log2(RSEM+1) 原尺度
X_lin <- t(2^X[TRIA21, , drop = FALSE] - 1)      # 线性 RSEM（手稿主分析）
X_lin[X_lin < 0 & X_lin > -1e-6] <- 0

# ---------- 临床协变量（grade/stage，与 13_multivariable_cox 同口径） ----------
grade_f <- file.path(P2, "data/external_validation/kirc_tcga_grade.csv")
stage_f <- file.path(P2, "data/external_validation/kirc_tcga_pan_can_atlas_2018_clin_patient.csv")
grade <- if (file.exists(grade_f)) read.csv(grade_f, stringsAsFactors = FALSE) else NULL
pat  <- if (file.exists(stage_f)) read.csv(stage_f, stringsAsFactors = FALSE) else NULL
get_stage <- function(patient) {
  if (is.null(pat)) return(NA)
  r <- pat[pat$Patient.ID %in% patient | pat[,1] %in% patient, ]
  if (!nrow(r)) return(NA)
  st <- r$AJCC.Pathologic.Tumor.Stage[1] %||% NA
  as.integer(gsub("\\D", "", substr(st %||% "", 6, 6)))
}

# ---------- 单管线：尺度 -> 结果 ----------
run_pipeline <- function(Xs, scale_label) {
  Z <- scale(Xs); Z[!is.finite(Z)] <- 0
  cvf <- cv.glmnet(Z, Surv(kclin$time, kclin$status), family = "cox",
                   alpha = 0.5, nfolds = 10)
  beta <- as.numeric(coef(cvf, s = "lambda.min")); names(beta) <- colnames(Z)
  sel <- names(beta)[beta != 0]
  score <- as.numeric(Z %*% beta)

  d <- data.frame(time = kclin$time, status = kclin$status, score = score)
  d$score_sd <- as.numeric(scale(d$score))
  # 单变量
  f1 <- coxph(Surv(time, status) ~ score_sd, data = d)
  # apparent C
  c_app <- as.numeric(concordance(Surv(time, status) ~ score, data = d, reverse = TRUE)$concordance)
  # KM 中位分组
  d$risk <- factor(ifelse(d$score > median(d$score), "High", "Low"), levels = c("Low","High"))
  km_hr <- exp(coef(coxph(Surv(time, status) ~ risk, data = d))["riskHigh"])

  # 多变量（年龄/性别/分期/分级；列名以 CSV 实际为准）
  hr_adj <- NA
  if (!is.null(pat) && !is.null(grade)) {
    m <- match(common_p, pat$patientId)
    age <- suppressWarnings(as.numeric(pat$AGE[m]))
    sex <- as.character(pat$SEX[m])
    roman <- c("I"=1L,"II"=2L,"III"=3L,"IV"=4L)
    stage_str <- trimws(gsub("STAGE", "", toupper(pat$AJCC_PATHOLOGIC_TUMOR_STAGE[m])))
    stage <- suppressWarnings(as.integer(roman[stage_str]))
    gm <- match(common_p, grade$patientId)
    gr <- suppressWarnings(as.integer(gsub("^G", "", toupper(grade$GRADE[gm]))))
    d2 <- data.frame(d, age = age, sex = sex, stage = stage, grade = gr)
    d2$grade <- pmax(d2$grade, 2)  # grade 1-2 合并（13 脚本口径）
    d2 <- d2[!is.na(d2$stage) & !is.na(d2$age) & !is.na(d2$grade), ]
    f2 <- tryCatch(coxph(Surv(time, status) ~ score_sd + age + sex + stage + grade, data = d2),
                   error = function(e) NULL)
    if (!is.null(f2)) hr_adj <- exp(coef(f2)["score_sd"])
  }

  # OOF 5 折（折内重选 lambda；concordance 口径，与手稿 IPCW 口径不同处已在报告注明）
  K <- 5; folds <- sample(rep(1:K, length.out = nrow(Z)))
  pred <- rep(NA_real_, nrow(Z))
  for (k in 1:K) {
    tr <- folds != k; te <- folds == k
    cvk <- cv.glmnet(Z[tr, , drop = FALSE], Surv(kclin$time[tr], kclin$status[tr]),
                     family = "cox", alpha = 0.5, nfolds = 10)
    pred[te] <- as.numeric(predict(cvk, newx = Z[te, , drop = FALSE], s = "lambda.min"))
  }
  c_oof <- as.numeric(concordance(Surv(kclin$time, kclin$status) ~ pred, reverse = TRUE)$concordance)

  list(label = scale_label, beta = beta, sel = sel, score = score,
       hr_uni = exp(coef(f1)["score_sd"]), c_app = c_app, km_hr = km_hr,
       hr_adj = hr_adj, c_oof = c_oof)
}

cat("\n>>> linear (手稿主尺度)\n"); r_lin <- run_pipeline(X_lin, "linear")
cat("   选中", length(r_lin$sel), "基因:", paste(r_lin$sel, collapse = ","), "\n")
cat(">>> log2(RSEM+1)\n");          r_log <- run_pipeline(X_log, "log")
cat("   选中", length(r_log$sel), "基因:", paste(r_log$sel, collapse = ","), "\n")

# ---------- 外部验证：冻结系数 -> E-MTAB-1980 / CPTAC-3 ----------
MIN_GENES <- 8
score_frozen <- function(X, beta, sel) {
  g <- intersect(sel, colnames(X))
  if (length(g) < MIN_GENES) return(NULL)
  Z <- scale(X[, g, drop = FALSE]); Z[!is.finite(Z)] <- 0
  as.numeric(Z %*% beta[g])
}
eval_os <- function(d, cohort) {
  d <- d[!is.na(d$os_t) & !is.na(d$os_s) & d$os_t > 0, ]
  cd <- survival::concordance(Surv(os_t, os_s) ~ TAS, data = d, reverse = TRUE)
  list(n = nrow(d), ev = sum(d$os_s), c = as.numeric(cd$concordance))
}

# E-MTAB-1980（log 尺度芯片；REFSEQ->SYMBOL 用 org.Hs.eg.db）
suppressPackageStartupMessages({ library(org.Hs.eg.db); library(AnnotationDbi) })
ext_em <- function(r) {
  expr_f <- file.path(MP, "data/external/E-MTAB-1980_processed/ccRCC_exp_log_quantile_normalized.txt")
  if (!file.exists(expr_f)) return(list(n = NA, ev = NA, c = NA))
  m <- read.delim(expr_f, stringsAsFactors = FALSE, check.names = FALSE)
  uniq_ref <- unique(m$SystematicName)
  ann <- AnnotationDbi::select(org.Hs.eg.db, keys = uniq_ref, columns = "SYMBOL", keytype = "REFSEQ")
  ann <- ann[!is.na(ann$SYMBOL) & ann$SYMBOL != "", ]
  sym <- toupper(setNames(ann$SYMBOL, ann$REFSEQ)[m$SystematicName])
  hit <- which(sym %in% r$sel)
  Xp <- as.matrix(m[hit, -(1:3), drop = FALSE]); storage.mode(Xp) <- "numeric"
  spl <- split(seq_len(nrow(Xp)), sym[hit])
  X <- t(vapply(spl, function(ii) colMeans(Xp[ii, , drop = FALSE], na.rm = TRUE),
                FUN.VALUE = numeric(ncol(Xp))))
  X <- t(X); colnames(X) <- names(spl); rownames(X) <- colnames(m)[-(1:3)]
  sc <- score_frozen(X, r$beta, r$sel)
  if (is.null(sc)) return(list(n = NA, ev = NA, c = NA))
  sc_em <- data.frame(sample = rownames(X), TAS = sc)
  sc_em$patientId <- gsub("-", "_", sc_em$sample)   # 与 12_external_os_validation.R 同口径
  sato <- read.csv(file.path(MP, "data/external/cbioportal/ccrcc_utokyo_2013_patient_clin_full.csv"),
                   stringsAsFactors = FALSE)
  d <- merge(sc_em, sato, by = "patientId", all.x = TRUE)
  d$os_t <- as.numeric(d$FOLLOW_UP_TIME_MONTHS)
  d$os_s <- as.integer(toupper(d$VITAL_STATUS) == "DECEASED")
  eval_os(d, "E-MTAB-1980")
}
ext_cp <- function(r) {
  expr <- read.csv(file.path(P2, "data/external_validation/cptac_tas_expr21.csv"), stringsAsFactors = FALSE)
  expr$patientId <- sub("-[0-9]+$", "", expr$sampleId)
  ord <- ifelse(grepl("-01$", expr$sampleId), 0, 1)
  expr <- expr[order(expr$patientId, ord, expr$sampleId), ]
  expr <- expr[!duplicated(expr$patientId), ]
  genes21 <- colnames(expr)[-1]
  X <- as.matrix(expr[, genes21]); storage.mode(X) <- "numeric"
  X[is.na(X)] <- 0   # 与 12 脚本同口径：z=0 (队列均值)
  rownames(X) <- expr$patientId
  sc <- score_frozen(X, r$beta, r$sel)
  if (is.null(sc)) return(list(n = NA, ev = NA, c = NA))
  surv <- read.delim(file.path(MP, "data/external/cptac/CCRCC_survival.txt"), stringsAsFactors = FALSE)
  d <- merge(data.frame(patientId = rownames(X), TAS = sc), surv, by.x = "patientId", by.y = "case_id")
  d$os_t <- as.numeric(d$OS_days) / 30.4375; d$os_s <- as.integer(d$OS_event)
  eval_os(d, "CPTAC-3")
}

cat("\n>>> 外部验证（linear 冻结系数）\n"); e_lin_em <- ext_em(r_lin); e_lin_cp <- ext_cp(r_lin)
cat(">>> 外部验证（log 冻结系数）\n");    e_log_em <- ext_em(r_log); e_log_cp <- ext_cp(r_log)

# ---------- 汇总对比 ----------
summ <- data.frame(
  metric = c("n_genes_selected", "gene_overlap_with_linear", "HR_perSD_univariate",
             "HR_perSD_adjusted(age/sex/stage/grade)", "KM_HR_median_split",
             "C_apparent", "C_OOF_5fold(concordance)",
             "C_EMTAB1980", "C_CPTAC3"),
  linear = c(length(r_lin$sel), NA, r_lin$hr_uni, r_lin$hr_adj, r_lin$km_hr,
             r_lin$c_app, r_lin$c_oof, e_lin_em$c, e_lin_cp$c),
  log    = c(length(r_log$sel), length(intersect(r_log$sel, r_lin$sel)),
             r_log$hr_uni, r_log$hr_adj, r_log$km_hr,
             r_log$c_app, r_log$c_oof, e_log_em$c, e_log_cp$c)
)
print(summ, digits = 4)
fwrite(summ, file.path(OUT, "scale_sensitivity_summary.csv"), row.names = FALSE)

genes_tab <- data.frame(
  gene = TRIA21,
  in_linear = TRIA21 %in% r_lin$sel, beta_linear = as.numeric(r_lin$beta[TRIA21]),
  in_log = TRIA21 %in% r_log$sel, beta_log = as.numeric(r_log$beta[TRIA21])
)
fwrite(genes_tab, file.path(OUT, "scale_sensitivity_genes.csv"), row.names = FALSE)
cat("\n输出:", OUT, "\n")
cat("判读: 若选中基因高度重合且 HR/C 差异小, 则结论对表达尺度稳健;\n")
cat("若明显变动, 需在手稿中如实报告两种尺度的结果区间.\n")
cat("注: OOF 此处为 concordance 口径; 手稿主分析为 IPCW Uno 口径, 两者数值不直接可比,\n")
cat("    但 linear vs log 在同口径下的差值是有意义的.\n")
cat("=== done ===\n")
