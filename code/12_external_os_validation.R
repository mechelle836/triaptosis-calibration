#!/usr/bin/env Rscript
# =============================================================================
# 12_external_os_validation.R — Paper2-Triaptosis 独立 OS 队列验证 (2026-09-02)
#
# 冻结 TAS (results/02_model/KIRC/model.rds 的弹性网 beta, 不重训练) 打分两个
# 独立 ccRCC OS 队列 (均与 TCGA-KIRC 训练集无关):
#   A. E-MTAB-1980 (Sato 2013, n=101, 23 events) — ArrayExpress 微阵列
#   B. CPTAC-3 ccRCC (LinkedOmics OS, n~100, ~20 events) — cBioPortal TPM
#
# 协议 (与 03_focus_model.R 平行癌种验证、MitoScore 42_external_os_validation.R 一致):
#   - 队列内 z-score 标准化 10 个选中基因 -> 乘冻结 beta
#   - 命中基因 < 8/10 则中止 (覆盖度保护)
#   - C-index: survival::concordance(reverse=TRUE) (高 TAS = 高风险)
#   - 分层: 队列内中位数 High/Low; KM log-rank; Cox HR per 1-SD (连续)
#   - 事件数少 (~20), 结论表述为 "方向一致性/探索性", 不称 multi-cohort confirmation
# =============================================================================
suppressPackageStartupMessages({
  library(survival); library(org.Hs.eg.db); library(AnnotationDbi)
})
P2  <- "/Users/xiaolin/WorkBuddy/2026-09-04-20-00-08/paper2_triaptosis"
MP  <- "/Users/xiaolin/CodeBuddy/Claw/pan-cancer_disulfidptosis_project"
OUT <- file.path(P2, "results", "07_external_os")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ---------- 冻结模型 ----------
mod <- readRDS(file.path(P2, "results/02_model/KIRC/model.rds"))
beta <- mod$beta
names(beta) <- mod$genes  # coef(cvf) 经 as.numeric 后丢名; genes 顺序 = glmnet 输入列序
stopifnot(!any(is.na(beta)))
sel <- mod$genes_sel
cat("冻结 TAS: ", length(sel), "基因:", paste(sel, collapse = ","), "\n")
MIN_GENES <- 8

score_frozen <- function(X) {  # X: samples x genes (仅 21 基因列或其子集)
  g <- intersect(sel, colnames(X))
  if (length(g) < MIN_GENES) stop("基因覆盖不足: ", length(g), "/", length(sel))
  Z <- scale(X[, g, drop = FALSE]); Z[!is.finite(Z)] <- 0
  as.numeric(Z %*% beta[g])
}

eval_os <- function(d, cohort, n_genes_used) {
  d <- d[!is.na(d$os_t) & !is.na(d$os_s) & d$os_t > 0, ]
  n <- nrow(d); ev <- sum(d$os_s)
  cd <- survival::concordance(Surv(os_t, os_s) ~ TAS, data = d, reverse = TRUE)
  cidx <- as.numeric(cd$concordance); cse <- sqrt(as.numeric(cd$var))
  d$TAS_sd <- as.numeric(scale(d$TAS))
  fit <- coxph(Surv(os_t, os_s) ~ TAS_sd, data = d)
  s <- summary(fit)
  hr <- exp(coef(fit)["TAS_sd"])
  hr_lo <- s$conf.int["TAS_sd", "lower .95"]
  hr_hi <- s$conf.int["TAS_sd", "upper .95"]
  d$risk <- ifelse(d$TAS >= median(d$TAS), "High", "Low")
  d$risk <- factor(d$risk, levels = c("Low", "High"))
  sd1 <- survdiff(Surv(os_t, os_s) ~ risk, data = d)
  km_p <- 1 - pchisq(sd1$chisq, df = max(1, length(sd1$n) - 1))
  list(d = d, tab = data.frame(
    cohort = cohort, n = n, events = ev,
    genes_used = n_genes_used,
    Cindex = cidx, Cindex_se = cse,
    Cindex_lo = cidx - 1.96 * cse, Cindex_hi = cidx + 1.96 * cse,
    cox_HR_perSD = hr, cox_HR_lo = hr_lo, cox_HR_hi = hr_hi, cox_p = s$waldtest["pvalue"], km_p = km_p))
}

# ---------- A. E-MTAB-1980 ----------
cat("\n=== A. E-MTAB-1980 (Sato 2013) ===\n")
expr_f <- file.path(MP, "data/external/E-MTAB-1980_processed/ccRCC_exp_log_quantile_normalized.txt")
m <- read.delim(expr_f, stringsAsFactors = FALSE, check.names = FALSE)
cat("  probe rows:", nrow(m), " samples:", ncol(m) - 3, "\n")
uniq_ref <- unique(m$SystematicName)
ann <- tryCatch(
  AnnotationDbi::select(org.Hs.eg.db, keys = uniq_ref, columns = "SYMBOL",
                        keytype = "REFSEQ"),
  error = function(e) NULL)
if (is.null(ann)) stop("REFSEQ->SYMBOL 注释失败")
ann <- ann[!is.na(ann$SYMBOL) & ann$SYMBOL != "", ]
map <- setNames(ann$SYMBOL, ann$REFSEQ)
sym <- toupper(map[m$SystematicName])
hit <- which(sym %in% sel)
cat("  probes mapping to TAS genes:", length(hit),
    " unique genes:", length(unique(sym[hit])), "\n")
Xp <- as.matrix(m[hit, -(1:3), drop = FALSE])
storage.mode(Xp) <- "numeric"
spl <- split(seq_len(nrow(Xp)), sym[hit])
X <- t(vapply(spl, function(ii) colMeans(Xp[ii, , drop = FALSE], na.rm = TRUE),
              FUN.VALUE = numeric(ncol(Xp))))
X <- t(X); colnames(X) <- names(spl); rownames(X) <- colnames(m)[-(1:3)]
tas_em <- score_frozen(X)
sc_em <- data.frame(sample = rownames(X), TAS = tas_em)
sc_em$patientId <- gsub("-", "_", sc_em$sample)

sato <- read.csv(file.path(MP, "data/external/cbioportal/ccrcc_utokyo_2013_patient_clin_full.csv"),
                 stringsAsFactors = FALSE)
d_em <- merge(sc_em, sato, by = "patientId", all.x = TRUE)
d_em$os_t <- as.numeric(d_em$FOLLOW_UP_TIME_MONTHS)
d_em$os_s <- as.integer(toupper(d_em$VITAL_STATUS) == "DECEASED")
res_em <- eval_os(d_em, "E-MTAB-1980", ncol(X))
print(res_em$tab)
write.csv(res_em$d[, c("sample", "patientId", "TAS", "os_t", "os_s", "risk")],
          file.path(OUT, "EMTAB1980_TAS_scored.csv"), row.names = FALSE)

# ---------- B. CPTAC-3 ----------
cat("\n=== B. CPTAC-3 (LinkedOmics OS) ===\n")
expr <- read.csv(file.path(P2, "data/external_validation/cptac_tas_expr21.csv"),
                 stringsAsFactors = FALSE)
# 一个患者一份切片: -01 优先 (42 号脚本口径)
expr$patientId <- sub("-[0-9]+$", "", expr$sampleId)
ord <- ifelse(grepl("-01$", expr$sampleId), 0, 1)
expr <- expr[order(expr$patientId, ord, expr$sampleId), ]
expr <- expr[!duplicated(expr$patientId), ]
genes21 <- colnames(expr)[-1]
X <- as.matrix(expr[, genes21]); storage.mode(X) <- "numeric"
n_na <- sum(is.na(X))
cat("  表达矩阵非数值单元格:", n_na, "/", length(X),
    sprintf("(%.2f%%), 按 z=0 (队列均值) 处理\n", 100 * n_na / length(X)))
rownames(X) <- expr$patientId
tas_cp <- score_frozen(X)
sc_cp <- data.frame(patientId = rownames(X), TAS = tas_cp)

surv <- read.delim(file.path(MP, "data/external/cptac/CCRCC_survival.txt"),
                   stringsAsFactors = FALSE)
d_cp <- merge(sc_cp, surv, by.x = "patientId", by.y = "case_id")
d_cp$os_t <- as.numeric(d_cp$OS_days) / 30.4375
d_cp$os_s <- as.integer(d_cp$OS_event)
res_cp <- eval_os(d_cp, "CPTAC-3", length(intersect(sel, genes21)))
print(res_cp$tab)
write.csv(res_cp$d[, c("patientId", "TAS", "os_t", "os_s", "risk")],
          file.path(OUT, "CPTAC3_TAS_scored.csv"), row.names = FALSE)

# ---------- 汇总 ----------
summ <- rbind(res_em$tab, res_cp$tab)
summ$genes_used <- NULL
write.csv(summ, file.path(OUT, "external_OS_TAS_summary.csv"), row.names = FALSE)
cat("\n", paste(capture.output(print(summ)), collapse = "\n"), "\n")

# ---------- KM 图 ----------
plot_km <- function(res, title, fname) {
  d <- res$d
  fit <- survfit(Surv(os_t, os_s) ~ risk, data = d)
  dat <- NULL
  for (s in names(fit$strata)) {
    st <- strsplit(s, "=")[[1]][2]; n <- fit$strata[[s]]
    tt <- fit$time[cumsum(fit$strata)[[s]] - n + seq_len(n)]
    dat <- rbind(dat, data.frame(time = tt, surv = fit$surv[cumsum(fit$strata)[[s]] - n + seq_len(n)],
                                 group = st))
  }
  p <- sprintf("%.3g", res$tab$km_p)
  hr <- sprintf("%.2f", res$tab$cox_HR_perSD)
  png(file.path(OUT, paste0(fname, ".png")), width = 1400, height = 1300, res = 200)
  par(mar = c(4.2, 4.2, 3, 1))
  cols <- c("High" = "#C0392B", "Low" = "#2471A3")
  plot(0, 0, type = "n", xlim = c(0, max(dat$time, na.rm = TRUE)),
       ylim = c(0, 1), xlab = "Time (months)", ylab = "Overall survival",
       main = title, cex.lab = 1.05)
  for (g in c("Low", "High")) {
    dd <- dat[dat$group == g, ]
    dd <- rbind(dd, dd[nrow(dd), ]); dd$time[nrow(dd)] <- max(dat$time)
    lines(dd$time, dd$surv, col = cols[[g]], lwd = 2.4)
  }
  legend("bottomleft", c(sprintf("High (n=%d)", sum(d$risk == "High")),
                         sprintf("Low (n=%d)", sum(d$risk == "Low"))),
         col = cols, lwd = 2.4, bty = "n", cex = 0.9)
  legend("topright", c(sprintf("log-rank P = %s", p),
                       sprintf("HR(per SD) = %s", hr),
                       sprintf("C = %.3f [%.3f, %.3f]", res$tab$Cindex,
                               res$tab$Cindex_lo, res$tab$Cindex_hi)),
         bty = "n", cex = 0.9)
  dev.off()
  cat("  figure:", file.path(OUT, paste0(fname, ".png")), "\n")
}
plot_km(res_em, sprintf("E-MTAB-1980 (Sato 2013)  n=%d, events=%d",
                        res_em$tab$n, res_em$tab$events), "Fig_ext_EMTAB1980_TAS_KM")
plot_km(res_cp, sprintf("CPTAC-3 ccRCC  n=%d, events=%d",
                        res_cp$tab$n, res_cp$tab$events), "Fig_ext_CPTAC3_TAS_KM")

note <- c(
  "Independent OS validation of frozen TAS (elastic-net beta from TCGA-KIRC, not retrained).",
  sprintf("E-MTAB-1980: n=%d events=%d C=%.3f [%.3f,%.3f] KM p=%.3g HR/perSD=%.3f (p=%.3g)",
          summ$n[1], summ$events[1], summ$Cindex[1], summ$Cindex_lo[1], summ$Cindex_hi[1],
          summ$km_p[1], summ$cox_HR_perSD[1], summ$cox_p[1]),
  sprintf("CPTAC-3: n=%d events=%d C=%.3f [%.3f,%.3f] KM p=%.3g HR/perSD=%.3f (p=%.3g)",
          summ$n[2], summ$events[2], summ$Cindex[2], summ$Cindex_lo[2], summ$Cindex_hi[2],
          summ$km_p[2], summ$cox_HR_perSD[2], summ$cox_p[2]),
  "Scoring: cohort-internal z-score of the 10 selected genes x frozen beta (same protocol as parallel-cancer validation).",
  "Event counts modest (~23/~20): report as directional/exploratory unless KM significant; do not overclaim multi-cohort confirmation.",
  "Reference (MitoScore same cohorts): E-MTAB-1980 C=0.666 KM p=0.093; CPTAC-3 C=0.749 KM p=0.0215."
)
writeLines(note, file.path(OUT, "external_OS_TAS_README.txt"))
cat(paste(note, collapse = "\n"), "\n>>> DONE", OUT, "\n")
