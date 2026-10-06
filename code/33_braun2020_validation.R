#!/usr/bin/env Rscript
# Third external OS cohort: Braun et al. Nat Med 2020 (CheckMate 009/010/025), frozen TAS.
# Previously treated advanced ccRCC receiving nivolumab or everolimus. OS_CNSR = 1 is death;
# CM-025 medians with this coding match the published trial.
#
# Analysis plan, fixed before fitting:
#   Primary   : same protocol as 12_external_os_validation.R (cohort z-score x frozen beta;
#               per-SD Cox HR; Harrell C with reverse=TRUE; median-split KM).
#   Secondary : Cox stratified by arm; plus MSKCC risk group; within-arm HR; TAS x arm interaction.
#   Sensitivity: primary-tumour samples only.

suppressPackageStartupMessages(library(survival))

P2  <- "/Users/xiaolin/WorkBuddy/2026-09-04-20-00-08/paper2_triaptosis"
D   <- file.path(P2, "data/external_validation/braun2020")
OUT <- file.path(P2, "results", "07_external_os")

mod <- readRDS(file.path(P2, "results/02_model/KIRC/model.rds"))
beta <- mod$beta
names(beta) <- mod$genes
sel <- mod$genes_sel
MIN_GENES <- 8

score_frozen <- function(X) {
  g <- intersect(sel, colnames(X))
  if (length(g) < MIN_GENES) stop("gene coverage ", length(g), "/", length(sel))
  Z <- scale(X[, g, drop = FALSE]); Z[!is.finite(Z)] <- 0
  list(score = as.numeric(Z %*% beta[g]), n_genes = length(g))
}

ex <- read.csv(file.path(D, "braun2020_tria21_expr.csv"), stringsAsFactors = FALSE, check.names = FALSE)
X <- t(as.matrix(ex[, -1])); colnames(X) <- ex$gene_name
storage.mode(X) <- "numeric"

cl <- read.csv(file.path(D, "braun2020_clinical.csv"), stringsAsFactors = FALSE,
               na.strings = c("NA", ""))
cl <- cl[!is.na(cl$RNA_ID), ]
cl$os_t <- as.numeric(cl$OS)
cl$os_s <- as.integer(cl$OS_CNSR == 1)
cl$Arm <- factor(cl$Arm, levels = c("EVEROLIMUS", "NIVOLUMAB"))
cl$MSKCC <- factor(cl$MSKCC, levels = c("FAVORABLE", "INTERMEDIATE", "POOR"))

sc <- score_frozen(X)
d0 <- data.frame(RNA_ID = rownames(X), TAS = sc$score)
d0 <- merge(d0, cl, by = "RNA_ID")
d0 <- d0[!is.na(d0$os_t) & !is.na(d0$os_s) & d0$os_t > 0, ]
cat("genes used", sc$n_genes, "| patients", nrow(d0), "| deaths", sum(d0$os_s), "\n")

eval_os <- function(d, label) {
  d$TASz <- as.numeric(scale(d$TAS))
  cd <- concordance(Surv(os_t, os_s) ~ TAS, data = d, reverse = TRUE)
  fit <- coxph(Surv(os_t, os_s) ~ TASz, data = d, ties = "efron")
  s <- summary(fit)
  d$risk <- factor(ifelse(d$TAS >= median(d$TAS), "High", "Low"), levels = c("Low", "High"))
  lr <- survdiff(Surv(os_t, os_s) ~ risk, data = d)
  fb <- coxph(Surv(os_t, os_s) ~ risk, data = d, ties = "efron")
  sb <- summary(fb)
  data.frame(
    analysis = label, n = nrow(d), events = sum(d$os_s),
    Cindex = cd$concordance, Cindex_lo = cd$concordance - 1.96 * sqrt(cd$var),
    Cindex_hi = cd$concordance + 1.96 * sqrt(cd$var),
    HR_perSD = s$conf.int["TASz", "exp(coef)"], HR_lo = s$conf.int["TASz", "lower .95"],
    HR_hi = s$conf.int["TASz", "upper .95"], cox_p = s$coefficients["TASz", "Pr(>|z|)"],
    HR_high_vs_low = sb$conf.int[1, "exp(coef)"], HRb_lo = sb$conf.int[1, "lower .95"],
    HRb_hi = sb$conf.int[1, "upper .95"], km_p = 1 - pchisq(lr$chisq, df = 1)
  )
}

d0$TASz <- as.numeric(scale(d0$TAS))
d0$risk <- factor(ifelse(d0$TAS >= median(d0$TAS), "High", "Low"), levels = c("Low", "High"))

primary <- eval_os(d0, "Primary: all RNA patients")
arm_rows <- do.call(rbind, lapply(levels(d0$Arm), function(a) eval_os(d0[d0$Arm == a, ], paste("Within arm:", a))))
prim_only <- eval_os(d0[d0$Tumor_Sample_Primary_or_Metastasis %in% "PRIMARY", ], "Sensitivity: primary-tumour samples")
tab <- rbind(primary, arm_rows, prim_only)
print(tab, digits = 3)

adj_row <- function(fit, label, term = "TASz") {
  s <- summary(fit)
  data.frame(model = label, n = fit$n, events = fit$nevent,
             HR = s$conf.int[term, "exp(coef)"], lo = s$conf.int[term, "lower .95"],
             hi = s$conf.int[term, "upper .95"], p = s$coefficients[term, "Pr(>|z|)"])
}
f_str <- coxph(Surv(os_t, os_s) ~ TASz + strata(Arm), data = d0, ties = "efron")
dm <- d0[!is.na(d0$MSKCC), ]
f_msk0 <- coxph(Surv(os_t, os_s) ~ MSKCC + strata(Arm), data = dm, ties = "efron")
f_msk <- coxph(Surv(os_t, os_s) ~ TASz + MSKCC + strata(Arm), data = dm, ties = "efron")
f_int0 <- coxph(Surv(os_t, os_s) ~ TASz + Arm, data = d0, ties = "efron")
f_int <- coxph(Surv(os_t, os_s) ~ TASz * Arm, data = d0, ties = "efron")
adj <- rbind(
  adj_row(f_str, "TAS, stratified by arm"),
  adj_row(f_msk, "TAS + MSKCC risk, stratified by arm")
)
a_int <- anova(f_int0, f_int)
s_int <- summary(f_int)
int_term <- "TASz:ArmNIVOLUMAB"
inter <- data.frame(
  term = "TAS x nivolumab (reference everolimus)",
  ratio_of_HR = s_int$conf.int[int_term, "exp(coef)"],
  lo = s_int$conf.int[int_term, "lower .95"],
  hi = s_int$conf.int[int_term, "upper .95"],
  lrt_chisq = a_int$Chisq[2], lrt_p = a_int$`Pr(>|Chi|)`[2]
)
a_msk <- anova(f_msk0, f_msk)
inc <- data.frame(
  contrast = "TAS added to MSKCC risk, stratified by arm",
  n = f_msk$n, events = f_msk$nevent,
  C_MSKCC = unname(summary(f_msk0)$concordance["C"]),
  C_MSKCC_TAS = unname(summary(f_msk)$concordance["C"]),
  lrt_chisq = a_msk$Chisq[2], lrt_p = a_msk$`Pr(>|Chi|)`[2]
)
inc$delta_C <- inc$C_MSKCC_TAS - inc$C_MSKCC
print(adj, digits = 3); print(inter, digits = 3); print(inc, digits = 3)

zph <- cox.zph(coxph(Surv(os_t, os_s) ~ TASz + MSKCC + strata(Arm), data = dm, ties = "efron"),
               transform = "km")
zph_df <- data.frame(term = rownames(zph$table), zph$table, row.names = NULL)
colnames(zph_df) <- c("term", "chisq", "df", "p")
print(zph_df, digits = 3)

write.csv(tab, file.path(OUT, "braun2020_TAS_summary.csv"), row.names = FALSE)
write.csv(adj, file.path(OUT, "braun2020_TAS_adjusted.csv"), row.names = FALSE)
write.csv(inter, file.path(OUT, "braun2020_TAS_arm_interaction.csv"), row.names = FALSE)
write.csv(inc, file.path(OUT, "braun2020_TAS_incremental_mskcc.csv"), row.names = FALSE)
write.csv(zph_df, file.path(OUT, "braun2020_TAS_ph.csv"), row.names = FALSE)
write.csv(d0[, c("SUBJID", "RNA_ID", "Cohort", "Arm", "MSKCC", "Tumor_Sample_Primary_or_Metastasis",
                 "TAS", "os_t", "os_s", "risk")],
          file.path(OUT, "braun2020_TAS_scored.csv"), row.names = FALSE)

fit_km <- survfit(Surv(os_t, os_s) ~ risk, data = d0)
png(file.path(OUT, "Fig_ext_Braun2020_TAS_KM.png"), width = 1400, height = 1300, res = 200)
par(mar = c(4.2, 4.2, 3, 1))
plot(fit_km, col = c("#3C5488", "#E64B35"), lwd = 2.4, xlab = "Time (months)",
     ylab = "Overall survival", main = sprintf("Braun 2020 CheckMate  n=%d, deaths=%d",
                                               primary$n, primary$events))
legend("bottomleft", c(sprintf("Low (n=%d)", sum(d0$risk == "Low")),
                       sprintf("High (n=%d)", sum(d0$risk == "High"))),
       col = c("#3C5488", "#E64B35"), lwd = 2.4, bty = "n", cex = 0.9)
legend("topright", c(sprintf("log-rank P = %.2g", primary$km_p),
                     sprintf("HR per SD = %.2f", primary$HR_perSD),
                     sprintf("C = %.3f [%.3f, %.3f]", primary$Cindex, primary$Cindex_lo, primary$Cindex_hi)),
       bty = "n", cex = 0.9)
dev.off()
cat("=== done ===\n")
