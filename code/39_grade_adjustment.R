#!/usr/bin/env Rscript
# Does TAS add to stage AND histologic grade?
# Grade from cBioPortal PanCanAtlas (code/39a_fetch_kirc_grade.py). G1 (n=13) is pooled with G2;
# GX is missing. Complete cases for age, sex, stage, and grade.
# Same machinery as code/31_incremental_cindex.R: apparent Harrell C, LRT, 1,000-resample bootstrap
# of the C difference, separate 5-fold out-of-fold Harrell C (seed 20260923), Schoenfeld tests,
# and a stage-stratified model. TAS is per SD of the full TCGA-KIRC cohort.

suppressPackageStartupMessages(library(survival))
P2  <- "/Users/xiaolin/WorkBuddy/2026-09-04-20-00-08/paper2_triaptosis"
OUT <- file.path(P2, "results", "08_multivariable")
SEED <- 20260923; B <- 1000; K <- 5

score <- read.csv(file.path(P2, "results/02_model/KIRC/score_table.csv"), row.names = 1)
score$patientId <- substr(rownames(score), 1, 12)
score$TAS <- as.numeric(scale(score$TAS))
clin <- read.csv(file.path(P2, "data/external_validation/kirc_tcga_pan_can_atlas_2018_clin_patient.csv"))
st <- gsub("^STAGE\\s*", "", toupper(trimws(clin$AJCC_PATHOLOGIC_TUMOR_STAGE)))
st[!st %in% c("I", "II", "III", "IV")] <- NA
clin$STAGE <- factor(st, levels = c("I", "II", "III", "IV"))
clin$AGE <- as.numeric(clin$AGE)
clin$SEX <- factor(toupper(trimws(clin$SEX)))
gr <- read.csv(file.path(P2, "data/external_validation/kirc_tcga_grade.csv"))
gr$GRADE <- factor(ifelse(gr$GRADE %in% c("G1", "G2"), "G1-2",
                   ifelse(gr$GRADE %in% c("G3", "G4"), gr$GRADE, NA)), levels = c("G1-2", "G3", "G4"))

d0 <- merge(score[, c("patientId", "time", "status", "TAS", "group")],
            clin[, c("patientId", "AGE", "SEX", "STAGE")], by = "patientId")
d0 <- merge(d0, gr[, c("patientId", "GRADE")], by = "patientId", all.x = TRUE)
d0 <- d0[complete.cases(d0[, c("time", "status", "TAS", "AGE", "SEX", "STAGE")]) & d0$time > 0, ]
d <- d0[!is.na(d0$GRADE), ]
cat(sprintf("stage-complete %d (%d deaths); with grade %d (%d deaths); grade missing %d\n",
            nrow(d0), sum(d0$status), nrow(d), sum(d$status), sum(is.na(d0$GRADE))))
print(table(d$GRADE))
ct <- cor.test(d$TAS, as.numeric(d$GRADE), method = "spearman", exact = FALSE)
cat(sprintf("TAS vs grade Spearman rho %.3f, P %.2g\n", ct$estimate, ct$p.value))

forms <- list(
  stage_grade = Surv(time, status) ~ STAGE + GRADE,
  stage_grade_tas = Surv(time, status) ~ TAS + STAGE + GRADE,
  clinical = Surv(time, status) ~ AGE + SEX + STAGE + GRADE,
  full = Surv(time, status) ~ TAS + AGE + SEX + STAGE + GRADE
)
fits <- lapply(forms, function(f) coxph(f, data = d, ties = "efron"))
hc <- function(f) unname(summary(f)$concordance["C"])
app <- data.frame(model = names(fits), harrell_C = sapply(fits, hc), row.names = NULL)
lrt <- function(a, b) { x <- anova(fits[[a]], fits[[b]])
  data.frame(contrast = paste(b, "vs", a), chisq = x$Chisq[2], df = x$Df[2], p = x$`Pr(>|Chi|)`[2]) }
lrts <- rbind(lrt("stage_grade", "stage_grade_tas"), lrt("clinical", "full"))
# Converse question: does grade still add once TAS is in the model?
g0 <- anova(coxph(Surv(time, status) ~ AGE + SEX + STAGE, data = d, ties = "efron"), fits$clinical)
g1 <- anova(coxph(Surv(time, status) ~ TAS + AGE + SEX + STAGE, data = d, ties = "efron"), fits$full)
lrts <- rbind(lrts,
  data.frame(contrast = "grade added to age + sex + stage", chisq = g0$Chisq[2], df = g0$Df[2], p = g0$`Pr(>|Chi|)`[2]),
  data.frame(contrast = "grade added to TAS + age + sex + stage", chisq = g1$Chisq[2], df = g1$Df[2], p = g1$`Pr(>|Chi|)`[2]))

s <- summary(fits$full)
hr <- data.frame(term = rownames(s$conf.int), HR = s$conf.int[, 1], lo = s$conf.int[, 3],
                 hi = s$conf.int[, 4], p = s$coefficients[, 5], row.names = NULL)

set.seed(SEED)
boot <- t(replicate(B, {
  db <- d[sample(nrow(d), replace = TRUE), ]
  cc <- function(f) unname(summary(coxph(f, data = db, ties = "efron"))$concordance["C"])
  c(cc(forms$stage_grade_tas) - cc(forms$stage_grade), cc(forms$full) - cc(forms$clinical))
}))
bsum <- data.frame(contrast = c("TAS added to stage + grade", "TAS added to age + sex + stage + grade"),
                   delta_C = c(app$harrell_C[2] - app$harrell_C[1], app$harrell_C[4] - app$harrell_C[3]),
                   lo = apply(boot, 2, quantile, 0.025), hi = apply(boot, 2, quantile, 0.975))

set.seed(SEED)
folds <- sample(rep(seq_len(K), length.out = nrow(d)))
oof <- function(form) {
  lp <- rep(NA_real_, nrow(d))
  for (k in seq_len(K)) {
    fit <- coxph(form, data = d[folds != k, ], ties = "efron")
    lp[folds == k] <- predict(fit, newdata = d[folds == k, ], type = "lp")
  }
  unname(concordance(Surv(time, status) ~ lp, data = d, reverse = TRUE)$concordance)
}
app$oof_C <- sapply(forms, oof)

zph <- cox.zph(coxph(Surv(time, status) ~ TAS + AGE + SEX + STAGE + GRADE, data = d, ties = "efron"),
               transform = "km")
zph_df <- data.frame(term = rownames(zph$table), zph$table, row.names = NULL)
f_str0 <- coxph(Surv(time, status) ~ AGE + SEX + GRADE + strata(STAGE), data = d, ties = "efron")
f_str <- coxph(Surv(time, status) ~ TAS + AGE + SEX + GRADE + strata(STAGE), data = d, ties = "efron")
ss <- summary(f_str)
zs <- cox.zph(coxph(Surv(time, status) ~ TAS + AGE + SEX + GRADE + strata(STAGE), data = d, ties = "efron"),
              transform = "km")
strat <- data.frame(model = "TAS + age + sex + grade, stratified by stage", n = f_str$n, events = f_str$nevent,
                    HR = ss$conf.int["TAS", 1], lo = ss$conf.int["TAS", 3], hi = ss$conf.int["TAS", 4],
                    p = ss$coefficients["TAS", 5], lrt_p = anova(f_str0, f_str)$`Pr(>|Chi|)`[2],
                    ph_TAS_p = zs$table["TAS", "p"], ph_global_p = zs$table["GLOBAL", "p"])

sub <- do.call(rbind, lapply(list(`G1-2` = "G1-2", `G3-4` = c("G3", "G4")), function(g) {
  dd <- d[d$GRADE %in% g, ]
  dd$group <- factor(dd$group, c("Low", "High"))
  fb <- summary(coxph(Surv(time, status) ~ group, data = dd))
  fc <- summary(coxph(Surv(time, status) ~ TAS + AGE + SEX + STAGE, data = dd))
  data.frame(grade = paste(g, collapse = "/"), n = nrow(dd), events = sum(dd$status),
             HR_high_vs_low = fb$conf.int[1, 1], lo_b = fb$conf.int[1, 3], hi_b = fb$conf.int[1, 4],
             km_p = 1 - pchisq(survdiff(Surv(time, status) ~ group, data = dd)$chisq, 1),
             HR_perSD_adj = fc$conf.int["TAS", 1], lo_adj = fc$conf.int["TAS", 3],
             hi_adj = fc$conf.int["TAS", 4], p_adj = fc$coefficients["TAS", 5])
}))

print(app, digits = 4); print(lrts, digits = 4); print(bsum, digits = 3); print(hr, digits = 3)
print(zph_df, digits = 3); print(strat, digits = 3); print(sub, digits = 3)
write.csv(app, file.path(OUT, "grade_harrell_c.csv"), row.names = FALSE)
write.csv(lrts, file.path(OUT, "grade_lrt.csv"), row.names = FALSE)
write.csv(bsum, file.path(OUT, "grade_bootstrap_delta_c.csv"), row.names = FALSE)
write.csv(hr, file.path(OUT, "grade_multivariable_cox.csv"), row.names = FALSE)
write.csv(zph_df, file.path(OUT, "grade_ph_schoenfeld.csv"), row.names = FALSE)
write.csv(strat, file.path(OUT, "grade_stage_stratified.csv"), row.names = FALSE)
write.csv(sub, file.path(OUT, "grade_subgroups.csv"), row.names = FALSE)
