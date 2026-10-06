#!/usr/bin/env Rscript
# Incremental discrimination of TAS over stage and clinical covariates.
# Apparent Harrell C, likelihood-ratio tests, Schoenfeld PH check,
# bootstrap CI for the C difference, and a separate 5-fold out-of-fold Harrell C.
# Does not reuse the TAS-only IPCW folds from 07_model_validation.R.

suppressPackageStartupMessages(library(survival))

P2  <- "/Users/xiaolin/WorkBuddy/2026-09-04-20-00-08/paper2_triaptosis"
OUT <- file.path(P2, "results", "08_multivariable")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

score <- read.csv(file.path(P2, "results/02_model/KIRC/score_table.csv"),
                  stringsAsFactors = FALSE, row.names = 1)
score$patientId <- substr(rownames(score), 1, 12)
# Continuous TAS HRs are per SD of the full TCGA-KIRC cohort.
score$TAS <- as.numeric(scale(score$TAS))

clin <- read.csv(file.path(P2, "data/external_validation/kirc_tcga_pan_can_atlas_2018_clin_patient.csv"),
                 stringsAsFactors = FALSE)
st <- toupper(trimws(clin$AJCC_PATHOLOGIC_TUMOR_STAGE))
st <- gsub("^STAGE\\s*", "", st)
st[!st %in% c("I", "II", "III", "IV")] <- NA
clin$STAGE <- factor(st, levels = c("I", "II", "III", "IV"))
clin$AGE <- as.numeric(clin$AGE)
clin$SEX <- factor(toupper(trimws(clin$SEX)))

d <- merge(score[, c("patientId", "time", "status", "TAS", "group")],
           clin[, c("patientId", "AGE", "SEX", "STAGE")],
           by = "patientId")
d <- d[complete.cases(d[, c("time", "status", "TAS", "AGE", "SEX", "STAGE")]), ]
d <- d[d$time > 0, ]
d$group <- factor(d$group, levels = c("Low", "High"))
cat("complete cases:", nrow(d), " events:", sum(d$status), "\n")

# Binary high vs low on this complete-case set, and on the full score table.
bin_hr <- function(dd) {
  dd$group <- factor(dd$group, levels = c("Low", "High"))
  fx <- coxph(Surv(time, status) ~ group, data = dd, ties = "efron")
  ci <- exp(confint(fx))
  lr <- survdiff(Surv(time, status) ~ group, data = dd)
  data.frame(
    n = nrow(dd), events = sum(dd$status),
    HR = unname(exp(coef(fx))),
    lo = unname(ci[1]), hi = unname(ci[2]),
    cox_p = summary(fx)$coefficients[1, "Pr(>|z|)"],
    km_p = 1 - pchisq(lr$chisq, df = 1),
    C = summary(fx)$concordance["C"]
  )
}
score_bin <- score[score$time > 0 & !is.na(score$status) & !is.na(score$group), ]
bin <- rbind(
  full_score_table = bin_hr(score_bin),
  complete_case = bin_hr(d)
)
print(bin)
write.csv(bin, file.path(OUT, "binary_hr_audit.csv"))

forms <- list(
  stage = Surv(time, status) ~ STAGE,
  clinical = Surv(time, status) ~ AGE + SEX + STAGE,
  tas = Surv(time, status) ~ TAS,
  tas_stage = Surv(time, status) ~ TAS + STAGE,
  full = Surv(time, status) ~ TAS + AGE + SEX + STAGE
)
fits <- lapply(forms, function(f) coxph(f, data = d, ties = "efron"))

c_row <- function(name, fit) {
  sm <- summary(fit)
  data.frame(
    model = name,
    n = fit$n,
    events = fit$nevent,
    harrell_C = unname(sm$concordance["C"]),
    C_se = unname(sm$concordance["se(C)"]),
    loglik = fit$loglik[2],
    df = length(coef(fit))
  )
}
tab <- do.call(rbind, Map(c_row, names(fits), fits))
rownames(tab) <- NULL

lrt_pair <- function(reduced, full, label) {
  a <- anova(fits[[reduced]], fits[[full]])
  data.frame(
    contrast = label,
    df = a$Df[2],
    chisq = a$Chisq[2],
    p = a$`Pr(>|Chi|)`[2],
    delta_apparent_C = tab$harrell_C[tab$model == full] - tab$harrell_C[tab$model == reduced]
  )
}
lrt <- rbind(
  lrt_pair("stage", "tas_stage", "TAS added to stage"),
  lrt_pair("clinical", "full", "TAS added to age, sex, and stage"),
  lrt_pair("tas", "full", "age, sex, and stage added to TAS")
)

# Schoenfeld residuals on the full model.
zph <- cox.zph(fits$full, transform = "km")
zph_df <- data.frame(term = rownames(zph$table), zph$table, row.names = NULL)
colnames(zph_df) <- c("term", "chisq", "df", "p")

# Bootstrap percentile CI for apparent Harrell C difference. Patients resampled.
set.seed(20260923)
B <- 1000
boot_delta <- matrix(NA_real_, B, 2, dimnames = list(NULL, c("stage", "clinical")))
for (b in seq_len(B)) {
  ii <- sample.int(nrow(d), replace = TRUE)
  db <- d[ii, ]
  if (length(unique(db$STAGE)) < 4 || sum(db$status) < 30) next
  ok <- tryCatch({
    fs <- coxph(forms$stage, data = db, ties = "efron")
    ft <- coxph(forms$tas_stage, data = db, ties = "efron")
    fc <- coxph(forms$clinical, data = db, ties = "efron")
    ff <- coxph(forms$full, data = db, ties = "efron")
    boot_delta[b, 1] <- summary(ft)$concordance["C"] - summary(fs)$concordance["C"]
    boot_delta[b, 2] <- summary(ff)$concordance["C"] - summary(fc)$concordance["C"]
    TRUE
  }, error = function(e) FALSE)
}
boot_sum <- function(x, name) {
  x <- x[is.finite(x)]
  data.frame(
    contrast = name,
    n_boot = length(x),
    mean_delta_C = mean(x),
    lo = unname(quantile(x, 0.025)),
    hi = unname(quantile(x, 0.975))
  )
}
boot <- rbind(
  boot_sum(boot_delta[, 1], "TAS added to stage"),
  boot_sum(boot_delta[, 2], "TAS added to age, sex, and stage")
)

# 5-fold out-of-fold Harrell C. New seed; not the IPCW folds.
set.seed(20260923)
K <- 5
folds <- sample(rep(seq_len(K), length.out = nrow(d)))
oof <- function(form) {
  lp <- rep(NA_real_, nrow(d))
  for (k in seq_len(K)) {
    tr <- folds != k
    te <- folds == k
    fit <- coxph(form, data = d[tr, ], ties = "efron")
    lp[te] <- predict(fit, newdata = d[te, ], type = "lp")
  }
  conc <- concordance(Surv(time, status) ~ lp, data = d, reverse = TRUE)
  data.frame(C = conc$concordance, se = sqrt(conc$var))
}
oof_tab <- do.call(rbind, lapply(names(forms), function(nm) {
  r <- oof(forms[[nm]])
  data.frame(model = nm, oof_harrell_C = r$C, oof_se = r$se)
}))
oof_delta <- data.frame(
  contrast = c("TAS added to stage", "TAS added to age, sex, and stage"),
  delta_oof_C = c(
    oof_tab$oof_harrell_C[oof_tab$model == "tas_stage"] - oof_tab$oof_harrell_C[oof_tab$model == "stage"],
    oof_tab$oof_harrell_C[oof_tab$model == "full"] - oof_tab$oof_harrell_C[oof_tab$model == "clinical"]
  )
)

cat("\n--- Harrell C ---\n"); print(tab, digits = 4)
cat("\n--- LRT ---\n"); print(lrt, digits = 4)
cat("\n--- PH ---\n"); print(zph_df, digits = 4)
cat("\n--- bootstrap delta C ---\n"); print(boot, digits = 4)
cat("\n--- OOF C ---\n"); print(oof_tab, digits = 4)
cat("\n--- OOF delta ---\n"); print(oof_delta, digits = 4)

write.csv(tab, file.path(OUT, "incremental_harrell_c.csv"), row.names = FALSE)
write.csv(lrt, file.path(OUT, "incremental_lrt.csv"), row.names = FALSE)
write.csv(zph_df, file.path(OUT, "ph_schoenfeld_full.csv"), row.names = FALSE)
write.csv(boot, file.path(OUT, "incremental_bootstrap_delta_c.csv"), row.names = FALSE)
write.csv(oof_tab, file.path(OUT, "incremental_oof_harrell_c.csv"), row.names = FALSE)
write.csv(oof_delta, file.path(OUT, "incremental_oof_delta_c.csv"), row.names = FALSE)

# Stage violates proportional hazards, so also fit TAS with stage as a stratum.
fit_str0 <- coxph(Surv(time, status) ~ AGE + SEX + strata(STAGE), data = d, ties = "efron")
fit_str1 <- coxph(Surv(time, status) ~ TAS + AGE + SEX + strata(STAGE), data = d, ties = "efron")
sm <- summary(fit_str1)
strat <- data.frame(
  model = "TAS + age + sex, stratified by stage",
  n = fit_str1$n,
  events = fit_str1$nevent,
  HR = sm$conf.int["TAS", "exp(coef)"],
  lo = sm$conf.int["TAS", "lower .95"],
  hi = sm$conf.int["TAS", "upper .95"],
  p = sm$coefficients["TAS", "Pr(>|z|)"],
  harrell_C = unname(sm$concordance["C"]),
  lrt_chisq = anova(fit_str0, fit_str1)$Chisq[2],
  lrt_df = anova(fit_str0, fit_str1)$Df[2],
  lrt_p = anova(fit_str0, fit_str1)$`Pr(>|Chi|)`[2]
)
zph_s <- cox.zph(fit_str1, transform = "km")
zph_s_df <- data.frame(term = rownames(zph_s$table), zph_s$table, row.names = NULL)
colnames(zph_s_df) <- c("term", "chisq", "df", "p")
cat("\n--- stratified ---\n"); print(strat, digits = 4)
cat("\n--- stratified PH ---\n"); print(zph_s_df, digits = 4)
write.csv(strat, file.path(OUT, "tas_stratified_by_stage.csv"), row.names = FALSE)
write.csv(zph_s_df, file.path(OUT, "ph_schoenfeld_stratified.csv"), row.names = FALSE)
cat("\n=== done ===\n")
