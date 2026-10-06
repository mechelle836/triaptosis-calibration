#!/usr/bin/env Rscript
# Does TAS stand in for VHL, PBRM1, SETD2, or BAP1 mutation?
# Restricted to patients in the cBioPortal sequenced list; unsequenced patients are not wild type.
# Inputs: score table, clinical table, data/external_validation/kirc_tcga_driver_mutations.csv
# (from code/32a_fetch_kirc_mutations.py). Output: results/08_multivariable/mutation_*.csv

suppressPackageStartupMessages(library(survival))

P2  <- "/Users/xiaolin/WorkBuddy/2026-09-04-20-00-08/paper2_triaptosis"
OUT <- file.path(P2, "results", "08_multivariable")
GENES <- c("VHL", "PBRM1", "SETD2", "BAP1")

score <- read.csv(file.path(P2, "results/02_model/KIRC/score_table.csv"),
                  stringsAsFactors = FALSE, row.names = 1)
score$patientId <- substr(rownames(score), 1, 12)
clin <- read.csv(file.path(P2, "data/external_validation/kirc_tcga_pan_can_atlas_2018_clin_patient.csv"),
                 stringsAsFactors = FALSE)
st <- gsub("^STAGE\\s*", "", toupper(trimws(clin$AJCC_PATHOLOGIC_TUMOR_STAGE)))
st[!st %in% c("I", "II", "III", "IV")] <- NA
clin$STAGE <- factor(st, levels = c("I", "II", "III", "IV"))
clin$AGE <- as.numeric(clin$AGE)
clin$SEX <- factor(toupper(trimws(clin$SEX)))
mut <- read.csv(file.path(P2, "data/external_validation/kirc_tcga_driver_mutations.csv"),
                stringsAsFactors = FALSE)

d_all <- merge(score[, c("patientId", "time", "status", "TAS")],
               clin[, c("patientId", "AGE", "SEX", "STAGE")], by = "patientId")
d_all <- d_all[complete.cases(d_all) & d_all$time > 0, ]
d_all$sequenced <- d_all$patientId %in% mut$patientId

# Selection check: are sequenced patients different from the rest?
sel <- data.frame(
  group = c("sequenced", "not_sequenced"),
  n = c(sum(d_all$sequenced), sum(!d_all$sequenced)),
  events = c(sum(d_all$status[d_all$sequenced]), sum(d_all$status[!d_all$sequenced])),
  median_TAS = c(median(d_all$TAS[d_all$sequenced]), median(d_all$TAS[!d_all$sequenced])),
  stage_III_IV = c(mean(d_all$STAGE[d_all$sequenced] %in% c("III", "IV")),
                   mean(d_all$STAGE[!d_all$sequenced] %in% c("III", "IV")))
)
sel$wilcox_TAS_p <- wilcox.test(TAS ~ sequenced, data = d_all)$p.value
print(sel)

d <- merge(d_all[d_all$sequenced, ], mut, by = "patientId")
cat("analysis set:", nrow(d), "patients,", sum(d$status), "deaths\n")

# TAS distribution by mutation status.
assoc <- do.call(rbind, lapply(GENES, function(g) {
  data.frame(
    gene = g,
    n_mut = sum(d[[g]] == 1),
    n_wt = sum(d[[g]] == 0),
    median_TAS_mut = median(d$TAS[d[[g]] == 1]),
    median_TAS_wt = median(d$TAS[d[[g]] == 0]),
    wilcox_p = wilcox.test(d$TAS[d[[g]] == 1], d$TAS[d[[g]] == 0])$p.value
  )
}))
assoc$fdr <- p.adjust(assoc$wilcox_p, method = "BH")
print(assoc)

# Descriptives above stay on the raw TAS scale; Cox HRs are per SD of the full TCGA-KIRC cohort,
# so they are comparable with Table 1.
d$TAS <- d$TAS / sd(score$TAS)

mut_terms <- paste(GENES, collapse = " + ")
forms <- list(
  clinical = Surv(time, status) ~ AGE + SEX + STAGE,
  clinical_tas = Surv(time, status) ~ TAS + AGE + SEX + STAGE,
  clinical_mut = as.formula(paste("Surv(time, status) ~ AGE + SEX + STAGE +", mut_terms)),
  clinical_mut_tas = as.formula(paste("Surv(time, status) ~ TAS + AGE + SEX + STAGE +", mut_terms)),
  strat_tas = Surv(time, status) ~ TAS + AGE + SEX + strata(STAGE),
  strat_mut_tas = as.formula(paste("Surv(time, status) ~ TAS + AGE + SEX +", mut_terms, "+ strata(STAGE)"))
)
fits <- lapply(forms, coxph, data = d, ties = "efron")

hr_row <- function(name, fit, term) {
  sm <- summary(fit)
  data.frame(model = name, term = term,
             HR = sm$conf.int[term, "exp(coef)"],
             lo = sm$conf.int[term, "lower .95"],
             hi = sm$conf.int[term, "upper .95"],
             p = sm$coefficients[term, "Pr(>|z|)"],
             harrell_C = unname(sm$concordance["C"]))
}
hr <- rbind(
  hr_row("clinical + TAS", fits$clinical_tas, "TAS"),
  hr_row("clinical + mutations + TAS", fits$clinical_mut_tas, "TAS"),
  hr_row("stage-stratified + TAS", fits$strat_tas, "TAS"),
  hr_row("stage-stratified + mutations + TAS", fits$strat_mut_tas, "TAS"),
  do.call(rbind, lapply(GENES, function(g) hr_row("clinical + mutations + TAS", fits$clinical_mut_tas, g)))
)
print(hr, digits = 4)

lrt <- function(reduced, full, label) {
  a <- anova(fits[[reduced]], fits[[full]])
  data.frame(contrast = label, df = a$Df[2], chisq = a$Chisq[2], p = a$`Pr(>|Chi|)`[2],
             delta_C = summary(fits[[full]])$concordance["C"] - summary(fits[[reduced]])$concordance["C"])
}
lrts <- rbind(
  lrt("clinical", "clinical_tas", "TAS added to clinical"),
  lrt("clinical_mut", "clinical_mut_tas", "TAS added to clinical + mutations"),
  lrt("clinical", "clinical_mut", "mutations added to clinical"),
  lrt("clinical_tas", "clinical_mut_tas", "mutations added to clinical + TAS")
)
rownames(lrts) <- NULL
print(lrts, digits = 4)

set.seed(20260923)
B <- 1000
bd <- matrix(NA_real_, B, 2, dimnames = list(NULL, c("clinical", "clinical_mut")))
for (b in seq_len(B)) {
  db <- d[sample.int(nrow(d), replace = TRUE), ]
  if (length(unique(db$STAGE)) < 4 || any(sapply(GENES, function(g) sum(db[[g]])) < 5)) next
  r <- tryCatch({
    cc <- function(f) summary(coxph(f, data = db, ties = "efron"))$concordance["C"]
    c(cc(forms$clinical_tas) - cc(forms$clinical),
      cc(forms$clinical_mut_tas) - cc(forms$clinical_mut))
  }, error = function(e) c(NA, NA))
  bd[b, ] <- r
}
boot <- data.frame(
  contrast = c("TAS added to clinical", "TAS added to clinical + mutations"),
  n_boot = colSums(is.finite(bd)),
  mean_delta_C = colMeans(bd, na.rm = TRUE),
  lo = apply(bd, 2, quantile, 0.025, na.rm = TRUE),
  hi = apply(bd, 2, quantile, 0.975, na.rm = TRUE),
  row.names = NULL
)
print(boot, digits = 4)

fit_zph <- coxph(Surv(time, status) ~ TAS + AGE + SEX + VHL + PBRM1 + SETD2 + BAP1 + strata(STAGE),
                 data = d, ties = "efron")
zph <- cox.zph(fit_zph, transform = "km")
zph_df <- data.frame(term = rownames(zph$table), zph$table, row.names = NULL)
colnames(zph_df) <- c("term", "chisq", "df", "p")
print(zph_df, digits = 4)

write.csv(sel, file.path(OUT, "mutation_selection_check.csv"), row.names = FALSE)
write.csv(assoc, file.path(OUT, "mutation_tas_association.csv"), row.names = FALSE)
write.csv(hr, file.path(OUT, "mutation_adjusted_hr.csv"), row.names = FALSE)
write.csv(lrts, file.path(OUT, "mutation_lrt.csv"), row.names = FALSE)
write.csv(boot, file.path(OUT, "mutation_bootstrap_delta_c.csv"), row.names = FALSE)
write.csv(zph_df, file.path(OUT, "mutation_ph_schoenfeld.csv"), row.names = FALSE)
cat("=== done ===\n")
