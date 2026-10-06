#!/usr/bin/env Rscript
# TAS versus stage + grade and two published, externally trained ccRCC signatures on the same patients.
#   ClearCode34 (Brooks et al., Eur Urol 2014): PAM nearest-centroid ccA/ccB classifier, retrained from the
#     authors' reference data (github.com/aedin/KidneyCancerSubtypes) and applied with their log2 +
#     TCGA-median centring; continuous score = posterior probability of ccB.
#   16-gene recurrence score (Rini et al., Lancet Oncol 2015): unscaled RS = -0.45 VN - 0.31 IR + 0.27 CGD
#     + 0.04 IL6 (US patent 10181008; Rini et al., Clin Cancer Res 2018, Supplementary Fig. S1), computed
#     from log2(RSEM+1) normalised to the mean of the five reference genes. RT-PCR compression of
#     low-expressing genes is not reproduced.
# TAS is trained on TCGA-KIRC, the two signatures are not; TAS is therefore reported both as the
# apparent score and as the 5-fold out-of-fold score (elastic net refitted within folds), and the
# out-of-fold score is the fair comparator.
# Expression for the comparators: PanCanAtlas EB++ matrix (UCSC Xena), same samples as the TAS cohort.

suppressPackageStartupMessages({
  library(survival); library(data.table); library(jsonlite); library(pamr)
})
P2  <- "/Users/xiaolin/WorkBuddy/2026-09-04-20-00-08/paper2_triaptosis"
OUT <- file.path(P2, "results", "08_multivariable")
KCS <- file.path(P2, "data/external/KidneyCancerSubtypes/data")
SEED <- 20260923; B <- 1000

# ---- cohort: stage + grade complete (as code/39_grade_adjustment.R) --------------------
src <- readLines(file.path(P2, "code/39_grade_adjustment.R"))
eval(parse(text = src[seq_len(grep("^d <- d0\\[!is.na\\(d0\\$GRADE\\), \\]", src))]))
oof <- read.csv(file.path(P2, "results/05_validation/oof_scores.csv"), row.names = 1)
d$TAS_oof <- oof[d$patientId, "risk"]
stopifnot(!anyNA(d$TAS_oof))

# ---- expression for the comparator genes ------------------------------------------------
kexpr <- fromJSON(file.path(P2, "data/KIRC_expr.json"))
sid <- names(kexpr); pid <- substr(sid, 1, 12)
d$sample <- sid[match(d$patientId, pid)]
load(file.path(KCS, "xn.rda")); load(file.path(KCS, "classes.rda"))
cc_genes <- rownames(xn)
rs_genes <- c("APOLD1", "EDNRB", "NOS3", "PPAP2B", "EIF4EBP1", "TUBB2A", "LMNB1",
              "CEACAM1", "CX3CL1", "CCL5", "IL6", "AAMP", "ARF1", "ATP5E", "GPX1", "RPLP1")
xf <- file.path(P2, "data/pancan/EBpp_geneExp.xena.gz")
hdr <- names(fread(cmd = paste("gzip -dc", shQuote(xf), "| head -1"), header = TRUE))
cols <- c(1, which(hdr %in% d$sample & !duplicated(hdr)))
X <- fread(cmd = paste("gzip -dc", shQuote(xf)), select = cols, header = TRUE, showProgress = FALSE)
X <- X[X[[1]] %in% c(cc_genes, rs_genes)]
L <- as.matrix(X[, -1]); rownames(L) <- X[[1]]           # log2(RSEM + 1)
stopifnot(all(c(cc_genes, rs_genes) %in% rownames(L)))
d <- d[d$sample %in% colnames(L), ]
L <- L[, d$sample]
cat("patients with comparator expression:", nrow(d), "deaths:", sum(d$status), "\n")

# ---- ClearCode34 -------------------------------------------------------------------------
med533 <- scan(file.path(KCS, "Median533_scalingFactor.txt"), quiet = TRUE)
stopifnot(length(med533) == length(cc_genes))
ours <- apply(L[cc_genes, ], 1, median)
cat(sprintf("ClearCode34 centring check: r(TCGA-533 medians, our medians) = %.3f\n", cor(med533, ours)))
set.seed(SEED)
cc_model <- pamr.train(list(x = scale(xn), y = classes, geneid = cc_genes, genenames = cc_genes))
newx <- sweep(L[cc_genes, ], 1, med533, "-")
post <- pamr.predict(cc_model, scale(newx), threshold = 0, type = "posterior")
d$CC34_ccB <- post[, "ccB"]
d$CC34_class <- factor(ifelse(d$CC34_ccB > 0.5, "ccB", "ccA"), c("ccA", "ccB"))
cat("ClearCode34 classes:"); print(table(d$CC34_class))

# ---- 16-gene recurrence score --------------------------------------------------------------
ref <- colMeans(L[c("AAMP", "ARF1", "ATP5E", "GPX1", "RPLP1"), ])
N <- sweep(L[rs_genes[1:11], ], 2, ref, "-")
VN  <- (0.5 * N["APOLD1", ] + 0.5 * N["EDNRB", ] + N["NOS3", ] + N["PPAP2B", ]) / 4
CGD <- (N["EIF4EBP1", ] + 1.3 * N["LMNB1", ] + N["TUBB2A", ]) / 3
IR  <- (0.5 * N["CCL5", ] + N["CEACAM1", ] + N["CX3CL1", ]) / 3
d$RS16 <- as.numeric(-0.45 * VN - 0.31 * IR + 0.27 * CGD + 0.04 * N["IL6", ])

# ---- same-scale comparison -------------------------------------------------------------------
scores <- c(TAS_apparent = "TAS", TAS_oof = "TAS_oof", ClearCode34 = "CC34_ccB", RS16 = "RS16")
for (s in scores) d[[paste0(s, "_z")]] <- as.numeric(scale(d[[s]]))
cstat <- function(v, dd = d) {
  cc <- concordance(Surv(time, status) ~ v, data = cbind(dd, v = v), reverse = TRUE)
  c(C = cc$concordance, lo = cc$concordance - 1.96 * sqrt(cc$var), hi = cc$concordance + 1.96 * sqrt(cc$var))
}
base <- Surv(time, status) ~ AGE + SEX + STAGE + GRADE
f_base <- coxph(base, data = d, ties = "efron")
C_base <- unname(summary(f_base)$concordance["C"])
uni <- do.call(rbind, lapply(names(scores), function(n) {
  z <- paste0(scores[[n]], "_z")
  s1 <- summary(coxph(as.formula(paste("Surv(time, status) ~", z)), data = d))
  fa <- coxph(as.formula(paste("Surv(time, status) ~ AGE + SEX + STAGE + GRADE +", z)), data = d, ties = "efron")
  sa <- summary(fa)
  data.frame(signature = n, t(cstat(d[[scores[[n]]]])),
             HR_perSD = s1$conf.int[1, 1], HR_lo = s1$conf.int[1, 3], HR_hi = s1$conf.int[1, 4],
             HR_adj = sa$conf.int[z, 1], HR_adj_lo = sa$conf.int[z, 3], HR_adj_hi = sa$conf.int[z, 4],
             p_adj = sa$coefficients[z, 5],
             C_clinical = C_base, C_clinical_plus = unname(sa$concordance["C"]),
             lrt_p = anova(f_base, fa)$`Pr(>|Chi|)`[2])
}))
uni$delta_C <- uni$C_clinical_plus - uni$C_clinical
cat("stage + grade only: C =", round(unname(summary(coxph(Surv(time, status) ~ STAGE + GRADE, data = d))$concordance["C"]), 4), "\n")

set.seed(SEED)
bt <- t(replicate(B, {
  db <- d[sample(nrow(d), replace = TRUE), ]
  cb <- unname(summary(coxph(base, data = db, ties = "efron"))$concordance["C"])
  vapply(scores, function(s) {
    f <- as.formula(paste("Surv(time, status) ~ AGE + SEX + STAGE + GRADE +", paste0(s, "_z")))
    unname(summary(coxph(f, data = db, ties = "efron"))$concordance["C"]) - cb
  }, numeric(1))
}))
uni$delta_lo <- apply(bt, 2, quantile, 0.025)
uni$delta_hi <- apply(bt, 2, quantile, 0.975)
pair <- data.frame(
  contrast = c("TAS_oof - ClearCode34", "TAS_oof - RS16", "TAS_apparent - ClearCode34", "TAS_apparent - RS16"),
  diff = c(uni$delta_C[2] - uni$delta_C[3], uni$delta_C[2] - uni$delta_C[4],
           uni$delta_C[1] - uni$delta_C[3], uni$delta_C[1] - uni$delta_C[4]),
  lo = c(quantile(bt[, 2] - bt[, 3], 0.025), quantile(bt[, 2] - bt[, 4], 0.025),
         quantile(bt[, 1] - bt[, 3], 0.025), quantile(bt[, 1] - bt[, 4], 0.025)),
  hi = c(quantile(bt[, 2] - bt[, 3], 0.975), quantile(bt[, 2] - bt[, 4], 0.975),
         quantile(bt[, 1] - bt[, 3], 0.975), quantile(bt[, 1] - bt[, 4], 0.975)))

# Head to head: clinical + all three (out-of-fold TAS)
f_all <- coxph(Surv(time, status) ~ AGE + SEX + STAGE + GRADE + TAS_oof_z + CC34_ccB_z + RS16_z,
               data = d, ties = "efron")
drop1_p <- sapply(c("TAS_oof_z", "CC34_ccB_z", "RS16_z"), function(v) {
  f0 <- update(f_all, as.formula(paste(". ~ . -", v)))
  anova(f0, f_all)$`Pr(>|Chi|)`[2]
})
sa <- summary(f_all)
h2h <- data.frame(term = c("TAS_oof_z", "CC34_ccB_z", "RS16_z"),
                  HR = sa$conf.int[c("TAS_oof_z", "CC34_ccB_z", "RS16_z"), 1],
                  lo = sa$conf.int[c("TAS_oof_z", "CC34_ccB_z", "RS16_z"), 3],
                  hi = sa$conf.int[c("TAS_oof_z", "CC34_ccB_z", "RS16_z"), 4],
                  lrt_drop_p = drop1_p, C_model = unname(sa$concordance["C"]))
cors <- cor(d[, c("TAS", "TAS_oof", "CC34_ccB", "RS16")], method = "spearman")
cc_km <- summary(coxph(Surv(time, status) ~ CC34_class, data = d))

print(uni, digits = 3); print(pair, digits = 3); print(h2h, digits = 3); print(round(cors, 3))
cat(sprintf("ClearCode34 ccB vs ccA HR %.2f (%.2f-%.2f)\n", cc_km$conf.int[1, 1], cc_km$conf.int[1, 3], cc_km$conf.int[1, 4]))
write.csv(uni, file.path(OUT, "signature_comparison.csv"), row.names = FALSE)
write.csv(pair, file.path(OUT, "signature_comparison_paired.csv"), row.names = FALSE)
write.csv(h2h, file.path(OUT, "signature_head_to_head.csv"), row.names = FALSE)
write.csv(cors, file.path(OUT, "signature_spearman.csv"))
write.csv(d[, c("patientId", "time", "status", "STAGE", "GRADE", "TAS", "TAS_oof", "CC34_ccB", "CC34_class", "RS16")],
          file.path(OUT, "signature_scores.csv"), row.names = FALSE)
