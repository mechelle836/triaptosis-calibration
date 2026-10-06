#!/usr/bin/env Rscript
# Frozen KIRC TAS in every soft-pass cancer (KIRC, STAD, SARC, LUSC).
# Same scoring as the other parallel cohorts: z-score the 21 genes within the cohort
# (missing values set to 0 after scaling), multiply by the frozen elastic-net beta,
# split at the cohort median. KIRC row is the TCGA training estimate from score_table.csv.

suppressPackageStartupMessages({ library(jsonlite); library(survival) })
ROOT <- "/Users/xiaolin/WorkBuddy/2026-09-04-20-00-08/paper2_triaptosis"
RES  <- file.path(ROOT, "results")
mod <- readRDS(file.path(RES, "02_model/KIRC/model.rds"))
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || (length(a) == 1 && is.na(a))) b else a

hr_row <- function(cancer, time, status, tas) {
  g <- factor(ifelse(tas > median(tas), "High", "Low"), c("Low", "High"))
  fx <- coxph(Surv(time, status) ~ g)
  ci <- exp(confint(fx))
  lr <- survdiff(Surv(time, status) ~ g)
  data.frame(cancer = cancer, n = length(time), events = sum(status),
             HR = unname(exp(coef(fx))), lo = ci[1], hi = ci[2],
             km_p = 1 - pchisq(lr$chisq, 1), row.names = NULL)
}

parallel <- function(cancer) {
  expr <- fromJSON(file.path(ROOT, "data", paste0(cancer, "_expr.json")))
  clin <- fromJSON(file.path(ROOT, "data", paste0(cancer, "_clin.json")))
  mat <- t(vapply(expr, function(x) {
    v <- rep(NA_real_, length(mod$genes)); names(v) <- mod$genes
    for (g in mod$genes) if (!is.null(x[[g]])) v[g] <- as.numeric(x[[g]])
    v
  }, FUN.VALUE = numeric(length(mod$genes))))
  rownames(mat) <- substr(names(expr), 1, 12)
  mat <- mat[!duplicated(rownames(mat)), , drop = FALSE]
  os_status <- sapply(clin, function(x) as.character(x$OS_STATUS %||% ""))
  os_months <- suppressWarnings(as.numeric(sapply(clin, function(x) x$OS_MONTHS %||% NA)))
  dead <- grepl("DECEASED|1", toupper(os_status))
  keep <- !is.na(os_months) & os_months > 0 & os_status != ""
  common <- intersect(rownames(mat), names(clin)[keep])
  idx <- match(common, names(clin)[keep])
  z <- scale(mat[common, , drop = FALSE]); z[!is.finite(z)] <- 0
  cat(sprintf("%s: n=%d, samples with any missing gene=%d\n", cancer, length(common),
              sum(!complete.cases(mat[common, , drop = FALSE]))))
  hr_row(cancer, os_months[keep][idx], as.numeric(dead[keep][idx]), as.numeric(z %*% mod$beta))
}

sc <- read.csv(file.path(RES, "02_model/KIRC/score_table.csv"), row.names = 1)
out <- rbind(hr_row("KIRC", sc$time, sc$status, sc$TAS),
             do.call(rbind, lapply(c("STAD", "SARC", "LUSC"), parallel)))
print(out, digits = 3)
write.csv(out, file.path(RES, "02_model/parallel_cancer_hr.csv"), row.names = FALSE)
