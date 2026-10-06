#!/usr/bin/env Rscript
# Export 21-gene x hard-pass-cancer univariate Cox (HR, P).
# Same Cox protocol as 02_pancancer_screen.R; does not re-rank cancers.
suppressMessages({
  library(jsonlite); library(survival)
})

ROOT <- "/Users/xiaolin/WorkBuddy/2026-09-04-20-00-08/paper2_triaptosis"
DATA <- file.path(ROOT, "data")
OUT  <- file.path(ROOT, "results/01_screen")
screen <- read.csv(file.path(OUT, "pancancer_screen_matrix.csv"),
                   stringsAsFactors = FALSE)
hard <- screen$cancer[isTRUE(screen$pass_hard) | screen$pass_hard %in% c(TRUE, "TRUE")]
# pass_hard may be logical or character after read.csv
hard <- screen$cancer[as.logical(screen$pass_hard)]
hard <- hard[!is.na(hard)]
cat("Hard-pass cancers:", length(hard), "\n")

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || (length(a) == 1 && is.na(a))) b else a

score_genes <- function(cancer) {
  expr <- fromJSON(file.path(DATA, paste0(cancer, "_expr.json")))
  clin <- fromJSON(file.path(DATA, paste0(cancer, "_clin.json")))
  if (is.null(clin) || length(clin) == 0) return(NULL)
  genes_all <- sort(unique(unlist(lapply(expr, names))))
  mat <- t(vapply(expr, function(x) {
    v <- rep(NA_real_, length(genes_all)); names(v) <- genes_all
    for (g in genes_all) if (!is.null(x[[g]])) v[g] <- as.numeric(x[[g]])
    v
  }, FUN.VALUE = numeric(length(genes_all))))
  colnames(mat) <- genes_all
  os_status <- sapply(clin, function(x) as.character(x$OS_STATUS %||% ""))
  os_months <- suppressWarnings(as.numeric(sapply(clin, function(x) x$OS_MONTHS %||% NA)))
  dead <- grepl("DECEASED|1", toupper(os_status))
  keep <- !is.na(os_months) & os_months > 0 & !is.na(dead) & (os_status != "")
  rownames(mat) <- substr(rownames(mat), 1, 12)
  common <- intersect(rownames(mat), names(clin)[keep])
  mat2 <- mat[common, , drop = FALSE]
  idx <- match(common, names(clin)[keep])
  sdf <- data.frame(time = os_months[keep][idx], status = as.numeric(dead[keep][idx]))
  rows <- lapply(colnames(mat2), function(g) {
    v <- mat2[, g]
    if (sum(!is.na(v)) < 2 || isTRUE(sd(v, na.rm = TRUE) == 0)) return(NULL)
    fit <- tryCatch(coxph(Surv(time, status) ~ v, data = sdf, na.action = na.omit),
                    error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    sm <- summary(fit)
    data.frame(cancer = cancer, gene = g,
               HR = as.numeric(sm$coefficients[2]),
               p = as.numeric(sm$coefficients[5]),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows[!sapply(rows, is.null)])
}

res <- lapply(hard, function(c) {
  cat("  ", c, "\n")
  score_genes(c)
})
gcox <- do.call(rbind, res[!sapply(res, is.null)])
write.csv(gcox, file.path(OUT, "gene_cancer_cox.csv"), row.names = FALSE)
cat("Wrote", nrow(gcox), "rows ->", file.path(OUT, "gene_cancer_cox.csv"), "\n")
