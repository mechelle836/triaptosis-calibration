#!/usr/bin/env Rscript
# Sensitivity of the RNA-seq 16-gene score to low-expression compression.
# The published assay floors IL6 at 4 CT when expression is lower (US patent 10181008, Table 3 note).
# That cutoff is on the RT-PCR CT scale and has no published RSEM equivalent, so it is not applied here.
# As a sensitivity analysis, reference-centered log2(RSEM+1) values below a gene's own 10th percentile
# are set to that percentile, first for IL6 alone and then for all 11 cancer genes.

suppressPackageStartupMessages({ library(survival); library(data.table); library(jsonlite) })
P2 <- "/Users/xiaolin/WorkBuddy/2026-09-04-20-00-08/paper2_triaptosis"
sc <- read.csv(file.path(P2, "results/08_multivariable/signature_scores.csv"))
kexpr <- fromJSON(file.path(P2, "data/KIRC_expr.json"))
sid <- names(kexpr); pid <- substr(sid, 1, 12)
sc$sample <- sid[match(sc$patientId, pid)]
stopifnot(!anyNA(sc$sample))

genes <- c("APOLD1","EDNRB","NOS3","PPAP2B","EIF4EBP1","TUBB2A","LMNB1",
           "CEACAM1","CX3CL1","CCL5","IL6","AAMP","ARF1","ATP5E","GPX1","RPLP1")
xf <- file.path(P2, "data/pancan/EBpp_geneExp.xena.gz")
hdr <- names(fread(cmd = paste("gzip -dc", shQuote(xf), "| head -1"), header = TRUE))
cols <- c(1, which(hdr %in% sc$sample & !duplicated(hdr)))
X <- fread(cmd = paste("gzip -dc", shQuote(xf)), select = cols, header = TRUE, showProgress = FALSE)
X <- X[X[[1]] %in% genes]
L <- as.matrix(X[, -1]); rownames(L) <- X[[1]]
L <- L[, sc$sample]

ref <- colMeans(L[c("AAMP","ARF1","ATP5E","GPX1","RPLP1"), ])
N <- sweep(L[genes[1:11], ], 2, ref, "-")
rs <- function(M) {
  VN  <- (0.5 * M["APOLD1", ] + 0.5 * M["EDNRB", ] + M["NOS3", ] + M["PPAP2B", ]) / 4
  CGD <- (M["EIF4EBP1", ] + 1.3 * M["LMNB1", ] + M["TUBB2A", ]) / 3
  IR  <- (0.5 * M["CCL5", ] + M["CEACAM1", ] + M["CX3CL1", ]) / 3
  as.numeric(-0.45 * VN - 0.31 * IR + 0.27 * CGD + 0.04 * M["IL6", ])
}
base <- rs(N)
stopifnot(max(abs(base - sc$RS16)) < 1e-6)

floor_at <- function(M, which) {
  Z <- M
  for (g in which) {
    thr <- as.numeric(quantile(Z[g, ], 0.10))
    Z[g, Z[g, ] < thr] <- thr
  }
  Z
}
cmp <- function(alt, label) {
  rho <- cor(base, alt, method = "spearman")
  cd <- concordance(Surv(time, status) ~ alt, data = sc, reverse = TRUE)
  data.frame(version = label, spearman_vs_unfloored = rho,
             C = cd$concordance, C_lo = cd$concordance - 1.96 * sqrt(cd$var),
             C_hi = cd$concordance + 1.96 * sqrt(cd$var))
}
# score-only C of the original, for reference
out <- rbind(
  cmp(base, "unfloored"),
  cmp(rs(floor_at(N, "IL6")), "IL6 floored at its 10th percentile"),
  cmp(rs(floor_at(N, rownames(N))), "all 11 genes floored at their 10th percentile")
)
print(out, digits = 4)
write.csv(out, file.path(P2, "results/08_multivariable/rini_compression_sensitivity.csv"), row.names = FALSE)
