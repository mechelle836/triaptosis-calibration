#!/usr/bin/env Rscript
# 48_gene_drivers_by_cancer.R — 分解 OV/COADREAD/LUSC 高 Stouffer 的驱动基因
# 目的：triaptosis 21 基因在 KIRC/OV/COADREAD/LUSC 的逐基因 Cox P 对比，
#   判断高 Stouffer 癌种是否由同一批（应激类）基因驱动。
# 数据与口径与 45/46/47 完全一致（同矩阵、同线性尺度、同样本筛选）。
# 输出: results/01_screen/gene_drivers_by_cancer.csv

suppressPackageStartupMessages({ library(data.table); library(jsonlite); library(survival) })
ROOT <- "/Users/xiaolin/WorkBuddy/2026-09-04-20-00-08/paper2_triaptosis"
DATA <- file.path(ROOT, "data"); OUT <- file.path(ROOT, "results/01_screen")
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a)) b else a

TRIA21 <- c("PIK3C3","PIK3R4","MTM1","KEAP1","NFE2L2","SH3GL3","ELMO2","NCKAP1","ACTN2",
            "KXD1","WDR91","FYCO1","RAB9A","SCARB2","ANXA8","CLCN4","ATP13A2",
            "UVRAG","ATG14","NRBF2","BECN1")
FOCUS <- c("KIRC","OV","COADREAD","LUSC")

cohorts <- lapply(FOCUS, function(cc) {
  expr <- fromJSON(file.path(DATA, paste0(cc, "_expr.json")))
  clin <- fromJSON(file.path(DATA, paste0(cc, "_clin.json")))
  os_status <- sapply(clin, function(x) as.character(x$OS_STATUS %||% ""))
  os_months <- suppressWarnings(as.numeric(sapply(clin, function(x) x$OS_MONTHS %||% NA)))
  dead <- grepl("DECEASED|1", toupper(os_status))
  keep <- !is.na(os_months) & os_months > 0 & (os_status != "")
  sid <- names(expr); pid <- substr(sid, 1, 12)
  common <- intersect(pid, names(clin)[keep])
  s_for_p <- sid[match(common, pid)]
  idx <- match(common, names(clin)[keep])
  list(cancer = cc, samples = s_for_p,
       time = os_months[keep][idx], status = as.numeric(dead[keep][idx]))
})
names(cohorts) <- FOCUS

xf <- file.path(DATA, "pancan/EBpp_geneExp.xena.gz")
all_samples <- unique(unlist(lapply(cohorts, `[[`, "samples")))
hdr <- names(fread(cmd = paste("gzip -dc", shQuote(xf), "| head -1"), header = TRUE))
keep_cols <- c(1, which(hdr %in% all_samples))
X <- fread(cmd = paste("gzip -dc", shQuote(xf)), select = keep_cols, header = TRUE, showProgress = FALSE)
genes <- X[[1]]; X <- as.matrix(X[, -1]); rownames(X) <- genes
X <- X[, !duplicated(colnames(X)), drop = FALSE]
alias <- c(C19orf50 = "KXD1", KIAA0831 = "ATG14")
rownames(X)[rownames(X) %in% names(alias)] <- alias[rownames(X)[rownames(X) %in% names(alias)]]
stopifnot(all(TRIA21 %in% rownames(X)))
X <- 2^X - 1
X[X < 0 & X > -1e-6] <- 0

res <- do.call(rbind, lapply(cohorts, function(co) {
  s <- co$samples; inx <- s %in% colnames(X)
  mat2 <- t(X[TRIA21, s[inx], drop = FALSE])
  sdf <- data.frame(time = co$time[inx], status = co$status[inx])
  out <- lapply(TRIA21, function(g) {
    v <- mat2[, g]
    fit <- tryCatch(coxph(Surv(time, status) ~ v, data = sdf), error = function(e) NULL)
    p <- if (is.null(fit)) NA_real_ else summary(fit)$coefficients[5]
    z <- if (is.null(fit)) NA_real_ else summary(fit)$coefficients[4]
    data.frame(cancer = co$cancer, gene = g, cox_z = round(z, 2), cox_p = signif(p, 3),
               neglog10p = round(-log10(max(p, 1e-300)), 2))
  })
  do.call(rbind, out)
}))
fwrite(res, file.path(OUT, "gene_drivers_by_cancer.csv"))

# 每癌种 Stouffer 贡献最大的 5 个基因（按 -log10P）
cat("\n=== 各癌种 -log10P 最高的 5 个 triaptosis 基因 ===\n")
for (cc in FOCUS) {
  rc <- res[res$cancer == cc, ]; rc <- rc[order(-rc$neglog10p), ]
  cat(sprintf("%-9s: %s\n", cc,
      paste(sprintf("%s(%.1f)", rc$gene[1:5], rc$neglog10p[1:5]), collapse = ", ")))
}
# 应激类 vs VPS34 核心类分组比较
stress <- c("KEAP1","NFE2L2"); vps34core <- c("PIK3C3","PIK3R4","UVRAG","ATG14","BECN1","NRBF2")
cat("\n=== 分组 mean(-log10P)：应激类(KEAP1/NFE2L2) vs VPS34 核心复合体 ===\n")
for (cc in FOCUS) {
  rc <- res[res$cancer == cc, ]
  cat(sprintf("%-9s 应激类 %.2f | VPS34核心 %.2f | 其余 %.2f\n", cc,
      mean(rc$neglog10p[rc$gene %in% stress]),
      mean(rc$neglog10p[rc$gene %in% vps34core]),
      mean(rc$neglog10p[!rc$gene %in% c(stress, vps34core)])))
}
cat("=== done ===\n")
