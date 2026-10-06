# =============================================================================
# 07_triaptosis_null.R — Triaptosis 签名 (TAS) 的随机基因零分布检验
#
# 要回答的问题:
#   TAS (弹性网从 21 个 triaptosis 基因选出 10 个) 在 TCGA-KIRC 上的
#   预后能力, 是 triaptosis 基因集的功劳, 还是"弹性网从任意 21 个基因里
#   都能榨出类似信号"的普遍现象?
#
# 设计 (与 mitoxyperilysis 的检验同构, 保证可比):
#   参照: 在 21 个 triaptosis 基因上跑弹性网 (alpha=0.5, cv.glmnet) -> TAS
#   零分布: 从全转录组随机抽 21 个基因 x N, 跑完全相同的弹性网流程
#   指标: C-index (固定分数 10 折, concordance reverse=TRUE) + 嵌套 CV
#   基因重叠: 检查 triaptosis 与铁死亡/凋亡核心机器的重叠
# =============================================================================

suppressPackageStartupMessages({
  library(survival); library(jsonlite); library(glmnet)
})
PROJ <- "/Users/xiaolin/CodeBuddy/Claw/paper2_triaptosis"
KIRC <- "/Users/xiaolin/CodeBuddy/Claw/pan-cancer_disulfidptosis_project"
OUT  <- file.path(PROJ, "results/05_null")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

N_RANDOM <- as.integer(Sys.getenv("N_RANDOM", "20"))

TRIA21 <- c("ACTN2","ANXA8","ATG14","ATP13A2","BECN1","CLCN4","ELMO2","FYCO1",
            "KEAP1","KXD1","MTM1","NCKAP1","NFE2L2","NRBF2","PIK3C3","PIK3R4",
            "RAB9A","SCARB2","SH3GL3","UVRAG","WDR91")

# ------------------------------------------------------------- 数据加载 ---
cat(">>> load KIRC full transcriptome + locked cohort\n")
locked <- readRDS(file.path(KIRC, "results/04_mitoxy/KIRC_model_results.rds"))
d0 <- locked$d[!duplicated(locked$d$patient), ]
mat_list <- fromJSON(file.path(KIRC, "data/KIRC_full_transcriptome.json"))
common <- intersect(names(mat_list), d0$sample)
d <- d0[match(common, d0$sample), ]
allg <- sort(unique(unlist(lapply(mat_list[common], names))))
mat <- matrix(NA_real_, length(allg), length(common), dimnames = list(allg, common))
for (i in seq_along(common)) {
  v <- mat_list[[common[i]]]; mat[names(v), i] <- unlist(v)
}
X <- t(mat)
keep <- colnames(X)[colSums(is.na(X)) == 0 & apply(X, 2, sd, na.rm = TRUE) > 0]
X <- X[, keep, drop = FALSE]
dat <- data.frame(os_t = d$os_t, os_s = d$os_s); rownames(dat) <- d$sample
dat <- dat[!is.na(dat$os_t) & !is.na(dat$os_s) & dat$os_t > 0, ]
X <- X[rownames(dat), , drop = FALSE]
cat("   cohort:", nrow(X), "x", ncol(X), "| events =", sum(dat$os_s), "\n")

# ------------------------------------------------------------- 基因重叠 ---
FERRO <- c("GPX4","SLC7A11","SLC3A2","GCLC","GCLM","GSR","GSS","TFRC","FTH1",
           "FTL","NCOA4","HMOX1","IREB2","ACSL4","LPCAT3","ALOX15","NFE2L2","AIFM2","SAT1")
APOP  <- c("BAX","BAK1","BID","BCL2","BCL2L11","BCL2L1","BCL2L2","MCL1","BAD",
           "BBC3","PMAIP1","APAF1","CASP9","CASP3","CASP7","CYCS","DIABLO","XIAP","TP53")
DISUL <- c("SLC7A11","NCKAP1","NCF1","RAC1","ACTB","MYH9","ACTN4","FLNA","FLNB",
           "IQGAP1","CD2AP","DSTN","LIMK1","PAK1","WASF2","NCKAP1L","CYBA","NCF2","RAC2")
cat("\n>>> 基因重叠:\n")
cat("   Triaptosis ∩ Ferroptosis:", paste(intersect(TRIA21, FERRO), collapse = ", "), "\n")
cat("   Triaptosis ∩ Apoptosis  :", paste(intersect(TRIA21, APOP),  collapse = ", "), "\n")
cat("   Triaptosis ∩ Disulfidptosis:", paste(intersect(TRIA21, DISUL), collapse = ", "), "\n")

# ------------------------------------------------------------- 弹性网流程 ---
fit_enet <- function(Z, time, status, alpha = 0.5) {
  cvf <- cv.glmnet(Z, Surv(time, status), family = "cox", alpha = alpha,
                   nfolds = 10, cox.ties = "efron")
  beta <- as.numeric(coef(cvf, s = "lambda.min"))
  names(beta) <- colnames(Z)
  beta
}
score_from_beta <- function(Z, beta) as.numeric(Z %*% beta)

K <- 10
set.seed(2026)
FOLDS <- sample(rep(1:K, length.out = nrow(dat)))
cv_cindex <- function(score) {
  pred <- rep(NA_real_, nrow(dat))
  for (k in 1:K) {
    tr <- FOLDS != k; te <- FOLDS == k
    if (sum(dat$os_s[tr] == 1) < 5) next
    tmp <- data.frame(dat, .sc = score)
    fit <- try(coxph(Surv(os_t, os_s) ~ .sc, data = tmp[tr, ]), silent = TRUE)
    if (inherits(fit, "try-error")) next
    pred[te] <- predict(fit, newdata = tmp[te, ], type = "lp")
  }
  ok <- !is.na(pred)
  as.numeric(concordance(Surv(os_t, os_s) ~ pred, data = dat[ok, ], reverse = TRUE)$concordance)
}
# 嵌套 CV: 弹性网选择 + 权重都折内重做 (最诚实)
cv_nested_enet <- function(genes) {
  pred <- rep(NA_real_, nrow(dat))
  for (k in 1:K) {
    tr <- FOLDS != k; te <- FOLDS == k
    if (sum(dat$os_s[tr] == 1) < 5) next
    Ztr <- X[tr, genes, drop = FALSE]; Ztr <- scale(Ztr)
    Ztr[!is.finite(Ztr)] <- 0
    beta <- try(fit_enet(Ztr, dat$os_t[tr], dat$os_s[tr]), silent = TRUE)
    if (inherits(beta, "try-error")) next
    Zte <- X[te, genes, drop = FALSE]
    mu <- attr(Ztr, "scaled:center"); sdv <- attr(Ztr, "scaled:scale"); sdv[sdv == 0] <- 1
    Zte <- scale(Zte, center = mu, scale = sdv); Zte[!is.finite(Zte)] <- 0
    pred[te] <- as.numeric(Zte %*% beta)
  }
  ok <- !is.na(pred)
  as.numeric(concordance(Surv(os_t, os_s) ~ pred, data = dat[ok, ], reverse = TRUE)$concordance)
}

run_set <- function(genes) {
  genes <- intersect(genes, colnames(X))
  Z <- scale(X[, genes, drop = FALSE]); Z[!is.finite(Z)] <- 0
  beta <- try(fit_enet(Z, dat$os_t, dat$os_s), silent = TRUE)
  if (inherits(beta, "try-error")) return(NULL)
  sc <- score_from_beta(Z, beta)
  nsel <- sum(beta != 0)
  if (nsel == 0) return(NULL)
  # 单变量 cox 与 KM
  tmp <- data.frame(dat, .sc = sc)
  cox <- summary(coxph(Surv(os_t, os_s) ~ .sc, data = tmp))
  hr <- cox$conf.int[1, "exp(coef)"]; pv <- cox$coefficients[1, "Pr(>|z|)"]
  med <- ifelse(sc > median(sc), "High", "Low")
  km_p <- 1 - pchisq(survdiff(Surv(os_t, os_s) ~ med, data = tmp)$chisq, 1)
  list(cv_fixed = cv_cindex(sc), cv_nested = cv_nested_enet(genes),
       HR = hr, p = pv, km_p = km_p, nsel = nsel)
}

# ------------------------------------------------------------- 参照: TAS ---
cat("\n>>> [1/2] TAS 参照\n")
r_tas <- run_set(TRIA21)
cat(sprintf("   选中基因 %d 个 | cv_fixed=%.4f  cv_nested=%.4f  HR=%.3f  p=%.3g  KM p=%.3g\n",
            r_tas$nsel, r_tas$cv_fixed, r_tas$cv_nested, r_tas$HR, r_tas$p, r_tas$km_p))
cat("   (主项目记录: HR=3.52, KM p=2.90e-14)\n")

# ------------------------------------------------------------- 随机零分布 ---
cat("\n>>> [2/2] 随机 21 基因零分布 (N =", N_RANDOM, ")\n")
pool <- colnames(X)
set.seed(2026)
rand_sets <- lapply(seq_len(N_RANDOM), function(i) sample(pool, length(TRIA21)))
res <- data.frame()
t0 <- Sys.time()
for (i in seq_len(N_RANDOM)) {
  r <- run_set(rand_sets[[i]])
  if (!is.null(r)) res <- rbind(res, data.frame(i = i, cv_fixed = r$cv_fixed,
                                                cv_nested = r$cv_nested, HR = r$HR,
                                                p = r$p, km_p = r$km_p, nsel = r$nsel))
  if (i %% 20 == 0) cat(sprintf("   %4d/%d  elapsed %.1f min\n", i, N_RANDOM,
                                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
cat("\n=== 随机 21 基因零分布 (固定分数 / 嵌套) ===\n")
qs <- function(v) round(quantile(v, c(0, .05, .25, .5, .75, .95, 1)), 4)
tabq <- rbind(fixed = qs(res$cv_fixed), nested = qs(res$cv_nested))
colnames(tabq) <- c("min", "5%", "25%", "50%", "75%", "95%", "max")
print(tabq)
cat("\n均值±SD: fixed", round(mean(res$cv_fixed), 4), "±", round(sd(res$cv_fixed), 4),
    " | nested", round(mean(res$cv_nested), 4), "±", round(sd(res$cv_nested), 4), "\n")

cmp <- data.frame(
  metric = c("cv_fixed", "cv_nested"),
  TAS = round(c(r_tas$cv_fixed, r_tas$cv_nested), 4),
  null_median = round(c(median(res$cv_fixed), median(res$cv_nested)), 4),
  null_mean = round(c(mean(res$cv_fixed), mean(res$cv_nested)), 4),
  null_sd = round(c(sd(res$cv_fixed), sd(res$cv_nested)), 4),
  pct_vs_null = round(c(mean(res$cv_fixed >= r_tas$cv_fixed),
                        mean(res$cv_nested >= r_tas$cv_nested)), 4),
  z_vs_null = round(c((r_tas$cv_fixed - mean(res$cv_fixed)) / sd(res$cv_fixed),
                      (r_tas$cv_nested - mean(res$cv_nested)) / sd(res$cv_nested)), 2)
)
cat("\n=== TAS vs 随机零分布 ===\n"); print(cmp, row.names = FALSE)
cat("\n随机签名达到 Cox p<0.05 的比例:", round(mean(res$p < 0.05), 3),
    " | 达到 KM p<0.05 的比例:", round(mean(res$km_p < 0.05), 3), "\n")

write.csv(res, file.path(OUT, "triaptosis_random_null.csv"), row.names = FALSE)
write.csv(cmp, file.path(OUT, "triaptosis_vs_null.csv"), row.names = FALSE)

cat("\n=== 判读 ===\n")
cat("  若 TAS 在零分布中位于高百分位 (如 >95%), 则 triaptosis 基因集\n")
cat("  确实携带超出随机的预后信号, 可支撑生物学论文; 反之则与 mitoxyperilysis 同命运。\n")
cat("\n=== done ===\n")
