#!/usr/bin/env Rscript
# Paper2-Triaptosis: ICI 队列真实验证 (IMmotion150, 肾癌 anti-PD-L1, n=263)
# 用主线癌种模型权重对 ICI 队列打分 -> 响应/生存验证
suppressMessages({ library(jsonlite); library(survival) })

DATA <- "/Users/xiaolin/CodeBuddy/Claw/paper2_triaptosis/data"
MODEL_DIR <- "/Users/xiaolin/CodeBuddy/Claw/paper2_triaptosis/results/02_model/KIRC"
OUT <- "/Users/xiaolin/CodeBuddy/Claw/paper2_triaptosis/results/03_ici"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
set.seed(2026)

model <- readRDS(file.path(MODEL_DIR, "model.rds"))
beta <- model$beta
genes_all <- model$genes  # 与 beta 位置对齐的基因顺序
cat("模型基因数:", length(beta), "\n")

expr <- fromJSON(file.path(DATA, "IMmotion150_expr.json"))
clin <- fromJSON(file.path(DATA, "IMmotion150_clin.json"))

genes_all <- model$genes  # 与 beta 位置对齐的基因顺序
mat <- t(vapply(expr, function(x) {
  v <- rep(NA_real_, length(genes_all)); names(v) <- genes_all
  for (g in genes_all) if (!is.null(x[[g]])) v[g] <- as.numeric(x[[g]])
  v
}, FUN.VALUE = numeric(length(genes_all))))
colnames(mat) <- genes_all
rownames(mat) <- names(expr)
cat("IMmotion150 表达样本:", nrow(mat), "\n")

# 样本对齐 (IMmotion150 样本ID = p000xx)
common <- intersect(rownames(mat), names(clin))
cat("临床匹配:", length(common), "\n")
mat2 <- mat[common, , drop = FALSE]

# 用主线模型的 z 变换参数: 直接用模型 beta 打分 (表达已可比)
# 注意: 主线模型对 TCGA 做了 z-scale; ICI 队列单独 z-scale 后打分
z <- scale(mat2)
tas <- as.numeric(z %*% beta)
names(tas) <- common

med <- median(tas)
grp <- ifelse(tas > med, "High", "Low")
tab <- table(grp)
cat("分组: High=", tab["High"], " Low=", tab["Low"], "\n")

# 1) RESPONDER 差异 (核心结果)
resp <- sapply(clin[common], function(x) as.character(x$RESPONDER %||% x$RESPONSE %||% ""))
resp_bin <- ifelse(toupper(resp) %in% c("TRUE","1","YES","RESPONDER","CR","PR","COMPLETE RESPONSE","PARTIAL RESPONSE","R"), "R", "NR")
resp_bin[resp == ""] <- NA
tbl <- table(grp, resp_bin)
cat("\n=== RESPONDER 交叉表 ===\n"); print(tbl)
ft <- tryCatch(fisher.test(tbl), error = function(e) NULL)
if (!is.null(ft)) cat(sprintf("Fisher P = %.4f\n", ft$p.value))
cat(sprintf("High 组响应率 = %.1f%%, Low 组响应率 = %.1f%%\n",
            100*mean(resp_bin[grp=="High"]=="R", na.rm=TRUE),
            100*mean(resp_bin[grp=="Low"]=="R", na.rm=TRUE)))

# 2) PFS KM (免疫治疗生存)
pfs <- suppressWarnings(as.numeric(sapply(clin[common], function(x) x$PFS_MONTHS %||% NA)))
pfs_status <- sapply(clin[common], function(x) as.character(x$PFS_STATUS %||% ""))
pfs_dead <- grepl("1|Progressed|TRUE", toupper(pfs_status))
keep <- !is.na(pfs) & pfs > 0
if (sum(keep) > 20) {
  sdf <- data.frame(time = pfs[keep], status = as.numeric(pfs_dead[keep]), grp = grp[keep])
  lr <- survdiff(Surv(time, status) ~ grp, data = sdf)
  p_pfs <- 1 - pchisq(lr$chisq, df = 1)
  cat(sprintf("\n=== PFS KM: n=%d, log-rank P = %.4f ===\n", nrow(sdf), p_pfs))
  mcox <- tryCatch(coxph(Surv(time, status) ~ grp, data = sdf), error = function(e) NULL)
  if (!is.null(mcox)) cat(sprintf("Cox HR(High vs Low) = %.3f, P = %.4f\n",
      exp(coef(mcox)), summary(mcox)$coefficients[5]))
}

# 3) CLINICAL_BENEFIT
cb <- sapply(clin[common], function(x) as.character(x$CLINICAL_BENEFIT %||% ""))
if (any(cb != "")) {
  cb_bin <- ifelse(toupper(cb) %in% c("TRUE","1","YES"), "Benefit", "No")
  tbl2 <- table(grp, cb_bin)
  cat("\n=== CLINICAL_BENEFIT 交叉表 ===\n"); print(tbl2)
  ft2 <- tryCatch(fisher.test(tbl2), error = function(e) NULL)
  if (!is.null(ft2)) cat(sprintf("Fisher P = %.4f\n", ft2$p.value))
}

# 保存
out <- data.frame(sample = common, TAS = tas, group = grp,
                  responder = resp_bin, pfs = pfs, pfs_status = pfs_status)
write.csv(out, file.path(OUT, "immotion150_validation.csv"), row.names = FALSE)
cat("\n结果已保存:", file.path(OUT, "immotion150_validation.csv"), "\n")

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a)) b else a
