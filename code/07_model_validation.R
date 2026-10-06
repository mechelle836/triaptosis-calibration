#!/usr/bin/env Rscript
# Paper2-Triaptosis: 预测级模型验证 (OOF 5折) — 依据 bio-machine-learning-survival-analysis 技能规范
# 指标: Uno IPCW C(tau) / 时间依赖 AUC(IPCW) / Brier & IBS vs KM / 校准曲线(分组)
# 原则: 全部 out-of-fold 预测, 无样本内泄漏
suppressMessages({
  library(jsonlite); library(survival); library(glmnet); library(pec)
})

DATA <- "/Users/xiaolin/CodeBuddy/Claw/paper2_triaptosis/data"
MODEL_DIR <- "/Users/xiaolin/CodeBuddy/Claw/paper2_triaptosis/results/02_model/KIRC"
OUT <- "/Users/xiaolin/CodeBuddy/Claw/paper2_triaptosis/results/05_validation"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
set.seed(2026)

expr <- fromJSON(file.path(DATA, "KIRC_expr.json"))
clin <- fromJSON(file.path(DATA, "KIRC_clin.json"))
genes_all <- sort(unique(unlist(lapply(expr, names))))
mat <- t(vapply(expr, function(x) {
  v <- rep(NA_real_, length(genes_all)); names(v) <- genes_all
  for (g in genes_all) if (!is.null(x[[g]])) v[g] <- as.numeric(x[[g]])
  v
}, FUN.VALUE = numeric(length(genes_all))))
colnames(mat) <- genes_all
rownames(mat) <- substr(names(expr), 1, 12)
mat <- mat[!duplicated(rownames(mat)), , drop = FALSE]

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a)) b else a
os_status <- sapply(clin, function(x) as.character(x$OS_STATUS %||% ""))
os_months <- suppressWarnings(as.numeric(sapply(clin, function(x) x$OS_MONTHS %||% NA)))
dead <- grepl("DECEASED|1", toupper(os_status))
keep <- !is.na(os_months) & os_months > 0 & os_status != ""
common <- intersect(rownames(mat), names(clin)[keep])
idx <- match(common, names(clin)[keep])
sdf <- data.frame(time = os_months[keep][idx], status = as.numeric(dead[keep][idx]))
rownames(sdf) <- common
z <- scale(mat[common, , drop = FALSE])
cat("样本:", nrow(sdf), " 死亡:", sum(sdf$status), "\n")

# ---------- 5 折 OOF: 弹性网 Cox (alpha=0.5) ----------
K <- 5
folds <- sample(rep(1:K, length.out = nrow(sdf)))
oof_risk <- rep(NA_real_, nrow(sdf))
oof_surv <- matrix(NA_real_, nrow(sdf), 3, dimnames = list(NULL, c("t12", "t36", "t60")))
times_eval <- c(12, 36, 60)

for (k in 1:K) {
  tr <- folds != k; te <- folds == k
  cvf <- cv.glmnet(z[tr, , drop = FALSE], Surv(sdf$time[tr], sdf$status[tr]),
                   family = "cox", alpha = 0.5, nfolds = 10)
  eta_te <- as.numeric(predict(cvf, newx = z[te, , drop = FALSE], s = "lambda.min"))
  oof_risk[te] <- eta_te
  # 校准模型: 训练折内 risk->Cox, 预测验证折生存概率
  cal <- tryCatch(coxph(Surv(time, status) ~ eta, data = data.frame(
    time = sdf$time[tr], status = sdf$status[tr],
    eta = as.numeric(predict(cvf, newx = z[tr, , drop = FALSE], s = "lambda.min")))),
    error = function(e) NULL)
  if (!is.null(cal)) {
    nd <- data.frame(eta = eta_te, time = 0, status = 0)
    sf <- survfit(cal, newdata = nd)          # predict(type=survival) 在此版本有 bug, 用 survfit
    ss <- summary(sf, times = times_eval)
    sp <- if (!is.null(ss$surv)) t(ss$surv) else NULL  # 样本 x times
    if (!is.null(sp) && nrow(sp) == sum(te) && ncol(sp) == length(times_eval)) {
      oof_surv[te, ] <- sp
    } else {
      cat(sprintf("  折 %d: survfit 提取失败 (dim %s)\n", k, paste(dim(sp), collapse="x")))
    }
  }
  cat(sprintf("折 %d/%d 完成 (训练 n=%d)\n", k, K, sum(tr)))
}
sdf$risk <- oof_risk
sdf$S12 <- oof_surv[, 1]; sdf$S36 <- oof_surv[, 2]; sdf$S60 <- oof_surv[, 3]

# ---------- 指标 1: IPCW C-index (Uno-type, 12/36/60 月截断) ----------
uno_c <- NA
ipcw_c <- NULL
cat("预测概率 NA 数:", sum(is.na(oof_surv)), "\n")
if (requireNamespace("pec", quietly = TRUE) && sum(!is.na(sdf$S12)) > 100) {
  prob_mat_c <- as.matrix(sdf[, c("S12", "S36", "S60")])
  colnames(prob_mat_c) <- as.character(times_eval)
  cc <- tryCatch(pec::cindex(object = prob_mat_c,
                             formula = Surv(time, status) ~ 1, data = sdf,
                             eval.times = times_eval, pred.times = times_eval,
                             cens.model = "marginal", na.action = na.omit),
                 error = function(e) {cat("  cindex ERR:", conditionMessage(e), "\n"); NULL})
  if (!is.null(cc)) {
    mm <- cc$AppCindex$matrix  # 单预测器时可能是向量
    vals <- if (is.matrix(mm)) as.numeric(mm[1, ]) else as.numeric(mm)
    if (length(vals) == length(times_eval)) {
      ipcw_c <- data.frame(times = times_eval, C = vals)
      uno_c <- ipcw_c$C[ipcw_c$times == 60]
      cat("IPCW C (Uno-type):"); print(ipcw_c)
    } else {
      cat("  cindex 返回维度异常:", length(vals), "\n")
    }
  }
}

# ---------- 指标 2: 时间依赖 AUC (IPCW, 12/36/60 月) ----------
auc_res <- NULL
if (requireNamespace("riskRegression", quietly = TRUE)) {
  sc <- tryCatch(riskRegression::Score(
    list("TAS" = sdf$risk), formula = Surv(time, status) ~ 1, data = sdf,
    metrics = "auc", times = times_eval, cens.method = "ipcw", seed = 1),
    error = function(e) NULL)
  if (!is.null(sc)) {
    auc_res <- as.data.frame(sc$AUC$score)[, c("times", "AUC")]
    print(auc_res)
  }
}

# ---------- 指标 3: IPCW Brier + IBS vs KM (手写, Graf 1999; 预测概率来自 OOF) ----------
Gfit <- survfit(Surv(time, 1 - status) ~ 1, data = sdf)  # 删失分布 KM
Gfun <- function(t) {
  ss <- summary(Gfit, times = t)
  if (length(ss$surv) == 0) 1 else ss$surv[1]
}
brier_t <- function(t) {
  Gt <- Gfun(t); n <- nrow(sdf)
  num <- den <- 0
  for (i in seq_len(n)) {
    if (is.na(oof_surv[i, which(times_eval == t)])) next
    obs <- ifelse(sdf$time[i] > t, 1, 0)
    pred <- oof_surv[i, which(times_eval == t)]
    w <- if (sdf$time[i] <= t && sdf$status[i] == 1) 1 / max(Gfun(sdf$time[i]), 0.01)
         else if (sdf$time[i] > t) 1 / max(Gt, 0.01)
         else 0
    num <- num + (obs - pred)^2 * w; den <- den + w
  }
  num / den
}
b_vec <- sapply(times_eval, brier_t)
# KM 参照 Brier
km_b <- rep(NA, 3)
km_surv <- summary(survfit(Surv(time, status) ~ 1, data = sdf), times = times_eval)$surv
for (j in seq_along(times_eval)) {
  t <- times_eval[j]
  Gt <- Gfun(t); n <- nrow(sdf)
  num <- den <- 0
  for (i in seq_len(n)) {
    obs <- ifelse(sdf$time[i] > t, 1, 0)
    w <- if (sdf$time[i] <= t && sdf$status[i] == 1) 1 / max(Gfun(sdf$time[i]), 0.01)
         else if (sdf$time[i] > t) 1 / max(Gt, 0.01)
         else 0
    num <- num + (obs - km_surv[j])^2 * w; den <- den + w
  }
  km_b[j] <- num / den
}
# IBS = 梯形积分 (12-60 月)
ibs_tas <- sum(diff(c(0, times_eval)) * (c(0, b_vec[-length(b_vec)]) + b_vec) / 2) / max(times_eval)
ibs_km  <- sum(diff(c(0, times_eval)) * (c(0, km_b[-length(km_b)]) + km_b) / 2) / max(times_eval)
brier_res <- data.frame(model = "TAS", times = times_eval, Brier = b_vec)
brier_res <- rbind(brier_res, data.frame(model = "KM", times = times_eval, Brier = km_b))
cat("IPCW Brier (TAS vs KM):"); print(brier_res)
cat(sprintf("IBS: TAS=%.4f  KM=%.4f   (相对增益 %.1f%%)\n",
            ibs_tas, ibs_km, 100 * (ibs_km - ibs_tas) / ibs_km))

# ---------- 指标 4: 校准 (风险四分位组: 预测 vs KM 观察, 1/3/5 年) ----------
cal_tab <- NULL
q <- quantile(sdf$risk, probs = c(0, .25, .5, .75, 1), na.rm = TRUE)
grp4 <- cut(sdf$risk, breaks = q, include.lowest = TRUE, labels = paste0("Q", 1:4))
if (!any(is.na(grp4))) {
  cal_tab <- do.call(rbind, lapply(levels(grp4), function(g) {
    sel <- grp4 == g
    km <- survfit(Surv(time, status) ~ 1, data = sdf[sel, ])
    obs <- summary(km, times = times_eval)$surv
    pred <- colMeans(oof_surv[sel, , drop = FALSE], na.rm = TRUE)
    data.frame(group = g, n = sum(sel),
               pred_12 = pred[1], obs_12 = obs[1],
               pred_36 = pred[2], obs_36 = obs[2],
               pred_60 = pred[3], obs_60 = obs[3])
  }))
  cal_tab$absdiff <- rowMeans(abs(cal_tab[, c("pred_12","pred_36","pred_60")] -
                                    cal_tab[, c("obs_12","obs_36","obs_60")]))
  print(cal_tab)
}

# ---------- 保存 ----------
res_list <- list(uno_c = uno_c, ipcw_c = ipcw_c, auc = auc_res, brier = brier_res,
                 ibs_tas = ibs_tas, ibs_km = ibs_km, calibration = cal_tab)
saveRDS(res_list, file.path(OUT, "validation_metrics.rds"))
write.csv(sdf, file.path(OUT, "oof_scores.csv"), row.names = TRUE)
if (!is.null(auc_res))  write.csv(auc_res, file.path(OUT, "auc_ipcw.csv"), row.names = FALSE)
if (!is.null(brier_res)) write.csv(brier_res, file.path(OUT, "brier.csv"), row.names = FALSE)
if (!is.null(ipcw_c))   write.csv(ipcw_c, file.path(OUT, "ipcw_cindex.csv"), row.names = FALSE)
if (!is.null(cal_tab))  write.csv(cal_tab, file.path(OUT, "calibration.csv"), row.names = FALSE)
cat(sprintf("\n完成. Uno/IPCW C(60mo)=%.3f | 结果已保存: %s\n", uno_c, OUT))
