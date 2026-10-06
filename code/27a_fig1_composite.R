# =====================================================================
# 27a_fig1_composite.R —— Fig1 泛癌筛选（投稿版）
# A 计数条  B KIRC 加框气泡  C 分型 KM + at risk  D per-SD Cox 森林
# 输出: figures_pub/Fig1_screen_framework.{png,pdf}
# HR per SD: 与 02_pancancer_screen.R 同一队列，基因内 z-score 后重拟合。
# 单变量线性缩放不改变 P；原始每单位 HR 不再上图。
# =====================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(survival)
  library(jsonlite); library(showtext); library(patchwork)
  library(ggsurvfit)
})
ROOT <- "/Users/xiaolin/WorkBuddy/2026-09-04-20-00-08/paper2_triaptosis"
OUT  <- file.path(ROOT, "figures_pub")
RES  <- file.path(ROOT, "results")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
source(file.path(ROOT, "code", "lib", "viz_base.R"))

gc <- read.csv(file.path(RES, "01_screen/gene_cancer_cox.csv"))
sc <- read.csv(file.path(RES, "02_model/KIRC/score_table.csv"), row.names = 1)
screen <- read.csv(file.path(RES, "01_screen/pancancer_screen_matrix.csv"))
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a[1])) b else a

# ---- per-SD Cox（KIRC，对齐筛选脚本的样本） ---------------------------
expr <- fromJSON(file.path(ROOT, "data/KIRC_expr.json"))
clin <- fromJSON(file.path(ROOT, "data/KIRC_clin.json"))
genes_all <- sort(unique(unlist(lapply(expr, names))))
mat <- t(vapply(expr, function(x) {
  v <- rep(NA_real_, length(genes_all)); names(v) <- genes_all
  for (g in genes_all) if (!is.null(x[[g]])) v[g] <- as.numeric(x[[g]])
  v
}, FUN.VALUE = numeric(length(genes_all))))
colnames(mat) <- genes_all
rownames(mat) <- substr(names(expr), 1, 12)

os_status <- vapply(clin, function(x) as.character(x$OS_STATUS %||% ""), character(1))
os_months <- suppressWarnings(as.numeric(vapply(clin, function(x) {
  m <- x$OS_MONTHS %||% NA
  as.character(m)
}, character(1))))
dead <- grepl("DECEASED|1", toupper(os_status))
keep <- !is.na(os_months) & os_months > 0 & os_status != ""
names(os_months) <- names(clin); names(dead) <- names(clin)
common <- intersect(rownames(mat), names(clin)[keep])
mat2 <- mat[common, , drop = FALSE]
sdf <- data.frame(time = os_months[common], status = as.numeric(dead[common]))

hr_rows <- lapply(colnames(mat2), function(g) {
  vraw <- mat2[, g]
  if (sum(!is.na(vraw)) < 2 || isTRUE(sd(vraw, na.rm = TRUE) == 0)) return(NULL)
  vz <- as.numeric(scale(vraw))
  fit <- tryCatch(coxph(Surv(time, status) ~ vz, data = sdf, na.action = na.omit),
                  error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  sm <- summary(fit)
  data.frame(gene = g,
             HR = unname(sm$conf.int[1, 1]),
             lo = unname(sm$conf.int[1, 3]),
             hi = unname(sm$conf.int[1, 4]),
             p = unname(sm$coefficients[1, 5]),
             stringsAsFactors = FALSE)
})
hrsd <- bind_rows(hr_rows)
orig <- gc %>% filter(cancer == "KIRC") %>% select(gene, p_orig = p)
chk <- hrsd %>% left_join(orig, by = "gene")
dp <- max(abs(log10(chk$p) - log10(chk$p_orig)), na.rm = TRUE)
message(sprintf("KIRC per-SD vs original P: max |d log10 P| = %.3g (n genes = %d)",
                dp, nrow(chk)))
if (dp > 0.05) warning("per-SD P diverges from the screen table; check cohort alignment")
write.csv(hrsd, file.path(OUT, "fig1_kirc_hr_per_sd.csv"), row.names = FALSE)

# ---- A: 筛选计数 -------------------------------------------------------
n_all  <- nrow(screen)
n_hard <- sum(screen$pass_hard)
n_soft <- sum(screen$pass_soft)
soft_names <- paste(screen$cancer[screen$pass_soft], collapse = ", ")
# 条形右侧写出通过的癌种, 这一格才不只是 4 个数字
steps <- c("All TCGA types", "Hard pass", "Soft pass", "Selected")
cnt <- data.frame(
  step = factor(steps, levels = rev(steps)),
  n = c(n_all, n_hard, n_soft, 1),
  lab = c(as.character(n_all), as.character(n_hard),
          sprintf("%d", n_soft), "1  KIRC"),
  hit = c(FALSE, FALSE, FALSE, TRUE)
)
pA <- ggplot(cnt, aes(n, step, fill = hit)) +
  geom_col(width = 0.58, show.legend = FALSE,
           colour = "black", linewidth = 0.2) +
  geom_text(aes(label = lab), hjust = -0.15, size = 2.4, family = "Arial") +
  annotate("text", x = 7.2, y = "Soft pass", hjust = 0, vjust = 0.5,
           size = 2.0, colour = "grey30", family = "Arial", label = soft_names) +
  scale_fill_manual(values = c(`TRUE` = "#E64B35", `FALSE` = "#8491B4")) +
  scale_x_continuous(limits = c(0, 32), breaks = c(0, 10, 20, 30),
                     expand = expansion(mult = c(0, 0.30))) +
  labs(x = "Cancer types", y = NULL, title = "Screen counts", tag = "A") +
  theme_cns(base_size = 8) +
  theme(axis.text.y = element_text(size = 6.5),
        plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
        plot.margin = margin(2, 4, 2, 2, "mm"))

# ---- B: 气泡，KIRC 列最左并加框 ----------------------------------------
kirc_p <- gc %>% filter(cancer == "KIRC") %>% arrange(p)
gene_order <- kirc_p$gene
cancer_n <- gc %>% filter(p < 0.05) %>% count(cancer)
cancer_order <- c("KIRC", setdiff(cancer_n$cancer[order(-cancer_n$n)], "KIRC"))
gc2 <- gc %>%
  mutate(cancer = factor(cancer, levels = cancer_order),
         gene = factor(gene, levels = rev(gene_order)),
         log10p = pmin(-log10(p), 8),
         direction = case_when(p >= 0.05 ~ "ns",
                               HR > 1 ~ "risk",
                               TRUE ~ "protective"))
n_gene <- nlevels(gc2$gene)
pB <- ggplot(gc2, aes(cancer, gene)) +
  annotate("rect", xmin = 0.55, xmax = 1.45, ymin = 0.5, ymax = n_gene + 0.5,
           fill = "#E64B35", alpha = 0.08, colour = "#E64B35", linewidth = 0.4) +
  geom_point(data = function(d) filter(d, direction == "ns"),
             aes(colour = direction), size = 0.6) +
  geom_point(data = function(d) filter(d, direction != "ns"),
             aes(size = log10p, colour = direction)) +
  scale_colour_manual(values = c(risk = "#E64B35", protective = "#3C5488",
                                 ns = "grey80"),
                      breaks = c("protective", "risk", "ns"),
                      labels = c(protective = "protective", risk = "risk",
                                 ns = "P \u2265 0.05"),
                      name = NULL, drop = FALSE) +
  guides(colour = guide_legend(override.aes = list(size = c(1.9, 1.9, 0.9)))) +
  scale_size_continuous(range = c(1.2, 5.2), name = expression(-log[10](P)),
                        breaks = c(2, 4, 6, 8)) +
  labs(x = NULL, y = NULL, title = "Univariate Cox across hard-pass cancers", tag = "B") +
  theme_cns(base_size = 7.5) +
  theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 6),
        axis.text.y = element_text(face = "italic", size = 6.5),
        legend.position = "bottom", legend.box = "vertical",
        legend.text = element_text(size = 6),
        legend.spacing.y = unit(0.2, "mm"),
        legend.margin = margin(0.5, 0, 0, 0, "mm"),
        legend.box.margin = margin(0, 0, 0, 0, "mm"),
        plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
        plot.margin = margin(2, 6, 2, 2, "mm"))

# ---- C: k-means KM + at risk ------------------------------------------
# at-risk 表的行标签直接取分层名, 所以分层名必须是最终要显示的 C1/C2
sc$cluster <- factor(paste0("C", sc$cluster))
ev <- tapply(sc$status, sc$cluster, sum)
nn <- tapply(!is.na(sc$time), sc$cluster, sum)
med <- tapply(sc$time[sc$status == 1], sc$cluster[sc$status == 1], median)
worse <- names(which.min(med))
col_clu <- setNames(rep("#3C5488", nlevels(sc$cluster)), levels(sc$cluster))
col_clu[worse] <- "#E64B35"
lab_clu <- sprintf("%s (n=%d, events=%d)", levels(sc$cluster),
                   nn[levels(sc$cluster)], ev[levels(sc$cluster)])
names(lab_clu) <- levels(sc$cluster)
lr <- survdiff(Surv(time, status) ~ cluster, data = sc)
lrp <- 1 - pchisq(lr$chisq, 1)
pC <- survfit2(Surv(time, status) ~ cluster, data = sc) |>
  ggsurvfit(linewidth = 0.6) +
  add_confidence_interval(alpha = 0.12) +
  add_risktable(risktable_stats = "n.risk",
                stats_label = list(n.risk = "At risk"),
                size = 2.2,
                theme = km_risktable_theme()) +
  km_x_scale(breaks = seq(0, 140, 20)) +
  km_y_scale() +
  scale_colour_manual(values = col_clu, labels = lab_clu, name = NULL) +
  scale_fill_manual(values = col_clu, guide = "none") +
  annotate("text", x = 2, y = 0.06, hjust = 0, vjust = 0,
           label = sprintf("log-rank P = %.1e", lrp), size = 2.4, family = "Arial") +
  labs(x = "Time (months)", y = "Overall survival",
       title = sprintf("KIRC k-means clusters (n=%d)", nrow(sc)), tag = "C") +
  theme_cns(base_size = 8) +
  theme(legend.position = "top",
        legend.key.width = unit(5, "mm"),
        legend.text = element_text(size = 6.5),
        plot.title = element_text(size = 9, face = "bold", hjust = 0.5))
pC <- wrap_elements(full = ggsurvfit_build(pC))

# ---- D: per-SD 森林 ----------------------------------------------------
fd <- hrsd %>%
  mutate(sig = p < 0.05,
         direction = case_when(!sig ~ "ns", HR > 1 ~ "risk", TRUE ~ "protective"),
         gene = factor(gene, hrsd$gene[order(abs(log(hrsd$HR)))]))
n_sig_cbio <- sum(fd$sig)
pD <- ggplot(fd, aes(HR, gene, colour = direction)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.3, colour = "grey40") +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y",
                width = 0.15, linewidth = 0.35, show.legend = FALSE) +
  geom_point(size = 1.8) +
  scale_colour_manual(values = c(risk = "#E64B35", protective = "#3C5488", ns = "#B0B0B0"),
                      breaks = c("protective", "risk", "ns"),
                      labels = c(protective = "protective", risk = "risk",
                                 ns = "P \u2265 0.05"),
                      name = NULL) +
  scale_x_log10(expand = expansion(mult = c(0.04, 0.06))) +
  labs(x = "HR per 1 SD (95% CI)", y = NULL,
       title = sprintf("21 genes in KIRC (%d/21 P < 0.05, cBioPortal)", n_sig_cbio),
       tag = "D") +
  theme_cns(base_size = 8) +
  theme(axis.text.y = element_text(face = "italic", size = 6.5),
        legend.position = "bottom",
        plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
        plot.margin = margin(2, 3, 2, 2, "mm"))

# ---- E/F: random-gene calibration of the screen (code/38_screen_null.R) ------------
nul <- read.csv(file.path(RES, "01_screen/screen_null_sets.csv"))
rep_x <- read.csv(file.path(RES, "01_screen/screen_null_reproduction.csv"))
nsum <- read.csv(file.path(RES, "01_screen/screen_null_summary.csv"))
obs_k <- rep_x$n_sig[rep_x$cancer == "KIRC"]
kn <- nul[nul$cancer == "KIRC", ]
# 两段注释分层; 红字直接给经验 P, 灰字点明 KIRC 功效最大的 caveat (衔接 F 的逐癌种校准)
pE <- ggplot(kn, aes(n_sig)) +
  geom_histogram(binwidth = 1, fill = "#C8C8C8", colour = "white", linewidth = 0.2) +
  geom_vline(xintercept = obs_k, colour = "#E64B35", linewidth = 0.7) +
  annotate("text", x = obs_k - 0.5, y = 62, vjust = 1, hjust = 1, size = 2.2,
           colour = "#E64B35", family = "Arial", lineheight = 1.05,
           label = sprintf("Triaptosis set: %d/21\nEmpirical P = 0/%d\n(< 0.005)",
                           obs_k, nrow(kn))) +
  annotate("text", x = min(kn$n_sig) - 0.5, y = 62, vjust = 1, hjust = 0, size = 2.0,
           colour = "grey25", family = "Arial", lineheight = 1.05,
           label = sprintf("KIRC has the most\nevents: %d%% of random\nsets rank it first",
                           round(100 * nsum$kirc_top_frac))) +
  scale_y_continuous(limits = c(0, 64), breaks = seq(0, 40, 10), expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(4, 18, 2),
                     expand = expansion(mult = c(0.06, 0.06))) +
  labs(x = "Cox-significant genes in KIRC (random 21-gene sets)", y = "Random sets",
       title = sprintf("Null distribution in KIRC (%d random sets)", nrow(kn)), tag = "E") +
  theme_cns(base_size = 8) +
  theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5))

# F 改为 z 分数棒棒糖: 每个癌种用自己的随机零分布标准化, (obs - median)/IQR
# 只突出 KIRC 超过自身 95 分位, 其余灰蓝; 图例压到零
fsum <- nul %>% group_by(cancer) %>%
  summarise(med = median(n_sig), lo = quantile(n_sig, 0.05), hi = quantile(n_sig, 0.95),
            q1 = quantile(n_sig, 0.25), q3 = quantile(n_sig, 0.75), .groups = "drop") %>%
  left_join(rep_x %>% select(cancer, obs = n_sig), by = "cancer") %>%
  mutate(z = (obs - med) / (q3 - q1),
         above = obs > hi,
         cancer = factor(cancer, cancer[order(z)]))
z_kirc <- fsum$z[fsum$cancer == "KIRC"]
pF <- ggplot(fsum, aes(z, cancer)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3, colour = "grey40") +
  geom_segment(aes(x = 0, xend = z, yend = cancer), colour = "grey75",
               linewidth = 0.35, show.legend = FALSE) +
  geom_point(aes(colour = above), size = 2.0, show.legend = FALSE) +
  scale_colour_manual(values = c(`TRUE` = "#E64B35", `FALSE` = "#3C5488")) +
  annotate("text", x = z_kirc + 0.15, y = "KIRC", hjust = 0, vjust = 0.5, size = 2.1,
           colour = "#E64B35", family = "Arial",
           label = "above 95th pct") +
  scale_x_continuous(expand = expansion(mult = c(0.08, 0.42))) +
  labs(x = "Standardized excess, (observed \u2212 random median) / IQR", y = NULL,
       title = "Excess is specific to KIRC", tag = "F") +
  theme_cns(base_size = 8) +
  theme(axis.text.y = element_text(size = 6),
        plot.title = element_text(size = 9, face = "bold", hjust = 0.5))

fig1 <- ((pA | pB) + plot_layout(widths = c(0.72, 1.35))) /
  ((pC | pD) + plot_layout(widths = c(0.72, 1.35))) /
  (pE | pF) +
  plot_layout(heights = c(1.15, 1, 0.75)) &
  theme(plot.tag = element_text(face = "bold", size = 11))

save_pub(fig1, "Fig1_screen_framework", OUT, w = 180, h = 245)
message("soft-pass cancers: ", soft_names)
