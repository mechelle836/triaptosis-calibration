#!/usr/bin/env Rscript
# Elastic-net coefficient path for TAS construction (supplementary).
# Refits the full glmnet path (alpha = 0.5) on the same z-scored 21-gene TCGA-KIRC matrix
# as code/03_focus_model.R and marks the frozen lambda.min from model.rds.
# Stops if the refitted coefficients at that lambda do not match the frozen model.

suppressPackageStartupMessages({
  library(jsonlite); library(survival); library(glmnet)
  library(ggplot2); library(dplyr); library(showtext)
})
ROOT <- "/Users/xiaolin/WorkBuddy/2026-09-04-20-00-08/paper2_triaptosis"
OUT  <- file.path(ROOT, "figures_pub")
source(file.path(ROOT, "code", "lib", "viz_base.R"))
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a)) b else a

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
mat <- mat[!duplicated(rownames(mat)), , drop = FALSE]
os_status <- sapply(clin, function(x) as.character(x$OS_STATUS %||% ""))
os_months <- suppressWarnings(as.numeric(sapply(clin, function(x) x$OS_MONTHS %||% NA)))
dead <- grepl("DECEASED|1", toupper(os_status))
keep <- !is.na(os_months) & os_months > 0 & os_status != ""
common <- intersect(rownames(mat), names(clin)[keep])
idx <- match(common, names(clin)[keep])
z <- scale(mat[common, , drop = FALSE])
y <- Surv(os_months[keep][idx], as.numeric(dead[keep][idx]))

mod <- readRDS(file.path(ROOT, "results/02_model/KIRC/model.rds"))
stopifnot(identical(colnames(z), mod$genes))
fit <- glmnet(z, y, family = "cox", alpha = 0.5)
b_ref <- as.numeric(coef(fit, s = mod$lambda, exact = TRUE, x = z, y = y,
                         family = "cox", alpha = 0.5))
dev <- max(abs(b_ref - mod$beta))
cat(sprintf("n=%d; max |refit - frozen beta| at lambda.min = %.2e\n", nrow(z), dev))
if (dev > 1e-3) stop("refitted path does not reproduce the frozen coefficients")

B <- as.matrix(fit$beta)
lab <- data.frame(gene = mod$genes, coef = mod$beta) %>% filter(gene %in% mod$genes_sel)
path <- data.frame(gene = rep(rownames(B), ncol(B)),
                   loglam = rep(log(fit$lambda), each = nrow(B)),
                   coef = as.vector(B)) %>%
  left_join(transmute(lab, gene, sign = ifelse(coef > 0, "risk", "protective")), by = "gene") %>%
  mutate(sign = ifelse(is.na(sign), "dropped", sign))
x0 <- log(mod$lambda)
# 标签放在路径最左端: 系数已经完全分开, 不会像钉在 lambda.min 上那样叠成一摞
lab_left <- path %>%
  filter(sign != "dropped") %>%
  group_by(gene, sign) %>%
  slice_min(loglam, n = 1, with_ties = FALSE) %>%
  ungroup()

p <- ggplot(path, aes(loglam, coef, group = gene, colour = sign)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60") +
  geom_line(data = filter(path, sign == "dropped"), linewidth = 0.35) +
  geom_line(data = filter(path, sign != "dropped"), linewidth = 0.55) +
  geom_vline(xintercept = x0, linetype = "dashed", linewidth = 0.35) +
  ggrepel::geom_text_repel(
    data = lab_left, aes(label = gene),
    direction = "y", hjust = 1, nudge_x = -0.25,
    xlim = c(-Inf, min(path$loglam) - 0.08),
    size = 2.3, fontface = "italic", family = "Arial",
    segment.size = 0.2, segment.colour = "grey55",
    min.segment.length = 0.12, box.padding = 0.22, point.padding = 0.15,
    seed = 1, show.legend = FALSE, max.overlaps = Inf) +
  scale_colour_manual(values = c(risk = "#E64B35", protective = "#3C5488", dropped = "grey75"),
                      breaks = c("risk", "protective", "dropped"),
                      labels = c("Retained, HR > 1", "Retained, HR < 1", "Dropped (11 genes)"),
                      name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.22, 0.03))) +
  labs(x = expression(log(lambda)), y = "Cox coefficient (per SD of gene expression)",
       title = sprintf("Elastic-net path (alpha = 0.5); dashed line: lambda.min, %d genes retained",
                       length(mod$genes_sel))) +
  theme_cns(base_size = 8) +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 8.5, face = "bold", hjust = 0.5))
save_pub(p, "FigS_enet_path", OUT, w = 150, h = 100)
