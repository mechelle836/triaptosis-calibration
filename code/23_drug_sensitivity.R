# =====================================================================
# 23_drug_sensitivity.R —— GDSC 药敏预测 (oncoPredict 等价: ComBat + ridge)
# 数据: GDSC2 FPKM 表达 + LN_IC50 (cog.sanger.ac.uk)
# 方法: ComBat 批次校正(GDSC细胞系 vs TCGA组织) + glmnet ridge 回归
#       每个药物单独建模 → 预测 KIRC 患者 IC50 → TAS 分组差异
# 输出: Fig11_drug_sensitivity.png + results/14_drug/
# =====================================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readxl); library(httr)
  library(org.Hs.eg.db); library(glmnet); library(sva)
  library(ggplot2); library(showtext); library(patchwork)
  library(jsonlite); library(clusterProfiler)
})
ROOT <- "/Users/xiaolin/WorkBuddy/2026-09-04-20-00-08/paper2_triaptosis"
GDSC <- file.path(ROOT, "data/gdsc")
OUT  <- file.path(ROOT, "figures_v2")
RES  <- file.path(ROOT, "results", "14_drug")
dir.create(RES, showWarnings = FALSE, recursive = TRUE)
source(file.path(ROOT, "code", "lib", "viz_base.R"))
showtext_auto(FALSE)
select <- dplyr::select   # 防 select 被覆盖 (sva/MASS 冲突)
font_add("Arial", regular = "/System/Library/Fonts/Supplemental/Arial.ttf")
showtext_opts(dpi = 300); showtext_auto()

score <- read.csv(file.path(ROOT, "results/02_model/KIRC/score_table.csv"),
                  row.names = 1)

# ---- 1. GDSC 表达 (基因 symbol 行 × SIDM 细胞系列) --------------------
# GDSC 表达: 读 Python 导出的干净 tsv.gz (37279 symbol x 1047 SIDM)
gdsc_fpkm <- read.delim(gzfile(file.path(GDSC, "gdsc_expr.tsv.gz")),
                        check.names = FALSE, stringsAsFactors = FALSE)
expr_mat <- as.matrix(gdsc_fpkm[, -1])
rownames(expr_mat) <- gdsc_fpkm[, 1]
cat("GDSC expr:", dim(expr_mat), "\n")

# ---- 2. IC50 (药物 × 细胞系 LN_IC50) ---------------------------------
ic_wide <- read.delim(gzfile(file.path(GDSC, "ic50.tsv.gz")),
                       check.names = FALSE, stringsAsFactors = FALSE,
                       na.strings = c("", "NA"))
rownames(ic_wide) <- ic_wide[, 1]
ic_wide <- as.matrix(ic_wide[, -1])
storage.mode(ic_wide) <- "double"
cat("IC50 matrix:", dim(ic_wide), " drugs x cell lines\n")
# 过滤: 至少 50% 细胞系有 IC50 的药物
ic_wide <- ic_wide[rowSums(!is.na(ic_wide)) >= 0.5 * ncol(ic_wide), ]
cat("drugs after filter:", nrow(ic_wide), "\n")

# ---- 3. 交集细胞系 + 交集基因 ----------------------------------------
cells_common <- intersect(colnames(expr_mat), colnames(ic_wide))
genes_use <- rownames(expr_mat)[!duplicated(rownames(expr_mat))]
expr_mat <- expr_mat[genes_use, cells_common]
ic_wide <- ic_wide[, cells_common]
cat("common cell lines:", length(cells_common), "\n")

# ---- 4. 抓 TCGA KIRC 表达 (缓存) --------------------------------------
kirc_rds <- file.path(RES, "kirc_expr_full.rds")
if (file.exists(kirc_rds)) {
  kirc_expr <- readRDS(kirc_rds)
  cat("loaded cached TCGA KIRC expr:", dim(kirc_expr), "\n")
} else {
  eg <- bitr(genes_use, fromType = "SYMBOL", toType = "ENTREZID",
             OrgDb = org.Hs.eg.db)
  eg <- eg[!duplicated(eg$SYMBOL) & !duplicated(eg$ENTREZID), ]
  cat("mapped genes:", nrow(eg), "/", length(genes_use), "\n")
  # 分批 POST fetch (每批 500 基因, body 传 entrezGeneIds 数组)
  prof <- "kirc_tcga_pan_can_atlas_2018_rna_seq_v2_mrna"
  slist <- "kirc_tcga_pan_can_atlas_2018_rna_seq_v2_mrna"
  url <- paste0("https://www.cbioportal.org/api/molecular-profiles/", prof,
                "/molecular-data/fetch")
  res <- list()
  chunks <- split(eg$ENTREZID, ceiling(seq_along(eg$ENTREZID) / 500))
  for (i in seq_along(chunks)) {
    body <- list(entrezGeneIds = as.list(chunks[[i]]), sampleListId = slist)
    r <- tryCatch(httr::POST(url,
                             body = jsonlite::toJSON(body, auto_unbox = TRUE),
                             httr::content_type_json(), httr::timeout(120),
                             httr::user_agent("Mozilla/5.0")),
                  error = function(e) NULL)
    if (is.null(r) || httr::status_code(r) != 200) {
      cat("  batch", i, "failed\n"); next
    }
    d <- jsonlite::fromJSON(httr::content(r, as = "text", encoding = "UTF-8"))
    if (is.data.frame(d) && nrow(d) > 0 &&
        all(c("entrezGeneId", "sampleId", "value") %in% names(d))) res[[i]] <- d
    if (i %% 20 == 0) cat("  batch", i, "/", length(chunks), "\n")
  }
  cat("batches collected:", length(res), "\n")
  res <- Filter(function(x) is.data.frame(x) && nrow(x) > 0, res)
  all_md <- do.call(rbind, lapply(res, function(x)
    data.frame(ENTREZID = as.character(x$entrezGeneId),
               sample = x$sampleId,
               value = as.numeric(x$value), stringsAsFactors = FALSE)))
  # entrez -> symbol 映射
  all_md <- all_md %>%
    left_join(eg %>% select(ENTREZID, SYMBOL) %>%
                mutate(ENTREZID = as.character(ENTREZID)), by = "ENTREZID") %>%
    filter(!is.na(SYMBOL))
  kirc_expr <- all_md %>%
    select(SYMBOL, sample, value) %>%
    group_by(SYMBOL, sample) %>% summarise(value = mean(value), .groups = "drop") %>%
    pivot_wider(names_from = sample, values_from = value) %>%
    tibble::column_to_rownames("SYMBOL")
  saveRDS(kirc_expr, kirc_rds)
  cat("TCGA KIRC expr:", dim(kirc_expr), "\n")
}

# ---- 5. 对齐基因 + ComBat + ridge ------------------------------------
genes_final <- intersect(rownames(expr_mat), rownames(kirc_expr))
cat("final genes:", length(genes_final), "\n")
gdsc_g <- expr_mat[genes_final, , drop = FALSE]
kirc_g <- kirc_expr[genes_final, , drop = FALSE]

# log2 变换 + 合并 + ComBat
gdsc_l <- log2(gdsc_g + 1)
kirc_l <- log2(kirc_g + 1)
M <- cbind(gdsc_l, kirc_l)
batch <- c(rep("GDSC", ncol(gdsc_l)), rep("TCGA", ncol(kirc_l)))
combat <- tryCatch(ComBat(dat = M, batch = batch), error = function(e) M)
gdsc_b <- combat[, 1:ncol(gdsc_l)]
kirc_b <- combat[, (ncol(gdsc_l) + 1):ncol(combat)]
# 过滤含 NA 的基因 (GDSC FPKM 有缺失, glmnet 不接受 NA)
keep_g <- rownames(gdsc_b)[rowSums(is.na(gdsc_b)) == 0 &
                           rowSums(is.na(kirc_b)) == 0]
gdsc_b <- gdsc_b[keep_g, , drop = FALSE]
kirc_b <- kirc_b[keep_g, , drop = FALSE]
cat("genes after NA filter:", length(keep_g), "\n")
# top 5000 高变基因 (标准基因选择, 加速 ridge 且降噪)
gdsc_var <- apply(gdsc_b, 1, var, na.rm = TRUE)
top_genes <- names(sort(gdsc_var, decreasing = TRUE))[1:min(5000, length(gdsc_var))]
gdsc_b <- gdsc_b[top_genes, , drop = FALSE]
kirc_b <- kirc_b[top_genes, , drop = FALSE]
cat("training genes:", length(top_genes), "\n")

# 对齐 KIRC 样本 (barcode 12 位) + TAS
kirc_samp <- substr(colnames(kirc_b), 1, 12)
common_s <- intersect(kirc_samp, rownames(score))
kirc_idx <- which(kirc_samp %in% common_s)
kirc_b2 <- kirc_b[, kirc_idx]
colnames(kirc_b2) <- kirc_samp[kirc_idx]
kirc_b2 <- kirc_b2[, common_s]
tas <- score[common_s, "TAS"]; names(tas) <- common_s

# ---- 6. ridge 训练 + 预测 --------------------------------------------
library(parallel)
train_one <- function(d) {
  y <- as.numeric(ic_wide[d, ])
  ok <- !is.na(y)
  if (sum(ok) < 30) return(NULL)
  cv <- cv.glmnet(t(gdsc_b[, ok]), y[ok], alpha = 0, nfolds = 5,
                  standardize = TRUE, family = "gaussian")
  predict(cv, newx = t(kirc_b2), s = "lambda.min")[, 1]
}
pred_rds <- file.path(RES, "predicted_ic50.rds")
if (file.exists(pred_rds)) {
  pred_ic50 <- readRDS(pred_rds)
  cat("loaded cached predicted IC50:", dim(pred_ic50), "\n")
} else {
  pred_list <- mclapply(rownames(ic_wide), train_one, mc.cores = 4)
  names(pred_list) <- rownames(ic_wide)
  pred_ic50 <- do.call(rbind, pred_list[!sapply(pred_list, is.null)])
  colnames(pred_ic50) <- common_s            # mclapply 向量无名, 补样本列名
  rownames(pred_ic50) <- names(pred_list)[!sapply(pred_list, is.null)]
  saveRDS(pred_ic50, pred_rds)
}
pred_df <- as.data.frame(t(pred_ic50))
pred_df$sample <- rownames(pred_df)
pred_df$TAS <- tas[pred_df$sample]
pred_df$group <- score[pred_df$sample, "group"]
cat("predicted IC50 for", nrow(pred_df), "samples\n")

# ---- 7. 药敏差异 (TAS High vs Low) -----------------------------------
drug_res <- data.frame(drug = rownames(ic_wide), p = NA, diff = NA,
                       high_mean = NA, low_mean = NA)
for (d in rownames(ic_wide)) {
  v <- pred_df[[d]]
  hi <- v[pred_df$group == "High"]; lo <- v[pred_df$group == "Low"]
  if (sum(!is.na(hi)) < 5 || sum(!is.na(lo)) < 5) next
  drug_res$p[drug_res$drug == d] <- wilcox.test(hi, lo)$p.value
  drug_res$diff[drug_res$drug == d] <- mean(hi, na.rm = TRUE) -
    mean(lo, na.rm = TRUE)
  drug_res$high_mean[drug_res$drug == d] <- mean(hi, na.rm = TRUE)
  drug_res$low_mean[drug_res$drug == d] <- mean(lo, na.rm = TRUE)
}
drug_res <- drug_res %>% filter(!is.na(p)) %>%
  mutate(fdr = p.adjust(p, "BH")) %>%
  arrange(p)
write.csv(drug_res, file.path(RES, "drug_sensitivity_diff.csv"),
          row.names = FALSE)
cat("drugs tested:", nrow(drug_res),
    " significant (FDR<0.05):", sum(drug_res$fdr < 0.05), "\n")

# ---- 8. 出图: A 差异排序 + B 关键药物箱线 -----------------------------
top_drugs <- drug_res %>% arrange(p) %>% slice(1:min(20, nrow(.)))
pA <- ggplot(top_drugs, aes(diff, reorder(drug, diff))) +
  geom_col(aes(fill = diff < 0), width = 0.6, show.legend = FALSE) +
  scale_fill_manual(values = c(`TRUE` = "#3C5488", `FALSE` = "#E64B35")) +
  geom_vline(xintercept = 0, colour = "grey30", linewidth = 0.4) +
  labs(x = "IC50 difference (High - Low)", y = NULL,
       title = "Drugs with differential sensitivity by TAS") +
  theme_prism(base_size = 8, base_family = "Arial") +
  theme(axis.text.y = element_text(size = 6),
        plot.title = element_text(size = 9.5, face = "bold", hjust = 0.5),
        plot.margin = margin(4, 2, 2, 2, "mm"))

top5 <- drug_res %>% slice(1:min(4, nrow(.))) %>% pull(drug)
d5 <- pred_df %>%
  select(sample, group, all_of(top5)) %>%
  pivot_longer(-c(sample, group), names_to = "drug", values_to = "ic50") %>%
  mutate(group = factor(group, c("Low", "High")),
         drug = factor(drug, top5))
pB <- ggplot(d5, aes(group, ic50, fill = group)) +
  geom_boxplot(width = 0.55, outlier.size = 0.4, show.legend = FALSE) +
  scale_fill_manual(values = cols_risk) +
  facet_wrap(~drug, scales = "free_y", nrow = 2) +
  labs(x = NULL, y = "Predicted IC50",
       title = "Top differential drugs (sensitivity)") +
  theme_prism(base_size = 8, base_family = "Arial") +
  theme(axis.text.x = element_text(size = 7),
        strip.text = element_text(size = 6.5, face = "bold"),
        axis.title.y = element_text(size = 7.5),
        plot.title = element_text(size = 9.5, face = "bold", hjust = 0.5),
        plot.margin = margin(4, 2, 2, 2, "mm"))

fig11 <- pA | pB + plot_layout(widths = c(1, 1.5))
ggsave(file.path(OUT, "Fig11_drug_sensitivity.png"), fig11,
       width = 210, height = 155, units = "mm", dpi = 300)
ggsave(file.path(OUT, "Fig11_drug_sensitivity.pdf"), fig11,
       width = 210, height = 155, units = "mm", device = cairo_pdf,
       )
message("saved: Fig11_drug_sensitivity")