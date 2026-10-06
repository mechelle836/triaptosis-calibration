# =====================================================================
# 20_fig9_mechanism.R —— Fig9「机制自洽」富集组图
# A: GO BP 富集气泡 (top10)  B: KEGG 条形 (msigdbr 离线)
# C: 基因-通路二分网络 (ggraph, 基因按模型 beta 符号着色)
# 输入: results/02_model/KIRC/model.rds (21 基因 + beta)
# =====================================================================
suppressPackageStartupMessages({
  library(clusterProfiler); library(org.Hs.eg.db); library(msigdbr)
  library(ggplot2); library(dplyr); library(tidyr); library(tibble)
  library(tidygraph); library(ggraph); library(showtext); library(patchwork); library(ggprism)
})
font_add("Arial", regular = "/System/Library/Fonts/Supplemental/Arial.ttf")
showtext_opts(dpi = 300)
showtext_auto(FALSE)   # macOS 系统 Arial 直渲

ROOT <- "/Users/xiaolin/WorkBuddy/2026-09-04-20-00-08/paper2_triaptosis"
OUT  <- file.path(ROOT, "figures_v2")
RES  <- file.path(ROOT, "results", "12_mechanism")
dir.create(RES, showWarnings = FALSE, recursive = TRUE)

mod <- readRDS(file.path(ROOT, "results/02_model/KIRC/model.rds"))
genes21 <- mod$genes
beta_map <- setNames(mod$beta, genes21)
tas10 <- mod$genes_sel

# ---- 1. SYMBOL -> ENTREZ -------------------------------------------
eg <- bitr(genes21, fromType = "SYMBOL", toType = "ENTREZID",
           OrgDb = org.Hs.eg.db)
cat("mapped:", nrow(eg), "/", length(genes21), "\n")

# ---- 2. GO BP 富集 --------------------------------------------------
ego <- enrichGO(gene = eg$ENTREZID, OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH",
                pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
go_df <- as.data.frame(ego)
cat("GO BP terms:", nrow(go_df), "\n")
go_top <- go_df %>% arrange(p.adjust) %>% slice(1:min(10, nrow(.))) %>%
  mutate(GeneRatio_n = sapply(strsplit(GeneRatio, "/"),
                              function(x) as.numeric(x[1]) / as.numeric(x[2])))
write.csv(go_df, file.path(RES, "go_bp_enrichment.csv"), row.names = FALSE)

# ---- 3. GO CC 富集 (离线; VPS34/PI3P = 内体定位, CC 比 KEGG 更贴题) ----
ecc <- enrichGO(gene = eg$ENTREZID, OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID", ont = "CC", pAdjustMethod = "BH",
                pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
cc_df <- as.data.frame(ecc) %>%
  mutate(Description = gsub("\\(GO:[0-9]+\\)", "", Description) %>%
           trimws())
cat("GO CC terms:", nrow(cc_df), "\n")
kegg_top <- cc_df %>% arrange(p.adjust) %>% slice(1:min(10, nrow(.))) %>%
  mutate(GeneRatio_n = sapply(strsplit(GeneRatio, "/"),
                              function(x) as.numeric(x[1]) / as.numeric(x[2])))
write.csv(cc_df, file.path(RES, "go_cc_enrichment.csv"), row.names = FALSE)

# ---- 4. Panel A/B: 气泡 + 条形 --------------------------------------
pA <- ggplot(go_top, aes(GeneRatio_n, reorder(Description, GeneRatio_n))) +
  geom_point(aes(size = Count, colour = -log10(p.adjust))) +
  scale_colour_gradient(low = "#BDD7E7", high = "#08519C",
                        name = expression(-log[10](p[adj])),
                        guide = "none") +
  scale_size(range = c(2.5, 6), name = "Count", guide = "none") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
  scale_y_discrete(labels = function(x) sapply(x, function(t) {
    t <- gsub("\\(GO:[0-9]+\\)", "", t)
    ifelse(nchar(t) > 52, paste0(substr(t, 1, 50), "..."), t)
  })) +
  labs(x = "Gene ratio", y = NULL,
       title = "GO biological process (21 triaptosis genes)") +
  theme_prism(base_size = 8.5, base_family = "Arial") +
  theme(axis.text.y = element_text(size = 7.5),
        axis.text.x = element_text(size = 6.5),
        legend.position = "none",
        plot.title = element_text(size = 10, face = "bold", hjust = 0))

pB <- ggplot(kegg_top, aes(reorder(Description, -log10(p.adjust)),
                           -log10(p.adjust))) +
  geom_col(aes(fill = GeneRatio_n), width = 0.62) +
  coord_flip() +
  scale_fill_gradient(low = "#CDEEE3", high = "#0F6E56",
                      name = "Gene ratio", guide = "none") +
  labs(x = NULL, y = expression(-log[10](p[adj])),
       title = "GO cellular component") +
  theme_prism(base_size = 8.5, base_family = "Arial") +
  theme(axis.text.y = element_text(size = 7.5),
        axis.text.x = element_text(size = 6.5),
        plot.title = element_text(size = 10, face = "bold", hjust = 0))

# ---- 5. Panel C: 基因-通路二分网络 -----------------------------------
# 用 GO+KEGG top 通路的 geneID 建边 (symbol 化)
edge_go <- go_top %>% select(Description, geneID) %>%
  mutate(Description = paste0("GO: ", Description))
edge_kegg <- kegg_top %>% select(Description, geneID) %>%
  mutate(Description = paste0("CC: ", Description))
trunc46 <- function(x) ifelse(nchar(x) > 46, paste0(substr(x, 1, 44), "..."), x)
edges <- bind_rows(edge_go, edge_kegg) %>%
  mutate(geneID = gsub("/", "; ", geneID),
         Description = trunc46(trimws(Description))) %>%
  separate_rows(geneID, sep = "; ") %>%
  transmute(from = trimws(geneID), to = trimws(Description)) %>%
  filter(from %in% genes21)
edge_w <- edges %>% count(from, to, name = "w")

nodes <- data.frame(
  name = c(unique(edges$from), unique(edges$to)),
  type = c(rep("gene", length(unique(edges$from))),
           rep("pathway", length(unique(edges$to))))) %>%
  distinct(name, .keep_all = TRUE) %>%
  mutate(beta = beta_map[name]) %>%
  mutate(role = case_when(type == "pathway" ~ "pathway",
                          beta > 0 ~ "risk",
                          beta < 0 ~ "protective",
                          TRUE ~ "neutral"))

# 手动两列布局
gn <- nodes %>% filter(type == "gene") %>%
  arrange(-abs(beta)) %>%
  mutate(x = 0, y = rev(seq_len(n())) - n()/2)
pn <- nodes %>% filter(type == "pathway") %>%
  arrange(name) %>%
  mutate(x = 2.2, y = rev(seq_len(n())) - n()/2)
node_pos <- bind_rows(gn, pn)

gr <- tbl_graph(nodes = node_pos, edges = edge_w, directed = FALSE)
pC <- ggraph(gr, layout = "manual", x = x, y = y) +
  geom_edge_diagonal(aes(width = w), colour = "grey72", alpha = 0.75,
                     edge_linewidth = 0.5) +
  scale_edge_width(range = c(0.2, 1.2), guide = "none") +
  geom_node_point(aes(colour = role, size = type == "gene")) +
  geom_node_text(aes(label = name, x = x - ifelse(type == "gene", 0.12, 0),
                     y = y),
                 hjust = ifelse(node_pos$type == "gene", 1, 0), size = 2.4,
                 fontface = ifelse(node_pos$type == "gene", "italic", "plain"),
                 family = "Arial") +
  scale_colour_manual(values = c(risk = "#E64B35", protective = "#3C5488",
                                 neutral = "grey55", pathway = "#00A087"),
                      name = NULL,
                      labels = c(risk = "Risk gene (beta>0)",
                                 protective = "Protective (beta<0)",
                                 neutral = "Beta = 0",
                                 pathway = "Enriched pathway")) +
  scale_size_manual(values = c(`TRUE` = 4.2, `FALSE` = 2.6), guide = "none") +
  expand_limits(x = c(-2.6, 7.2)) +
  labs(title = "Gene-pathway bipartite network",
       x = NULL, y = NULL) +
  theme_void(base_family = "Arial") +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 10, face = "bold", hjust = 0),
        legend.text = element_text(size = 7),
        plot.margin = margin(4, 4, 4, 4, "mm"))

# ---- 6. 组图输出 ------------------------------------------------------
fig9 <- (pA | pB) / pC + plot_annotation(tag_levels = "A") +
  plot_layout(heights = c(1.1, 1), widths = c(1.25, 1))
ggsave(file.path(OUT, "Fig9_mechanism.png"), fig9,
       width = 200, height = 210, units = "mm", dpi = 300)
showtext_auto()
ggsave(file.path(OUT, "Fig9_mechanism.pdf"), fig9,
       width = 200, height = 210, units = "mm", device = cairo_pdf,
       )
message("saved: Fig9_mechanism")