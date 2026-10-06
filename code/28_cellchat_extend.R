# =====================================================================
# 28_cellchat_extend.R
# 在已有圆环总图之外，补三张同类网络图（不重跑 CellChat）
#   Fig12C  TAS-high / TAS-low 肿瘤细胞发出的通讯强度
#   Fig12E  同一布局下的配体-受体条数
#   Fig12D  差异最大的 6 条通路（只保留两条肿瘤边）
# 输出: figures_v2/ 与 figures_pub/
# =====================================================================
suppressPackageStartupMessages({
  library(CellChat)
})
ROOT <- "/Users/xiaolin/WorkBuddy/2026-09-04-20-00-08/paper2_triaptosis"
DIRS <- c(file.path(ROOT, "figures_v2"), file.path(ROOT, "figures_pub"))
RES  <- file.path(ROOT, "results", "15_cellchat")
cc <- readRDS(file.path(RES, "cellchat_TAS.rds"))

label_map <- c(
  "B/Plasma/Mast" = "B/Mast",
  "Endothelial"   = "Endothelial",
  "Macrophages"   = "Macrophage",
  "T cells"       = "T cell",
  "Tumor_TAShigh" = "TAS-high",
  "Tumor_TASlow"  = "TAS-low",
  "vSMC/Pericyte" = "vSMC"
)
node_col <- c(
  "B/Mast"      = "#F39B7F",
  "Endothelial" = "#4DBBD5",
  "Macrophage"  = "#00A087",
  "T cell"      = "#8491B4",
  "TAS-high"    = "#E64B35",
  "TAS-low"     = "#3C5488",
  "vSMC"        = "#7E6148"
)

relabel <- function(m) {
  rn <- unname(label_map[rownames(m)])
  cn <- unname(label_map[colnames(m)])
  if (anyNA(rn) || anyNA(cn)) stop("unmapped cell group")
  dimnames(m) <- list(rn, cn)
  m
}

# 只保留一个肿瘤来源的发出边，去掉自环，其余节点留在圆上
outgoing_only <- function(m, source) {
  out <- m * 0
  out[source, ] <- m[source, ]
  diag(out) <- 0
  out
}

# 两条肿瘤来源都保留，去掉自环
tumor_outgoing <- function(m) {
  out <- m * 0
  src <- c("Tumor_TAShigh", "Tumor_TASlow")
  out[src, ] <- m[src, ]
  diag(out) <- 0
  out
}

draw_circle <- function(m, title, edge_max, label_cex = 0.85, edge_max_width = 8) {
  mm <- relabel(m)
  netVisual_circle(
    mm,
    color.use = node_col[rownames(mm)],
    weight.scale = TRUE,
    edge.weight.max = edge_max,
    edge.width.max = edge_max_width,
    vertex.weight = 1,
    vertex.size.max = 8,
    vertex.label.cex = label_cex,
    arrow.size = 0.18,
    arrow.width = 1,
    alpha.edge = 0.75,
    edge.curved = 0.25,
    margin = 0.28,
    title.name = title
  )
}

save_base <- function(draw, file, w_in, h_in) {
  for (dir in DIRS) {
    dir.create(dir, showWarnings = FALSE, recursive = TRUE)
    grDevices::png(file.path(dir, paste0(file, ".png")),
                   width = w_in, height = h_in, units = "in", res = 300,
                   family = "Arial", bg = "white")
    draw()
    grDevices::dev.off()
    if (exists("quartz", mode = "function")) {
      grDevices::quartz(file = file.path(dir, paste0(file, ".pdf")),
                        type = "pdf", width = w_in, height = h_in,
                        family = "Arial")
      draw()
      grDevices::dev.off()
    }
    message("saved: ", file.path(dir, file))
  }
}

# ---- Fig12C: 发出强度，左右同一标尺 ------------------------------------
w <- cc@net$weight
hi <- outgoing_only(w, "Tumor_TAShigh")
lo <- outgoing_only(w, "Tumor_TASlow")
edge_max_w <- max(hi, lo)
save_base(function() {
  par(mfrow = c(1, 2), mar = c(0.2, 0.2, 1.6, 0.2), family = "Arial", xpd = NA)
  draw_circle(hi, "TAS-high tumor outgoing", edge_max_w, label_cex = 0.95)
  draw_circle(lo, "TAS-low tumor outgoing", edge_max_w, label_cex = 0.95)
}, "Fig12C_outgoing_circles", w_in = 10.6, h_in = 5.6)

# ---- Fig12E: 配体-受体条数，同一布局 ------------------------------------
n <- cc@net$count
hi_n <- outgoing_only(n, "Tumor_TAShigh")
lo_n <- outgoing_only(n, "Tumor_TASlow")
edge_max_n <- max(hi_n, lo_n)
save_base(function() {
  par(mfrow = c(1, 2), mar = c(0.2, 0.2, 1.6, 0.2), family = "Arial", xpd = NA)
  draw_circle(hi_n, "TAS-high tumor  |  interaction count", edge_max_n,
              label_cex = 0.95, edge_max_width = 8)
  draw_circle(lo_n, "TAS-low tumor  |  interaction count", edge_max_n,
              label_cex = 0.95, edge_max_width = 8)
}, "Fig12E_outgoing_counts", w_in = 10.6, h_in = 5.6)

# ---- Fig12D: 通路圆环，只画肿瘤发出 ------------------------------------
# 按 |high - low| 发出总量排序后取前 6，MIF 是唯一明显反向的通路
pw <- cc@netP$pathways
pw_tab <- do.call(rbind, lapply(pw, function(p) {
  m <- cc@netP$prob[, , p]
  data.frame(pathway = p,
             high = sum(m["Tumor_TAShigh", ], na.rm = TRUE),
             low  = sum(m["Tumor_TASlow", ], na.rm = TRUE),
             stringsAsFactors = FALSE)
}))
pw_tab$delta <- pw_tab$high - pw_tab$low
pw_tab <- pw_tab[order(-abs(pw_tab$delta)), ]
write.csv(pw_tab, file.path(RES, "pathway_outgoing_delta.csv"), row.names = FALSE)
keep <- pw_tab$pathway[seq_len(6)]
message("pathways: ", paste(sprintf("%s (%.3f vs %.3f)", keep,
                                    pw_tab$high[match(keep, pw_tab$pathway)],
                                    pw_tab$low[match(keep, pw_tab$pathway)]),
                            collapse = "; "))

save_base(function() {
  par(mfrow = c(2, 3), mar = c(0.15, 0.15, 1.5, 0.15), family = "Arial", xpd = NA)
  for (p in keep) {
    raw <- cc@netP$prob[, , p]
    mm <- tumor_outgoing(raw)
    row <- pw_tab[pw_tab$pathway == p, ]
    ttl <- sprintf("%s    high %.3f   low %.3f", p, row$high, row$low)
    # 每条通路自己的最大值，避免弱通路被 APP 的标尺压成看不见
    draw_circle(mm, ttl, edge_max = max(mm), label_cex = 0.72, edge_max_width = 6)
  }
}, "Fig12D_pathway_circles", w_in = 11.2, h_in = 7.6)

message("done")
