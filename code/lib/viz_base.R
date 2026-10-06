# =====================================================================
# viz_base.R —— 生信出版级绘图基建（CNS 极简风格 v2）
# 覆盖 ggprism::theme_prism 为 CNS 极简主题, 所有脚本 source 后自动升级
# 配色: Nature Publishing Group (npg) 色板, 色盲安全, 避红绿对比
# 用法: source("<项目>/code/lib/viz_base.R") 后直接使用:
#   viz_palette / cols_risk / theme_publication() / theme_prism() /
#   theme_cns() / save_plot() / fmt_p()
# =====================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(showtext)
})

# ---------------------------------------------------------------------
# 1) 集中式语义调色板（NPG 顶刊色, 色盲安全）
# ---------------------------------------------------------------------
viz_palette <- list(
  risk_high = "#E64B35",   # 高风险 / 上调 / 不良预后 (npg red)
  risk_low  = "#3C5488",   # 低风险 / 下调 / 良好预后 (npg navy)
  accent    = "#00A087",   # 强调 (模型/显著, npg green)
  secondary = "#4DBBD5",   # 次要系列 (npg cyan)
  highlight = "#F39B7F",   # 高亮 (npg salmon)
  neutral   = "#8491B4",   # 中性 / 对照 / 零模型 (npg grey-purple)
  sig       = "#E64B35",   # 显著
  ns        = "#B4B2A9"    # 不显著
)
cols_risk <- c(High = viz_palette$risk_high, Low = viz_palette$risk_low)

# 风险组双色快捷 (兼容部分脚本直接引用单值)
COL_HIGH <- viz_palette$risk_high
COL_LOW  <- viz_palette$risk_low
COL_MID  <- viz_palette$secondary
COL_GREEN <- viz_palette$accent
COL_ORANGE <- viz_palette$highlight
COL_GREY <- viz_palette$neutral

# ---------------------------------------------------------------------
# 字体注册 (macOS Arial; 系统无该字体时静默回退 sans)
# ---------------------------------------------------------------------
.viz_register_font <- function(family = "Arial") {
  candidates <- c(
    "/System/Library/Fonts/Supplemental/Arial.ttf",
    "/Library/Fonts/Arial.ttf",
    "/usr/share/fonts/truetype/msttcorefonts/Arial.ttf")
  f <- candidates[file.exists(candidates)]
  if (length(f) > 0) {
    tryCatch({
      font_add(family, regular = f[1])
      showtext_opts(dpi = 300); showtext_auto()
    }, error = function(e) invisible(NULL))
  }
  invisible(NULL)
}
.viz_register_font()

# ---------------------------------------------------------------------
# 2) CNS 极简主题 (覆盖 ggprism::theme_prism)
#    特征: 纯白背景 / 无网格 / 细线 0.3 / 刻度朝外 / 统一字重 / 标题居中
# ---------------------------------------------------------------------
theme_cns <- function(base_size = 8, base_family = "Arial",
                      base_line_size = 0.3) {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      # 轴: 仅左/下, 细线, 刻度朝外
      axis.line = element_line(colour = "black", linewidth = base_line_size),
      axis.ticks = element_line(colour = "black", linewidth = base_line_size),
      axis.ticks.length = unit(1.4, "mm"),
      axis.text = element_text(size = base_size - 1, colour = "black"),
      axis.title = element_text(size = base_size, colour = "black"),
      # 背景: 纯白, 无网格
      panel.background = element_rect(fill = "white", colour = NA),
      panel.grid = element_blank(),
      plot.background = element_rect(fill = "white", colour = NA),
      # 图例: 无框无背景
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.title = element_text(size = base_size - 1),
      legend.text = element_text(size = base_size - 1),
      # 标题: 居中加粗
      plot.title = element_text(size = base_size + 0.5, face = "bold",
                                hjust = 0.5, colour = "black"),
      # facet strip: 无背景加粗
      strip.background = element_blank(),
      strip.text = element_text(size = base_size, face = "bold", colour = "black")
    )
}

# 覆盖 ggprism::theme_prism (兼容脚本原有调用签名)
theme_prism <- function(base_size = 8, base_family = "Arial",
                        base_line_size = 0.3, axis_text_angle = 0) {
  theme_cns(base_size = base_size, base_family = base_family,
            base_line_size = base_line_size) +
    theme(axis.text.x = element_text(size = base_size - 1, colour = "black",
                                     angle = axis_text_angle))
}

# theme_publication 保持向后兼容 (指向 CNS 主题)
theme_publication <- function(base_size = 8, base_family = "Arial",
                              base_line_size = 0.3) {
  theme_cns(base_size = base_size, base_family = base_family,
            base_line_size = base_line_size) +
    theme(plot.title = element_text(size = base_size, face = "bold",
                                    hjust = 0, colour = "black"),
          legend.position = "top")
}

# ---------------------------------------------------------------------
# 3) save_plot: 双输出 (PNG 300dpi + 矢量 PDF), 统一 & theme 兜底
# ---------------------------------------------------------------------
save_pub <- function(p, file, dir = ".", w = 180, h = 130, dpi = 300) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  png_path <- file.path(dir, paste0(file, ".png"))
  pdf_path <- file.path(dir, paste0(file, ".pdf"))
  showtext_auto(FALSE)
  ggsave(png_path, p, width = w, height = h, units = "mm", dpi = dpi)
  # cairo DLL is missing (no X11). Quartz PDF embeds Arial on macOS.
  showtext_auto(FALSE)
  wrote_pdf <- tryCatch({
    if (exists("quartz")) {
      ggsave(pdf_path, p, width = w, height = h, units = "mm", device = function(filename, width, height, ...) {
        grDevices::quartz(file = filename, type = "pdf", width = width, height = height)
      })
    } else {
      ggsave(pdf_path, p, width = w, height = h, units = "mm", device = grDevices::pdf)
    }
    file.exists(pdf_path) && file.info(pdf_path)$size > 20000
  }, error = function(e) {
    message("PDF skipped (", conditionMessage(e), ")")
    FALSE
  })
  message("saved: ", file, if (isTRUE(wrote_pdf)) " (png+pdf)" else " (png only)")
}

save_plot <- function(p, file, dir = ".", w = 180, h = 130, dpi = 300,
                      apply_theme = TRUE) {
  if (apply_theme) p <- p & theme_cns()
  showtext_auto(FALSE)   # PNG: 系统 Arial 直渲 (宽度准确)
  ggsave(file.path(dir, paste0(file, ".png")), p,
         width = w, height = h, units = "mm", dpi = dpi)
  showtext_auto()        # PDF: showtext 嵌入字体 (cairo_pdf 需要)
  ggsave(file.path(dir, paste0(file, ".pdf")), p,
         width = w, height = h, units = "mm",
         device = cairo_pdf)
  showtext_auto(FALSE)
  message("saved: ", file)
}

# ---------------------------------------------------------------------
# 4) 连续色条必须用 raster = FALSE
#    save_pub 的 PDF 走 grDevices::quartz(type = "pdf"), 该设备把色条的
#    栅格渐变上下翻转, 于是 PDF 里色条的方向和散点的取色相反 (PNG 正常)。
#    raster = FALSE 让色条改用一叠实心矩形绘制, 不经过栅格路径。
#    所有 scale_*_gradient* 都要传 guide = guide_cbar()。
# ---------------------------------------------------------------------
guide_cbar <- function(...) {
  # ggplot2 >= 3.5 用 display = "rectangles", 旧版只有 raster = FALSE
  if ("display" %in% names(formals(ggplot2::guide_colourbar))) {
    ggplot2::guide_colourbar(display = "rectangles", ...)
  } else {
    ggplot2::guide_colourbar(raster = FALSE, ...)
  }
}

# ---------------------------------------------------------------------
# 5) KM 图公共设置
#    at-risk 表首列(t=0)的数字以刻度为中心, 面板左侧必须留出半个字宽,
#    否则数字会压到分层标签上。x/y 尺度与 at-risk 字号在此统一。
# ---------------------------------------------------------------------
km_x_scale <- function(breaks = waiver(), ...) {
  scale_x_continuous(breaks = breaks,
                     expand = expansion(mult = c(0.055, 0.03)), ...)
}

# 生存概率统一用百分数, 全文各图一致
km_y_scale <- function() {
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25),
                     labels = function(x) paste0(x * 100, "%"),
                     expand = expansion(mult = c(0.02, 0.03)))
}

km_risktable_theme <- function(base_size = 7) {
  list(
    ggsurvfit::theme_risktable_default(
      axis.text.y.size = base_size - 0.5,
      plot.title.size = base_size),
    theme(plot.title = element_text(face = "plain", hjust = 0),
          axis.text.y = element_text(colour = "black", hjust = 1,
                                     margin = margin(r = 1.2, unit = "mm")))
  )
}

# ---------------------------------------------------------------------
# 常用小工具
# ---------------------------------------------------------------------
fmt_p <- function(p) {
  ifelse(p < 1e-4, formatC(p, format = "e", digits = 1), sprintf("%.3f", p))
}

# 期刊体例: P < 0.001 不给具体数值, 其余两位有效数字
fmt_p_star <- function(p) {
  ifelse(p < 0.001, "P < 0.001", sprintf("P = %.3f", p))
}

invisible(NULL)