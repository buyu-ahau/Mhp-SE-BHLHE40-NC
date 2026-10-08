# ============================================================================
# _palette_MHP.R — 全部作图脚本共用的细胞类型 / 分组配色、ggplot2 主题与存图函数（由其他脚本 source）
# 运行位置：本地；路径变量见 config/paths.R
# ============================================================================

source(file.path(Sys.getenv("MHP_CODE", "."), "config", "paths.R"))

# ---------- 1. scRNA 11 类：数据值 → 颜色 / 图中显示名（顺序 = 髓系 → 淋巴 → 结构细胞）----------
ct_rna_col <- c("Macrophages"       = "#f69886", "Neutrophils"       = "#d984b7", "DC cells"    = "#00a18a",
                "CD4 cells"         = "#2f5885", "CD8 cells"         = "#ef3f33", "NK cells"    = "#3c9dd7",
                "B cells"           = "#26a870", "Plasma cells"      = "#895da7",
                "Epithelial cells"  = "#d4c1ee", "Endothelial cells" = "#a8d490", "Fibroblasts" = "#fdc182")
ct_rna_order <- names(ct_rna_col)
# 图中显示名（全 ASCII；数据里的 celltype_new 值不改）
ct_rna_lab <- c("Macrophages" = "Macrophages", "Neutrophils" = "Neutrophils", "DC cells" = "Dendritic cells",
                "CD4 cells" = "CD4+ T cells", "CD8 cells" = "CD8+ T cells", "NK cells" = "\u03b3\u03b4 T cells",
                "B cells" = "B cells", "Plasma cells" = "Plasma cells", "Epithelial cells" = "Epithelial cells",
                "Endothelial cells" = "Endothelial cells", "Fibroblasts" = "Fibroblasts")

# ---------- 2. snATAC 8 类：颜色取自代表性 RNA 类型，保证同型同色 ----------
atac_color_src <- c("Macrophages" = "Macrophages", "Neutrophils" = "Neutrophils", "T_NK_cells" = "CD8 cells",
                    "B_cells" = "B cells", "Plasma_cells" = "Plasma cells", "Epithelial_cells" = "Epithelial cells",
                    "Endothelial_cells" = "Endothelial cells", "Fibroblasts" = "Fibroblasts")
ct_atac_col   <- setNames(unname(ct_rna_col[atac_color_src]), names(atac_color_src))
ct_atac_order <- names(ct_atac_col)
ct_atac_lab   <- c("Macrophages" = "Macrophages", "Neutrophils" = "Neutrophils", "T_NK_cells" = "T/NK cells",
                   "B_cells" = "B cells", "Plasma_cells" = "Plasma cells", "Epithelial_cells" = "Epithelial cells",
                   "Endothelial_cells" = "Endothelial cells", "Fibroblasts" = "Fibroblasts")
ct_atac_lab_legend <- replace(ct_atac_lab, "T_NK_cells", "T/NK cells (CD4+ T, CD8+ T, \u03b3\u03b4 T)")   # 图例用长名，图内标签用短名
atac_fine2main <- c("AT2_cells" = "Epithelial_cells", "Macrophages_Inf" = "Macrophages")

# ---------- 3. 两模态类型对应（11 → 8；DC 在 snATAC 无对应）----------
ct_map <- data.frame(
  rna  = ct_rna_order,
  atac = c("Macrophages", "Neutrophils", NA, "T_NK_cells", "T_NK_cells", "T_NK_cells",
           "B_cells", "Plasma_cells", "Epithelial_cells", "Endothelial_cells", "Fibroblasts"),
  stringsAsFactors = FALSE)

# ---------- 4. 分组与样本 ----------
grp_levels <- c("Con", "Mhp")
grp_col    <- c("Con" = "#9FB8CD", "Mhp" = "#D6604D")
grp_raw2std <- c("Day0" = "Con", "Day28" = "Mhp", "Mock" = "Con", "Mhp" = "Mhp", "Con" = "Con")
to_grp <- function(x) {                                   # RNA 的 Day0/Day28 与 ATAC 的 Mock/Mhp 统一成 Con/Mhp
  out <- unname(grp_raw2std[as.character(x)])
  if (anyNA(out)) stop("未知分组值: ", paste(unique(x[is.na(out)]), collapse = ", "))
  factor(out, levels = grp_levels) }
to_smp <- function(x) {                                   # Day0_1 / Mock_1 → Con_1；Day28_1 → Mhp_1
  out <- sub("^(Day0|Mock)_", "Con_", sub("^Day28_", "Mhp_", as.character(x)))
  if (!all(grepl("^(Con|Mhp)_[0-9]+$", out))) stop("未知样本名: ", paste(unique(x[!grepl("^(Con|Mhp)_[0-9]+$", out)]), collapse = ", "))
  factor(out, levels = paste0(rep(grp_levels, each = 3), "_", 1:3)) }
# ★两模态的样本编号是否为同一头猪未核实，不得按编号跨模态配对

# 样本色：Con 三只为蓝色系深浅、Mhp 三只为红色系深浅（Supp Fig 1 按样本着色、逐动物柱用）
smp_col <- c("Con_1" = "#6F93B3", "Con_2" = "#9FB8CD", "Con_3" = "#C9D8E5",
             "Mhp_1" = "#B03A2E", "Mhp_2" = "#D6604D", "Mhp_3" = "#EFA08F")

# ---------- 5. 显示名换算 ----------
lab_rna  <- function(x) unname(ct_rna_lab[as.character(x)])
lab_atac <- function(x, legend = FALSE) unname((if (legend) ct_atac_lab_legend else ct_atac_lab)[as.character(x)])

# ---------- 6. 版面与字号（Nature 系：单栏 88 mm、双栏 183 mm；图内字号 5–7 pt）----------
W1 <- 88; W2 <- 183                                       # mm
FS_BASE <- 7; FS_SMALL <- 6; FS_TITLE <- 8                # pt
pt2size <- function(pt) pt / ggplot2::.pt                 # geom_text 的 size 以 mm 计，pt 需换算

# Nature 风格：Arial、正文 7 pt / 刻度 6 pt、线宽以 pt 计（Nature 要求 ≥ 0.25 pt）、无网格、黑字
lw <- function(pt) pt / (ggplot2::.pt * 0.75)             # pt → ggplot linewidth（lwd 1 = 0.75 pt）
FONT <- "Arial"
theme_mhp <- function(base = FS_BASE) ggplot2::theme_classic(base_size = base, base_family = FONT) +
  ggplot2::theme(text = ggplot2::element_text(colour = "black"), axis.text = ggplot2::element_text(size = FS_SMALL, colour = "black"),
                 axis.line = ggplot2::element_line(linewidth = lw(0.5)), axis.ticks = ggplot2::element_line(linewidth = lw(0.5)),
                 axis.ticks.length = grid::unit(1, "mm"), plot.title = ggplot2::element_text(size = FS_TITLE, face = "bold", hjust = 0),
                 legend.title = ggplot2::element_text(size = FS_BASE), legend.text = ggplot2::element_text(size = FS_SMALL),
                 legend.key.size = grid::unit(3, "mm"), strip.background = ggplot2::element_blank(),
                 strip.text = ggplot2::element_text(size = FS_BASE), plot.background = ggplot2::element_rect(fill = "white", colour = NA))
theme_mhp_void <- function(base = FS_BASE) ggplot2::theme_void(base_size = base, base_family = FONT) +
  ggplot2::theme(text = ggplot2::element_text(colour = "black"), plot.title = ggplot2::element_text(size = FS_TITLE, face = "bold", hjust = 0),
                 legend.text = ggplot2::element_text(size = FS_SMALL), legend.key.size = grid::unit(3, "mm"),
                 plot.background = ggplot2::element_rect(fill = "white", colour = NA))

# ---------- 7. 路径 ----------
MHP_DIR <- list(
  rna_tf    = file.path(MHP_RES, "8、scRNA-seq/5.转录因子分析_v2_2026-09"),
  rna_prop  = file.path(MHP_RES, "8、scRNA-seq/Supp_Fig1_CellProportion"),
  atac      = file.path(MHP_RES, "9.scatac 补充分析"),
  fig1      = file.path(MHP_RES, "Figure/Fig1_v2_双模态图谱_2026-09"))   # 其下 Plots/<编号_名称>/、Tables/、_panels/

# ---------- 8. 存图 ----------
# PDF 需支持 Unicode（κ/ρ/Δ/−）：本机 Mac 未装 XQuartz，capabilities("cairo") 报 TRUE 但 cairo_pdf 只警告不出文件
# → macOS 用自带 quartz PDF 设备，其他系统用 cairo_pdf
.pdf_dev <- if (Sys.info()[["sysname"]] == "Darwin") {
  function(filename, width, height, ...) grDevices::quartz(type = "pdf", file = filename, width = width, height = height)
} else grDevices::cairo_pdf

.svg_dev <- function(filename, width, height, ...) svglite::svglite(filename, width = width, height = height, fix_text_size = FALSE)
# ggplot 对象：final = FALSE 只出 PNG 300 dpi 自检；final = TRUE 出 SVG（矢量，供 Illustrator）+ TIFF 900 dpi 终稿
# panel_dir 非空时同时 saveRDS，供组装脚本只摆放不重画
sav <- function(p, stem, dir, w_mm, h_mm, final = FALSE, panel_dir = NULL) {
  if (identical(Sys.getenv("MHP_FINAL"), "1")) final <- TRUE     # 全局开关：MHP_FINAL=1 Rscript xx.R 时各脚本不改代码即出 SVG+TIFF 终稿
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  f <- file.path(dir, stem)
  if (!final) {
    out <- paste0(f, ".png")
    ggplot2::ggsave(out, p, width = w_mm, height = h_mm, units = "mm", dpi = 300, bg = "white", device = ragg::agg_png)
  } else {
    out <- paste0(f, c(".svg", ".tiff"))
    ggplot2::ggsave(out[1], p, width = w_mm, height = h_mm, units = "mm", device = .svg_dev)
    ggplot2::ggsave(out[2], p, width = w_mm, height = h_mm, units = "mm", dpi = 900, bg = "white",
                    device = ragg::agg_tiff, compression = "lzw")
  }
  stopifnot(all(file.exists(out)))
  if (!is.null(panel_dir)) {
    dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)
    saveRDS(p, file.path(panel_dir, paste0(stem, ".rds")))
  }
  cat("✓", paste(basename(out), collapse = " + "), sprintf("(%g × %g mm)\n", w_mm, h_mm))
  invisible(out) }

# 非 ggplot（ComplexHeatmap 等）：传入一个无参绘图函数，如 function() draw(ht)
sav_draw <- function(draw_fun, stem, dir, w_mm, h_mm, final = FALSE) {
  if (identical(Sys.getenv("MHP_FINAL"), "1")) final <- TRUE
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  f <- file.path(dir, stem)
  run <- function(open_dev) { open_dev(); on.exit(grDevices::dev.off()); draw_fun() }
  if (!final) {
    out <- paste0(f, ".png")
    run(function() ragg::agg_png(out, width = w_mm, height = h_mm, units = "mm", res = 300, background = "white"))
  } else {
    out <- paste0(f, c(".svg", ".tiff"))
    run(function() .svg_dev(out[1], width = w_mm / 25.4, height = h_mm / 25.4))
    run(function() ragg::agg_tiff(out[2], width = w_mm, height = h_mm, units = "mm", res = 900,
                                  background = "white", compression = "lzw"))
  }
  stopifnot(all(file.exists(out)))
  cat("✓", paste(basename(out), collapse = " + "), sprintf("(%g × %g mm)\n", w_mm, h_mm))
  invisible(out) }

# ---------- 9. 自洽检查（只检查本文件内部定义，不读数据）----------
stopifnot(
  length(ct_rna_col) == 11, identical(names(ct_rna_lab), ct_rna_order),
  length(ct_atac_col) == 8,  identical(names(ct_atac_lab), ct_atac_order),
  all(na.omit(ct_map$atac) %in% ct_atac_order), all(ct_atac_order %in% ct_map$atac),
  all(atac_fine2main %in% ct_atac_order),
  ct_atac_col[["T_NK_cells"]] == "#ef3f33", length(smp_col) == 6)
