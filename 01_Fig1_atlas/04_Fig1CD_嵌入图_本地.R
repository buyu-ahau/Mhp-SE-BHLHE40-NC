# ============================================================================
# 04_Fig1CD_嵌入图_本地.R — Fig. 1b、1d：scRNA-seq 与 snATAC-seq 的 t-SNE 嵌入图（文件名中的面板字母为旧编号）
# 运行位置：本地；路径变量见 config/paths.R
# ============================================================================

source(file.path(Sys.getenv("MHP_CODE", "."), "config", "paths.R"))

source(file.path(MHP_CODE, "shared", "_palette_MHP.R"))
suppressPackageStartupMessages({ library(ggplot2); library(data.table); library(patchwork); library(ggrepel); library(grid) })

# ---------- 路径 ----------
T30 <- file.path(MHP_DIR$fig1, "Tables/30_Fig1v2_RNA导出")
T31 <- file.path(MHP_DIR$fig1, "Tables/31_Fig1v2_ATAC导出")
OUT <- file.path(MHP_DIR$fig1, "Plots/01_嵌入图")
PAN <- file.path(MHP_DIR$fig1, "_panels")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE); dir.create(PAN, recursive = TRUE, showWarnings = FALSE)

# ---------- 读数据（fread 读 csv.gz 需 R.utils，缺则退回 read.csv(gzfile())）----------
rd <- function(f) if (requireNamespace("R.utils", quietly = TRUE)) fread(f) else as.data.table(read.csv(gzfile(f), stringsAsFactors = FALSE))
rna0 <- rd(file.path(T30, "RNA_embed.csv.gz"))
rna <- data.table(x = rna0$umap_1, y = rna0$umap_2, grp = to_grp(rna0$Groups), smp = to_smp(rna0$orig.ident),
                  ct = factor(rna0$celltype_new, levels = ct_rna_order))
stopifnot(isTRUE(all.equal(rna0$tsne_1, rna0$tsneH_1)), isTRUE(all.equal(rna0$tsne_2, rna0$tsneH_2)))
rna_t <- copy(rna)[, `:=`(x = rna0$tsne_1, y = rna0$tsne_2)]
rm(rna0)
atac <- rd(file.path(T31, "ATAC_embeddings_all.csv.gz"))
atac <- data.table(x = atac$TSNE_1, y = atac$TSNE_2, grp = to_grp(atac$group),
                   ct = factor(atac$celltype, levels = ct_atac_order))
stopifnot(!anyNA(rna$ct), !anyNA(atac$ct), nrow(rna) == 35330, nrow(atac) == 9610)
cat("RNA:", nrow(rna), "细胞 |"); print(table(rna$grp))
cat("ATAC:", nrow(atac), "核 |"); print(table(atac$grp))

# ---------- 本脚本内的小工具（不动 _palette_MHP.R）----------
BG_COL  <- "grey88"                                     # 两联图的灰底（全部细胞）
PT_MULT_MAIN <- 0.8; ALPHA_MAIN <- 1                     # 单图（scRNA / snATAC t-SNE 细胞类型图）点径倍数与不透明度
PT_MULT_PAIR <- 0.8; ALPHA_PAIR <- 1                     # 两联图同上（点变小时须不透明，否则颜色被冲淡）
FS_LAB  <- 5.5; FS_AXIS <- 5                             # 类型标签 / L 轴轴名 字号（pt）；6 pt 在 60 mm 内 11 个标签互相挤，5.5 更合适
LAB_FACE <- "plain"; LAB_COL <- "black"                  # 标签字重/颜色（图谱 01 用 bold + grey15；5.5 pt 下 bold 偏重，审图对照后定 plain 黑）
LAB_HALO <- FALSE
AVOID_N <- 0;   LAB_PAD <- 0.15                          # 标签避让障碍点数 / 标签框留白（lines）
# ★AVOID_N 取 500/2000 实测过（标签被推到团旁 + 引线）：引线交叉、标签落到别的团上，比直接压在密度峰上更乱 → 定 0
TITLE_H <- 3.5                                          # 标题行占高（mm，7 pt × 1.2 行距 + 1 pt 下边距），只用于估算面板实际宽度
# 手动偏移（单位 = 各轴范围的比例；加在密度峰上）：只给"标签压在团上会把小团遮没"的类型，其余仍在峰上
LAB_ADJ <- list(
  rna  = list("NK cells"         = c(-0.10,  0.00),
              "Epithelial cells" = c( 0.00, -0.07),     # 624 细胞细长条：标签放到条下方（上方贴中性粒）
              "CD4 cells"        = c( 0.00,  0.06)),    # 1,838 细胞小团：标签上移到团顶，不盖住深蓝团
  atac = list("B_cells"          = c( 0.13, -0.01),     # 135 核紧贴浆细胞：标签放到团右侧空白处
              "T_NK_cells"       = c( 0.06,  0.01)),    # 密度峰偏团左，右移居中（纯美观）
  rna_tsne = list("Endothelial cells" = c( 0.00,  0.15),   # 团贴左缘：标签放团上方（与上皮碎片之间的空隙），压点 0
                  "Epithelial cells"  = c( 0.00,  0.08),   # 主碎片在左上：标签放碎片上方，压点 0
                  "NK cells"          = c( 0.00, -0.10),   # 左右分别贴内皮/CD8：标签放团正下方（与 B 团之间），压点 15（★键 = 数据值，显示 "γδ T cells"）
                  "B cells"           = c(-0.18, -0.02),   # 标签放团左侧空白（左缘与团之间），压点 0
                  "Fibroblasts"       = c( 0.16,  0.06),   # 团贴上缘：标签放团右上空白，压点 0
                  "DC cells"          = c( 0.00, -0.02),   # 上有巨噬、下有 CD4、右是边缘：只能压自身团中部（自身 257/其他 25）
                  "CD4 cells"         = c( 0.10,  0.00)))  # 上 CD8、下 CD8 亚团：标签放团右端并伸入右侧空白（其他类型 1 点）
pt_of   <- function(n) round(45 / sqrt(n), 2)            # 点径随细胞数自适应：3.5 万 → 0.24，1 万 → 0.46
fmt_n   <- function(n) format(n, big.mark = ",")
peak2d  <- function(x, y, n = 128) {                     # 二维核密度峰（标签锚点）
  k <- MASS::kde2d(x, y, n = n); i <- which(k$z == max(k$z), arr.ind = TRUE)[1, ]
  c(k$x[i[1]], k$y[i[2]]) }
# 实测字符串宽度（mm）：临时开一个 ragg 设备用 Arial 度量，供标签夹紧与图例列宽用
str_w_mm <- function(s, pt, face = "plain") {
  tf <- tempfile(fileext = ".png"); ragg::agg_png(tf, width = 50, height = 50, units = "mm", res = 72)
  on.exit({ grDevices::dev.off(); unlink(tf) })
  vapply(s, function(z) grid::convertWidth(grid::grobWidth(grid::textGrob(z, gp = grid::gpar(
    fontsize = pt, fontfamily = FONT, fontface = face))), "mm", valueOnly = TRUE), numeric(1)) }
# 面板范围：左、下各留 9% 放 L 轴，右、上留 2%
lim_of <- function(xr, yr) list(x = c(xr[1] - diff(xr) * 0.09, xr[2] + diff(xr) * 0.02),
                                y = c(yr[1] - diff(yr) * 0.09, yr[2] + diff(yr) * 0.02))
# L 型小箭头坐标轴：原点在点云左下角外 2%，箭头长 15% 范围；位置量用 annotate 固化
add_L_axis <- function(p, xr, yr, axis, frac = 0.15) {
  xs <- diff(xr); ys <- diff(yr); ox <- xr[1] - xs * 0.02; oy <- yr[1] - ys * 0.02
  ar <- arrow(length = unit(0.8, "mm"), type = "closed")
  p + annotate("segment", x = ox, y = oy, xend = ox + xs * frac, yend = oy, arrow = ar, linewidth = lw(0.5), colour = "black") +
      annotate("segment", x = ox, y = oy, xend = ox, yend = oy + ys * frac, arrow = ar, linewidth = lw(0.5), colour = "black") +
      annotate("text", x = ox + xs * frac / 2, y = oy - ys * 0.025, label = axis[1], vjust = 1,
               size = pt2size(FS_AXIS), family = FONT, colour = "black") +
      annotate("text", x = ox - xs * 0.025, y = oy + ys * frac / 2, label = axis[2], angle = 90, vjust = 0,
               size = pt2size(FS_AXIS), family = FONT, colour = "black") }
th_emb <- theme_mhp_void() + theme(plot.title = element_text(size = FS_BASE, face = "bold", hjust = 0, margin = margin(0, 0, 1, 0)),
                                   legend.position = "none", plot.margin = margin(1, 1, 1, 1, "mm"))

# ---------- 通用嵌入图（C 与 D 共用；两联图的每一联也走这里）----------
# d: 含 x, y, ct；col/lab: 命名色向量与显示名；axis: 轴名对；bg: 灰底用的全体坐标（NULL 则不画）；
# n_pt: 点径依据的细胞数（两联图传全体 n）；lim: 坐标范围（两联图传全体范围）；label: 是否标类型名；show_axis: 是否画 L 轴
# avoid_n > 0 时抽 avoid_n 个细胞作空标签"障碍点"，让标签避开点云落到团旁空白处（ggrepel 官方技巧）；0 = 标签直接压在密度峰上
# adj: 手动偏移表（LAB_ADJ$rna / $atac）；w_mm/h_mm: 该面板的存图尺寸，只用来估算"1 mm = 多少数据单位"以夹紧标签
mk_emb <- function(d, col, lab, title, axis, bg = NULL, n_pt = nrow(d), lim = NULL, label = TRUE, show_axis = TRUE, seed = 1,
                   avoid_n = AVOID_N, fs_lab = FS_LAB, pad = LAB_PAD, adj = NULL, w_mm = 60, h_mm = 60, pt_mult = PT_MULT_MAIN, bg_col = BG_COL, alpha_pt = ALPHA_MAIN) {
  xr <- if (is.null(lim)) range(d$x) else lim$x; yr <- if (is.null(lim)) range(d$y) else lim$y
  ps <- pt_of(n_pt) * pt_mult                                   # pt_mult：点径倍数（两联图用更小的点）
  set.seed(seed); ds <- d[sample(nrow(d)), ]                     # ★打乱绘制顺序
  L <- lim_of(xr, yr)
  p <- ggplot()
  if (!is.null(bg)) p <- p + geom_point(data = bg, aes(x, y), colour = bg_col, size = ps, shape = 16, stroke = 0)
  p <- p + geom_point(data = ds, aes(x, y, colour = ct), size = ps, shape = 16, stroke = 0, alpha = alpha_pt) +
    scale_colour_manual(values = col, drop = FALSE)
  if (label) {                                                   # 类型名标在各自密度峰，白色描边防压色
    lb <- d[, { pk <- peak2d(x, y); .(x = pk[1], y = pk[2]) }, by = ct][order(ct)]
    lb[, lab := unname(lab[as.character(ct)])]
    for (k in names(adj)) if (k %in% lb$ct)                      # 手动偏移（范围比例 → 数据单位）
      lb[ct == k, `:=`(x = x + adj[[k]][1] * diff(xr), y = y + adj[[k]][2] * diff(yr))]
    # 夹紧到面板内：coord_fixed 下面板实际宽 = min(可用宽, 可用高 / 纵横比)，由此换算 1 mm 对应的数据单位
    asp <- diff(L$y) / diff(L$x); pw <- min(w_mm - 2, (h_mm - 2 - TITLE_H) / asp); upm <- diff(L$x) / pw
    hw <- str_w_mm(lb$lab, fs_lab, LAB_FACE) / 2 * upm; hh <- fs_lab / 72 * 25.4 / 2 * upm   # 半宽 / 半高（数据单位）
    lb[, `:=`(x = pmin(pmax(x, L$x[1] + hw), L$x[2] - hw), y = pmin(pmax(y, L$y[1] + hh), L$y[2] - hh))]
    if (avoid_n > 0) { set.seed(seed); lb <- rbind(lb, d[sample(nrow(d), min(avoid_n, nrow(d))), .(ct, x, y, lab = "")]) }
    # ★point.size = NA：关闭 ggrepel 对"被标注点"的排斥（默认 point.size = 1 会把标签推离密度峰）
    p <- p + geom_text_repel(data = lb, aes(x, y, label = lab), inherit.aes = FALSE,
                             size = pt2size(fs_lab), family = FONT, fontface = LAB_FACE, colour = LAB_COL,
                             bg.color = if (LAB_HALO) "white" else NA, bg.r = 0.12, seed = 1, max.overlaps = Inf, max.iter = 20000, max.time = 2,
                             box.padding = pad, point.padding = 0, point.size = NA, min.segment.length = 0.3,
                             segment.size = 0.2, segment.colour = "grey40") }
  if (show_axis) p <- add_L_axis(p, xr, yr, axis)
  p + coord_fixed(xlim = L$x, ylim = L$y, expand = FALSE, clip = "off") + labs(title = title) + th_emb }

# 分组两联：灰底 = 全体细胞，彩色 = 本组；两联同坐标范围、同点径；L 轴只画左联
mk_pair <- function(d, col, axis, unit, pt_mult = PT_MULT_PAIR, bg_col = BG_COL, alpha_pt = ALPHA_PAIR) {
  lim <- list(x = range(d$x), y = range(d$y)); bg <- d[, .(x, y)]
  pl <- lapply(grp_levels, function(g) {
    s <- d[grp == g]
    mk_emb(s, col, lab = NULL, title = sprintf("%s (n = %s %s)", g, fmt_n(nrow(s)), unit), axis = axis,
           bg = bg, n_pt = nrow(d), lim = lim, label = FALSE, show_axis = (g == grp_levels[1]), pt_mult = pt_mult, bg_col = bg_col, alpha_pt = alpha_pt) })
  wrap_plots(pl, nrow = 1) }

mk_legend <- function(ncol = 1, UNIT = 2.6, byrow = ncol > 1, gap = 1.6) {
  ent <- data.table(ct = ct_rna_order, lab = lab_rna(ct_rna_order), col = unname(ct_rna_col))
  ent[ct == "DC cells", lab := paste(lab, "(scRNA-seq only)")]
  ent <- rbind(ent[1:6], data.table(ct = "T_NK_cells", lab = paste0(sub("\\)$", "", lab_atac("T_NK_cells", legend = TRUE)), "; snATAC-seq)"),
                                    col = ct_atac_col[["T_NK_cells"]]), ent[7:.N])
  n <- nrow(ent); nrow_ <- ceiling(n / ncol); ent[, i := .I]
  if (byrow) ent[, `:=`(row = (i - 1) %/% ncol + 1, cid = (i - 1) %% ncol + 1)] else ent[, `:=`(row = (i - 1) %% nrow_ + 1, cid = (i - 1) %/% nrow_ + 1)]
  ent[, w := str_w_mm(lab, FS_SMALL) / UNIT]                     # 文字宽度（数据单位）
  cw <- ent[, .(w = max(w) + 0.6 + gap), by = cid][order(cid)]   # 每列宽 = 该列最长文字 + 点到文字间距 + 列间距
  ent[, colx := c(0, cumsum(cw$w))[cid]]
  ent[, ct := factor(ct, levels = ct)]
  xmax <- ent[, max(colx + 0.6 + w)] + 0.2; ymin <- -(nrow_ + 0.6)
  p <- ggplot(ent) +
    geom_point(aes(colx, -row, colour = ct), size = 1.8, shape = 16, stroke = 0) +
    geom_text(aes(colx + 0.6, -row, label = lab), hjust = 0, size = pt2size(FS_SMALL), family = FONT, colour = "black") +
    scale_colour_manual(values = setNames(ent$col, levels(ent$ct)), guide = "none") +
    coord_fixed(xlim = c(-0.4, xmax), ylim = c(ymin, -0.4), expand = FALSE, clip = "off") +
    theme_mhp_void() + theme(plot.margin = margin(0, 0, 0, 0))
  attr(p, "wh") <- c((xmax + 0.4) * UNIT, (-0.4 - ymin) * UNIT)
  p }

EXTRA <- FALSE
AX_U <- c("UMAP 1", "UMAP 2"); AX_T <- c("t-SNE 1", "t-SNE 2")
pC  <- mk_emb(rna,  ct_rna_col,  ct_rna_lab,  sprintf("scRNA-seq (n = %s cells)",   fmt_n(nrow(rna))),  AX_U, adj = LAB_ADJ$rna)
pD  <- mk_emb(atac, ct_atac_col, ct_atac_lab, sprintf("snATAC-seq (n = %s nuclei)", fmt_n(nrow(atac))), AX_T, adj = LAB_ADJ$atac)
pC2 <- mk_pair(rna,  ct_rna_col,  AX_U, "cells")
pS  <- mk_pair(atac, ct_atac_col, AX_T, "nuclei")
lgV <- mk_legend(1); lgH <- mk_legend(4)

if (EXTRA) sav(pC,  "Fig1C_scRNA_UMAP_细胞类型",       OUT, 60, 60, panel_dir = PAN)
sav(pD,  "Fig1D_snATAC_tSNE_细胞类型",      OUT, 60, 60, panel_dir = PAN)
if (EXTRA) sav(pC2, "Fig1C2_scRNA_UMAP_分组Con_Mhp",   OUT, 90, 45, panel_dir = PAN)
sav(pS,  "SuppFig_snATAC_tSNE_分组Con_Mhp", OUT, 90, 45, panel_dir = PAN)
whV <- attr(lgV, "wh"); whH <- attr(lgH, "wh")
if (EXTRA) sav(lgV, "Fig1_共用图例_细胞类型_竖排", OUT, whV[1], whV[2], panel_dir = PAN)
if (EXTRA) sav(lgH, "Fig1_共用图例_细胞类型_横排", OUT, whH[1], whH[2], panel_dir = PAN)

# 自检拼版预览（不入 _panels）：上行 C | D | 竖排图例，下行 C' | S-ATAC
if (EXTRA) {
  prev <- (pC | pD | lgV) + plot_layout(widths = c(60, 60, whV[1]))
  prev <- wrap_elements(full = prev) / wrap_elements(full = (pC2 | pS)) + plot_layout(heights = c(60, 45))
  sav(prev, "_预览_Fig1CD_嵌入图拼版", OUT, W2, 108)
}

pCt  <- mk_emb(rna_t, ct_rna_col, ct_rna_lab, sprintf("scRNA-seq (n = %s cells)", fmt_n(nrow(rna_t))), AX_T, adj = LAB_ADJ$rna_tsne)
# 2) 分组两联的 t-SNE 版
pC2t <- mk_pair(rna_t, ct_rna_col, AX_T, "cells")
sav(pCt,  "Fig1C_scRNA_tSNE_细胞类型",     OUT, 60, 60, panel_dir = PAN)
sav(pC2t, "Fig1C2_scRNA_tSNE_分组Con_Mhp", OUT, 90, 45, panel_dir = PAN)

# 3) 诊断：按样本着色，看 t-SNE 是否算在 Harmony 空间
#    判据：若同组三只动物在同一细胞类型内各自成团（尤其巨噬/中性粒等大团出现三条平行"带"），说明是在未校正 PCA 上算的；
#    若六个样本在每个团内充分混杂、仅 Con/Mhp 有生物学偏移，则为 Harmony 后。下行放 UMAP（已知 Harmony 后）作同版式对照
mk_smp_row <- function(d, axis, tag) {
  ds  <- copy(d)[, ct := smp]                                  # 借 mk_emb 的 colour = ct 通道画样本色
  lim <- list(x = range(d$x), y = range(d$y)); bg <- d[, .(x, y)]
  ov  <- mk_emb(ds, smp_col, lab = NULL, title = sprintf("%s: 6 samples overlaid", tag), axis = axis, label = FALSE)
  hl  <- lapply(levels(d$smp), function(s) {
    x <- ds[smp == s]
    mk_emb(x, smp_col, lab = NULL, title = sprintf("%s (n = %s)", s, fmt_n(nrow(x))), axis = axis,
           bg = bg, n_pt = nrow(d), lim = lim, label = FALSE, show_axis = FALSE, w_mm = 39, h_mm = 32) })
  (ov | wrap_elements(full = wrap_plots(hl, nrow = 2, byrow = TRUE))) + plot_layout(widths = c(65, 117)) }
if (EXTRA) {
  diag <- wrap_elements(full = mk_smp_row(rna_t, AX_T, "t-SNE")) / wrap_elements(full = mk_smp_row(rna, AX_U, "UMAP (Harmony)")) +
    plot_layout(heights = c(65, 65))
  sav(diag, "_诊断_scRNA_tSNE_按样本", OUT, W2, 132, panel_dir = PAN)
}

# 4) 对比：一行三联 RNA UMAP | RNA t-SNE | ATAC t-SNE（标题写明降维方法，其余与 C/D 面板完全一致）
if (EXTRA) {
cmp <- (mk_emb(rna,   ct_rna_col,  ct_rna_lab,  "scRNA-seq UMAP",   AX_U, adj = LAB_ADJ$rna) |
        mk_emb(rna_t, ct_rna_col,  ct_rna_lab,  "scRNA-seq t-SNE",  AX_T, adj = LAB_ADJ$rna_tsne) |
        mk_emb(atac,  ct_atac_col, ct_atac_lab, "snATAC-seq t-SNE", AX_T, adj = LAB_ADJ$atac)) + plot_layout(widths = c(60, 60, 60))
sav(cmp, "_对比_RNA_UMAP_vs_tSNE_与ATAC", OUT, 180, 60, panel_dir = PAN)
}

# ---------- 关键数字 ----------
cat("\n实际点径(mm): RNA", pt_of(nrow(rna)) * PT_MULT_MAIN, "| ATAC", pt_of(nrow(atac)) * PT_MULT_MAIN, "| 两联图", pt_of(nrow(rna)) * PT_MULT_PAIR, "\n")
cat("RNA 各类型:\n"); print(table(rna$ct))
cat("ATAC 各类型 × 分组:\n"); print(table(atac$ct, atac$grp))
mix_ratio <- function(d) { cen <- d[, .(cx = median(x), cy = median(y)), by = .(ct, smp)]
  sp <- d[, .(s = sqrt(var(x) + var(y))), by = ct]
  r <- cen[, .(mp = mean(as.numeric(dist(cbind(cx, cy))))), by = ct]; merge(r, sp, by = "ct")[, .(ct, ratio = round(mp / s, 2))] }
#   基线：Harmony 校正前 UMAP（03 脚本 RunUMAP(reduction = "pca")，07 脚本补充图 1C 所用）——若 t-SNE 的比值接近 Harmony 后 UMAP
#   而远小于校正前基线，则 t-SNE 是在 Harmony 空间上算的
pre <- rd(file.path(T30, "RNA_cellQC_umapPre.csv.gz"))
pre <- data.table(x = pre$umapPre_1, y = pre$umapPre_2, ct = factor(pre$celltype_new, levels = ct_rna_order), smp = to_smp(pre$orig.ident))
mr <- Reduce(function(a, b) merge(a, b, by = "ct"), list(setnames(mix_ratio(pre), "ratio", "ratio_umapPre"),
             setnames(mix_ratio(rna), "ratio", "ratio_umapHarmony"), setnames(mix_ratio(rna_t), "ratio", "ratio_tsne")))
cat("样本质心分离比（校正前 UMAP | Harmony UMAP | t-SNE）:\n"); print(mr[order(ct)])
cat(sprintf("中位数: 校正前 %.2f | Harmony UMAP %.2f | t-SNE %.2f\n", median(mr$ratio_umapPre), median(mr$ratio_umapHarmony), median(mr$ratio_tsne)))
cat("\n✅ 04 完成 →", OUT, "\n")
