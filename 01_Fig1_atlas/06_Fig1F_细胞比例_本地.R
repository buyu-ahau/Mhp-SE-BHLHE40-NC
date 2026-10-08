# ============================================================================
# 06_Fig1F_细胞比例_本地.R — Fig. 1g：Con 与 Mhp 组的细胞类型比例（文件名中的面板字母为旧编号）
# 运行位置：本地；路径变量见 config/paths.R
# ============================================================================

source(file.path(Sys.getenv("MHP_CODE", "."), "config", "paths.R"))

source(file.path(MHP_CODE, "shared", "_palette_MHP.R"))
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(ggrepel) })

T_RNA      <- file.path(MHP_DIR$fig1, "Tables/30_Fig1v2_RNA导出/RNA_embed.csv.gz")
T_ATAC     <- file.path(MHP_DIR$fig1, "Tables/31_Fig1v2_ATAC导出/ATAC_embeddings_all.csv.gz")
T_SUPP     <- file.path(MHP_DIR$rna_prop, "Supp_cellcount_table.csv")
T_ATAC_OLD <- file.path(MHP_DIR$atac, "Plots/10_CellProportion/CellProportion_byGroup.csv")
OUT   <- file.path(MHP_DIR$fig1, "Plots/03_细胞比例")
PANEL <- file.path(MHP_DIR$fig1, "_panels")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. 计数：按组（合并口径）与按样本（逐动物口径）
# ============================================================
rna <- fread(T_RNA, select = c("cell", "orig.ident", "Groups", "celltype_new"))
rna[, `:=`(grp = to_grp(Groups), smp = to_smp(orig.ident), ct = factor(celltype_new, levels = ct_rna_order))]
stopifnot(nrow(rna) == 35330, !anyNA(rna$ct))
atac <- fread(T_ATAC, select = c("cell", "Sample", "group", "celltype"))
atac[, `:=`(grp = to_grp(group), smp = to_smp(Sample), ct = factor(celltype, levels = ct_atac_order))]
stopifnot(nrow(atac) == 9610, !anyNA(atac$ct))

# 按 by（"grp" 或 "smp"）统计每类细胞数与占比；缺失组合补 0，保证每组都有全部类型
prop_by <- function(d, by) {
  cnt  <- d[, .(N = .N), by = c(by, "ct")]
  full <- CJ(k = levels(d[[by]]), ct = levels(d$ct)); setnames(full, "k", by)
  cnt  <- merge(full, cnt, by = c(by, "ct"), all.x = TRUE)
  cnt[is.na(N), N := 0L]
  cnt[, n_total := sum(N), by = by][, perc := 100 * N / n_total]
  cnt[, (by) := factor(get(by), levels = levels(d[[by]]))][, ct := factor(ct, levels = levels(d$ct))]
  setorderv(cnt, c(by, "ct")); cnt[] }
rna_g <- prop_by(rna, "grp");  rna_s <- prop_by(rna, "smp")
atac_g <- prop_by(atac, "grp"); atac_s <- prop_by(atac, "smp")

supp <- fread(T_SUPP)
chk  <- merge(rna_g[, .(ct = as.character(ct), grp = as.character(grp), N)],
              supp[, .(ct = CellType, grp = Group, Total)], by = c("ct", "grp"))
stopifnot(nrow(chk) == 22, all(chk$N == chk$Total))
old  <- fread(T_ATAC_OLD)[, .(ct = celltype, grp = as.character(to_grp(group)), count)]
chk2 <- merge(atac_g[, .(ct = as.character(ct), grp = as.character(grp), N)], old, by = c("ct", "grp"))
stopifnot(nrow(chk2) == 16, all(chk2$N == chk2$count))
stopifnot(identical(unique(rna_g[, .(grp, n_total)])$n_total, c(18695L, 16635L)),
          identical(unique(atac_g[, .(grp, n_total)])$n_total, c(5044L, 4566L)))
cat("计数核对通过：RNA Con/Mhp =", unique(rna_g$n_total), "| ATAC Con/Mhp =", unique(atac_g$n_total), "\n")

# 落表（中文文件名；显示名一并给出）
fwrite(rna_g[, .(group = grp, celltype = ct, label = lab_rna(ct), n = N, n_total, perc = round(perc, 2))],  file.path(OUT, "细胞比例_RNA按组.csv"))
fwrite(rna_s[, .(sample = smp, celltype = ct, label = lab_rna(ct), n = N, n_total, perc = round(perc, 2))], file.path(OUT, "细胞比例_RNA按样本.csv"))
fwrite(atac_g[, .(group = grp, celltype = ct, label = lab_atac(ct), n = N, n_total, perc = round(perc, 2))],  file.path(OUT, "细胞比例_ATAC按组.csv"))
fwrite(atac_s[, .(sample = smp, celltype = ct, label = lab_atac(ct), n = N, n_total, perc = round(perc, 2))], file.path(OUT, "细胞比例_ATAC按样本.csv"))

# ============================================================
# 2. 流图（两柱 + 彩带 + 表格）—— 坐标单位 = mm，coord_fixed 保证 1 数据单位 = 1 mm
# ============================================================
rib_col <- function(col) colorspace::lighten(col, 0.15)          # 彩带用类型色的浅色版（本脚本新增）
FLOW <- list(row_h = 4.8, W = 88, bar = list(Con = c(1, 11), Mhp = c(23, 33)),   # 行高 / 总宽 / 两柱 x 范围（柱宽 10，彩带间距 12）；左右各留 1 mm
             xt = c(37, 47.5, 58, 87))                                            # 表格 4 条竖线：Con% | Mhp% | 色块+名称（10.5 / 10.5 / 29 mm）

mk_flow <- function(pg, cols, labs, title_bar, title_tab, total_lab) {
  # pg: grp(因子 Con/Mhp) / ct(因子，顺序即行序) / N / n_total / perc
  cts <- levels(pg$ct); N <- length(cts); rh <- FLOW$row_h; xt <- FLOW$xt; bx <- FLOW$bar
  y0 <- 1 + rh; y1 <- y0 + rh * N                      # 类型行区底（= 总数行顶）与顶（= 柱顶）
  yh <- y1 + 2.3; yT <- y1 + 5.9; H <- y1 + 8.0        # 表头行、标题行中心 y；总高
  xc_tab <- c(Con = mean(xt[1:2]), Mhp = mean(xt[2:3]))
  # 柱段：从顶向下按 ct 顺序堆叠
  seg <- rbindlist(lapply(grp_levels, function(g) {
    s <- pg[grp == g][order(ct)]; h <- s$perc / 100 * (y1 - y0)
    data.table(grp = g, ct = as.character(s$ct), xmin = bx[[g]][1], xmax = bx[[g]][2],
               ymax = y1 - c(0, cumsum(h)[-N]), ymin = y1 - cumsum(h)) }))
  seg[, fill := unname(cols[ct])]
  # S 形彩带：同一类型左右柱段相连，smoothstep 插值（沿用 10_CellProportion.R 段 C）
  sfun <- function(t) t^2 * (3 - 2 * t); tt <- seq(0, 1, length.out = 80)
  rib <- rbindlist(lapply(cts, function(k) {
    a <- seg[grp == "Con" & ct == k]; b <- seg[grp == "Mhp" & ct == k]
    x <- a$xmax + (b$xmin - a$xmax) * tt
    data.table(ct = k, x = c(x, rev(x)),
               y = c(a$ymax + (b$ymax - a$ymax) * sfun(tt), rev(a$ymin + (b$ymin - a$ymin) * sfun(tt)))) }))
  rib[, fill := rib_col(unname(cols[ct]))]
  # 表格行
  rows <- data.table(ct = cts, i = seq_len(N))[, yc := y1 - (i - 0.5) * rh][, fill := unname(cols[ct])][, lab := unname(labs[ct])]
  tc  <- merge(pg[, .(grp = as.character(grp), ct = as.character(ct), perc)], rows[, .(ct, yc)], by = "ct")
  tc[, `:=`(xc = xc_tab[grp], lab = sprintf("%.2f", perc))]
  tot <- unique(pg[, .(grp = as.character(grp), n_total)])[, `:=`(xc = xc_tab[grp], lab = prettyNum(n_total, big.mark = ","))]
  hl  <- data.table(y = c(y1 - (0:N) * rh, 1))          # 水平线：每行边界 + 总数行底
  hdr <- data.table(x = c(mean(bx$Con), mean(bx$Mhp), xc_tab, mean(xt[3:4])), y = yh, lab = c("Con", "Mhp", "Con", "Mhp", "Cell type"))
  ttl <- data.table(x = c(mean(unlist(bx)), mean(xt[c(1, 3)])), y = yT, lab = c(title_bar, title_tab))
  p <- ggplot() +
    geom_polygon(data = rib, aes(x, y, group = ct, fill = fill), alpha = 0.5, colour = NA) +
    geom_rect(data = seg, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill), colour = "black", linewidth = lw(0.3)) +
    geom_segment(data = hl, aes(x = xt[1], xend = xt[4], y = y, yend = y), linewidth = lw(0.3), colour = "black") +
    annotate("segment", x = xt, xend = xt, y = 1, yend = y1, linewidth = lw(0.3), colour = "black") +
    geom_text(data = tc,  aes(xc, yc, label = lab), size = pt2size(6), family = FONT) +
    geom_text(data = tot, aes(xc, 1 + rh / 2, label = lab), size = pt2size(6), family = FONT) +
    annotate("text", x = mean(xt[3:4]), y = 1 + rh / 2, label = total_lab, size = pt2size(6), family = FONT) +
    geom_rect(data = rows, aes(xmin = xt[3] + 1.6, xmax = xt[3] + 4.6, ymin = yc - 1.5, ymax = yc + 1.5, fill = fill), colour = "black", linewidth = lw(0.3)) +
    geom_text(data = rows, aes(xt[3] + 6, yc, label = lab), hjust = 0, size = pt2size(6), family = FONT) +
    geom_text(data = hdr, aes(x, y, label = lab), size = pt2size(6), family = FONT, fontface = "bold") +
    geom_text(data = ttl, aes(x, y, label = lab), size = pt2size(7), family = FONT, fontface = "bold") +
    scale_fill_identity() +
    coord_fixed(ratio = 1, xlim = c(0, FLOW$W), ylim = c(0, H), expand = FALSE, clip = "off") +
    theme_mhp_void() + theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0))
  attr(p, "size_mm") <- c(FLOW$W, H); p }

pF <- mk_flow(rna_g, ct_rna_col, ct_rna_lab, "scRNA-seq", "Cells (%)", "Total no. of cells")
sz <- attr(pF, "size_mm")
sav(pF, "Fig1F_scRNA细胞比例_流图", OUT, sz[1], sz[2], final = FALSE, panel_dir = PANEL)

pS2a <- mk_flow(atac_g, ct_atac_col, ct_atac_lab, "snATAC-seq", "Nuclei (%)", "Total no. of nuclei")
sz <- attr(pS2a, "size_mm")
sav(pS2a, "SuppFig_snATAC细胞比例_流图", OUT, sz[1], sz[2], final = FALSE, panel_dir = PANEL)

# ============================================================
# 3. 逐样本 100% 堆叠柱（每头猪各自归一；柱下写细胞数；coord_cartesian 只裁视野）
#    面板宽度用 theme(panel.widths) 固定为 BAR_PW mm（两个 facet 合计），S1 / S2b 柱宽因此完全一致；
#    图例 legend.justification 左贴面板，剩余空白留在右侧（88 mm 画布内 snATAC 长图例名也放得下）
# ============================================================
BAR_PW <- 52                                                    # mm（两 facet 合计）；柱间距 ≈ 7.8 mm，6 pt 的 "Mhp_1" 两侧各余 ~1 mm
# 图例 key：固定 2.6 mm 方块（本脚本新增）。ggplot2 ≥ 3.5 的 key 会随多行图例文字拉高，snATAC 长名折成两行时默认方块会变长条
draw_key_sq <- function(data, params, size) {
  g <- grid::rectGrob(width = grid::unit(2.6, "mm"), height = grid::unit(2.6, "mm"),
                      gp = grid::gpar(fill = data$fill, col = NA))
  attr(g, "width") <- 0.26; attr(g, "height") <- 0.26           # cm，向 ggplot2 申报 glyph 尺寸
  g }
wrap_leg <- function(s, max_ch = 16) vapply(s, function(z) {
  if (!grepl(" (", z, fixed = TRUE)) return(z)
  head <- sub(" \\(.*$", "", z); parts <- strsplit(sub("^.*\\(", "(", z), ", ", fixed = TRUE)[[1]]
  lines <- parts[1]
  for (p in parts[-1]) {
    cand <- paste(lines[length(lines)], p, sep = ", ")
    if (nchar(cand) <= max_ch) lines[length(lines)] <- cand
    else { lines[length(lines)] <- paste0(lines[length(lines)], ","); lines <- c(lines, p) } }   # 换行处补回被 strsplit 吃掉的逗号
  paste(c(head, lines), collapse = "\n") }, character(1))            # 保留 names（= 数据值），scale_fill_manual 按名匹配
mk_sample_bar <- function(ps, cols, leg_labs, title, ylab) {
  ps <- copy(ps)[, grp := to_grp(sub("_.*$", "", as.character(smp)))]
  xl <- unique(ps[, .(smp, n_total)])[, lab := paste0(smp, "\n", prettyNum(n_total, big.mark = ","))]
  leg_labs <- wrap_leg(leg_labs)                                # 长名折行（见 wrap_leg），图例宽度与 scRNA 版一致；长名取自 lab_atac(legend = TRUE)，不写死
  ggplot(ps, aes(smp, perc, fill = ct)) +
    geom_col(width = 0.8, colour = "white", linewidth = lw(0.25), key_glyph = draw_key_sq) +
    scale_fill_manual(values = cols, labels = leg_labs, name = NULL) +
    scale_x_discrete(labels = setNames(xl$lab, as.character(xl$smp))) +
    scale_y_continuous(breaks = seq(0, 100, 25), expand = expansion(0)) +
    coord_cartesian(ylim = c(0, 100), clip = "off") +
    facet_grid(~grp, scales = "free_x", space = "free_x") +
    labs(x = NULL, y = ylab, title = title) +
    theme_mhp() +
    theme(strip.text = element_text(size = FS_BASE, face = "bold"), plot.title = element_text(size = FS_BASE),
          axis.text.x = element_text(lineheight = 0.95), axis.line.x = element_blank(), axis.ticks.x = element_blank(),
          legend.key.size = unit(2.6, "mm"), legend.key.spacing.y = unit(0.4, "mm"), legend.margin = margin(0, 0, 0, 0),
          legend.key.justification = "top",
          legend.justification = c(0, 0.5), panel.widths = unit(BAR_PW, "mm"),
          panel.spacing.x = unit(2, "mm"), plot.margin = margin(1, 1, 1, 1)) }

pS1 <- mk_sample_bar(rna_s, ct_rna_col, ct_rna_lab, "scRNA-seq", "Cells (%)")
sav(pS1, "SuppFig_scRNA逐动物细胞比例", OUT, W1, 52, final = FALSE, panel_dir = PANEL)
pS2b <- mk_sample_bar(atac_s, ct_atac_col, ct_atac_lab_legend, "snATAC-seq", "Nuclei (%)")
sav(pS2b, "SuppFig_snATAC逐样本细胞比例", OUT, W1, 52, final = FALSE, panel_dir = PANEL)

# ============================================================
# 4. RNA vs ATAC 比例 log2FC(Mhp/Con) 散点：8 个可对应类型（CD4+ T、CD8+ T、γδ T（数据值 "NK cells"）→ T/NK，按 ct_map 合并；DC 无对应，不参与）
# ============================================================
rna_m <- merge(rna_g[, .(ct = as.character(ct), grp, N, n_total)], ct_map, by.x = "ct", by.y = "rna")[!is.na(atac)]
rna_m <- rna_m[, .(N = sum(N), n_total = n_total[1]), by = .(atac, grp)][, perc := 100 * N / n_total]
fc <- merge(dcast(rna_m, atac ~ grp, value.var = "perc")[, .(ct = atac, rna_con = Con, rna_mhp = Mhp)],
            dcast(atac_g[, .(ct = as.character(ct), grp, perc)], ct ~ grp, value.var = "perc")[, .(ct, atac_con = Con, atac_mhp = Mhp)], by = "ct")
fc[, `:=`(rna_log2fc = log2(rna_mhp / rna_con), atac_log2fc = log2(atac_mhp / atac_con))]
fc[, `:=`(same_dir = sign(rna_log2fc) == sign(atac_log2fc), lab = lab_atac(ct), col = unname(ct_atac_col[ct]))]
fc[ct == "B_cells", lab := "B cells (3 nuclei in Mhp)"]         # ATAC Mhp 组 B 细胞仅 3 核，log2FC 极端
fc[, ct := factor(ct, levels = ct_atac_order)]; setorder(fc, ct)
sp <- suppressWarnings(cor.test(fc$rna_log2fc, fc$atac_log2fc, method = "spearman", exact = TRUE))   # n = 8，精确 P
rho <- unname(sp$estimate); p_sp <- sp$p.value; n_same <- sum(fc$same_dir)
fwrite(fc[, .(celltype = ct, label = lab_atac(ct), rna_con, rna_mhp, atac_con, atac_mhp, rna_log2fc, atac_log2fc, same_dir)],
       file.path(OUT, "细胞比例_RNA_vs_ATAC_log2FC.csv"))
cat(sprintf("RNA vs ATAC log2FC：Spearman rho = %.3f，P = %.4f，同向 %d/8\n", rho, p_sp, n_same)); print(fc[, .(ct, rna_log2fc, atac_log2fc, same_dir)])

# 标签逐类手工 nudge（按 ct_atac_order 顺序；单位 = 数据单位），避开 y = 0 / x = 0 虚线与相邻点；ggrepel 只做微调
NUDGE <- data.table(ct = ct_atac_order,
                    dx = c(-0.40, 0.00, 0.80, -2.50, -1.75, 1.60, 0.00, 0.65),
                    dy = c( 0.55, 0.50, -0.50, 0.00, 0.00, -0.12, -0.60, 0.55))
fc <- merge(fc, NUDGE, by = "ct", sort = FALSE); setorder(fc, ct)
lim <- c(-6.5, 3.2)
pS3 <- ggplot(fc, aes(rna_log2fc, atac_log2fc)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", linewidth = lw(0.3), colour = "grey60") +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = lw(0.3), colour = "grey40") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = lw(0.3), colour = "grey40") +
  geom_point(aes(colour = col), size = 2.2, stroke = 0) +
  geom_text_repel(aes(label = lab), nudge_x = fc$dx, nudge_y = fc$dy, force = 0.3, size = pt2size(6), family = FONT,
                  min.segment.length = 0.4, segment.size = lw(0.25), segment.colour = "grey50",
                  box.padding = 0.15, point.padding = 0.2, max.overlaps = Inf, seed = 3) +
  scale_colour_identity() +
  scale_x_continuous(limits = lim, breaks = seq(-6, 3, 3)) + scale_y_continuous(limits = lim, breaks = seq(-6, 3, 3)) +
  coord_fixed(ratio = 1, expand = FALSE) +
  labs(x = expression("scRNA-seq "*log[2]*"FC (Mhp/Con)"), y = expression("snATAC-seq "*log[2]*"FC (Mhp/Con)"),
       # 统计注释放面板上方（面板内四象限都有点/标签）；rho / P 用 plotmath 渲染，源码全 ASCII；n = 8 类，P 如实给出
       title = bquote("Spearman "*italic(rho)*" = "*.(sprintf("%.2f", rho))*", "*italic(P)*" = "*.(sprintf("%.3f", p_sp))*", "*italic(n)*" = 8"),
       subtitle = sprintf("Same direction in both assays: %d of 8", n_same)) +
  theme_mhp() + theme(plot.margin = margin(1, 2, 1, 1),
                      plot.title = element_text(size = FS_SMALL, face = "plain", hjust = 0, margin = margin(b = 0.5)),
                      plot.subtitle = element_text(size = FS_SMALL, hjust = 0, margin = margin(b = 1.5)))
sav(pS3, "SuppFig_RNA_vs_ATAC比例log2FC散点", OUT, 55, 57, final = FALSE, panel_dir = PANEL)

# ============================================================
# 5. 关键数字
# ============================================================
cat("\n---- RNA 合并口径（%）----\n"); print(dcast(rna_g[, .(ct = lab_rna(ct), grp, perc = round(perc, 2))], ct ~ grp, value.var = "perc"))
cat("---- ATAC 合并口径（%）----\n"); print(dcast(atac_g[, .(ct = lab_atac(ct), grp, perc = round(perc, 2))], ct ~ grp, value.var = "perc"))
cat("✅ 06 完成 →", OUT, "\n")
