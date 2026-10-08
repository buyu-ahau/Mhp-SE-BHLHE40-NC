# ============================================================================
# 11_补充图5_四因子mRNA_基因符号_本地.R — 补充图 5：四个炎症因子 mRNA 在 6、12、24 h 的敲低效应时间序列
# 运行位置：本地；路径变量见 config/paths.R
# ============================================================================

source(file.path(Sys.getenv("MHP_CODE", "."), "config", "paths.R"))

suppressPackageStartupMessages({ library(xml2); library(data.table); library(ggplot2); library(ggtext) })
source(file.path(MHP_CODE, "shared", "_palette_MHP.R"))   # theme_mhp、sav、lw、pt2size、FONT、.pdf_dev

D   <- file.path(MHP_RES, "10.验证结果/BHLHE40补充验证/20260618肺炎项目数据汇总")
OUT <- file.path(MHP_RES, "10.验证结果/BHLHE40补充验证/新Fig3_作图")
REF <- file.path(OUT, "Plots/面板")                                        # 01 现有输出（基准，只读）
PAN <- if (nzchar(Sys.getenv("SF5_TEST_DIR"))) file.path(Sys.getenv("SF5_TEST_DIR"), "Plots/面板") else REF
dir.create(PAN, recursive = TRUE, showWarnings = FALSE)
FINAL <- identical(Sys.getenv("MHP_FINAL"), "1")
STEM_OLD <- "补充_四因子mRNA_6-12-24h"; STEM <- "补充图5_四因子mRNA_6-12-24h_基因符号"; W_MM <- 150; H_MM <- 55
targets <- file.path(PAN, paste0(STEM, if (FINAL) c(".svg", ".tiff", ".pdf") else ".png"))
if (any(file.exists(targets))) stop("输出已存在，按规定不覆盖：", paste(basename(targets[file.exists(targets)]), collapse = ", "))
TMP <- file.path(tempdir(), "补充图5核对"); dir.create(TMP, showWarnings = FALSE)   # 复现与透明标题的核对图（R 退出即清理）

# ---------------------------------------------------------------- 读数与统计（逐字复制自 01）
GRP <- c("si-NC", "si-NC+Mhp", "si-BHLHE40", "si-BHLHE40+Mhp")
COL <- c("si-NC" = "#707070", "si-NC+Mhp" = "#C2507B", "si-BHLHE40" = "#E9A06B", "si-BHLHE40+Mhp" = "#7B3F8C")
CMP <- list(c("si-NC", "si-NC+Mhp"), c("si-NC+Mhp", "si-BHLHE40+Mhp"), c("si-NC", "si-BHLHE40"))
pz_tables <- function(f) { x <- read_xml(f); xml_ns_strip(x); xml_find_all(x, "//Table") }
two_way <- function(t) {                                   # TwoWay 表：行 = 时间点，子列 = 重复
  rows <- xml_text(xml_find_all(t, "./RowTitlesColumn/Subcolumn/d"))
  rbindlist(lapply(xml_find_all(t, "./YColumn"), function(yc) {
    subs <- xml_find_all(yc, "./Subcolumn")
    rbindlist(lapply(seq_along(subs), function(r) { v <- suppressWarnings(as.numeric(xml_text(xml_find_all(subs[[r]], "./d"))))
      data.table(group = trimws(xml_text(xml_find_first(yc, "./Title"))), time = rows[seq_along(v)], rep = r, value = v) })) }))
}
tb <- function(sub) pz_tables(file.path(D, sub))
tukey <- function(d) {
  d <- copy(d)[, g := factor(paste0("G", as.integer(group)))]
  tk <- TukeyHSD(aov(value ~ g, data = d))$g
  rbindlist(lapply(CMP, function(cp) { a <- paste0("G", match(cp, GRP)); nm <- c(paste(a[2], a[1], sep = "-"), paste(a[1], a[2], sep = "-"))
    k <- intersect(nm, rownames(tk)); data.table(g1 = cp[1], g2 = cp[2], p = tk[k, "p adj"]) }))
}
star <- function(p) fcase(p < 0.001, "***", p < 0.01, "**", p < 0.05, "*", default = "ns")

q <- tb("2.1 qPCR检测炎症因子/GraphPad.pzfx"); qn <- c("IL-1β", "IL-6", "IL-8", "TNF-α")   # 表序 = qPCR 的 IL-1β / IL-6 / IL-8 / TNF-α（01 已按均值核对）
qa <- rbindlist(lapply(1:4, function(i) two_way(q[[i]])[, marker := qn[i]]))[!is.na(value)]
qa[, `:=`(group = factor(group, levels = GRP), h = as.numeric(sub("h", "", time)), marker = factor(marker, levels = qn))]
qs <- qa[, .(mean = mean(value)), by = .(marker, h, group)]
qst <- qa[, tukey(.SD[, .(group, value)]), by = .(marker, h)][g1 == "si-NC+Mhp" & g2 == "si-BHLHE40+Mhp"][, star := star(p)]
qst <- merge(qst, qa[group %in% c("si-NC+Mhp", "si-BHLHE40+Mhp"), .(ytop = max(value)), by = .(marker, h)], by = c("marker", "h"))
qst[, y := ytop * 1.07]

# 核对 ①：与 01 写出的统计表一致（24 h si-NC+Mhp 组均值同 01 的 chk）
ref_st <- fread(file.path(OUT, "Tables/补充_四因子mRNA时间序列_统计.csv"), encoding = "UTF-8")
ck <- merge(qst[, .(marker = as.character(marker), time_h = h, p3 = signif(p, 3), star)], ref_st, by = c("marker", "time_h"))
stopifnot(nrow(ref_st) == 12, nrow(ck) == 12, all(ck$p3 == ck$p), all(ck$star.x == ck$star.y))
m24 <- qa[h == 24 & group == "si-NC+Mhp", mean(value), by = marker][, setNames(V1, marker)]
stopifnot(abs(m24 - c(2.63, 2.68, 2.35, 2.15)) < 0.02 * c(2.63, 2.68, 2.35, 2.15))
cat("✓ 核对①：12 个 P 值与星号 = 补充_四因子mRNA时间序列_统计.csv；24 h 组均值 = 01 核对值\n")

GENE_MD <- c("IL-1β" = "*IL1B*", "IL-6" = "*IL6*", "IL-8" = "*CXCL8*", "TNF-α" = "*TNF*")
stopifnot(identical(names(GENE_MD), qn))
M0 <- grid::convertUnit(calc_element("strip.text.x.top", theme_mhp() + theme(strip.text = element_text(size = FS_BASE)))$margin, "pt", valueOnly = TRUE)
stopifnot(length(M0) == 4)
mk_ps <- function(gene_symbol, strip_col = NULL, dh = 0, s = 0) {
  fw <- if (gene_symbol) facet_wrap(~ marker, nrow = 1, scales = "free_y", labeller = as_labeller(GENE_MD))
        else facet_wrap(~ marker, nrow = 1, scales = "free_y")
  st <- if (!gene_symbol) element_text(size = FS_BASE, colour = strip_col)
        else if (dh == 0 && s == 0) element_markdown(size = FS_BASE, colour = strip_col)
        else element_markdown(size = FS_BASE, colour = strip_col, margin = margin(M0[1] - dh + s, M0[2], M0[3] - s, M0[4], unit = "pt"))
  ggplot(qs, aes(h, mean, colour = group)) +
    geom_point(data = qa, aes(h, value, colour = group), size = 2.0, alpha = 0.30, stroke = 0, shape = 16, inherit.aes = FALSE) +
    geom_line(linewidth = lw(0.9)) + geom_point(size = 1.6, stroke = 0.7, shape = 21, fill = "white") +
    geom_text(data = qst, aes(h, y, label = star), inherit.aes = FALSE, vjust = ifelse(qst$star == "ns", -0.2, 0.3),
              size = pt2size(ifelse(qst$star == "ns", FS_SMALL, FS_BASE + 1)), family = FONT) +
    fw +
    scale_colour_manual(values = COL, name = NULL, breaks = GRP) +
    scale_x_continuous(breaks = c(6, 12, 24), expand = expansion(add = c(2.5, 2.5))) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
    labs(x = "Time after infection (h)", y = "Relative mRNA") +
    theme_mhp() + theme(strip.background = element_blank(), strip.text = st, legend.position = "top",
                        legend.key.size = unit(3, "mm"), axis.text = element_text(colour = "black"))
}
P_OLD <- mk_ps(FALSE)

# 输出设备（与 _palette_MHP.R 的 sav() / 01 的 pdf1() 参数完全一致）
DEV <- list(png  = function(f) ragg::agg_png(f, width = W_MM, height = H_MM, units = "mm", res = 300),
            tiff = function(f) ragg::agg_tiff(f, width = W_MM, height = H_MM, units = "mm", res = 900),
            svg  = function(f) .svg_dev(f, width = W_MM / 25.4, height = H_MM / 25.4),
            pdf  = function(f) .pdf_dev(f, width = W_MM / 25.4, height = H_MM / 25.4))
rend <- function(p, file) {
  ext <- tools::file_ext(file)
  if (ext == "png")  ggsave(file, p, width = W_MM, height = H_MM, units = "mm", dpi = 300, bg = "white", device = ragg::agg_png)
  if (ext == "svg")  ggsave(file, p, width = W_MM, height = H_MM, units = "mm", device = .svg_dev)
  if (ext == "tiff") ggsave(file, p, width = W_MM, height = H_MM, units = "mm", dpi = 900, bg = "white", device = ragg::agg_tiff, compression = "lzw")
  if (ext == "pdf")  ggsave(file, p, width = W_MM, height = H_MM, units = "mm", device = .pdf_dev)
  stopifnot(file.exists(file)); invisible(file)
}

fit_new <- function(ext, strip_col = NULL) {
  DEV[[ext]](file.path(TMP, paste0("measure.", ext))); on.exit(invisible(dev.off()))
  hpt <- function(g) grid::convertHeight(g$heights, "pt", valueOnly = TRUE); wpt <- function(g) grid::convertWidth(g$widths, "pt", valueOnly = TRUE)
  g0 <- ggplotGrob(mk_ps(FALSE, strip_col))
  r <- unique(g0$layout$t[grepl("^strip-t", g0$layout$name)]); stopifnot(length(r) == 1)
  dh <- hpt(ggplotGrob(mk_ps(TRUE, strip_col)))[r] - hpt(g0)[r]
  desc <- function(face) grid::convertHeight(grid::grobDescent(grid::textGrob("gjpqyQ", gp = grid::gpar(fontsize = FS_BASE, fontfamily = FONT, fontface = face))), "pt", valueOnly = TRUE)
  s <- desc("italic") - desc("plain") / 2
  stopifnot(s > 0, s < M0[3])                                                                 # 下边距够减
  p <- mk_ps(TRUE, strip_col, dh, s); g1 <- ggplotGrob(p)
  stopifnot(identical(dim(g0), dim(g1)), identical(g0$layout[, c("t", "l", "b", "r", "name")], g1$layout[, c("t", "l", "b", "r", "name")]),
            isTRUE(all.equal(hpt(g0), hpt(g1), tolerance = 1e-12)), isTRUE(all.equal(wpt(g0), wpt(g1), tolerance = 1e-12)),
            !any(grid::unitType(g1$heights[seq_len(r)]) == "null"))                            # 标题行以上均为绝对单位 → 可累加得像素带
  mm <- hpt(g1) / 72.27 * 25.4                                                                # grid 的 pt = 1/72.27 in
  list(p = p, dh = dh, s = s, band_mm = c(sum(mm[seq_len(r - 1)]), sum(mm[seq_len(r)])), row = r)
}
md5 <- function(f) unname(tools::md5sum(f))
rd <- function(f) round(255 * if (grepl("\\.png$", f)) png::readPNG(f, native = FALSE) else tiff::readTIFF(f))   # 8 位图读成 0–1 → 还原为 0–255 整数
px_check <- function(ext, dpi) {
  f_old <- rend(P_OLD, file.path(TMP, paste0("old.", ext)))
  stopifnot(md5(f_old) == md5(file.path(REF, paste0(STEM_OLD, ".", ext))))
  nt <- fit_new(ext, "transparent")
  stopifnot(md5(rend(mk_ps(FALSE, "transparent"), file.path(TMP, paste0("old_notitle.", ext)))) ==
            md5(rend(nt$p, file.path(TMP, paste0("new_notitle.", ext)))))
  fn <- fit_new(ext); stopifnot(fn$dh == nt$dh, fn$s == nt$s)
  f_new <- rend(fn$p, file.path(PAN, paste0(STEM, ".", ext)))
  a <- rd(f_old); b <- rd(f_new); stopifnot(identical(dim(a), dim(b)))
  dif <- Reduce(`|`, lapply(seq_len(dim(a)[3]), function(k) a[, , k] != b[, , k]))
  rows <- which(rowSums(dif) > 0); band_px <- fn$band_mm / 25.4 * dpi
  stopifnot(length(rows) > 0, min(rows) - 1 >= floor(band_px[1]), max(rows) <= ceiling(band_px[2]))
  cols <- which(colSums(dif) > 0); runs <- split(cols, cumsum(c(1, diff(cols) > 5 * dpi / 300)))
  r0 <- max(1, floor(band_px[1])); r1 <- min(dim(a)[1], ceiling(band_px[2]))
  ink_bottom <- function(img) {
    ink <- Reduce(`|`, lapply(1:3, function(k) img[r0:r1, , k] < 254.5))
    cc <- which(colSums(ink) > 0); cl <- split(cc, cumsum(c(1, diff(cc) > 20 * dpi / 300)))
    cl <- cl[vapply(cl, function(j) any(j %in% cols), logical(1))]; stopifnot(length(cl) == 4)
    sapply(cl, function(j) r0 - 1 + max(which(rowSums(ink[, j, drop = FALSE]) > 0)))
  }
  bo <- ink_bottom(a); bn <- ink_bottom(b); tol <- max(1, round(dpi / 450))
  stopifnot(all(abs(bo[2:4] - bn[2:4]) <= tol))
  cat(sprintf("✓ %s（%d dpi，dh = %.4g pt，s = %.4g pt）：旧版复现与 01 现有文件逐字节相同；标题透明后新旧逐字节相同；差异像素 %d 个，全部在第 %d–%d 行（标题行像素带 %.1f–%.1f），%d 簇，横向 %s；标题墨迹底行 旧 %s → 新 %s\n",
              toupper(ext), dpi, fn$dh, fn$s, sum(dif), min(rows), max(rows), band_px[1], band_px[2], length(runs),
              paste(sapply(runs, function(r) paste(range(r), collapse = "–")), collapse = " / "),
              paste(bo, collapse = "/"), paste(bn, collapse = "/")))
  invisible(f_new)
}

if (!FINAL) {
  px_check("png", 300)
} else {
  f_old <- rend(P_OLD, file.path(TMP, "old.svg")); stopifnot(md5(f_old) == md5(file.path(REF, paste0(STEM_OLD, ".svg"))))
  fs_ <- fit_new("svg")
  f_svg <- rend(fs_$p, file.path(PAN, paste0(STEM, ".svg")))
  so <- readLines(f_old, encoding = "UTF-8", warn = FALSE); sn <- readLines(f_svg, encoding = "UTF-8", warn = FALSE)
  stopifnot(length(sn) == length(so))
  k <- which(so != sn); txt <- function(x) sub(".*>([^<]*)</text>.*", "\\1", x)
  ya <- function(x) as.numeric(sub(".* y='([0-9.]+)'.*", "\\1", x))                            # 文字基线 y
  stopifnot(length(k) == 4, all(grepl("^<text", c(so[k], sn[k]))), identical(txt(so[k]), qn),
            identical(txt(sn[k]), c("IL1B", "IL6", "CXCL8", "TNF")), all(grepl("font-style: italic", sn[k])), !any(grepl("italic", so[-k])),
            all(abs(ya(so[k]) - ya(sn[k])) <= 0.011))                                         # 基线原位（SVG 坐标保留 2 位小数）
  # 终稿规范：无 <image>、无 Symbol 字体、最小字号 ≥ 5 pt（svglite 以 1 px = 1 pt 记字号）
  fs <- as.numeric(regmatches(sn, regexpr("(?<=font-size: )[0-9.]+(?=px)", sn, perl = TRUE)))
  stopifnot(!any(grepl("<image", sn)), !any(grepl("Symbol", sn)), min(fs) >= 5)
  cat(sprintf("✓ SVG（dh = %.4g pt，s = %.4g pt）：旧版复现与 01 现有文件逐字节相同；共 %d 行，新旧只差第 %s 行 = 四个分面标题（%s → %s，italic；基线 y 旧 %s / 新 %s）；无 <image>、无 Symbol，最小字号 %.2f pt\n",
              fs_$dh, fs_$s, length(sn), paste(k, collapse = "、"), paste(txt(so[k]), collapse = " / "), paste(txt(sn[k]), collapse = " / "),
              paste(ya(so[k]), collapse = "/"), paste(ya(sn[k]), collapse = "/"), min(fs)))
  # ---- TIFF 900 dpi：同 PNG
  px_check("tiff", 900)
  # ---- PDF（quartz，纯矢量）：栅格数 = 0、字体只含 Arial 系、pdftoppm 回渲染
  fp <- fit_new("pdf"); f_pdf <- rend(fp$p, file.path(PAN, paste0(STEM, ".pdf")))
  stopifnot(sum(grepl("/Subtype /Image", readLines(f_pdf, warn = FALSE, skipNul = TRUE), useBytes = TRUE)) == 0)
  pf <- system2("pdffonts", shQuote(f_pdf), stdout = TRUE)
  fonts <- unique(sub("^[A-Z]{6}\\+", "", sapply(strsplit(trimws(pf[-(1:2)]), "\\s+"), `[`, 1)))
  stopifnot(length(fonts) >= 2, all(grepl("^Arial", fonts)), any(grepl("Italic", fonts)), !any(grepl("Symbol", fonts)))
  prev <- file.path(Sys.getenv("SF5_PDF_PREVIEW", TMP), "补充图5_PDF回渲染")
  system2("pdftoppm", c("-png", "-r", "300", "-singlefile", shQuote(f_pdf), shQuote(prev))); stopifnot(file.exists(paste0(prev, ".png")))
  # 旧 PDF（01 现有文件）同法回渲染，逐像素比：差异只在标题行像素带内（PDF 含时间戳，不比字节）
  system2("pdftoppm", c("-png", "-r", "300", "-singlefile", shQuote(file.path(REF, paste0(STEM_OLD, ".pdf"))), shQuote(file.path(TMP, "pdf_old"))))
  a <- rd(file.path(TMP, "pdf_old.png")); b <- rd(paste0(prev, ".png")); stopifnot(identical(dim(a), dim(b)))
  dif <- Reduce(`|`, lapply(1:3, function(k) a[, , k] != b[, , k])); rows <- which(rowSums(dif) > 0); band_px <- fp$band_mm / 25.4 * 300
  stopifnot(length(rows) > 0, min(rows) - 1 >= floor(band_px[1]), max(rows) <= ceiling(band_px[2]))
  cat(sprintf("✓ PDF 回渲染（300 dpi）新旧差异像素 %d 个，全部在第 %d–%d 行（标题行像素带 %.1f–%.1f）\n", sum(dif), min(rows), max(rows), band_px[1], band_px[2]))
  cat(sprintf("✓ PDF（dh = %.4g pt，s = %.4g pt）：/Subtype /Image = 0；字体 %s；pdftoppm 回渲染 → %s.png\n", fp$dh, fp$s, paste(fonts, collapse = "、"), prev))
}
cat("✓ 补充图 5 基因符号版完成 →", paste(basename(targets), collapse = "、"), "\n")
