# ============================================================================
# 09_g_E2_BHLHE40结合_本地.R — Fig. 6a 右：E2 区域的 BHLHE40 CUT&Tag 结合轨迹与 Peak_32739 归一化计数（本脚本按早期 Fig. 5 版式输出，文件名为 Fig5h_*；Fig. 6a 终版只调整了版式与字号，数据相同）
# 运行位置：本地；路径变量见 config/paths.R
# ============================================================================

source(file.path(Sys.getenv("MHP_CODE", "."), "config", "paths.R"))

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(grid); library(patchwork)
  library(rtracklayer); library(GenomicRanges); library(GenomicFeatures); library(Biostrings); library(BSgenome.Sscrofa.UCSC.susScr11) })
source(file.path(MHP_CODE, "shared", "_palette_MHP.R"))   # theme_mhp、sav、lw、pt2size、FONT、grp_col、MHP_RES
BH_COL <- "#B2182B"                                                                # BHLHE40 高亮红（新 Fig 4 / Fig 5 各脚本共用常量）
TRACK_SCALE <- "RPKM"                                                              # ★"RPKM"（默认，保守）或 "DESeq2"（乘 B/A size factor）
stopifnot(TRACK_SCALE %in% c("RPKM", "DESeq2"))
BIN <- 50L                                                                         # 轨迹 bin（bp）：68 mm 画 23.4 kb 时 50 bp ≈ 0.145 mm ≈ 0.41 pt（≥ Nature 0.25 pt）

# ---------------------------------------------------------------- 路径
CT   <- file.path(MHP_RES, "10.验证结果/7.转录因子CUTtag")
BWD  <- file.path(CT, "09_BW文件")
F7T  <- file.path(CT, "新Fig7_作图/Tables")
OUT  <- file.path(MHP_RES, "7、HI-C结果/新Fig5_作图")
PAN  <- file.path(OUT, "Plots/面板"); TAB <- file.path(OUT, "Tables"); PNL <- file.path(OUT, "_panels")
for (d in c(PAN, TAB, PNL)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
F_SF   <- file.path(F7T, "02_sizeFactors_按归一化方法.csv")
F_LIB  <- file.path(F7T, "01_样本文库统计_FRiP.csv")
F_PEAK <- file.path(F7T, "05_BHLHE40位点_chr13_61.0-61.8Mb_峰表.csv")
F_SHAPE<- file.path(F7T, "05_Peak32739_逐重复信号形状.csv")
F_TXDB <- file.path(CT, "Sus_scrofa.Sscrofa11.1.114.sqlite")
stopifnot(file.exists(F_SF), file.exists(F_LIB), file.exists(F_PEAK), file.exists(F_SHAPE), file.exists(F_TXDB))

# ---------------------------------------------------------------- 已定坐标（任务给定）
CHR_ENS <- "13"; CHR_UCSC <- "chr13"
E2      <- c(61329246L, 61350682L); E2_LEN <- E2[2] - E2[1] + 1L; stopifnot(E2_LEN == 21437L)
FLANK   <- 1000L
WIN     <- c(E2[1] - FLANK, E2[2] + FLANK)                                          # 61,328,246–61,351,682（23,437 bp）
LUC     <- c(61341000L, 61344001L); stopifnot(LUC[2] - LUC[1] + 1L == 3002L)
SG      <- data.table(name = c("sg1", "sg2", "sg3"), start = c(61331074L, 61339468L, 61347062L),
                      end = c(61331093L, 61339487L, 61347081L), strand = c("+", "-", "-"))
P2G_PK  <- c(61335772L, 61336272L)                                                  # snATAC P2G 核心 peak（ArchR 定宽 501 bp；05_e_E2_P2G8_peak清单.csv 第 4 行，r = 0.72）
stopifnot(P2G_PK[2] - P2G_PK[1] + 1L == 501L)
PEAK    <- c(61335782L, 61336248L); PEAK_ID <- "Peak_32739"; stopifnot(PEAK[2] - PEAK[1] + 1L == 467L)
EXON_GIVEN <- data.table(start = c(61334154L, 61335906L, 61338231L, 61341226L, 61342342L, 61345810L, 61350446L),
                         end   = c(61334354L, 61336016L, 61338326L, 61341348L, 61342446L, 61346002L, 61350568L))
stopifnot(all(c(LUC, SG$start, SG$end, P2G_PK, PEAK) >= E2[1]), all(c(LUC, SG$start, SG$end, P2G_PK, PEAK) <= E2[2]),
          P2G_PK[1] <= PEAK[1], P2G_PK[2] >= PEAK[2])                                   # P2G peak 完整包含 Peak_32739（重叠 467/467 bp）
SMP <- data.table(sample = c("ZF23_BHL_1", "ZF23_BHL_2", "ZF15_BHL_1", "ZF15_BHL_2"),
                  group = c("Con", "Con", "Mhp", "Mhp"), rep = c(1L, 2L, 1L, 2L))              # 显示顺序：Con 1、Con 2、Mhp 1、Mhp 2（自上而下）
SMP[, `:=`(file = file.path(BWD, paste0(sample, ".bw")), lab = paste(group, rep))]
stopifnot(all(file.exists(SMP$file)))

# ---------------------------------------------------------------- ① 注释复核：ITPR1 外显子（TxDb）与 CACGTG（BSgenome）
txdb <- AnnotationDbi::loadDb(F_TXDB)
win_gr <- GRanges(CHR_ENS, IRanges(WIN[1], WIN[2]))
ex_itpr1 <- reduce(subsetByOverlaps(exonsBy(txdb, "gene")[["ENSSSCG00000023437"]], win_gr))   # ITPR1（Ensembl 114）在窗口内的外显子（各转录本合并）
EXON <- data.table(start = start(ex_itpr1), end = end(ex_itpr1))[order(start)]
stopifnot(nrow(EXON) == 7L, identical(EXON$start, EXON_GIVEN$start), identical(EXON$end, EXON_GIVEN$end),
          unique(as.character(strand(ex_itpr1))) == "+")
EXON[, `:=`(n = seq_len(.N), width = end - start + 1L, overlaps_peak = pmin(end, PEAK[2]) - pmax(start, PEAK[1]) + 1L > 0)]
peak_exon_ov <- EXON[overlaps_peak == TRUE, sum(pmin(end, PEAK[2]) - pmax(start, PEAK[1]) + 1L)]
stopifnot(sum(EXON$overlaps_peak) == 1L, peak_exon_ov == 111L)
g <- BSgenome.Sscrofa.UCSC.susScr11
seq_e2 <- getSeq(g, GRanges(CHR_UCSC, IRanges(E2[1], E2[2])))[[1]]
ebox_e2 <- start(matchPattern("CACGTG", seq_e2)) + E2[1] - 1L                                   # CACGTG 回文，正反链同一位置
seq_pk <- getSeq(g, GRanges(CHR_UCSC, IRanges(PEAK[1], PEAK[2])))[[1]]
ebox_pk <- start(matchPattern("CACGTG", seq_pk)) + PEAK[1] - 1L
stopifnot(length(ebox_pk) == 1L, ebox_pk >= PEAK[1], ebox_pk + 5L <= PEAK[2], all(ebox_pk %in% ebox_e2))
EBOX <- c(ebox_pk, ebox_pk + 5L)
cat(sprintf("ITPR1 窗口内外显子 7 个（与任务给定一致）；E2 21 kb 内 CACGTG %d 个（%s），Peak_32739 内 1 个：chr13:%s–%s（距外显子 61,336,016 末端 +%d bp）\n",
            length(ebox_e2), paste(format(ebox_e2, big.mark = ","), collapse = ", "), format(EBOX[1], big.mark = ","), format(EBOX[2], big.mark = ","), EBOX[1] - 61336016L))

# ---------------------------------------------------------------- ② size factor、Peak_32739 计数（Fig 7 校正表）
sf  <- fread(F_SF); lib <- fread(F_LIB); pk <- fread(F_PEAK)[peak_id == PEAK_ID]; shp <- fread(F_SHAPE)
stopifnot(nrow(pk) == 1L, pk$Start == PEAK[1], pk$End == PEAK[2], pk$Chr == 13L,
          setequal(sf$sample, SMP$sample), setequal(lib$sample, SMP$sample), setequal(shp$sample, SMP$sample))
lib[, B_check := total_alignments / exp(mean(log(total_alignments)))]
stopifnot(max(abs(merge(lib, sf, by = "sample")[, B_check - B_libSize])) < 1e-3)                  # B 法 = 总比对数 / 几何均值（核对口径）
D <- merge(SMP, sf[, .(sample, A_medianRatio, B_libSize)], by = "sample")
D <- merge(D, lib[, .(sample, total_alignments, FRiP)], by = "sample")
raw_vec <- as.integer(unlist(pk[1, SMP$sample, with = FALSE])); names(raw_vec) <- SMP$sample
D[, raw := unname(raw_vec[sample])]                                                                 # 原始计数（featureCounts）
D[, norm := raw / A_medianRatio]                                                                  # DESeq2 median-of-ratios 归一化计数
D[, track_factor := if (TRACK_SCALE == "DESeq2") B_libSize / A_medianRatio else 1]
D <- D[order(match(sample, SMP$sample))]
stopifnot(identical(D[group == "Mhp", raw], c(209L, 155L)), identical(D[group == "Con", raw], c(61L, 77L)))
LFC <- log2(D[group == "Mhp", mean(norm)] / D[group == "Con", mean(norm)])
LFC_tab <- pk$A_medianRatio_log2FC; PADJ <- pk$A_medianRatio_padj; LFC_B <- pk$B_libSize_log2FC
stopifnot(abs(LFC - LFC_tab) < 0.01, PADJ < 1e-4, abs(LFC_tab - 2.27) < 0.01)
P_LAB <- "P < 0.0001"                                                                              # 图上不写 6.8e-20 这类精确值
cat(sprintf("Peak_32739 归一化计数：%s；log2FC（复算）= %.3f，表 %.3f；padj = %.2g → 图上 \"%s\"\n",
            paste(D[, sprintf("%s %.1f", lab, norm)], collapse = "、"), LFC, LFC_tab, PADJ, P_LAB))

# ---------------------------------------------------------------- ③ bigwig → 50 bp bin（区间均值），共用 y 上限
read_win <- function(f) { x <- import(f, which = win_gr); x <- restrict(x, WIN[1], WIN[2])
  v <- rep(0, WIN[2] - WIN[1] + 1L)
  for (k in seq_along(x)) v[(start(x)[k] - WIN[1] + 1L):(end(x)[k] - WIN[1] + 1L)] <- x$score[k]
  xp <- restrict(x, PEAK[1], PEAK[2])                                              # 峰内 10 bp 区间（合并同值 run）的未加权均值 = Fig 7 表 05 的 mean_in_peak 口径
  list(v = v, peak_int_mean = mean(xp$score)) }
bin_start <- seq(WIN[1], WIN[2], by = BIN); bin_end <- pmin(bin_start + BIN - 1L, WIN[2])
bin_idx <- findInterval(seq(WIN[1], WIN[2]), bin_start)
TR <- rbindlist(lapply(seq_len(nrow(SMP)), function(i) {
  r <- read_win(SMP$file[i]); v <- r$v
  data.table(sample = SMP$sample[i], group = SMP$group[i], rep = SMP$rep[i], lab = SMP$lab[i],
             bin_start = bin_start, bin_end = bin_end, RPKM = as.numeric(tapply(v, bin_idx, mean)),
             peak_mean_RPKM = mean(v[(PEAK[1] - WIN[1] + 1L):(PEAK[2] - WIN[1] + 1L)]),        # bp 加权（本脚本口径）
             peak_int_mean_RPKM = r$peak_int_mean, win_max_RPKM = max(v)) }))
stopifnot(nrow(TR) == 4L * length(bin_start), !anyNA(TR$RPKM))
# 峰内均值与 Fig 7 表 05_Peak32739_逐重复信号形状.csv 逐样本核对：该表的 mean_in_peak 是 10 bp 区间的未加权均值（同一 bigwig），
# 用同口径复算应一致（< 3%）；图与关键数字表采用 bp 加权均值（对峰内宽 run 不欠权），两者差 ≤ 12%，方向与倍数同
chk <- merge(unique(TR[, .(sample, peak_mean_RPKM, peak_int_mean_RPKM)]), shp[, .(sample, mean_in_peak)], by = "sample")
stopifnot(max(abs(chk$peak_int_mean_RPKM / chk$mean_in_peak - 1)) < 0.03, max(abs(chk$peak_mean_RPKM / chk$mean_in_peak - 1)) < 0.12,
          all(chk[grepl("ZF15", sample), peak_mean_RPKM] > max(chk[grepl("ZF23", sample), peak_mean_RPKM]) * 3))
TR <- merge(TR, D[, .(sample, track_factor)], by = "sample")[, value := RPKM * track_factor][order(match(sample, SMP$sample), bin_start)]
YMAX <- ceiling(max(TR$value) / 50) * 50
UNIT <- if (TRACK_SCALE == "RPKM") "RPKM" else "norm. cov."
cat(sprintf("轨迹口径 %s；50 bp bin 窗口最大 %s；共用 y 上限 0–%g\n", TRACK_SCALE,
            paste(TR[, .(m = max(value)), by = lab][, sprintf("%s %.0f", lab, m)], collapse = "、"), YMAX))
peak_ratio_bw <- D[, .(m = mean(chk$peak_mean_RPKM[match(sample, chk$sample)])), by = group][, m[group == "Mhp"] / m[group == "Con"]]

# ---------------------------------------------------------------- ④ 版面（y = mm，自下而上）
TH <- 4.6; G <- 0.5                                                                # 单条轨迹高、轨迹间隙
base_tr <- setNames((3:0) * (TH + G), SMP$lab)                                     # Con 1 最上、Mhp 2 最下（贴 x 轴）
y_tr_top <- max(base_tr) + TH                                                      # 19.9
y_m0 <- y_tr_top + 0.9; H_M <- 3.0                                                 # Peak_32739 行
y_e0 <- y_m0 + H_M + 0.6; H_E <- 4.2                                               # 元件行（Luc + sgRNA；标签在上）
y_g0 <- y_e0 + H_E + 0.6; H_G <- 2.4                                               # ITPR1 外显子行
y_b  <- y_g0 + H_G + 0.9                                                           # E2 括号线
H_TOP <- y_b + 2.4
W_BP <- WIN[2] - WIN[1]
npx <- function(bp) (bp - WIN[1]) / W_BP; npy <- function(mm) mm / H_TOP
S5 <- 5; S55 <- 5.5
lwd_pt <- function(pt) pt / 0.75                                                   # grid gpar lwd：1 lwd = 0.75 pt
# 文字 grob：位置 = 面板 npc + mm 偏移（面板宽度变化时左右边距文字仍贴边）
tg <- function(label, x_npc, y_npc, size = FS_SMALL, face = "plain", col = "black", hjust = 0.5, vjust = 0.5, rot = 0, dx = 0, dy = 0)
  textGrob(label, x = unit(x_npc, "npc") + unit(dx, "mm"), y = unit(y_npc, "npc") + unit(dy, "mm"), hjust = hjust, vjust = vjust, rot = rot,
           gp = gpar(fontsize = size, fontfamily = FONT, fontface = face, col = col))
full <- function(grob) annotation_custom(grob, xmin = WIN[1], xmax = WIN[2], ymin = 0, ymax = H_TOP)
lab_left <- function(label, y_mm, ...) full(tg(label, 0, npy(y_mm), hjust = 1, dx = -1.5, ...))
lab_right<- function(label, y_mm, ...) full(tg(label, 1, npy(y_mm), hjust = 0, dx = 0.8, ...))
# 固定毫米尺寸的小图形（sgRNA 箭头、P2G 三角）：放在以 x_bp 为中心、±1 kb 的视口内，用 mm 单位绘制，不随 bin 缩放
sg_arrow <- function(x_bp, strand, y_mm, len = 1.6, head = 0.65) {
  s <- if (strand == "+") 1 else -1
  gr <- segmentsGrob(x0 = unit(0.5, "npc") - unit(s * len / 2, "mm"), x1 = unit(0.5, "npc") + unit(s * len / 2, "mm"),
                     y0 = unit(0.5, "npc"), y1 = unit(0.5, "npc"),
                     arrow = arrow(length = unit(head, "mm"), type = "closed", angle = 28), gp = gpar(col = "black", fill = "black", lwd = lwd_pt(0.6)))
  annotation_custom(gr, xmin = x_bp - 1000, xmax = x_bp + 1000, ymin = y_mm - 1, ymax = y_mm + 1) }
tri_down <- function(x_bp, y_top_mm, h = 1.0, w = 1.2, col = "black") {
  gr <- polygonGrob(x = unit(0.5, "npc") + unit(c(-w / 2, w / 2, 0), "mm"), y = unit(1, "npc") - unit(c(0, 0, h), "mm"), gp = gpar(col = NA, fill = col))
  annotation_custom(gr, xmin = x_bp - 1000, xmax = x_bp + 1000, ymin = y_top_mm - 2, ymax = y_top_mm) }
# 阶梯多边形（每条轨迹一个 path，SVG 无栅格、无 bin 接缝）
step_poly <- function(dt, base) { n <- nrow(dt)
  data.table(x = c(dt$bin_start[1], as.vector(rbind(dt$bin_start, dt$bin_end + 1L)), dt$bin_end[n] + 1L),
             y = base + c(0, as.vector(rbind(dt$value, dt$value)), 0) / YMAX * TH) }
POLY <- TR[, step_poly(.SD, base_tr[[lab[1]]]), by = .(lab, group)]
stopifnot(all(POLY$y <= y_tr_top + 1e-9))
gene_y <- y_g0 + H_G / 2
gaps <- data.table(s = c(WIN[1], EXON$end + 1L), e = c(EXON$start - 1L, WIN[2]))[, mid := (s + e) / 2][e - s >= 2000]   # 链向 ">" 放在 ≥ 2 kb 的外显子间隙中点
x_key1 <- WIN[1] + 0.52 * W_BP; x_key2 <- WIN[1] + 0.79 * W_BP                     # Peak_32739 行右侧图内小键位置（P2G 键文字约 12 mm，CACGTG 键须在其右侧 ≥ 1 mm）
P2G_COL <- "#4D4D4D"                                                               # P2G 501 bp peak 外框色（与 Fig 6a 的 P2G 方框同色）
KB <- seq(61330000L, 61350000L, by = 5000L)

p_tr <- ggplot() +
  # Peak_32739 淡红带贯穿四条轨迹（到峰盒底部）
  annotate("rect", xmin = PEAK[1], xmax = PEAK[2] + 1L, ymin = 0, ymax = y_m0, fill = BH_COL, alpha = 0.13) +
  # 四条轨迹：灰基线 + 阶梯多边形
  annotate("segment", x = WIN[1], xend = WIN[2], y = base_tr, yend = base_tr, colour = "grey72", linewidth = lw(0.3)) +
  geom_polygon(data = POLY, aes(x, y, group = lab, fill = group), colour = NA) +
  scale_fill_manual(values = grp_col, guide = "none") +
  # E2 括号 + 标签
  annotate("segment", x = E2[1], xend = E2[2] + 1L, y = y_b, yend = y_b, colour = BH_COL, linewidth = lw(0.6)) +
  annotate("segment", x = c(E2[1], E2[2] + 1L), xend = c(E2[1], E2[2] + 1L), y = y_b, yend = y_b - 0.6, colour = BH_COL, linewidth = lw(0.6)) +
  annotate("text", x = mean(E2), y = y_b + 1.3, label = "E2 (21.4 kb)", family = FONT, fontface = "bold", colour = BH_COL, size = pt2size(S55)) +
  # ITPR1：基因线 + 外显子灰盒 + 链向
  annotate("segment", x = WIN[1], xend = WIN[2], y = gene_y, yend = gene_y, colour = "black", linewidth = lw(0.6)) +
  geom_rect(data = EXON, aes(xmin = start, xmax = end + 1L, ymin = y_g0 + 0.3, ymax = y_g0 + H_G - 0.3), fill = "#7F7F7F", colour = NA) +
  annotate("text", x = gaps$mid, y = gene_y, label = ">", family = FONT, colour = "grey35", size = pt2size(S5)) +
  # 元件行：Luc 3.0 kb 深灰盒 + 标签；sgRNA 箭头 + 标签
  annotate("rect", xmin = LUC[1], xmax = LUC[2] + 1L, ymin = y_e0 + 0.4, ymax = y_e0 + 2.0, fill = "#595959", colour = NA) +
  annotate("text", x = mean(LUC), y = y_e0 + 3.2, label = "Luc (3.0 kb)", family = FONT, size = pt2size(S5)) +
  sg_arrow(SG$start[1], SG$strand[1], y_e0 + 1.2) + sg_arrow(SG$start[2], SG$strand[2], y_e0 + 1.2) + sg_arrow(SG$start[3], SG$strand[3], y_e0 + 1.2) +
  geom_text(data = SG, aes(x = (start + end) / 2, y = y_e0 + 3.2, label = name), family = FONT, size = pt2size(S5)) +
  # Peak_32739 行：红盒 + P2G 501 bp peak 灰色外框（完整包含红盒）+ CACGTG 竖线（外框上方、同 x）；右侧图内小键
  annotate("rect", xmin = PEAK[1], xmax = PEAK[2] + 1L, ymin = y_m0, ymax = y_m0 + 1.4, fill = BH_COL, colour = NA) +
  annotate("rect", xmin = P2G_PK[1], xmax = P2G_PK[2] + 1L, ymin = y_m0 - 0.3, ymax = y_m0 + 1.7, fill = NA, colour = P2G_COL, linewidth = lw(0.5)) +
  annotate("segment", x = mean(EBOX), xend = mean(EBOX), y = y_m0 + 1.95, yend = y_m0 + 2.95, colour = "black", linewidth = lw(0.6)) +
  annotate("rect", xmin = x_key1 - 230, xmax = x_key1 + 230, ymin = y_m0 + 0.5, ymax = y_m0 + 2.5, fill = NA, colour = P2G_COL, linewidth = lw(0.5)) +
  full(tg("P2G core peak", npx(x_key1), npy(y_m0 + 1.5), size = S5, hjust = 0, dx = 1.3)) +                    # 501 bp 写在图注
  annotate("segment", x = x_key2, xend = x_key2, y = y_m0 + 0.6, yend = y_m0 + 2.4, colour = "black", linewidth = lw(0.6)) +
  full(tg("CACGTG", npx(x_key2), npy(y_m0 + 1.5), size = S5, hjust = 0, dx = 1.0)) +
  # 左侧轨道标签、右侧 y 范围
  lab_left(expression(italic(ITPR1)), gene_y) +
  lab_left(PEAK_ID, y_m0 + H_M / 2) +
  lab_left("Con 1", base_tr[["Con 1"]] + TH / 2) + lab_left("Con 2", base_tr[["Con 2"]] + TH / 2) +
  lab_left("Mhp 1", base_tr[["Mhp 1"]] + TH / 2) + lab_left("Mhp 2", base_tr[["Mhp 2"]] + TH / 2) +
  full(tg("BHLHE40 CUT&Tag", 0, npy(y_tr_top / 2), rot = 90, dx = -11.3)) +                          # 正体：抗体靶标是蛋白（与 Fig 4h/i、Fig 6f 一致）
  lab_right(sprintf("0–%g", YMAX), base_tr[["Con 1"]] + TH - 0.2, size = S5, vjust = 1) +
  lab_right(sprintf("0–%g", YMAX), base_tr[["Con 2"]] + TH - 0.2, size = S5, vjust = 1) +
  lab_right(sprintf("0–%g", YMAX), base_tr[["Mhp 1"]] + TH - 0.2, size = S5, vjust = 1) +
  lab_right(sprintf("0–%g", YMAX), base_tr[["Mhp 2"]] + TH - 0.2, size = S5, vjust = 1) +
  lab_right(UNIT, base_tr[["Mhp 2"]] + TH - 2.3, size = S5, vjust = 1) +
  scale_x_continuous(breaks = KB, labels = function(x) format(x / 1000, big.mark = ","), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(xlim = WIN, ylim = c(0, H_TOP), clip = "off") +
  labs(x = "chr13 position (kb)", y = NULL) +
  theme_mhp() +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(),
        axis.text.x = element_text(colour = "black"), plot.margin = margin(0.8, 6.5, 0.3, 14, unit = "mm"))

# ---------------------------------------------------------------- ⑤ 右侧小图：Peak_32739 逐重复归一化计数（三边空心柱 = 均值；n = 2 不画误差线）
PT_COL <- c(Con = "#6F93B3", Mhp = "#B03A2E")                                       # 点色比柱色深一档（smp_col 的深色样本），与 03_c 同法
D[, group := factor(group, levels = grp_levels)]
sm_in <- D[, .(mean = mean(norm), n = .N), by = group]
box_in <- sm_in[, .(x = as.integer(group) + c(-0.34, -0.34, 0.34, 0.34), y = c(0, mean, mean, 0)), by = group]
y_top <- max(D$norm); y_br <- y_top * 1.10
p_in <- ggplot() +
  geom_path(data = box_in, aes(x, y, colour = group, group = group), linewidth = lw(1.0), lineend = "butt", linejoin = "mitre") +
  scale_colour_manual(values = grp_col, guide = "none") +
  ggnewscale::new_scale_colour() +
  geom_point(data = D, aes(as.integer(group) + ifelse(rep == 1L, -0.13, 0.13), norm, colour = group), shape = 21, fill = "white", size = 1.4, stroke = 0.55) +
  scale_colour_manual(values = PT_COL, guide = "none") +
  annotate("segment", x = 1, xend = 2, y = y_br, yend = y_br, linewidth = lw(0.5)) +
  annotate("text", x = 1.5, y = y_br * 1.07, label = sprintf("italic(P) < '%s'", sub("^P < ", "", P_LAB)), parse = TRUE, family = FONT, size = pt2size(FS_SMALL), vjust = 0) +   # 字符串避免 0.0001 被 deparse 成 1e-04
  annotate("text", x = 1.5, y = y_br * 1.22, label = sprintf("log[2]*FC == %.2f", LFC_tab), parse = TRUE, family = FONT, size = pt2size(FS_SMALL), vjust = 0) +
  scale_x_continuous(breaks = 1:2, labels = grp_levels, limits = c(0.45, 2.55), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, NA), breaks = scales::breaks_extended(4), expand = expansion(mult = c(0, 0.08))) +   # 文字注释已计入数据范围，只留 8%
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = sprintf("Normalized counts\nat %s", PEAK_ID)) +
  theme_mhp() + theme(axis.text = element_text(colour = "black"), axis.title.y = element_text(size = FS_SMALL, lineheight = 0.95),
                      plot.margin = margin(0.8, 1, 0.3, 1.5, unit = "mm"))
b <- ggplot_build(p_in)$layout$panel_params[[1]]
stopifnot(b$y.range[1] == 0, b$y.range[2] > y_br * 1.22 + 15, b$y.range[2] < y_top * 1.6)   # 文字在面板内（clip off 仍留余量）且纵轴不空得过多

# ---------------------------------------------------------------- ⑥ 拼版（左轨迹弹性宽、右小图面板 16 mm 定宽）与输出
p <- p_tr + p_in + plot_layout(widths = unit(c(1, 16), c("null", "mm")))
sav(p, "Fig5h_E2_BHLHE40结合", PAN, 116, 44, panel_dir = PNL)

# Source Data / 注释表 / 关键数字
fwrite(TR[, .(sample, group, rep, chr = CHR_UCSC, bin_start, bin_end, RPKM = signif(RPKM, 5), track_factor = signif(track_factor, 4), plotted_value = signif(value, 5))],
       file.path(TAB, "09_h_轨迹_50bp.csv"))
ANN <- rbind(
  data.table(feature = "window", name = "E2 ± 1 kb", start = WIN[1], end = WIN[2], strand = ".", note = "轨迹窗口"),
  data.table(feature = "E2", name = "E2 (CRISPRi frame)", start = E2[1], end = E2[2], strand = ".", note = sprintf("%s bp；位于 SE1 内（Mhp SE1 61,209,361–61,517,999；Con SE1 61,209,367–61,406,081）", format(E2_LEN, big.mark = ","))),
  EXON[, .(feature = "ITPR1_exon", name = sprintf("ITPR1 exon (window #%d)", n), start, end, strand = "+", note = ifelse(overlaps_peak, sprintf("与 Peak_32739 重叠 %d bp", peak_exon_ov), ""))],
  data.table(feature = "luciferase", name = "Luc fragment", start = LUC[1], end = LUC[2], strand = ".", note = "3,002 bp；Fig. 5d 双荧光素酶实验所用的 E2 克隆片段"),
  SG[, .(feature = "sgRNA", name, start, end, strand, note = "dCas9-KRAB CRISPRi（Fig 6）")],
  data.table(feature = "P2G_core_peak", name = "P2G core accessibility peak (ArchR 501 bp)", start = P2G_PK[1], end = P2G_PK[2], strand = ".", note = sprintf("501 bp，中心 %s；完整包含 Peak_32739（重叠 467/467 bp）；r = 0.72（05_e_E2_P2G8_peak清单.csv）；图上灰色外框", format(as.integer(mean(P2G_PK)), big.mark = ","))),
  data.table(feature = "peak", name = PEAK_ID, start = PEAK[1], end = PEAK[2], strand = ".", note = sprintf("467 bp；raw %s；A 法 log2FC %.2f padj %.1e", paste(D[, sprintf("%s %d", lab, raw)], collapse = " / "), LFC_tab, PADJ)),
  data.table(feature = "E-box", name = "CACGTG", start = EBOX[1], end = EBOX[2], strand = "+/-", note = sprintf("BSgenome.Sscrofa.UCSC.susScr11 chr13；E2 21 kb 内唯一；距外显子 61,336,016 末端 +%d bp", EBOX[1] - 61336016L)))
fwrite(ANN, file.path(TAB, "09_h_注释元件.csv"))
fwrite(D[, .(sample, group, rep, label = lab, raw_count = raw, sizeFactor_A_medianRatio = signif(A_medianRatio, 5), sizeFactor_B_libSize = signif(B_libSize, 5),
             total_alignments, FRiP = signif(FRiP, 4), normalized_count = round(norm, 2), peak_mean_RPKM_bigwig = round(chk$peak_mean_RPKM[match(sample, chk$sample)], 2),
             track_factor = signif(track_factor, 4))],
       file.path(TAB, "09_h_Peak32739_归一化计数.csv"))
con_rpkm <- chk[match(SMP$sample, sample), peak_mean_RPKM]                                       # Con 1 / Con 2 / Mhp 1 / Mhp 2 峰内 bp 加权均值
legend_txt <- sprintf(paste0(
  "h, BHLHE40 CUT&Tag in 3D4/21 cells (Con and Mhp, two biological replicates each) across the 21.4-kb E2 interval (chr13:%s\u2013%s, plus 1-kb flanks). ",
  "Tracks show deepTools RPKM (library-size normalized; 50-bp bins; shared y-axis 0\u2013%g). On this shared scale the Con tracks are low but not empty: ",
  "mean RPKM within Peak_32739 is %.0f and %.0f in Con 1 and Con 2 versus %.0f and %.0f in Mhp 1 and Mhp 2. ",
  "Annotation rows: ITPR1 exons (grey boxes; + strand), the 3.0-kb luciferase fragment (Luc; panel d), ",
  "the three dCas9-KRAB sgRNAs used for CRISPRi (sg1 (+), sg2 (\u2212), sg3 (\u2212); arrows indicate strand; Fig. 6), Peak_32739 (chr13:%s\u2013%s; red bar and shaded band), ",
  "the snATAC peak-to-gene core peak linked to BHLHE40 (grey outline; chr13:%s\u2013%s, 501 bp, which fully contains Peak_32739) ",
  "and the single CACGTG E-box in E2 (tick; chr13:%s\u2013%s). Right, DESeq2 median-of-ratios-normalized counts at Peak_32739 per replicate (open circles; open bar = mean; n = 2, no error bar); ",
  "log2 fold change and Benjamini\u2013Hochberg-adjusted P value from the DESeq2 Wald test (Mhp versus Con; %s)."),
  format(E2[1], big.mark = ","), format(E2[2], big.mark = ","), YMAX, con_rpkm[1], con_rpkm[2], con_rpkm[3], con_rpkm[4],
  format(PEAK[1], big.mark = ","), format(PEAK[2], big.mark = ","), format(P2G_PK[1], big.mark = ","), format(P2G_PK[2], big.mark = ","),
  format(EBOX[1], big.mark = ","), format(EBOX[2], big.mark = ","), P_LAB)
key <- data.table(
  item = c("面板 / 尺寸", "轨迹窗口", "E2 区间", "轨迹数据", "分组对应", "轨迹刻度口径（TRACK_SCALE）", "RPKM → DESeq2 换算因子（B_libSize / A_medianRatio）", "bin", "共用 y 上限",
           "窗口内 50 bp bin 最大值", "Peak_32739 坐标", "Peak_32739 原始计数（featureCounts）", "size factor A（median-of-ratios）", "Peak_32739 归一化计数（raw / A）",
           "log2FC（A 法，DESeq2）", "log2FC（B 文库大小法，对照）", "padj（A 法）", "图上统计标注", "Peak_32739 峰内 RPKM 均值（bigwig 直接算）", "峰内 RPKM 均值 Mhp/Con",
           "CACGTG（E-box）", "P2G 核心可及性 peak", "ITPR1 外显子（窗口内 7 个，TxDb Ensembl 114 复核）", "Peak_32739 与外显子重叠", "荧光素酶片段", "sgRNA（dCas9-KRAB）",
           "n / 统计", "对照文库", "文库统计", "Con 轨迹共用刻度说明", "拼版归属", "图注建议句"),
  value = c("Fig 5 h，116 × 44 mm；左轨迹图 + 右小图（patchwork，小图面板定宽 16 mm）",
            sprintf("chr13:%s–%s（E2 ± 1 kb，%s bp）", format(WIN[1], big.mark = ","), format(WIN[2], big.mark = ","), format(WIN[2] - WIN[1] + 1L, big.mark = ",")),
            sprintf("chr13:%s–%s（%s bp；SE1 内、ITPR1 3′ 区）", format(E2[1], big.mark = ","), format(E2[2], big.mark = ","), format(E2_LEN, big.mark = ",")),
            "ZF23_BHL_1/2（Con）、ZF15_BHL_1/2（Mhp）deepTools bamCoverage RPKM、10 bp bin、去重 BAM（09_TF_CUTTag_bw探针 核对）",
            "ZF15 = Mhp、ZF23 = Con（3D4/21 生物学重复各 2）",
            sprintf("%s（默认 RPKM：只校深度、两种口径中更保守；改 \"DESeq2\" 则轨迹乘 B/A、单位 norm. cov.）", TRACK_SCALE),
            paste(D[, sprintf("%s %.3f", lab, B_libSize / A_medianRatio)], collapse = "；"),
            sprintf("%d bp 区间均值（68 mm 画 23.4 kb ≈ 0.145 mm/bin ≈ 0.41 pt）", BIN),
            sprintf("0–%g %s（四条轨迹共用）", YMAX, UNIT),
            paste(TR[, .(m = max(value)), by = lab][, sprintf("%s %.1f", lab, m)], collapse = "；"),
            sprintf("chr13:%s–%s（467 bp）", format(PEAK[1], big.mark = ","), format(PEAK[2], big.mark = ",")),
            paste(D[, sprintf("%s %d", lab, raw)], collapse = "；"),
            paste(D[, sprintf("%s %.3f", lab, A_medianRatio)], collapse = "；"),
            paste(D[, sprintf("%s %.1f", lab, norm)], collapse = "；"),
            sprintf("%.2f（复算 %.3f；Mhp 均值 %.1f / Con 均值 %.1f）", LFC_tab, LFC, D[group == "Mhp", mean(norm)], D[group == "Con", mean(norm)]),
            sprintf("%.2f", LFC_B), sprintf("%.2e", PADJ), sprintf("log2FC = %.2f；\"%s\"（BH 校正 P = %.1e）", LFC_tab, P_LAB, PADJ),
            sprintf("%s（bp 加权；Fig 7 表 05 mean_in_peak 为 10 bp 区间未加权均值：%s，同口径复算差 < 3%%）",
                    paste(chk[match(SMP$sample, sample), sprintf("%s %.1f", SMP$lab, peak_mean_RPKM)], collapse = "；"),
                    paste(chk[match(SMP$sample, sample), sprintf("%.1f", mean_in_peak)], collapse = " / ")),
            sprintf("×%.2f", peak_ratio_bw),
            sprintf("E2 21 kb 内 %d 个：chr13:%s–%s（Peak_32739 内；BSgenome chr 前缀；回文，两链同位）", length(ebox_e2), format(EBOX[1], big.mark = ","), format(EBOX[2], big.mark = ",")),
            sprintf("chr13:%s–%s（ArchR 定宽 501 bp peak，中心 %s；完整包含 Peak_32739，重叠 467/467 bp；r = 0.72；图上灰色外框，与 Fig 6a 的 501 bp 方框一致）", format(P2G_PK[1], big.mark = ","), format(P2G_PK[2], big.mark = ","), format(as.integer(mean(P2G_PK)), big.mark = ",")),
            paste(EXON[, sprintf("%s–%s", format(start, big.mark = ","), format(end, big.mark = ","))], collapse = "; "),
            sprintf("外显子 61,335,906–61,336,016 重叠 %d bp（467 bp 的 %.0f%%）；CACGTG 在外显子末端下游 %d bp", peak_exon_ov, peak_exon_ov / 467 * 100, EBOX[1] - 61336016L),
            sprintf("chr13:%s–%s（3,002 bp，E2 的 14%%）", format(LUC[1], big.mark = ","), format(LUC[2], big.mark = ",")),
            paste(SG[, sprintf("%s %s–%s (%s)", name, format(start, big.mark = ","), format(end, big.mark = ","), strand)], collapse = "；"),
            "n = 2 生物学重复 / 组；小图只画均值柱 + 点，不画误差线；统计 = DESeq2 Wald（Fig 7 脚本 02，A 法）",
            "RPKM（按测序深度归一化）",
            paste(D[, sprintf("%s %.2f M 比对、FRiP %.1f%%", lab, total_alignments / 1e6, FRiP * 100)], collapse = "；"),
            sprintf("共用 0–%g 刻度下 Con 峰高仅约 %.2f–%.2f mm（轨迹高 %.1f mm）；图注写明 Con 峰内均值 %.0f / %.0f RPKM vs Mhp %.0f / %.0f，避免读成「Con 无结合」",
                    YMAX, min(con_rpkm[1:2]) / YMAX * TH, max(con_rpkm[1:2]) / YMAX * TH, TH, con_rpkm[1], con_rpkm[2], con_rpkm[3], con_rpkm[4]),
            "Fig 5 面板 h",
            legend_txt))
fwrite(key, file.path(TAB, "09_h_关键数字.csv"))
for (f in c("09_h_关键数字.csv", "09_h_轨迹_50bp.csv", "09_h_注释元件.csv", "09_h_Peak32739_归一化计数.csv")) stopifnot(file.size(file.path(TAB, f)) > 0)
cat("✓ Tables: 09_h_关键数字.csv / 09_h_轨迹_50bp.csv / 09_h_注释元件.csv / 09_h_Peak32739_归一化计数.csv\n")
print(D[, .(lab, raw, A = round(A_medianRatio, 3), norm = round(norm, 1), peak_RPKM = round(chk$peak_mean_RPKM[match(sample, chk$sample)], 1))])
