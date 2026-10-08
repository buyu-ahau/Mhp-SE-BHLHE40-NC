# ============================================================================
# 02_b_BHLHE40位点多组学_本地.R — Fig. 5c：BHLHE40 位点多组学轨迹（文件名中的面板字母为旧编号）
# 运行位置：本地；路径变量见 config/paths.R
# ============================================================================

source(file.path(Sys.getenv("MHP_CODE", "."), "config", "paths.R"))

suppressPackageStartupMessages({
  library(data.table); library(rhdf5); library(GenomicRanges); library(rtracklayer); library(ggplot2); library(scales)
})
source(file.path(MHP_CODE, "shared", "_palette_MHP.R"))

# ---------------------------------------------------------------- 路径
HIC   <- MHP_HIC
COOL  <- c(Con = file.path(HIC, "single_result/02.align/_Con.5000.cool"), Mhp = file.path(HIC, "single_result/02.align/Mhp.5000.cool"))
LOOPF <- c(Con = file.path(HIC, "single_result/05.loop/_Con.combine.bedpe"), Mhp = file.path(HIC, "single_result/05.loop/Mhp.combine.bedpe"))
BWF   <- list(H3K27ac = setNames(file.path(MHP_RES, "2、H3K27ac结果/4、bw文件", sprintf("27ac%s%d.CPM.bw", rep(grp_levels, each = 3), 1:3)),
                                 paste0(rep(grp_levels, each = 3), "_", 1:3)),
              H3K4me1 = setNames(file.path(MHP_RES, "3、H3K4me1结果/4、bw文件", sprintf("me1%s%d.CPM.bw", rep(grp_levels, each = 3), 1:3)),
                                 paste0(rep(grp_levels, each = 3), "_", 1:3)))
SEBED <- c(Con = file.path(MHP_ROSE, "Con_SuperEnhancers.bed"), Mhp = file.path(MHP_ROSE, "Mhp_SuperEnhancers.bed"))
GTF   <- file.path(MHP_RES, "10.验证结果/7.转录因子CUTtag/Sus_scrofa.Sscrofa11.1.114.gtf.gz")
P2GF  <- file.path(MHP_RES, "9.scatac 补充分析/Tables/24_SE/P2G_links.csv")
FIG5  <- file.path(MHP_RES, "7、HI-C结果/新Fig5_作图")
TAB_IN <- file.path(FIG5, "Tables")                          # 00 探针产物与服务器下载文件所在（只读）
TEST_OUT <- Sys.getenv("MHP_FIG5_TEST_OUT")                  # ★仅供开发自测：设了就把全部输出改写到该目录（scratchpad），项目目录不动
OUT_ROOT <- if (nzchar(TEST_OUT)) TEST_OUT else FIG5
PAN <- file.path(OUT_ROOT, "Plots/面板"); TAB <- file.path(OUT_ROOT, "Tables"); PNL <- file.path(OUT_ROOT, "_panels")
for (d in c(PAN, TAB, PNL)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
W00   <- file.path(TAB_IN, "00_HiC探针_chr13_ICE权重.csv")   # 00 探针产物（先跑 00）
PX00  <- file.path(TAB_IN, "00_HiC探针_窗口像素.csv")
FRAGD <- if (nzchar(Sys.getenv("MHP_FIG5_TEST_FRAGD"))) Sys.getenv("MHP_FIG5_TEST_FRAGD") else file.path(TAB_IN, "BHLHE40位点子集_chr13_61.0-61.8Mb")
RNAD  <- file.path(TAB_IN, "RNA_bigwig")                      # 公司 RNA bigwig（{Con,Mhp}{1,2,3}_RNA.bw；不存在 → 不画）
META  <- file.path(MHP_RES, "9.scatac 补充分析/Tables/21_QC_final/cellColData_final.csv")   # 9,610 核：Sample / group / celltype / nFrags（snATAC 轨迹归一化用）
P2G8F <- file.path(TAB_IN, "05_e_E2_P2G8_peak清单.csv")        # 05_e 的 8 个 E2 内 P2G peak 坐标（E2 核对用；先跑 05）
E05F  <- file.path(TAB_IN, "05_e_E2_P2G8_作图数据.csv")        # 05_e 的 P2G8 CPM（Con / Mhp × 7 类型；E2 方向 / 幅度核对的参照）
stopifnot(file.exists(c(COOL, LOOPF, unlist(BWF), SEBED, GTF, P2GF, W00, PX00, META, P2G8F, E05F)))
HAS_FRAG <- dir.exists(FRAGD) && length(list.files(FRAGD, pattern = "\\.tsv(\\.gz)?$")) >= 6
HAS_RNA  <- dir.exists(RNAD) && length(list.files(RNAD, pattern = "\\.bw$")) >= 2
if (nzchar(TEST_OUT)) cat("测试模式：输出改写到", OUT_ROOT, "；fragments 目录", FRAGD, "\n")
BH_COL <- "#B2182B"                                           # 与 Fig 4 c/d/e 一致的 BHLHE40 高亮红
GREY   <- "#8C8C8C"

# ---------------------------------------------------------------- 参数：窗口与元件（chr13，1-based 闭区间）
CHR <- "chr13"; CHR_ENS <- "13"; BIN <- 5000L
WIN <- c(61100001L, 61700000L)
TSS <- 61526116L
PROM <- c(61520001L, 61540000L)                               # 方案 §2.2 "连启动子"加粗口径
EL <- data.table(name = c("E1", "E2", "E3", "E4"), start = c(61301100L, 61341000L, 61367000L, 61383500L),
                 end = c(61304101L, 61344001L, 61370001L, 61386501L))
E2_21 <- c(61329246L, 61350682L)
SE_NAMED <- data.table(group = c("Con", "Mhp", "Con", "Mhp"), name = c("SE1", "SE1", "SE2", "SE2"),
                       bed_start = c(61209367L, 61209361L, 61531856L, 61531355L), bed_end = c(61406081L, 61517999L, 61581260L, 61583348L))
SIG_BIN <- 500L                                               # 信号轨迹 bin
HIC_DEPTH <- 350000L                                          # 三角热图显示的最大距离
LOOP_BOLD <- "promoter"                                       # 弧线粗细口径："promoter" = 任一锚点与 PROM 重叠的加粗（方案 §2.2）；"none" = 全部同粗（更稳妥的备选）
stopifnot(LOOP_BOLD %in% c("promoter", "none"))
ATAC_BIN <- 500L; ATAC_GAP <- 0.3
                                                              #   也窄于 300 dpi 的一个像素；500 bp = 0.15 mm，与 CUT&Tag 轨迹同 bin）、类型行之间的行距（mm）
ATAC_MODE <- "stacked"
ATAC_H <- switch(ATAC_MODE, stacked = 4.8, mirror = 4.4, overlay = 3.5)
ATAC_SUBGAP <- 0.3                                             # stacked 时两条子行之间的空隙（mm）
ATAC_YCAP <- 25
ATAC_NORM <- "nFrags"                                          # snATAC 伪 bulk 归一化分母："nFrags"（该组该类型全部核 nFrags 之和 × 2 = 全基因组插入位点数）或 "ReadsInTSS"（同法用 ReadsInTSS 之和）
B_CELL_MODE <- "omit"                                          # B 细胞（Mhp 仅 3 核）："omit" = 不画 + 右缘小注（与面板 f / g 一致）；"grey" = 画但两组都灰
stopifnot(ATAC_NORM %in% c("nFrags", "ReadsInTSS"), B_CELL_MODE %in% c("omit", "grey"), ATAC_MODE %in% c("stacked", "mirror", "overlay"), is.na(ATAC_YCAP) || ATAC_YCAP > 0)
ct_atac_draw <- if (B_CELL_MODE == "omit") setdiff(ct_atac_order, "B_cells") else ct_atac_order   # 实际画的类型行（顺序 = ct_atac_order）
E2_LEN <- E2_21[2] - E2_21[1] + 1L; stopifnot(E2_LEN == 21437L)
fmt <- function(x, d = 0) formatC(x, format = "f", digits = d, big.mark = ",")
ov  <- function(s1, e1, s2, e2) s1 <= e2 & s2 <= e1
save_csv <- function(dt, name) {                              # 中文文件名只由 R 写出并核验大小
  f <- file.path(TAB, name); fwrite(dt, f); sz <- file.info(f)$size
  stopifnot(!is.na(sz), sz > 0); cat(sprintf("  [写出] %s  %s 行  %s 字节\n", name, fmt(nrow(dt)), fmt(sz))) }
KEY <- list(); key <- function(item, value, note = "") KEY[[length(KEY) + 1]] <<- data.table(item = item, value = as.character(value), note = note)
key("window", sprintf("%s:%s-%s", CHR, fmt(WIN[1]), fmt(WIN[2])), "1-based；与 ArchR 轨迹 PDF 一致")
key("BHLHE40_TSS", fmt(TSS), "Ensembl 114 / P2G gene_tss")
key("promoter_bold_window", sprintf("%s-%s", fmt(PROM[1]), fmt(PROM[2])), "loop 任一锚点与之重叠即加粗")

rows <- list(data.table(id = "gene", h = 6.0, gap = if (HAS_RNA) 1.0 else 2.4))
if (HAS_RNA) rows <- c(rows, list(data.table(id = c("rna_Con", "rna_Mhp"), h = c(3.2, 3.2), gap = c(0.3, 2.4))))
rows <- c(rows, list(data.table(id = c("se", "elem"), h = c(4.2, 4.2), gap = c(0.3, 1.0))))
rows <- c(rows, list(data.table(id = c("k27_Con", "k27_Mhp", "k4_Con", "k4_Mhp", "loop"), h = c(4.2, 4.2, 3.2, 3.2, 10), gap = c(0.3, 0.9, 0.3, 0.9, 0.9)),
                     data.table(id = "hic", h = if (HAS_FRAG) 24 else 28, gap = 0.9)))
rows <- c(rows, list(if (HAS_FRAG) data.table(id = paste0("atac_", ct_atac_draw), h = ATAC_H, gap = c(rep(ATAC_GAP, length(ct_atac_draw) - 1), 0.9))
                     else data.table(id = "atac", h = 4.0, gap = 0.9)))
rows <- c(rows, list(data.table(id = "p2g", h = 7.5, gap = 0)))
TR <- rbindlist(rows)
stopifnot(!anyDuplicated(TR$id))
TR[, top := rev(cumsum(rev(h + gap))) ][, bottom := top - h]
H_ALL <- TR$top[1]
tr <- function(tid, what = "top") { r <- TR[TR$id == tid, ]; switch(what, top = r$top, bottom = r$bottom, mid = (r$top + r$bottom) / 2) }
ATAC_IDS <- TR[grepl("^atac", id), id]                          # 占位版 "atac"；补轨版 atac_<8 类>
X_LAB <- WIN[1] - 6000                                        # 左边距标签的 x（hjust = 1）
W_MM <- 180; LEFT_MM <- 19; RIGHT_MM <- 5
BP_MM <- (WIN[2] - WIN[1] + 1) / (W_MM - LEFT_MM - RIGHT_MM)   # 数据区 156 mm 画 600 kb ≈ 3,846 bp / mm（与下方 plot.margin 一致）
RNG_DX <- 0.4; RNG_COL <- "grey40"
S65 <- pt2size(6.5); S6 <- pt2size(6); S55 <- pt2size(5.5); S5 <- pt2size(5)
ink <- function(s, pt) { g <- systemfonts::glyph_info(strsplit(gsub("\n", "", s), "")[[1]], family = FONT, size = pt, res = 7200)
  bb <- do.call(rbind, g$bbox); c(top = max(bb[, "ymax"]), bot = min(bb[, "ymin"])) * 25.4 / 7200 }
capM <- function(pt) ink("M", pt)[["top"]]                      # 'M' 字高（mm）
lsp  <- function(pt, lh = 0.9) lh * 1.2 * pt * 25.4 / 72        # 多行文字的实际行距（mm）
strw <- function(s, pt) systemfonts::string_width(s, family = FONT, size = pt, res = 7200) * 25.4 / 7200   # 字串宽（mm）

# ============================================================================ ① 基因模型（GTF，只读 chr13 窗口附近）
cmd <- paste0("gzcat '", GTF, "' | awk -F'\\t' '$1==\"", CHR_ENS, "\" && ($3==\"gene\" || $3==\"exon\") && $5>=", WIN[1] - 100000L, " && $4<=", WIN[2] + 100000L, "'")
gtf <- fread(cmd = cmd, header = FALSE, sep = "\t", quote = "", col.names = c("chr", "src", "type", "start", "end", "score", "strand", "frame", "attr"))
attr_get <- function(a, k) { m <- regmatches(a, regexpr(paste0(k, ' "[^"]+"'), a)); out <- rep(NA_character_, length(a)); hit <- grepl(paste0(k, ' "'), a); out[hit] <- sub(paste0(k, ' "([^"]+)"'), "\\1", m); out }
gtf[, `:=`(gene_id = attr_get(attr, "gene_id"), gene_name = attr_get(attr, "gene_name"), biotype = attr_get(attr, "gene_biotype"))]
genes <- gtf[type == "gene" & ov(start, end, WIN[1], WIN[2]), .(gene_id, gene_name, biotype, strand, start, end)]
genes[, label := fifelse(!is.na(gene_name) & biotype == "protein_coding", gene_name, NA_character_)]
genes[, lane := fifelse(!is.na(label), 1L, 2L)]
stopifnot(c("ITPR1", "BHLHE40", "ARL8B") %in% genes$label, genes[gene_name == "BHLHE40", start] == TSS)
ex <- gtf[type == "exon" & gene_id %in% genes$gene_id]
ex <- as.data.table(reduce(split(GRanges(CHR, IRanges(ex$start, ex$end)), ex$gene_id)))[, .(gene_id = as.character(group_name), start, end)]
ex <- merge(ex, genes[, .(gene_id, lane, label)], by = "gene_id")
genes[, `:=`(vis_start = pmax(start, WIN[1]), vis_end = pmin(end, WIN[2]))]
cat("① 基因模型：窗口内基因", nrow(genes), "个；标名", paste(genes[!is.na(label), label], collapse = "/"),
    "；EDEM1 起点 61,703,621 在窗口外 3.6 kb，不入图\n")
key("genes_in_window", nrow(genes), paste(genes[!is.na(label), sprintf("%s %s-%s(%s)", label, fmt(start), fmt(end), strand)], collapse = "; "))
save_csv(genes[, .(gene_id, gene_name, biotype, strand, start, end, lane, labelled = !is.na(label))], "02_b_基因模型.csv")

# 画图数据：基因线、外显子、名字、链方向 ">" 记号
LANE_Y <- c(tr("gene", "top") - 3.6, tr("gene", "bottom") + 0.8)     # 第 1 行（有名基因）中心、第 2 行（其他）中心
EXH <- c(1.5, 0.9)                                                     # 外显子盒高
genes[, y := LANE_Y[lane]]; ex[, y := LANE_Y[lane]]
gene_line <- genes[, .(x = vis_start, xend = vis_end, y, col = fifelse(lane == 1L, "black", GREY))]
exon_box  <- ex[, .(xmin = pmax(start, WIN[1]), xmax = pmin(end, WIN[2]), ymin = y - EXH[lane] / 2, ymax = y + EXH[lane] / 2,
                    col = fifelse(lane == 1L, fifelse(label == "BHLHE40", BH_COL, "black"), GREY))][xmax > xmin]
gene_lab  <- genes[!is.na(label), .(x = (vis_start + vis_end) / 2, y = y + 1.55, label,
                                    col = fifelse(label == "BHLHE40", BH_COL, "black"), face = fifelse(label == "BHLHE40", "bold.italic", "italic"))]
# 链方向记号：每 60 kb 一个 ">"（只对可见长度 ≥ 30 kb 的有名基因）
chev <- genes[!is.na(label) & vis_end - vis_start >= 30000, .(x = seq(vis_start + 15000, vis_end - 15000, by = 60000), y = y,
                                                                lab = fifelse(strand == "+", ">", "<")), by = gene_id]

# ============================================================================ ② ROSE SE 条 + E1–E4
se <- rbindlist(lapply(grp_levels, function(g) { x <- fread(SEBED[[g]], header = FALSE)[, 1:4]; setnames(x, c("chr", "bed_start", "bed_end", "id")); x[, group := g] }))
se <- se[chr == CHR & ov(bed_start + 1L, bed_end, WIN[1], WIN[2])]
se <- merge(se, SE_NAMED, by = c("group", "bed_start", "bed_end"), all.x = TRUE)
stopifnot(nrow(se) == 7, se[!is.na(name), .N] == 4)                    # Con 3 段、Mhp 4 段；SE1/SE2 两组各一
SE_LANE <- c(Con = tr("se", "top") - 0.75, Mhp = tr("se", "bottom") + 0.75)
# 左标签 "SE Con" / "SE Mhp" 放在各自 SE 条中心（vjust 0.5 → 'M' 字高居中）；墨迹间距（Con 行字底 − Mhp 行字顶）须 ≥ 0.5 mm
SE_LAB_INK <- rbind(Con = SE_LANE[["Con"]] - capM(6.5) / 2 + ink("SE Con", 6.5), Mhp = SE_LANE[["Mhp"]] - capM(6.5) / 2 + ink("SE Mhp", 6.5))
SE_LAB_GAP <- SE_LAB_INK["Con", "bot"] - SE_LAB_INK["Mhp", "top"]
stopifnot(SE_LAB_GAP >= 0.5)
cat(sprintf("② SE 左标签：两条 SE 条中心距 %.2f mm；\"SE Con\" 字底与 \"SE Mhp\" 字顶墨迹间距 %.2f mm（≥ 0.5）；标签大写字高中心 = SE 条中心\n",
            SE_LANE[["Con"]] - SE_LANE[["Mhp"]], SE_LAB_GAP))
key("layout_SE_label_gap_mm", sprintf("%.2f", SE_LAB_GAP), sprintf("SE 条中心距 %.2f mm；左标签 6.5 pt 按 'M' 字高居中于各自 SE 条；墨迹间距 ≥ 0.5 mm（stopifnot）", SE_LANE[["Con"]] - SE_LANE[["Mhp"]]))
se[, `:=`(xmin = pmax(bed_start + 1L, WIN[1]), xmax = pmin(bed_end, WIN[2]), y = SE_LANE[group], col = grp_col[group])]
se_lab <- se[group == "Mhp" & !is.na(name), .(x = (xmin + xmax) / 2, y = tr("se", "top") + 0.3, label = name)]
cat("② SE：窗口内 Con", se[group == "Con", .N], "段、Mhp", se[group == "Mhp", .N], "段；命名 SE1/SE2 按 Fig 4；旁邻 SE（61.04–61.19 Mb 两组、61.67–61.77 Mb 仅 Mhp）画出不标名\n")
save_csv(rbind(se[, .(track = "SE", group, name = fifelse(is.na(name), "(unnamed neighbour)", name), bed_start, bed_end, start_1based = bed_start + 1L, end = bed_end, id)],
               EL[, .(track = "E", group = "-", name, bed_start = start - 1L, bed_end = end, start_1based = start, end, id = "primer BSgenome (§28)")],
               data.table(track = "E2_21kb", group = "-", name = "E2 frame", bed_start = E2_21[1] - 1L, bed_end = E2_21[2], start_1based = E2_21[1], end = E2_21[2], id = sprintf("CRISPRi %s bp", fmt(E2_LEN)))),
         "02_b_SE与元件.csv")
key("SE1", paste(se[!is.na(name) & name == "SE1", sprintf("%s %s-%s", group, fmt(bed_start), fmt(bed_end))], collapse = "; "), "ROSE bed 0-based 起点")
key("SE2", paste(se[!is.na(name) & name == "SE2", sprintf("%s %s-%s", group, fmt(bed_start), fmt(bed_end))], collapse = "; "))
key("neighbour_SEs_drawn_unlabelled", paste(se[is.na(name), sprintf("%s %s-%s", group, fmt(bed_start), fmt(bed_end))], collapse = "; "))
key("E1-E4", paste(EL[, sprintf("%s %s-%s", name, fmt(start), fmt(end))], collapse = "; "), "3,002 bp 双荧光素酶片段")
key("E2_21kb_frame", sprintf("%s-%s", fmt(E2_21[1]), fmt(E2_21[2])), sprintf("%s bp（1-based 闭区间）", fmt(E2_LEN)))
EY <- c(tr("elem", "bottom") + 0.3, tr("elem", "bottom") + 1.7)            # E 盒 y 范围
EL[, `:=`(fill = fifelse(name == "E2", BH_COL, "#5A5A5A"))]
el_lab <- EL[, .(x = (start + end) / 2, y = EY[2] + 0.35, label = name, col = fill, face = fifelse(name == "E2", "bold", "plain"))]

# ============================================================================ ③ CUT&Tag 组均值轨迹（500 bp bin）
gr_win <- GRanges(CHR, IRanges(WIN[1], WIN[2]))
n_bin <- (WIN[2] - WIN[1] + 1L) %/% SIG_BIN; stopifnot(n_bin * SIG_BIN == WIN[2] - WIN[1] + 1L)
bin_bw <- function(f) {
  x <- import(f, which = gr_win); x <- restrict(x, start(gr_win), end(gr_win))
  v <- rep(0, WIN[2] - WIN[1] + 1L)
  for (k in seq_along(x)) v[(start(x)[k] - WIN[1] + 1L):(end(x)[k] - WIN[1] + 1L)] <- x$score[k]   # 公司 bigwig 区间不重叠，直接赋值
  colMeans(matrix(v, nrow = SIG_BIN)) }
SIG <- rbindlist(lapply(names(BWF), function(m) rbindlist(lapply(names(BWF[[m]]), function(s)
  data.table(mark = m, sample = s, group = sub("_.*", "", s), bin = seq_len(n_bin), cpm = bin_bw(BWF[[m]][[s]]))))))
SIGM <- SIG[, .(cpm = mean(cpm)), by = .(mark, group, bin)][, `:=`(start = WIN[1] + (bin - 1L) * SIG_BIN, end = WIN[1] + bin * SIG_BIN - 1L)]
SIGM[, ymax := ceiling(max(cpm) * 10) / 10, by = mark]                # 同 mark 两组共用 y 上限
YMAX <- SIGM[, .(ymax = ymax[1]), by = mark]
cat("③ 信号轨迹：", n_bin, "个 500 bp bin；y 上限 H3K27ac", YMAX[mark == "H3K27ac", ymax], " H3K4me1", YMAX[mark == "H3K4me1", ymax], "CPM\n")
# 方向自检：SE1（Mhp 区间）上 Mhp/Con 组均值比应 > 1（Fig 4 f：H3K27ac ×1.78、H3K4me1 ×1.31）
se1 <- SE_NAMED[group == "Mhp" & name == "SE1"]
chk <- SIGM[start >= se1$bed_start + 1L & end <= se1$bed_end, .(m = mean(cpm)), by = .(mark, group)]
ratio <- chk[, .(ratio = m[group == "Mhp"] / m[group == "Con"]), by = mark]
stopifnot(all(ratio$ratio > 1))
cat("   SE1 上 Mhp/Con：", paste(ratio[, sprintf("%s ×%.2f", mark, ratio)], collapse = "，"), "\n")
key("signal_bin_bp", SIG_BIN); key("signal_ymax_H3K27ac_CPM", YMAX[mark == "H3K27ac", ymax]); key("signal_ymax_H3K4me1_CPM", YMAX[mark == "H3K4me1", ymax])
key("SE1_Mhp_over_Con_H3K27ac", sprintf("%.3f", ratio[mark == "H3K27ac", ratio]), "500 bp bin 组均值在 Mhp SE1 上的均值比（自检，Fig 4 f 为 1.78）")
key("SE1_Mhp_over_Con_H3K4me1", sprintf("%.3f", ratio[mark == "H3K4me1", ratio]))
save_csv(dcast(SIGM, mark + bin + start + end ~ group, value.var = "cpm"), "02_b_信号轨迹_500bp组均值.csv")
# 画图多边形（每个 bin 一个矩形；面积图）
SIG_TR <- c(H3K27ac_Con = "k27_Con", H3K27ac_Mhp = "k27_Mhp", H3K4me1_Con = "k4_Con", H3K4me1_Mhp = "k4_Mhp")
SIGM[, tid := SIG_TR[paste0(mark, "_", group)]]
SIGM[, `:=`(y0 = TR$bottom[match(tid, TR$id)], hh = TR$h[match(tid, TR$id)])]
sig_rect <- SIGM[cpm > 0, .(xmin = start, xmax = end + 1L, ymin = y0, ymax = y0 + hh * pmin(cpm / ymax, 1), col = grp_col[group])]
sig_lab  <- unique(SIGM[, .(tid, mark, group, y0, hh, ymax)])[, .(tid, x = X_LAB, y = y0 + hh / 2, label = sprintf("%s %s", mark, group),
                                                                   label_r = sprintf("[0–%g]", ymax))]

# ============================================================================ ④ loop 弧（Con 上 / Mhp 下）
read_bedpe6 <- function(f) { x <- fread(f, header = FALSE)[, 1:6]; setnames(x, c("chr1", "s1", "e1", "chr2", "s2", "e2")); x }
LP <- rbindlist(lapply(grp_levels, function(g) read_bedpe6(LOOPF[[g]])[, group := g]))
stopifnot(LP[group == "Con", .N] == 22441, LP[group == "Mhp", .N] == 19203)
LP <- LP[chr1 == CHR & chr2 == CHR & (ov(s1 + 1L, e1, WIN[1], WIN[2]) | ov(s2 + 1L, e2, WIN[1], WIN[2]))]
LP[, `:=`(c1 = (s1 + 1L + e1) / 2, c2 = (s2 + 1L + e2) / 2)]
LP[, `:=`(span = c2 - c1, both_in = s1 + 1L >= WIN[1] & e2 <= WIN[2],
          prom = ov(s1 + 1L, e1, PROM[1], PROM[2]) | ov(s2 + 1L, e2, PROM[1], PROM[2]))]
stopifnot(LP[group == "Con", .N] == 11, LP[group == "Mhp", .N] == 20, LP[group == "Con", sum(prom)] == 2, LP[group == "Mhp", sum(prom)] == 4)
LP[, id := paste0(group, "_", .I)]
# E2 21 kb 内无锚点（§4.2 红线自检）
stopifnot(LP[ov(s1 + 1L, e1, E2_21[1], E2_21[2]) | ov(s2 + 1L, e2, E2_21[1], E2_21[2]), .N] == 0)
cat("④ loop：窗口内 Con", LP[group == "Con", .N], "（两端都在窗口内", LP[group == "Con", sum(both_in)], "；连启动子加粗", LP[group == "Con", sum(prom)], "）；Mhp",
    LP[group == "Mhp", .N], "（", LP[group == "Mhp", sum(both_in)], "；", LP[group == "Mhp", sum(prom)], "）；E2 21 kb 内锚点 0\n")
key("loops_drawn_Con", LP[group == "Con", .N], sprintf("两端都在窗口内 %d；一端出窗（弧在边缘截断）%d", LP[group == "Con", sum(both_in)], LP[group == "Con", sum(!both_in)]))
key("loops_drawn_Mhp", LP[group == "Mhp", .N], sprintf("两端都在窗口内 %d；一端出窗 %d", LP[group == "Mhp", sum(both_in)], LP[group == "Mhp", sum(!both_in)]))
key("loops_promoter_bold_Con", LP[group == "Con", sum(prom)], paste(LP[group == "Con" & prom, sprintf("%s-%s<->%s-%s", fmt(s1), fmt(e1), fmt(s2), fmt(e2))], collapse = "; "))
key("loops_promoter_bold_Mhp", LP[group == "Mhp", sum(prom)], paste(LP[group == "Mhp" & prom, sprintf("%s-%s<->%s-%s", fmt(s1), fmt(e1), fmt(s2), fmt(e2))], collapse = "; "))
key("loop_anchors_in_E2_21kb", 0, "两组均无 → 图注 / Results 只能写整段 SE1 与启动子成环（§4.2）")
save_csv(LP[, .(group, chr1, s1, e1, chr2, s2, e2, span_bp = span, both_in_window = both_in, promoter_bold = prom)], "02_b_loop弧.csv")
H_HALF <- TR[id == "loop", h] / 2; Y_LOOP <- tr("loop", "mid")
arc_pts <- function(c1, c2, h, ybase, sgn, n = 90) {                 # 半椭圆弧；出窗的点丢掉 → 弧在窗口边缘截断
  t <- seq(0, pi, length.out = n); r <- (c2 - c1) / 2
  data.table(x = (c1 + c2) / 2 - r * cos(t), y = ybase + sgn * h * sin(t))[x >= WIN[1] & x <= WIN[2]] }
LP[, h := (H_HALF - 0.4) * sqrt(span / max(span))]
LP[, bold := if (LOOP_BOLD == "promoter") prom else FALSE]                      # LOOP_BOLD = "none" 时全部同粗、同透明度
loop_arc <- LP[, arc_pts(c1, c2, h, Y_LOOP, if (group == "Con") 1 else -1), by = .(id, group, bold)]
loop_arc[, `:=`(col = grp_col[group], lwd = fifelse(bold, lw(1.1), lw(0.4)), alpha = fifelse(bold, 1, 0.75))]
key("loop_bold_rule", LOOP_BOLD, if (LOOP_BOLD == "promoter") sprintf("bold = either anchor overlaps chr13:%s-%s（每组一个合并文库，不做逐 loop 差异判定）", fmt(PROM[1]), fmt(PROM[2])) else "全部弧线同粗，高度只按跨度")
key("loop_E4_promoter_Mhp_only_bold", "61,385,000-61,390,000<->61,525,000-61,530,000",
    "SE1 内唯一只有 Mhp 有的粗弧；Hi-C 像素 O/E ×0.97（每组一个合并文库，仅作描述）")

# ============================================================================ ⑤ Hi-C 5 kb 镜像三角热图（cool 原始计数 × 00 探针 ICE 权重）
wl <- (WIN[1] - 1L) %/% BIN; wh <- WIN[2] %/% BIN - 1L                # 窗口 bin（chr13 内 0-based）：12220–12339，共 120 个
read_cool_window <- function(f, chrom, wl, wh) {
  nm <- as.character(h5read(f, "chroms/name")); co <- as.numeric(h5read(f, "indexes/chrom_offset", bit64conversion = "double"))
  lo <- co[match(chrom, nm)]
  b1o <- as.numeric(h5read(f, "indexes/bin1_offset", start = lo + wl + 1, count = wh - wl + 2, bit64conversion = "double"))
  plo <- b1o[1]; n <- b1o[length(b1o)] - plo
  bin2 <- as.numeric(h5read(f, "pixels/bin2_id", start = plo + 1, count = n, bit64conversion = "double"))
  cnt  <- as.numeric(h5read(f, "pixels/count",  start = plo + 1, count = n))
  bin1 <- rep.int(seq.int(lo + wl, lo + wh), diff(b1o))
  k <- bin2 <= lo + wh
  data.table(i = as.integer(bin1[k] - lo), j = as.integer(bin2[k] - lo), raw = cnt[k]) }
W <- fread(W00); stopifnot(nrow(W) == 41667, all(W$start == (seq_len(nrow(W)) - 1L) * BIN))
PX <- CJ(i = wl:wh, j = wl:wh)[i <= j]
for (g in grp_levels) {
  t0 <- Sys.time(); p <- read_cool_window(COOL[[g]], CHR, wl, wh)
  PX <- p[PX, on = c("i", "j")][is.na(raw), raw := 0]
  setnames(PX, "raw", paste0("raw_", g))
  wv <- W[[paste0("w_", g)]]; stopifnot(!anyNA(wv[(wl:wh) + 1L]))
  PX[, (paste0("bal_", g)) := get(paste0("raw_", g)) * wv[i + 1L] * wv[j + 1L]]
  cat(sprintf("⑤ Hi-C %s：读窗口 %.0f 秒；像素 %s，原始总和 %s\n", g, as.numeric(difftime(Sys.time(), t0, units = "secs")), fmt(nrow(PX)), fmt(sum(PX[[paste0("raw_", g)]]))))
}
# 与 00 探针窗口像素表逐像素核对
PX0 <- fread(PX00)[, `:=`(i = start1 %/% BIN, j = start2 %/% BIN)]
chk <- merge(PX, PX0[, .(i, j, raw_Con0 = raw_Con, raw_Mhp0 = raw_Mhp, bal_Con0 = bal_Con, bal_Mhp0 = bal_Mhp)], by = c("i", "j"))
stopifnot(nrow(chk) == 7260, all(chk$raw_Con == chk$raw_Con0), all(chk$raw_Mhp == chk$raw_Mhp0),
          isTRUE(all.equal(chk$bal_Con, chk$bal_Con0, tolerance = 1e-8)), isTRUE(all.equal(chk$bal_Mhp, chk$bal_Mhp0, tolerance = 1e-8)))
stopifnot(sum(PX$raw_Con) == 308115, sum(PX$raw_Mhp) == 282145)
PX[, `:=`(dist = (j - i) * BIN, start1 = i * BIN + 1L, start2 = j * BIN + 1L)]
PX[, `:=`(v_Con = log1p(bal_Con * 1e3), v_Mhp = log1p(bal_Mhp * 1e3))]
# 共用上限 = 两组距离 ≥ 2 bin（与 ICE ignore_diags = 2 同口径）像素 log1p 值的 99.5 百分位数中的较大者；更高的像素（主要是主对角线及其邻线）饱和
VMAX <- max(PX[dist >= 2L * BIN, quantile(v_Con, 0.995)], PX[dist >= 2L * BIN, quantile(v_Mhp, 0.995)])
cat(sprintf("   颜色 = log1p(ICE 校正值 × 1e3)，共用上限 %.3f（距离 ≥ 10 kb 像素 99.5 百分位；两组去主对角线最大 %.3f / %.3f，饱和）；显示距离 ≤ %s kb\n",
            VMAX, PX[dist >= BIN, max(v_Con)], PX[dist >= BIN, max(v_Mhp)], fmt(HIC_DEPTH / 1e3)))
cat(sprintf("   参考像素：SE1 锚点 61,275 kb ↔ 启动子 bin 61,525 kb  Con %.2f / Mhp %.2f；SE1–SE2 共有 loop 61,270 ↔ 61,575 kb  Con %.2f / Mhp %.2f；窗口中位 %.2f\n",
            PX[i == 61275000 %/% BIN & j == 61525000 %/% BIN, v_Con], PX[i == 61275000 %/% BIN & j == 61525000 %/% BIN, v_Mhp],
            PX[i == 61270000 %/% BIN & j == 61575000 %/% BIN, v_Con], PX[i == 61270000 %/% BIN & j == 61575000 %/% BIN, v_Mhp], PX[dist >= 2L * BIN, median(c(v_Con, v_Mhp))]))
key("hic_resolution_bp", BIN); key("hic_normalisation", "ICE (cis-only, cooler-equivalent, 00 探针权重)", "两组 ICE 后行和≈1，深度差已吸收；同色标")
key("hic_colour", sprintf("log1p(balanced x 1e3), shared max %.3f", VMAX), "上限 = 两组距离 ≥ 10 kb 像素 99.5 百分位数的较大者；主对角线及邻线饱和")
key("hic_pixel_log1p_SE1anchor_vs_promoter", sprintf("Con %.2f / Mhp %.2f (window median %.2f)", PX[i == 61275000 %/% BIN & j == 61525000 %/% BIN, v_Con],
    PX[i == 61275000 %/% BIN & j == 61525000 %/% BIN, v_Mhp], PX[dist >= 2L * BIN, median(c(v_Con, v_Mhp))]), "61,275–61,280 kb ↔ 61,525–61,530 kb 像素的颜色值")
key("hic_depth_shown_bp", HIC_DEPTH)
key("hic_raw_sum_window_Con", fmt(sum(PX$raw_Con))); key("hic_raw_sum_window_Mhp", fmt(sum(PX$raw_Mhp)))
# 探针报告中 SE1 锚点 61,275 kb ↔ 启动子 bin 的像素（图注可引）
px_se <- PX[i == 61275000 %/% BIN & j == TSS %/% BIN * 1L]
key("hic_pixel_SE1anchor_61275kb_vs_promoter_bin", sprintf("raw Con %d / Mhp %d; balanced Con %.4f / Mhp %.4f", px_se$raw_Con, px_se$raw_Mhp, px_se$bal_Con, px_se$bal_Mhp))
save_csv(PX[, .(chr = CHR, start1, end1 = start1 + BIN - 1L, start2, end2 = start2 + BIN - 1L, dist_bp = dist, raw_Con, raw_Mhp, bal_Con, bal_Mhp, v_Con, v_Mhp)],
         "02_b_HiC窗口像素_校正值.csv")
# 菱形多边形：u = 基因组坐标，v = 距离/2（bp）→ 版面 mm；Con 向上、Mhp 向下
H_TRI <- TR[id == "hic", h] / 2; Y_HIC <- tr("hic", "mid"); VSC <- H_TRI / (HIC_DEPTH / 2)
HCOL <- colorRampPalette(c("#FFFFFF", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"))(256)
val2col <- function(v) HCOL[pmin(255L, as.integer(floor(pmin(v, VMAX) / VMAX * 255))) + 1L]
mk_tri <- function(g, sgn) {
  d <- PX[dist <= HIC_DEPTH, .(id = paste0(g, "_", i, "_", j), u0 = (i + j) / 2 * BIN, v0 = (j - i) / 2 * BIN, v = get(paste0("v_", g)))]
  d[, col := val2col(v)]
  corners <- data.table(du = c(0, BIN / 2, BIN, BIN / 2), dv = c(0, -BIN / 2, 0, BIN / 2), k = 1:4)
  out <- d[, .(k = corners$k, x = u0 + corners$du, vv = v0 + corners$dv), by = .(id, col)]
  out[, vv := pmin(pmax(vv, 0), HIC_DEPTH / 2)]                          # 对角线菱形只画基线以上一半；顶端不超出显示深度
  out[, y := Y_HIC + sgn * vv * VSC][, .(id, col, k, x, y)] }
hic_poly <- rbind(mk_tri("Con", 1), mk_tri("Mhp", -1))
cat("   菱形多边形", fmt(uniqueN(hic_poly$id)), "个\n")
CB_Y0 <- tr("hic", "mid") - 4.5; CB_H <- 10; CB_TOP <- CB_Y0 + CB_H
CB <- data.table(k = 0:63)[, .(k, ymin = CB_Y0 + k * (CB_H / 64), ymax = CB_Y0 + (k + 1) * (CB_H / 64), col = HCOL[k * 4 + 1L])]
CB_X <- WIN[2] + c(0.8, 2.3) * BP_MM                                 # 色条 1.5 mm 宽
cb_tick <- data.table(v = seq(0, floor(VMAX), by = 1))[, `:=`(y = CB_Y0 + v / VMAX * CB_H, x = CB_X[2] + 0.35 * BP_MM)]
CB_TICK_R <- (CB_X[2] - WIN[2]) / BP_MM + 0.35 + max(strw(as.character(cb_tick$v), 5))   # 刻度字右缘（mm，自 x 轴右端起）
CB_TITLE <- "Contacts (log)"
X_CBTIT <- WIN[2] + CB_TICK_R * BP_MM                                # 标题右缘 = 刻度列右缘（hjust = 1）
Y_CBTIT <- CB_TOP + 0.6 - ink(CB_TITLE, 5)[["bot"]]                  # 标题基线（vjust = 0）：括号下伸部与色条顶隔 0.6 mm
CB_TIT_L <- X_CBTIT - strw(CB_TITLE, 5) * BP_MM                      # 标题左缘（bp）
HIC_TOP_AT <- function(x) Y_HIC + pmin(WIN[2] - x, HIC_DEPTH / 2) * VSC   # Con 半区热图在 x 处的上缘（mm；右端三角渐收）
CB_INK <- c(tick_top = max(cb_tick$y) - capM(5) / 2 + ink(paste(cb_tick$v, collapse = ""), 5)[["top"]],   # 最高刻度字顶
            tit_bot = Y_CBTIT + ink(CB_TITLE, 5)[["bot"]], tit_top = Y_CBTIT + ink(CB_TITLE, 5)[["top"]])
CB_CLEAR <- CB_INK[["tit_bot"]] - HIC_TOP_AT(CB_TIT_L)               # 标题左缘处：标题墨迹下沿 − 热图上缘
stopifnot(CB_INK[["tit_bot"]] - CB_TOP >= 0.4, CB_INK[["tit_bot"]] - CB_INK[["tick_top"]] >= 0.4, CB_INK[["tit_top"]] <= tr("hic", "top"),
          CB_Y0 - capM(5) / 2 - 0.1 >= tr("hic", "bottom"), CB_CLEAR >= 1,                       # 标题 / 刻度 / 色条互不相压、不越出 Hi-C 轨、不压热图
          CB_TICK_R <= RIGHT_MM - 0.3)                                                            # 色条 + 刻度放得下 5 mm 右边距（右缘留 ≥ 0.3 mm）
cat(sprintf("   色标：色条 x 轴右端外 0.8–2.3 mm；刻度字右缘 %.2f mm（右边距 %g）；标题 \"%s\" 宽 %.2f mm，左缘伸入数据区 %.2f mm，与热图上缘间距 %.2f mm；标题顶距 Hi-C 轨顶 %.2f mm\n",
            CB_TICK_R, RIGHT_MM, CB_TITLE, strw(CB_TITLE, 5), (WIN[2] - CB_TIT_L) / BP_MM, CB_CLEAR, tr("hic", "top") - CB_INK[["tit_top"]]))
key("layout_HiC_colourbar", sprintf("title '%s' right-aligned above bar; right margin %g mm", CB_TITLE, RIGHT_MM),
    sprintf("色标标题与单位写在图注（ICE 校正、5 kb bin、颜色 = log1p(校正值 × 10³)）；刻度字右缘 %.2f mm ≤ 右边距 %g mm；标题左缘伸入 Hi-C 右上空白三角 %.2f mm，与热图上缘间距 %.2f mm",
            CB_TICK_R, RIGHT_MM, (WIN[2] - CB_TIT_L) / BP_MM, CB_CLEAR))

# ============================================================================ ⑥ P2G 弧（BHLHE40 全部 49 条）
P2G <- fread(P2GF)[gene == "BHLHE40" & peak_chr == CHR_ENS]
stopifnot(nrow(P2G) == 49, all(P2G$gene_tss == TSS))
P2G[, `:=`(pc = (peak_start + peak_end) / 2, in_win = peak_end <= WIN[2] & peak_start >= WIN[1],
           in_E2 = ov(peak_start, peak_end, E2_21[1], E2_21[2]))]
stopifnot(sum(P2G$in_win) == 40, sum(P2G$in_E2) == 8)
P2G[, `:=`(id = paste0("p2g_", .I), span = abs(pc - TSS))]
cat("⑥ P2G：BHLHE40", nrow(P2G), "条；peak 在窗口内", sum(P2G$in_win), "；E2 21 kb 内", sum(P2G$in_E2), "（r ",
    paste(sprintf("%.2f", P2G[in_E2 == TRUE, cor]), collapse = "/"), "）；窗口外 9 条 peak 位于 61.70–61.77 Mb，弧在右缘截断\n")
key("p2g_links_total", nrow(P2G), "peak 61,291,843–61,770,311，r 0.45–0.75"); key("p2g_links_peak_in_window", sum(P2G$in_win), "其余 9 条弧在右缘截断")
key("p2g_links_in_E2_21kb", sum(P2G$in_E2), paste(P2G[in_E2 == TRUE, sprintf("%s(r=%.2f)", fmt(peak_start), cor)], collapse = "; "))
save_csv(P2G[, .(peak_chr, peak_start, peak_end, gene, gene_tss, cor, FDR, peak_in_window = in_win, peak_in_E2_21kb = in_E2)], "02_b_P2G弧.csv")
H_P2G <- TR[id == "p2g", h]; Y_P2G <- tr("p2g", "bottom") + 0.6
P2G[, h := (H_P2G - 1.0) * sqrt(span / max(span))]
p2g_arc <- P2G[, arc_pts(pmin(pc, TSS), pmax(pc, TSS), h, Y_P2G, 1), by = .(id, in_E2, cor)]
p2g_arc[, `:=`(col = fifelse(in_E2, BH_COL, GREY), lwd = fifelse(in_E2, lw(0.6), lw(0.35)), alpha = fifelse(in_E2, 0.95, 0.55))]
p2g_tick <- P2G[in_win == TRUE, .(x = pc, y = Y_P2G, col = fifelse(in_E2, BH_COL, GREY))]   # peak 位置小刻度只画窗口内的 40 条

# ============================================================================ ⑦ snATAC 8 类伪 bulk 轨迹 / RNA 轨迹：文件到位即画，否则占位 + 需服务器文件清单
NEED <- data.table(
  用途 = c(rep("b 面板 snATAC 伪 bulk 轨迹（必需）", 6), rep("b 面板 RNA 轨迹", 6)),
  服务器路径 = c(sprintf(file.path(MHP_FRAG, "%s.fragments.tsv.gz"), c("Mock_1", "Mock_2", "Mock_3", "Mhp_1", "Mhp_2", "Mhp_3")),
                 sprintf(file.path(MHP_RAW, "04.RNA/%s%d_RNA.bw"), rep(grp_levels, each = 3), 1:3)),
  处理 = c(rep("服务器上 tabix <文件> 13:61000000-61800000 > fig5_frag_<样本>.tsv；文件名保留样本名 Mock_1 … Mhp_3；列 = chr start end barcode count", 6),
           rep("直接下载", 6)),
  本地放置目录 = c(rep(FRAGD, 6), rep(RNAD, 6)),
  本地文件 = c(sprintf("fig5_frag_%s.tsv", c("Mock_1", "Mock_2", "Mock_3", "Mhp_1", "Mhp_2", "Mhp_3")), sprintf("%s%d_RNA.bw", rep(grp_levels, each = 3), 1:3)))
NEED[, 状态 := fifelse(file.exists(file.path(本地放置目录, 本地文件)), "已到位", "缺")]
save_csv(NEED, "02_b_需服务器文件.csv")
cat("⑦ snATAC fragments 子集本地存在：", HAS_FRAG, "；RNA bigwig 本地存在：", HAS_RNA, if (HAS_FRAG) "→ 画 8 类轨迹" else "→ 本版 snATAC 画占位行", if (HAS_RNA) "、画 RNA 轨" else "、不画 RNA 轨", "\n")
if (HAS_FRAG) {
  # ---- 核元数据：9,610 核的样本 / 组 / 类型 / nFrags（与面板 e 的 tSNE 表同一套核）
  meta <- fread(META)[, .(cell, Sample, group = as.character(to_grp(group)), celltype, nFrags, ReadsInTSS)]
  stopifnot(nrow(meta) == 9610, all(meta$celltype %in% ct_atac_order), all(meta$Sample %in% paste0(rep(c("Mock", "Mhp"), each = 3), "_", 1:3)), !anyNA(meta$nFrags), !anyNA(meta$ReadsInTSS))
  meta[, barcode := sub("^.*#", "", cell)]
  stopifnot(meta[group == "Mhp" & celltype == "B_cells", .N] == 3)                          # B 细胞 Mhp 仅 3 核（B_CELL_MODE 的依据）
  frag_files <- list.files(FRAGD, pattern = "\\.tsv(\\.gz)?$", full.names = TRUE)
  smp_of <- function(f) { m <- regmatches(basename(f), regexpr("(Mock|Mhp)_[123]", basename(f))); if (length(m) != 1) NA_character_ else m }
  fsmp <- vapply(frag_files, smp_of, character(1))
  stopifnot(!anyNA(fsmp), setequal(fsmp, unique(meta$Sample)), !anyDuplicated(fsmp))     # 6 个样本各一个文件，文件名含样本名
  n_abin <- (WIN[2] - WIN[1] + 1L) %/% ATAC_BIN; stopifnot(n_abin * ATAC_BIN == WIN[2] - WIN[1] + 1L)
  frag_stat <- list()
  ins_all <- rbindlist(lapply(seq_along(frag_files), function(k) {
    fr <- fread(frag_files[k], header = FALSE, sep = "\t")[, 1:4]; setnames(fr, c("chr", "start", "end", "barcode"))   # 第 5 列（读对支持数）不用：每行已是一个去重片段
    stopifnot(nrow(fr) > 0, all(sub("^chr", "", fr$chr) == CHR_ENS), all(fr$end > fr$start), !anyDuplicated(fr))
    mk <- meta[Sample == fsmp[k]]
    m <- merge(fr, mk[, .(barcode, group, celltype)], by = "barcode")   # 只保留通过 QC 且已注释的核；其余条码（非细胞 / 未过 QC）丢弃
    frag_stat[[fsmp[k]]] <<- data.table(sample = fsmp[k], group = as.character(to_grp(sub("_.*", "", fsmp[k]))), n_fragments_in_file = nrow(fr), n_fragments_in_QC_nuclei = nrow(m),
                                       pct_fragments_in_QC_nuclei = round(100 * nrow(m) / nrow(fr), 1), n_barcodes_in_file = uniqueN(fr$barcode),
                                       n_QC_nuclei = nrow(mk), n_QC_nuclei_hit = uniqueN(m$barcode), pct_QC_nuclei_hit = round(100 * uniqueN(m$barcode) / nrow(mk), 1))
    # Tn5 插入位点 = 片段两端（10x fragments 起点 0-based → 1-based 插入位点 start + 1；终点列即末位碱基的 1-based 坐标）
    ins <- rbind(m[, .(group, celltype, pos = start + 1L)], m[, .(group, celltype, pos = end)])[pos >= WIN[1] & pos <= WIN[2]]
    ins[, bin := (pos - WIN[1]) %/% ATAC_BIN + 1L]
    ins[, .(n_ins = .N), by = .(group, celltype, bin)] }))
  frag_stat <- rbindlist(frag_stat)
  stopifnot(all(frag_stat$pct_QC_nuclei_hit >= 90))                                      # 条码匹配自检：每个样本 ≥ 90% 的 QC 核在位点子集里有片段
  cat("   条码匹配：", paste(frag_stat[, sprintf("%s QC 核命中 %d/%d (%.1f%%)，QC 核内片段 %.1f%%", sample, n_QC_nuclei_hit, n_QC_nuclei, pct_QC_nuclei_hit, pct_fragments_in_QC_nuclei)], collapse = "；"), "\n")
  ins_all <- ins_all[, .(n_ins = sum(n_ins)), by = .(group, celltype, bin)]
  norm <- meta[, .(n_nuclei = .N, nFrags_total = sum(nFrags), ReadsInTSS_total = sum(ReadsInTSS)), by = .(group, celltype)]
  norm[, denom := if (ATAC_NORM == "nFrags") 2 * nFrags_total else ReadsInTSS_total]      # 分母：全基因组 Tn5 插入位点数（= 片段数 × 2）或 TSS 内读数
  ATAC <- merge(CJ(group = grp_levels, celltype = ct_atac_order, bin = seq_len(n_abin)), ins_all, by = c("group", "celltype", "bin"), all.x = TRUE)
  ATAC[is.na(n_ins), n_ins := 0L]
  ATAC <- merge(ATAC, norm, by = c("group", "celltype"))
  ATAC[, cpm := n_ins / denom * 1e6]                              # 每百万（分母按 ATAC_NORM；与面板 f 的 CPM 同为按组 × 类型归一）
  ATAC[, `:=`(start = WIN[1] + (bin - 1L) * ATAC_BIN, end = WIN[1] + bin * ATAC_BIN - 1L, tid = paste0("atac_", celltype), drawn = celltype %in% ct_atac_draw)]
  ATAC_MAXBIN <- ATAC[drawn == TRUE][which.max(cpm)]
  ATAC_YMAX_RAW <- ceiling(ATAC_MAXBIN$cpm * 10) / 10
  ATAC_YMAX <- if (is.na(ATAC_YCAP)) ATAC_YMAX_RAW else ATAC_YCAP
  n_capped <- ATAC[drawn == TRUE, sum(cpm > ATAC_YMAX)]                                          # 截断时顶到行高的 bin 数（不截断为 0）
  prom_bins <- ATAC[drawn == TRUE & start >= PROM[1] & end <= PROM[2]][order(-cpm)][1:2]        # 启动子口径内最高的两个 bin（供核对"由谁决定 y 上限"）
  maxbin_region <- if (ATAC_MAXBIN$start >= E2_21[1] & ATAC_MAXBIN$end <= E2_21[2]) "E2" else if (ATAC_MAXBIN$start >= PROM[1] & ATAC_MAXBIN$end <= PROM[2]) "promoter" else
    { hit <- EL[ov(start, end, ATAC_MAXBIN$start + 0, ATAC_MAXBIN$end + 0), name]; if (length(hit)) hit else "other" }
  stopifnot(ATAC_MAXBIN$cpm >= max(prom_bins$cpm))
  # 每类型行内两组的落位：stacked = Con 上 / Mhp 下两条子行（各自基线，同向）；mirror = 共用中线，Con 向上 / Mhp 向下；overlay = 同一基线叠画
  ATAC[, `:=`(row_y0 = TR$bottom[match(tid, TR$id)], row_h = TR$h[match(tid, TR$id)])]
  if (ATAC_MODE == "stacked") {
    sub_h <- (ATAC_H - ATAC_SUBGAP) / 2
    ATAC[, `:=`(y0 = fifelse(group == "Con", row_y0 + sub_h + ATAC_SUBGAP, row_y0), hh = sub_h, dir = 1)]
  } else if (ATAC_MODE == "mirror") {
    ATAC[, `:=`(y0 = row_y0 + row_h / 2, hh = row_h / 2, dir = fifelse(group == "Con", 1, -1))]
  } else ATAC[, `:=`(y0 = row_y0, hh = row_h, dir = 1)]
  atac_rect <- ATAC[drawn == TRUE & cpm > 0][order(group)][, .(xmin = start, xmax = end + 1L, ymin = pmin(y0, y0 + dir * hh * pmin(cpm / ATAC_YMAX, 1)),
                                                                ymax = pmax(y0, y0 + dir * hh * pmin(cpm / ATAC_YMAX, 1)),
                                                                col = if (B_CELL_MODE == "grey") fifelse(celltype == "B_cells", "grey60", grp_col[group]) else grp_col[group], group)]
  atac_base <- unique(ATAC[drawn == TRUE, .(y0, group)])[, `:=`(x = WIN[1], xend = WIN[2])]     # 每条子行（或每行）的基线；mirror 时两组共用一条中线
  atac_lab  <- data.table(x = X_LAB, y = TR[match(paste0("atac_", ct_atac_draw), id), (top + bottom) / 2], label = lab_atac(ct_atac_draw),
                          col = if (B_CELL_MODE == "grey") fifelse(ct_atac_draw == "B_cells", "grey50", "black") else "black")
  stopifnot(ct_atac_draw[1] == "Macrophages")
  ATAC_RNG <- ATAC[drawn == TRUE & celltype == ct_atac_draw[1] & group == "Con", .(top = y0[1] + hh[1], bottom = y0[1])]   # 巨噬 Con 子行上下缘（stacked：Con 在上）
  region_of <- function(s, e) vapply(seq_along(s), function(k) {
    if (s[k] >= E2_21[1] && e[k] <= E2_21[2]) return("E2")
    if (s[k] >= PROM[1] && e[k] <= PROM[2]) return("promoter")
    hit <- EL[ov(start, end, s[k], e[k]), name]; if (length(hit)) paste(hit, collapse = "/") else "other" }, character(1))
  ATAC_CAPPED <- ATAC[drawn == TRUE & cpm > ATAC_YMAX][order(group, match(celltype, ct_atac_order), bin)]
  stopifnot(nrow(ATAC_CAPPED) == n_capped, !is.na(ATAC_YCAP) || n_capped == 0)
  CAP_SEG <- data.table()
  if (n_capped > 0) {
    ATAC_CAPPED[, region := region_of(start, end)]
    ATAC_CAPPED[, run := cumsum(c(TRUE, diff(bin) != 1L)), by = .(group, celltype)]
    CAP_SEG <- ATAC_CAPPED[, .(start = min(start), end = max(end), n_bins = .N, max_cpm = max(cpm), region = paste(unique(region), collapse = "/")), by = .(group, celltype, run)]
    cat(sprintf("   截断：y 上限 %g CPM，超限 bin %d 个（相邻合并为 %d 段）：%s\n", ATAC_YMAX, n_capped, nrow(CAP_SEG),
                paste(CAP_SEG[, sprintf("%s %s %s–%s（%s，%d bin，最高 %.2f）", group, lab_atac(celltype), fmt(start), fmt(end), region, n_bins, max_cpm)], collapse = "；")))
  }
  key("snATAC_capped_bins", n_capped, if (n_capped > 0) paste(ATAC_CAPPED[, sprintf("%s %s chr13:%s-%s %.2f CPM (%s)", group, lab_atac(celltype), fmt(start), fmt(end), cpm, region)], collapse = "; ")
      else "不截断（ATAC_YCAP = NA）")
  key("snATAC_cap_markers", 0L, if (n_capped > 0) paste0(sprintf("截断只在图注写明；y 上限 %g CPM（7 类共用）；超限 bin 画到上限，相邻合并后 %d 段：", ATAC_YMAX, nrow(CAP_SEG)),
      paste(CAP_SEG[, sprintf("%s %s %s-%s (%s, %d bin%s)", group, lab_atac(celltype), fmt(start), fmt(end), region, n_bins, fifelse(n_bins > 1L, "s", ""))], collapse = "; ")) else "")
  # ---- E2 方向 / 幅度核对（与 05_e 的 P2G8 CPM 比）：① E2 21 kb 内 250 bp bin 均值；② 8 个 P2G peak 内插入位点合计 / 分母 × 1e6
  e2_atac <- ATAC[start >= E2_21[1] & end <= E2_21[2], .(E2_21kb_mean_CPM = mean(cpm)), by = .(group, celltype)]
  pk8 <- fread(P2G8F); stopifnot(nrow(pk8) == 8, all(pk8$gene == "BHLHE40"), all(pk8$peak_start >= E2_21[1] & pk8$peak_end <= E2_21[2]))
  in_pk8 <- ATAC[, sapply(seq_len(.N), function(r) any(ov(start[r], end[r], pk8$peak_start, pk8$peak_end)))]    # 与任一 peak 重叠的 250 bp bin（bin 与 501 bp peak 边界不齐，取重叠 bin）
  p8_atac <- ATAC[in_pk8, .(P2G8_bins_CPM = sum(n_ins) / denom[1] * 1e6, n_bins = .N), by = .(group, celltype)]
  e05 <- fread(E05F)[, .(celltype, Con_05e = Con, Mhp_05e = Mhp)]; stopifnot(nrow(e05) == 7, !"B_cells" %in% e05$celltype)
  E2CHK <- merge(dcast(e2_atac, celltype ~ group, value.var = "E2_21kb_mean_CPM"), dcast(p8_atac[, .(group, celltype, P2G8_bins_CPM)], celltype ~ group, value.var = "P2G8_bins_CPM"),
                 by = "celltype", suffixes = c("_E2mean", "_P2G8bins"))
  E2CHK <- merge(E2CHK, e05, by = "celltype", all.x = TRUE)
  E2CHK <- E2CHK[match(ct_atac_order, celltype)]
  E2CHK[, `:=`(ratio_E2mean = Mhp_E2mean / Con_E2mean, ratio_P2G8bins = Mhp_P2G8bins / Con_P2G8bins, ratio_05e = Mhp_05e / Con_05e)]
  E2CHK[, `:=`(rank_Con_E2mean = frank(-Con_E2mean), rank_Mhp_E2mean = frank(-Mhp_E2mean), rank_Con_05e = frank(-Con_05e, na.last = "keep"), rank_Mhp_05e = frank(-Mhp_05e, na.last = "keep"))]
  E2CHK[, direction_agrees_05e := fifelse(is.na(ratio_05e), NA, sign(ratio_E2mean - 1) == sign(ratio_05e - 1))]
  cmp <- E2CHK[celltype != "B_cells"]
  rho_E2 <- cor(c(cmp$Con_E2mean, cmp$Mhp_E2mean), c(cmp$Con_05e, cmp$Mhp_05e), method = "spearman")
  rho_P8 <- cor(c(cmp$Con_P2G8bins, cmp$Mhp_P2G8bins), c(cmp$Con_05e, cmp$Mhp_05e), method = "spearman")
  # 必须成立的方向：中性粒在两组都最高（7 类内）且感染后上升；巨噬第二且上升（面板 f：16.04 → 32.85、6.14 → 9.49）
  stopifnot(cmp[which.max(Con_E2mean), celltype] == "Neutrophils", cmp[which.max(Mhp_E2mean), celltype] == "Neutrophils",
            cmp[celltype == "Neutrophils", ratio_E2mean] > 1, cmp[celltype == "Macrophages", ratio_E2mean] > 1,
            cmp[order(-Mhp_E2mean), celltype][2] == "Macrophages", rho_E2 > 0.8)
  cat("   snATAC：", fmt(n_abin), "个", ATAC_BIN, "bp bin；分母 =", ATAC_NORM, "；共用 y 上限", ATAC_YMAX, "CPM", if (is.na(ATAC_YCAP)) "（不截断）" else sprintf("（截断；超限 bin %d 个）", n_capped),
      "；画", length(ct_atac_draw), "类（B 细胞", B_CELL_MODE, "；两组画法", ATAC_MODE, "，行高", ATAC_H, "mm）\n")
  cat(sprintf("   y 上限由最高 bin 决定：%s %s chr13:%s–%s（%s 区）%.2f CPM；启动子口径内最高两 bin %s / %s = %.1f / %.1f\n", ATAC_MAXBIN$group, lab_atac(ATAC_MAXBIN$celltype),
              fmt(ATAC_MAXBIN$start), fmt(ATAC_MAXBIN$end), maxbin_region, ATAC_MAXBIN$cpm, fmt(prom_bins$start[1]), fmt(prom_bins$start[2]), prom_bins$cpm[1], prom_bins$cpm[2]))
  cat("   E2 21 kb 均值 CPM（Con → Mhp，×）：", paste(E2CHK[, sprintf("%s %.2f → %.2f (×%.2f)", lab_atac(celltype), Con_E2mean, Mhp_E2mean, ratio_E2mean)], collapse = "；"), "\n")
  cat(sprintf("   与 05_e P2G8 CPM 的一致性：Spearman ρ（E2 均值）= %.3f，ρ（8 peak bin）= %.3f；方向一致 %d/7；中性粒两组最高、巨噬第二 ✓\n",
              rho_E2, rho_P8, sum(cmp$direction_agrees_05e)))
  save_csv(frag_stat, "02_b_snATAC_fragments读入统计.csv")
  save_csv(dcast(ATAC[, .(celltype, bin, start, end, group, cpm = signif(cpm, 5))], celltype + bin + start + end ~ group, value.var = "cpm"), sprintf("02_b_snATAC轨迹_%dbp组CPM.csv", ATAC_BIN))
  save_csv(E2CHK, "02_b_snATAC_E2核对_vs_05e.csv")
  key("snATAC_tracks", sprintf("drawn (%d cell types, %s; B cells: %s)", length(ct_atac_draw),
                               switch(ATAC_MODE, stacked = "Con upper / Mhp lower sub-track per cell type", mirror = "Con up / Mhp down mirrored per cell type", overlay = "Con / Mhp overlaid per row"), B_CELL_MODE),
      sprintf("Tn5 insertions per 10^6 (denominator = %s of the nuclei in each group x cell type), %d bp bins, shared y max %g; row height %g mm", if (ATAC_NORM == "nFrags") "2 x nFrags" else "ReadsInTSS", ATAC_BIN, ATAC_YMAX, ATAC_H))
  key("snATAC_ymax_determined_by", sprintf("%s %s chr13:%s-%s (%s region) %.2f CPM -> shared y max %g%s", ATAC_MAXBIN$group, lab_atac(ATAC_MAXBIN$celltype), fmt(ATAC_MAXBIN$start), fmt(ATAC_MAXBIN$end),
                                           maxbin_region, ATAC_MAXBIN$cpm, ATAC_YMAX, if (is.na(ATAC_YCAP)) " (no cap)" else sprintf(" (capped; %d bins clipped)", n_capped)),
      sprintf("最高 bin 在 E4 区，不是启动子 bin；启动子口径 %s-%s 内最高两 bin %s / %s = %.1f / %.1f CPM（均 Mhp 中性粒）", fmt(PROM[1]), fmt(PROM[2]),
              fmt(prom_bins$start[1]), fmt(prom_bins$start[2]), prom_bins$cpm[1], prom_bins$cpm[2]))
  key("snATAC_fragments_source", basename(FRAGD), "tabix 13:61000000-61800000 子集；每行 = 去重片段，第 5 列读对数不用；Mock = Con")
  key("snATAC_barcode_match", paste(frag_stat[, sprintf("%s %d/%d QC nuclei (%.1f%%), %.1f%% of file fragments", sample, n_QC_nuclei_hit, n_QC_nuclei, pct_QC_nuclei_hit, pct_fragments_in_QC_nuclei)], collapse = "; "),
      "未命中的 QC 核 = 该核在 800 kb 子集内无片段（预期）；文件其余片段属非细胞 / 未过 QC 条码")
  key("snATAC_nuclei_per_group_celltype", paste(norm[order(match(celltype, ct_atac_order), group), sprintf("%s %s n=%d nFrags=%s", celltype, group, n_nuclei, fmt(nFrags_total))], collapse = "; "))
  key("snATAC_E2_mean_CPM_Con_to_Mhp", paste(E2CHK[, sprintf("%s %.2f->%.2f (x%.2f)", celltype, Con_E2mean, Mhp_E2mean, ratio_E2mean)], collapse = "; "), sprintf("E2 21 kb 内 %d bp bin 均值；B 细胞 Mhp 仅 3 核，勿引用", ATAC_BIN))
  key("snATAC_E2_check_vs_05e", sprintf("Spearman rho E2-mean vs 05e P2G8 CPM = %.3f; 8-peak-bins vs 05e = %.3f; direction agrees %d/7; neutrophils top in both groups, macrophages second", rho_E2, rho_P8, sum(cmp$direction_agrees_05e)),
      "两量不同（本轨迹 = 区间内插入 / 全基因组插入；05_e = ArchR PeakMatrix CPM），只核方向与排序；Results 数字仍引 05_e")
} else {
  key("snATAC_tracks", "PLACEHOLDER (fragments not local)", "见 02_b_需服务器文件.csv；代码已备好，文件放到位即自动补画（Hi-C 半高同时压到 12 mm）")
}
if (HAS_RNA) {
  rna_files <- list.files(RNAD, pattern = "\\.bw$", full.names = TRUE)
  rna_grp <- ifelse(grepl("^Con", basename(rna_files)), "Con", ifelse(grepl("^Mhp", basename(rna_files)), "Mhp", NA_character_))
  stopifnot(!anyNA(rna_grp), all(grp_levels %in% rna_grp))
  bin_bw_any <- function(f) {                                     # RNA bigwig 染色体名可能是 "13" 或 "chr13"：按文件自身的 seqlevels 取窗口
    sl <- seqlevels(BigWigFile(f)); chr_use <- if (CHR %in% sl) CHR else if (CHR_ENS %in% sl) CHR_ENS else stop("RNA bigwig 无 chr13/13: ", f)
    x <- import(f, which = GRanges(chr_use, IRanges(WIN[1], WIN[2]))); x <- restrict(x, WIN[1], WIN[2])
    v <- rep(0, WIN[2] - WIN[1] + 1L)
    for (k in seq_along(x)) v[(start(x)[k] - WIN[1] + 1L):(end(x)[k] - WIN[1] + 1L)] <- x$score[k]
    colMeans(matrix(v, nrow = SIG_BIN)) }
  stopifnot(length(rna_files) == 6, sum(rna_grp == "Con") == 3, sum(rna_grp == "Mhp") == 3)
  RNA <- rbindlist(lapply(seq_along(rna_files), function(k) data.table(sample = basename(rna_files[k]), group = rna_grp[k], bin = seq_len(n_bin), v = bin_bw_any(rna_files[k]))))
  # ---- 归一化核对：每个 bigwig 的全基因组 sum(score × width)（按染色体 summary 均值 × 长度）与窗口总量；六库总量若近乎相等即已按深度归一化
  rna_norm <- rbindlist(lapply(seq_along(rna_files), function(k) {
    bw <- BigWigFile(rna_files[k]); s <- summary(bw, size = 1L, type = "mean")
    gsum <- sum(vapply(seq_along(s), function(i) mean(s[[i]]$score, na.rm = TRUE) * width(s[[i]]), numeric(1)), na.rm = TRUE)
    x <- import(rna_files[k], which = GRanges(if (CHR %in% seqlevels(bw)) CHR else CHR_ENS, IRanges(WIN[1], WIN[2]))); x <- restrict(x, WIN[1], WIN[2])
    data.table(sample = basename(rna_files[k]), group = rna_grp[k], n_seqlevels = length(seqlevels(bw)), genome_bp = sum(as.numeric(seqlengths(bw))),
               min_interval_bp = min(width(x)), genome_sum_score_x_width = gsum, window_sum_score_x_width = sum(x$score * width(x))) }))
  rna_norm[, genome_total_rel_to_mean := genome_sum_score_x_width / mean(genome_sum_score_x_width)]
  RNA_NORM_RANGE <- diff(range(rna_norm$genome_sum_score_x_width)) / mean(rna_norm$genome_sum_score_x_width)
  stopifnot(RNA_NORM_RANGE < 0.05)                                 # 六库全基因组总量极差 < 5% → 已归一化（原始覆盖度会随文库大小变化）
  cat(sprintf("   RNA bigwig：%d 个；染色体名 %s；最小区间 %d bp（= bin）；全基因组 sum(score×width) %.3g–%.3g（极差 %.1f%% → 已按深度归一化，与 deepTools RPKM 50 bp bin 一致）；窗口总量 Con %s / Mhp %s\n",
              nrow(rna_norm), if (CHR %in% seqlevels(BigWigFile(rna_files[1]))) CHR else CHR_ENS, min(rna_norm$min_interval_bp),
              min(rna_norm$genome_sum_score_x_width), max(rna_norm$genome_sum_score_x_width), 100 * RNA_NORM_RANGE,
              paste(fmt(rna_norm[group == "Con", window_sum_score_x_width]), collapse = "/"), paste(fmt(rna_norm[group == "Mhp", window_sum_score_x_width]), collapse = "/")))
  RNAM <- RNA[, .(v = mean(v)), by = .(group, bin)][, `:=`(start = WIN[1] + (bin - 1L) * SIG_BIN, end = WIN[1] + bin * SIG_BIN - 1L)]
  RNA_YMAX <- ceiling(max(RNAM$v))                                  # 两组共用 y 上限（= 组均值最高 bin，Mhp BHLHE40 3′ 外显子）
  bh_body <- RNAM[start >= TSS & end <= 61533674L, .(m = mean(v)), by = group]
  RNA_BH_RATIO <- bh_body[group == "Mhp", m] / bh_body[group == "Con", m]
  stopifnot(RNA_BH_RATIO > 1, RNAM[which.max(v), start] >= TSS & RNAM[which.max(v), start] <= 61533674L)   # BHLHE40 感染后上升（Fig 4 b 方向）；最高 bin 在 BHLHE40 基因体内
  cat(sprintf("   RNA 组均值 y 上限 %g；BHLHE40 基因体均值 Con %.1f / Mhp %.1f（×%.2f）；ARL8B 3′ 端 Con 最高 bin %.0f（%.0f%% 高度）\n", RNA_YMAX, bh_body[group == "Con", m], bh_body[group == "Mhp", m], RNA_BH_RATIO,
              RNAM[group == "Con" & start >= 61640671L & end <= 61698374L, max(v)], 100 * RNAM[group == "Con" & start >= 61640671L & end <= 61698374L, max(v)] / RNA_YMAX))
  RNAM[, tid := paste0("rna_", group)][, `:=`(y0 = TR$bottom[match(tid, TR$id)], hh = TR$h[match(tid, TR$id)])]
  rna_rect <- RNAM[v > 0, .(xmin = start, xmax = end + 1L, ymin = y0, ymax = y0 + hh * pmin(v / RNA_YMAX, 1), col = grp_col[group])]
  rna_lab  <- unique(RNAM[, .(tid, group, y0, hh)])[, .(tid, x = X_LAB, y = y0 + hh / 2, label = paste("RNA", group), label_r = sprintf("[0–%g]", RNA_YMAX))]
  save_csv(dcast(RNAM[, .(bin, start, end, group, v = signif(v, 5))], bin + start + end ~ group, value.var = "v"), "02_b_RNA轨迹_500bp组均值.csv")
  save_csv(rna_norm[, .(sample, group, n_seqlevels, genome_bp, min_interval_bp, genome_sum_score_x_width = signif(genome_sum_score_x_width, 5), genome_total_rel_to_mean = round(genome_total_rel_to_mean, 4),
                        window_sum_score_x_width = signif(window_sum_score_x_width, 5))], "02_b_RNA_bigwig归一化核对.csv")
  key("RNA_track", sprintf("drawn (%s; mean of 3 per group, %d bp bins, shared y max %g)", paste(sub("_RNA.bw", "", basename(rna_files)), collapse = ", "), SIG_BIN, RNA_YMAX),
      "两行放在基因模型正下方；BHLHE40 外显子信号与基因模型对齐（Mhp 3′ 外显子 = 全窗口最高 bin）")
  key("RNA_bigwig_normalisation", sprintf("depth-normalised: genome-wide sum(score x width) %.3g-%.3g across 6 libraries (range %.1f%% of mean); 50-bp bins; chr-prefixed names, %d seqlevels",
      min(rna_norm$genome_sum_score_x_width), max(rna_norm$genome_sum_score_x_width), 100 * RNA_NORM_RANGE, rna_norm$n_seqlevels[1]),
      "总量与 deepTools bamCoverage RPKM（50 bp bin、~150 bp 读长：≈1e9 × 每读跨的 bin 数 ≈ 4）一致；图注写 'depth-normalised coverage (RPKM-like)'")
  key("RNA_BHLHE40_body_mean_Mhp_over_Con", sprintf("%.2f (Con %.1f / Mhp %.1f)", RNA_BH_RATIO, bh_body[group == "Con", m], bh_body[group == "Mhp", m]), "500 bp bin 组均值在 BHLHE40 基因体上的均值比；与 Fig 4 b DESeq2 方向一致")
} else key("RNA_track", "not drawn (bigwig not local)")

# ============================================================================ ⑧ 组图
# 左侧标签（写在边距）
lab_left <- rbind(
  data.table(x = X_LAB, y = LANE_Y[1] + 0.5, label = "Genes", col = "black"),
  data.table(x = X_LAB, y = SE_LANE[["Con"]], label = "SE Con", col = "black"),
  data.table(x = X_LAB, y = SE_LANE[["Mhp"]], label = "SE Mhp", col = "black"),
  data.table(x = X_LAB, y = mean(EY), label = "E1–E4", col = "black"),
  if (HAS_RNA) rna_lab[, .(x, y, label, col = "black")] else NULL,
  sig_lab[, .(x, y, label, col = "black")],
  data.table(x = X_LAB, y = Y_LOOP + H_HALF / 2, label = "Loops Con", col = "black"),
  data.table(x = X_LAB, y = Y_LOOP - H_HALF / 2, label = "Loops Mhp", col = "black"),
  data.table(x = X_LAB, y = Y_HIC + H_TRI / 2, label = "Hi-C Con", col = "black"),
  data.table(x = X_LAB, y = Y_HIC - H_TRI / 2, label = "Hi-C Mhp", col = "black"),
  if (HAS_FRAG) atac_lab[, .(x, y, label, col)] else data.table(x = X_LAB, y = tr("atac", "mid"), label = "snATAC", col = GREY),
  data.table(x = X_LAB, y = tr("p2g", "mid"), label = "P2G links", col = "black"))
rng_lab <- rbind(
  if (HAS_RNA) rna_lab[, .(tid, label = label_r, top = TR$top[match(tid, TR$id)], bottom = TR$bottom[match(tid, TR$id)])] else NULL,
  sig_lab[, .(tid, label = label_r, top = TR$top[match(tid, TR$id)], bottom = TR$bottom[match(tid, TR$id)])],
  if (HAS_FRAG) data.table(tid = "atac_Macrophages_Con", label = sprintf("[0–%g]", ATAC_YMAX), top = ATAC_RNG$top, bottom = ATAC_RNG$bottom) else NULL)
rng_lab[, `:=`(x = WIN[1] + RNG_DX * BP_MM, y = top - 0.15 - vapply(label, function(s) ink(s, 5)[["top"]], numeric(1)))]
rng_lab[, `:=`(ink_bot = y + vapply(label, function(s) ink(s, 5)[["bot"]], numeric(1)), x_end = WIN[1] + (RNG_DX + strw(label, 5)) * BP_MM)]
# 碰撞核对：标签覆盖的 x 区间内，该轨柱顶（y0 + 行高 × 归一化值）须低于标签墨迹下沿 0.15 mm 以上
bar_top <- function(t, x0, x1) {
  d <- if (grepl("^rna_", t)) RNAM[tid == t, .(start, end, top = y0 + hh * pmin(v / RNA_YMAX, 1))]
       else if (grepl("^atac_", t)) ATAC[celltype == "Macrophages" & group == "Con", .(start, end, top = y0 + hh * pmin(cpm / ATAC_YMAX, 1))]
       else SIGM[tid == t, .(start, end, top = y0 + hh * pmin(cpm / ymax, 1))]
  d[start <= x1 & end >= x0, max(top)] }
rng_lab[, bar_max := mapply(bar_top, tid, x - 0.2 * BP_MM, x_end + 0.2 * BP_MM)]
rng_lab[, clear_mm := ink_bot - bar_max]
stopifnot(nrow(rng_lab) == 2 * HAS_RNA + 4 + HAS_FRAG, all(rng_lab$clear_mm >= 0.15), all(rng_lab$ink_bot - rng_lab$bottom >= 0.3))   # 字底在该轨基线之上 ≥ 0.3 mm
cat("⑧ 范围标签（轨内左上）：", paste(rng_lab[, sprintf("%s %s（柱顶距字底 %.2f mm）", tid, label, clear_mm)], collapse = "；"), "\n")
key("range_labels", paste(rng_lab[, sprintf("%s %s", tid, label)], collapse = "; "),
    sprintf("数据范围标在轨内左上灰色方括号（5 pt，%s，左缩进 %.1f mm）；snATAC 只在巨噬 Con 子行标共用范围；标签覆盖区内柱顶距字底 ≥ %.2f mm（stopifnot ≥ 0.15）",
            RNG_COL, RNG_DX, min(rng_lab$clear_mm)))
X_BRK <- c(WIN[1], seq(61200000L, WIN[2], by = 100000L))              # 首刻度 = 窗口起点（61,100,001，标 61,100），否则落在 xlim 外被丢
ATAC_TOP <- tr(ATAC_IDS[1], "top"); ATAC_BOT <- tr(ATAC_IDS[length(ATAC_IDS)], "bottom")   # snATAC 块（占位 1 行或 8 行）的上下缘
p <- ggplot() +
  # 竖向参考：E2 21 kb 淡红带（元件行 → loop 行；snATAC 块；P2G 行）、TSS 灰点线
  annotate("rect", xmin = E2_21[1], xmax = E2_21[2], ymin = tr("loop", "bottom"), ymax = tr("elem", "top"), fill = BH_COL, alpha = 0.08) +
  annotate("rect", xmin = E2_21[1], xmax = E2_21[2], ymin = tr("p2g", "bottom"), ymax = tr("p2g", "top"), fill = BH_COL, alpha = 0.08) +
  annotate("segment", x = TSS, xend = TSS, y = tr("loop", "bottom"), yend = tr("elem", "top"), colour = GREY, linewidth = lw(0.3), linetype = "22") +
  annotate("segment", x = TSS, xend = TSS, y = tr("p2g", "bottom"), yend = tr("p2g", "top"), colour = GREY, linewidth = lw(0.3), linetype = "22") +
  (if (HAS_FRAG) list(annotate("rect", xmin = E2_21[1], xmax = E2_21[2], ymin = ATAC_BOT, ymax = ATAC_TOP, fill = BH_COL, alpha = 0.08),
                      annotate("segment", x = TSS, xend = TSS, y = ATAC_BOT, yend = ATAC_TOP, colour = GREY, linewidth = lw(0.3), linetype = "22")) else NULL) +
  (if (HAS_RNA) list(annotate("rect", xmin = E2_21[1], xmax = E2_21[2], ymin = tr("rna_Mhp", "bottom"), ymax = tr("rna_Con", "top"), fill = BH_COL, alpha = 0.08),   # RNA 两行（基因模型下）也画 E2 带 + TSS 点线
                     annotate("segment", x = TSS, xend = TSS, y = tr("rna_Mhp", "bottom"), yend = tr("rna_Con", "top"), colour = GREY, linewidth = lw(0.3), linetype = "22")) else NULL) +
  # ① 基因
  geom_segment(data = gene_line, aes(x = x, xend = xend, y = y, yend = y, colour = col), linewidth = lw(0.4)) +
  geom_rect(data = exon_box, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = col), colour = NA) +
  geom_text(data = chev, aes(x = x, y = y, label = lab), colour = "black", size = S5, family = FONT, vjust = 0.45) +
  geom_text(data = gene_lab, aes(x = x, y = y, label = label, colour = col, fontface = face), size = S6, family = FONT, vjust = 0) +
  # ② SE 条 + SE1/SE2 名；E1–E4 + 21 kb 框
  geom_rect(data = se, aes(xmin = xmin, xmax = xmax, ymin = y - 0.65, ymax = y + 0.65, fill = col), colour = NA) +
  geom_text(data = se_lab, aes(x = x, y = y, label = label), size = S6, family = FONT, vjust = 0, fontface = "bold") +
  annotate("rect", xmin = E2_21[1], xmax = E2_21[2], ymin = EY[1] - 0.25, ymax = EY[2] + 0.25, fill = NA, colour = BH_COL, linewidth = lw(0.4), linetype = "22") +
  geom_rect(data = EL[, .(start, end, ymin = EY[1], ymax = EY[2], fill)], aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax, fill = fill), colour = NA) +
  geom_text(data = el_lab, aes(x = x, y = y, label = label, colour = col, fontface = face), size = S55, family = FONT, vjust = 0) +
  # ③ 信号轨迹（RNA 两行仅在 bigwig 到位时画）
  (if (HAS_RNA) list(geom_rect(data = rna_rect, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = col), colour = NA),
                     geom_segment(data = unique(RNAM[, .(y0)])[, `:=`(x = WIN[1], xend = WIN[2])], aes(x = x, xend = xend, y = y0, yend = y0), colour = "grey70", linewidth = lw(0.25))) else NULL) +
  geom_rect(data = sig_rect, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = col), colour = NA) +
  geom_segment(data = unique(SIGM[, .(y0)])[, `:=`(x = WIN[1], xend = WIN[2])], aes(x = x, xend = xend, y = y0, yend = y0), colour = "grey70", linewidth = lw(0.25)) +
  annotate("segment", x = WIN[1], xend = WIN[2], y = Y_LOOP, yend = Y_LOOP, colour = "grey70", linewidth = lw(0.25)) +
  geom_path(data = loop_arc[bold == FALSE], aes(x = x, y = y, group = id, colour = col, linewidth = lwd, alpha = alpha)) +
  (if (any(loop_arc$bold)) geom_path(data = loop_arc[bold == TRUE], aes(x = x, y = y, group = id, colour = col, linewidth = lwd, alpha = alpha)) else NULL) +
  # ⑤ Hi-C 镜像三角 + 色条
  geom_polygon(data = hic_poly, aes(x = x, y = y, group = id, fill = col), colour = NA) +
  geom_rect(data = CB[, `:=`(xmin = CB_X[1], xmax = CB_X[2])], aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = col), colour = NA) +
  geom_text(data = cb_tick, aes(x = x, y = y, label = v), size = S5, family = FONT, hjust = 0) +
  annotate("text", x = X_CBTIT, y = Y_CBTIT, label = CB_TITLE, size = S5, family = FONT, hjust = 1, vjust = 0) +
  (if (HAS_FRAG) list(
     geom_rect(data = atac_rect, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = col), colour = NA, alpha = if (ATAC_MODE == "overlay") 0.6 else 1),
     geom_segment(data = atac_base, aes(x = x, xend = xend, y = y0, yend = y0), colour = "grey70", linewidth = lw(0.25)))
   else list(
     annotate("rect", xmin = WIN[1], xmax = WIN[2], ymin = tr("atac", "bottom"), ymax = tr("atac", "top"), fill = "grey96", colour = "grey60", linewidth = lw(0.3), linetype = "22"),
     annotate("text", x = mean(WIN), y = tr("atac", "mid"), size = S6, family = FONT, colour = "grey40",
              label = "snATAC pseudo-bulk tracks (8 cell types × Con / Mhp) — placeholder: fragments subset to be exported from the server (see README)"))) +
  # ⑥ P2G
  annotate("segment", x = WIN[1], xend = WIN[2], y = Y_P2G, yend = Y_P2G, colour = "grey70", linewidth = lw(0.25)) +
  geom_path(data = p2g_arc[in_E2 == FALSE], aes(x = x, y = y, group = id, colour = col, linewidth = lwd, alpha = alpha)) +
  geom_path(data = p2g_arc[in_E2 == TRUE],  aes(x = x, y = y, group = id, colour = col, linewidth = lwd, alpha = alpha)) +
  geom_segment(data = p2g_tick, aes(x = x, xend = x, y = y - 0.5, yend = y, colour = col), linewidth = lw(0.4)) +
  geom_text(data = rng_lab, aes(x = x, y = y, label = label), size = S5, family = FONT, hjust = 0, vjust = 0, colour = RNG_COL) +
  # 左标签（6.5 pt；snATAC 8 类行只有 3 mm 高，用 5.5 pt）
  geom_text(data = lab_left[!label %in% lab_atac(ct_atac_order)], aes(x = x, y = y, label = label, colour = col), size = S65, family = FONT, hjust = 1) +
  (if (HAS_FRAG) geom_text(data = lab_left[label %in% lab_atac(ct_atac_order)], aes(x = x, y = y, label = label, colour = col), size = S55, family = FONT, hjust = 1) else NULL) +
  scale_fill_identity() + scale_colour_identity() + scale_linewidth_identity() + scale_alpha_identity() +
  scale_x_continuous(breaks = X_BRK, labels = fmt(round(X_BRK / 1e3)), name = "chr13 position (kb)") +
  coord_cartesian(xlim = WIN, ylim = c(0, H_ALL), expand = FALSE, clip = "off") +   # 边距里的标签 / 色条不被 scale 删掉
  theme_mhp() +
  theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(),
        axis.title.x = element_text(size = FS_BASE), plot.margin = margin(1.5, RIGHT_MM, 0.5, LEFT_MM, unit = "mm"))
H_MM <- H_ALL + 7.5
cat("⑧ 版面：子轨总高", round(H_ALL, 1), "mm；出图", W_MM, "×", round(H_MM, 1), "mm；整图预算 44 + 4 +", round(H_MM, 1), "+ 4 + 44 =", round(96 + H_MM, 1), "mm（上限 247）\n")
stopifnot(96 + H_MM <= 247)
key("panel_size_mm", sprintf("%g x %.1f", W_MM, H_MM), sprintf("左边距 %g mm 放标签，右边距 %g mm 只放 Hi-C 色条 + 刻度（与面板 b 同；x 轴右端与 b 对齐）；整图预算 44+4+%.1f+4+44 = %.1f mm（上限 247）", LEFT_MM, RIGHT_MM, H_MM, 96 + H_MM))
key("font_pt", "track labels 6.5; snATAC row labels 5.5; x-axis title 7 (FS_BASE); range labels / colour-bar title and ticks 5", "x 轴标题 7 pt，轨道标签 6.5 pt，范围标签 5 pt 灰")
key("图注建议句", sprintf(paste0(
  "The BHLHE40 locus (chr13:%s-%s, Sscrofa11.1). Tracks from top: Ensembl 114 genes; %sROSE super-enhancers of each group (SE1 and SE2 named as in Fig. 4; ",
  "the Mhp SE1 extends ~119 kb beyond the 3' end of ITPR1, so 'within an ITPR1 intron' holds for the Con SE1 and for E1-E4 only); E1-E4 (3,002-bp luciferase fragments) with the %s-bp E2 interval (dashed box; shading in all tracks); ",
  "H3K27ac and H3K4me1 CUT&Tag (mean CPM of three lungs per group, 500-bp bins, shared y-axis per mark); all Hi-C loops with at least one anchor in the window (Con above, Mhp below; arc height scales with span; ",
  "bold arcs, loops with an anchor within chr13:%s-%s (BHLHE40 promoter region); no per-loop differential call is made because each group is one merged library); Hi-C contact maps (ICE-balanced, cis-only, 5-kb bins; colour = log1p(balanced contacts x 10^3), shared scale; Con above / Mhp below); ",
  "%s; peak-to-gene links of BHLHE40 from snATAC-seq (49 links, 8 of them with the peak in E2, in red; 40 with the peak inside the window). ",
  "Bracketed numbers, data range of each track. No loop anchor falls within E2 in either group; the looping statement therefore refers to SE1 as a whole. Dotted line, BHLHE40 TSS (%s)."),
  fmt(WIN[1]), fmt(WIN[2]),
  if (HAS_RNA) "bulk lung RNA-seq, depth-normalised coverage (mean of three lungs per group, 500-bp bins, shared y-axis; the BHLHE40 exons carry the highest signal in the window); " else "",
  fmt(E2_LEN), fmt(PROM[1]), fmt(PROM[2]),
  if (HAS_FRAG) sprintf("snATAC-seq pseudo-bulk accessibility of %s cell types (Tn5 insertions per million %s of all nuclei of that cell type in each group, %d-bp bins; %s; shared scale%s; %s)",
                        c("seven", "eight")[1 + (B_CELL_MODE == "grey")], if (ATAC_NORM == "nFrags") "genome-wide insertions (2 x fragments)" else "reads in TSSs", ATAC_BIN,
                        switch(ATAC_MODE, stacked = "Con (upper, blue) and Mhp (lower, red) shown as separate sub-tracks for each cell type", mirror = "Con (upward, blue) and Mhp (downward, red) mirrored about a shared baseline for each cell type",
                               overlay = "Con and Mhp overlaid"),
                        if (is.na(ATAC_YCAP)) sprintf(" 0-%g CPM", ATAC_YMAX) else sprintf(" capped at %g CPM (the %d bins above %g CPM%s are drawn at the cap)",
                                                                                   ATAC_YMAX, n_capped, ATAC_YMAX,
                                                                                   if (n_capped > 0) sprintf(", all in %s and located in %s%s,", paste(unique(ATAC_CAPPED[, sprintf("%s %s", group, tolower(lab_atac(celltype)))]), collapse = " and "),
                                                                                                             paste(unique(ATAC_CAPPED$region), collapse = " and the "), if ("promoter" %in% ATAC_CAPPED$region) " region" else "") else ""),
                        if (B_CELL_MODE == "omit") "B cells not shown because the Mhp group contains only 3 B-cell nuclei" else "B cells in grey because the Mhp group contains only 3 B-cell nuclei")
  else "snATAC-seq pseudo-bulk tracks (to be added)",
  fmt(TSS)), "图注须说明：方括号 = 各轨数据范围；snATAC 25 CPM 截断与 3 个超限 bin（不再画 ▲）；B 细胞未画；RNA = 深度归一化覆盖度；粗弧口径；P2G 49 条 / E2 内 8 条；Hi-C ICE、5 kb、颜色 = log1p(校正值 × 10³)")
STEM <- "Fig5c_BHLHE40位点多组学"
sav(p, STEM, PAN, W_MM, H_MM, panel_dir = PNL)                                    # PNG 300 dpi 草稿（MHP_FINAL=1 时 SVG + TIFF）
svg_f <- file.path(PAN, paste0(STEM, ".svg"))                                     # 本面板另固定出 SVG（svglite，矢量；quartz PDF 仅 MHP_FINAL=1 时出，见下）
ggsave(svg_f, p, width = W_MM, height = H_MM, units = "mm", device = .svg_dev)
stopifnot(file.exists(svg_f), file.info(svg_f)$size > 0, !any(grepl("<image", readLines(svg_f, warn = FALSE), fixed = TRUE)))
cat("✓", basename(svg_f), fmt(file.info(svg_f)$size), "字节（<image> 0 个）\n")
if (identical(Sys.getenv("MHP_FINAL"), "1")) {
  pdf_f <- file.path(PAN, paste0(STEM, ".pdf")); ggsave(pdf_f, p, width = W_MM, height = H_MM, units = "mm", device = .pdf_dev)
  stopifnot(file.exists(pdf_f), length(grepRaw("/Subtype\\s*/Image", readBin(pdf_f, "raw", file.size(pdf_f)), all = TRUE)) == 0)
  cat("✓", basename(pdf_f), fmt(file.info(pdf_f)$size), "字节（quartz PDF，/Subtype /Image 0 个）\n") }
save_csv(rbindlist(KEY), "02_b_关键数字.csv")
cat("完成\n")
