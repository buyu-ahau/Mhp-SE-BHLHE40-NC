# ============================================================================
# 07_4C_虚拟4C_本地.R — Fig. 5b：以 BHLHE40 启动子为视点的 Hi-C 虚拟 4C 接触谱
# 运行位置：本地；路径变量见 config/paths.R
# ============================================================================

source(file.path(Sys.getenv("MHP_CODE", "."), "config", "paths.R"))

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork); library(systemfonts) })
source(file.path(MHP_CODE, "shared", "_palette_MHP.R"))
BH_COL <- "#B2182B"

OUT <- file.path(MHP_RES, "7、HI-C结果/新Fig5_作图")
PAN <- file.path(OUT, "Plots/面板"); TAB <- file.path(OUT, "Tables"); PNL <- file.path(OUT, "_panels")
for (d in c(PAN, TAB, PNL)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
F_V4C <- file.path(TAB, "00_HiC探针_虚拟4C_BHLHE40启动子.csv"); F_EL <- file.path(TAB, "00_HiC探针_元件坐标.csv")
F_SUM <- file.path(TAB, "00_HiC探针_虚拟4C_分区汇总.csv");    F_NUL <- file.path(TAB, "00_HiC探针_虚拟4C_chr13零分布.csv")
F_LOOP <- file.path(TAB, "00_HiC探针_窗口loop注释.csv")
stopifnot(file.exists(F_V4C), file.exists(F_EL), file.exists(F_SUM), file.exists(F_NUL), file.exists(F_LOOP))
RIGHT_MM <- 5
STEM_MAIN <- "Fig5b_虚拟4C_BHLHE40启动子"; W_MAIN <- 116; H_MAIN <- 44
STEM_SUPP <- "补充_虚拟4C_BHLHE40启动子";  W_SUPP <- 88;  H_SUPP <- 52               # 补充图版本（两行 + E1–E4 条带，设计不变）
MAIN_ROWS  <- "oe"                                                                 # 主图行数："oe" = 只画 O/E 一行（默认，理由见头注）；"both" = 与补充图一样两行
SHOW_RATIO <- FALSE
MAKE_SUPP  <- TRUE
stopifnot(MAIN_ROWS %in% c("oe", "both"))
BIN <- 5000L; WIN <- c(61100000L, 61700000L)                                       # 5 kb bin；窗口与面板 c 一致（0-based bin 起点）
TSS <- 61526116L; VP0 <- (TSS - 1L) %/% BIN * BIN                                  # 视点 bin 起点 61,525,000
D_MIN <- 2L                                                                        # 只画距视点 ≥ 2 bin 的位置

# ---------------------------------------------------------------- ① 读探针输出 + 核对
V <- fread(F_V4C)
stopifnot(nrow(V) == 400, all(c("start", "end", "d", "raw_Con", "bal_Con", "oe_Con", "rpm_Con", "raw_Mhp", "bal_Mhp", "oe_Mhp", "rpm_Mhp",
                                 "in_SE1_Mhp", "in_SE1_Con", "in_SE2_Mhp", "in_E1", "in_E2_3kb", "in_E3", "in_E4", "in_E2_21kb") %in% names(V)),
          all(V$end - V$start == BIN), V[d == 0, start] == VP0, all(V$d == abs(V$start - VP0) %/% BIN))
EL <- fread(F_EL); setkey(EL, element)
stopifnot(all(c("BHLHE40_TSS", "SE1_Mhp", "SE1_Con", "SE2_Mhp", "E1", "E2_3kb", "E3", "E4", "E2_21kb") %in% EL$element), EL["BHLHE40_TSS", start] == TSS,
          EL["SE1_Mhp", start] == 61209361, EL["SE1_Mhp", end] == 61517999, EL["SE2_Mhp", start] == 61531355, EL["E2_21kb", start] == 61329246)
SM <- fread(F_SUM)
W <- V[start >= WIN[1] & start < WIN[2]]
stopifnot(nrow(W) == 120, !anyNA(W$bal_Con), !anyNA(W$bal_Mhp), all(W$exp_Con > 0), all(W$exp_Mhp > 0))
W[d < D_MIN, c("oe_Con", "oe_Mhp", "rpm_Con", "rpm_Mhp") := NA_real_]             # 视点 ± 1 bin 留空
# 重算 SE1 汇总并与探针分区汇总表交叉核对（d ≥ 2 的 SE1 bin）
S1 <- V[in_SE1_Mhp == TRUE & d >= D_MIN]
se1 <- S1[, .(n_bins = .N, ratio_bal = sum(bal_Mhp) / sum(bal_Con), n_ge_bal = sum(bal_Mhp >= bal_Con),
              mean_oe_Con = mean(oe_Con), mean_oe_Mhp = mean(oe_Mhp), n_ge_oe = sum(oe_Mhp >= oe_Con), ratio_rpm = sum(rpm_Mhp) / sum(rpm_Con))]
ref <- SM[grepl("^SE1（Mhp", region)]
stopifnot(nrow(ref) == 1, se1$n_bins == ref$n_bins, abs(se1$ratio_bal - ref$ratio_bal) < 0.002, se1$n_ge_bal == ref$n_Mhp_ge_Con_bal,
          se1$n_ge_oe == ref$n_Mhp_ge_Con_oe, abs(se1$mean_oe_Mhp - ref$mean_oe_Mhp) < 0.002)
nse <- V[start >= WIN[1] & start < WIN[2] & in_SE1_Mhp == FALSE & in_SE2_Mhp == FALSE & d >= D_MIN,
         .(n_bins = .N, ratio_bal = sum(bal_Mhp) / sum(bal_Con), n_ge_bal = sum(bal_Mhp >= bal_Con))]
stopifnot(se1$ratio_bal > 1, se1$n_ge_bal > se1$n_bins / 2, abs(nse$ratio_bal - 0.989) < 0.005)   # 入图条件：SE1 上 Mhp ≥ Con；窗口内非 SE 区无此偏移
# 两组 O/E 峰位（距视点 > 10 kb、窗口内）：都应落在 SE1 内
pk <- W[d >= D_MIN, .(peak_Con = start[which.max(oe_Con)], peak_Mhp = start[which.max(oe_Mhp)],
                      oe_peak_Con = max(oe_Con, na.rm = TRUE), oe_peak_Mhp = max(oe_Mhp, na.rm = TRUE))]
stopifnot(V[start == pk$peak_Con, in_SE1_Mhp], V[start == pk$peak_Mhp, in_SE1_Mhp])
LP <- fread(F_LOOP)
stopifnot(all(c("set", "s1", "e1", "s2", "e2", "Prom_61520_61540", "far_in_SE1") %in% names(LP)))
anc_rel <- function(g, pk0) {
  a <- LP[set == paste0(g, "_all") & Prom_61520_61540 == TRUE & far_in_SE1 == TRUE, .(s1, e1)][order(s1)]   # 该组 SE1–启动子 loop 的 SE1 侧锚点
  stopifnot(nrow(a) >= 1, all(a$s1 < a$e1))
  ov <- a[s1 < pk0 + BIN & e1 > pk0]; adj <- a[e1 == pk0 | s1 == pk0 + BIN]
  hit <- if (nrow(ov)) ov[1] else if (nrow(adj)) adj[1] else data.table(s1 = NA_integer_, e1 = NA_integer_)
  list(rel = if (nrow(ov)) "重合" else if (nrow(adj)) "紧邻" else "无关", near_s1 = hit$s1, near_e1 = hit$e1,
       anchors = paste(a[, sprintf("%s–%s", format(s1 + 1, big.mark = ","), format(e1, big.mark = ","))], collapse = "、"))
}
rel_Con <- anc_rel("Con", pk$peak_Con); rel_Mhp <- anc_rel("Mhp", pk$peak_Mhp)
stopifnot(rel_Mhp$rel == "重合", rel_Con$rel == "紧邻")                          # 数据若变则此处报错提醒同改图注：Mhp 峰 bin = 其锚点；Con 峰 bin 紧邻其锚点
kb_rng <- function(s0, e0 = s0 + BIN) sprintf("%s–%s kb", format(s0 / 1e3, big.mark = ","), format(e0 / 1e3, big.mark = ","))   # 0-based 半开 → kb 区间文字
PEAK_TXT <- sprintf("Con %s（O/E %.2f）/ Mhp %s（O/E %.2f），均在 SE1 内；Mhp 峰 bin 即 Mhp 的 SE1–启动子 loop 锚点（%s），Con 峰 bin 不是 Con 的任何 loop 锚点、只是紧邻其锚点 %s 的下一个 bin（各组 SE1 侧锚点 1-based：Con %s；Mhp %s）",
                    format(pk$peak_Con + 1, big.mark = ","), pk$oe_peak_Con, format(pk$peak_Mhp + 1, big.mark = ","), pk$oe_peak_Mhp,
                    kb_rng(rel_Mhp$near_s1, rel_Mhp$near_e1), kb_rng(rel_Con$near_s1, rel_Con$near_e1), rel_Con$anchors, rel_Mhp$anchors)
PEAK_EN <- sprintf("In both groups the contact profile peaks within SE1, at or immediately adjacent to the SE1-promoter loop anchor of that group (Mhp, the anchor bin chr13:%s; Con, the bin chr13:%s adjacent to its anchor at %s)",
                   kb_rng(rel_Mhp$near_s1, rel_Mhp$near_e1), kb_rng(pk$peak_Con), kb_rng(rel_Con$near_s1, rel_Con$near_e1))
cat("峰位 vs 锚点：", PEAK_TXT, "\n")
# 经验零分布：视点上游 10–320 kb O/E 之和的 log2(Mhp/Con)，BHLHE40 启动子视点在 chr13 全部视点中的百分位
NB <- fread(F_NUL); vp <- NB[vp_start == VP0]; stopifnot(nrow(vp) == 1)
pct_all <- mean(NB$log2_Mhp_Con <= vp$log2_Mhp_Con) * 100
sim <- NB[(raw_Con + raw_Mhp) >= (vp$raw_Con + vp$raw_Mhp) * 0.5 & (raw_Con + raw_Mhp) <= (vp$raw_Con + vp$raw_Mhp) * 2]
pct_sim <- mean(sim$log2_Mhp_Con <= vp$log2_Mhp_Con) * 100
stopifnot(abs(vp$log2_Mhp_Con - 0.379) < 0.002, pct_all > 99)
# 单元件 bin
e2b <- V[in_E2_3kb == TRUE]; e4b <- V[start == 61385000]; e1b <- V[in_E1 == TRUE]
stopifnot(nrow(e2b) == 1, nrow(e4b) == 1, nrow(e1b) == 1, e2b$raw_Con == 18, e2b$raw_Mhp == 30)
RATIO_TXT <- sprintf("Mhp/Con ×%.2f", se1$ratio_bal)                        # 主图印的描述性比值（×1.22）
stopifnot(RATIO_TXT == "Mhp/Con ×1.22")                                       # 与 README / 图注句一致；数据若变则此处报错提醒同改图注
cat(sprintf("SE1（%d bin，d ≥ 2）：ICE 接触和 Mhp/Con ×%.3f（%d/%d bin Mhp ≥ Con）；O/E 均值 %.2f → %.2f；窗口内非 SE 区 ×%.3f\n",
            se1$n_bins, se1$ratio_bal, se1$n_ge_bal, se1$n_bins, se1$mean_oe_Con, se1$mean_oe_Mhp, nse$ratio_bal))
cat(sprintf("O/E 峰位 Con %s（%.2f）/ Mhp %s（%.2f），均在 SE1 内；上游块 log2 = %.3f，chr13 全体第 %.1f 百分位、接触密度相近视点中第 %.1f 百分位\n",
            format(pk$peak_Con + 1, big.mark = ","), pk$oe_peak_Con, format(pk$peak_Mhp + 1, big.mark = ","), pk$oe_peak_Mhp, vp$log2_Mhp_Con, pct_all, pct_sim))

# ---------------------------------------------------------------- ② 作图数据（长表；阶梯线以 bin 起点为 x，补一行终点使最后一个 bin 画满）
long <- melt(W[, .(start, end, d, oe_Con, oe_Mhp, rpm_Con, rpm_Mhp)], id.vars = c("start", "end", "d"))
long[, `:=`(metric = sub("_.*$", "", variable), group = factor(sub("^.*_", "", variable), levels = grp_levels))]
long[, panel := factor(fifelse(metric == "oe", "O/E", "CPM"), levels = c("O/E", "CPM"))]
long[, x := (start + 1) / 1e3]                                                     # 1-based 起点，kb（与面板 c 同单位）
tail_row <- long[start == max(start)][, `:=`(x = end / 1e3)]                       # 阶梯线最后一段
long <- rbind(long, tail_row)[order(panel, group, x)]
ymax <- long[, .(ymax = max(value, na.rm = TRUE)), by = panel]
WIN_KB <- WIN / 1e3
reg <- data.table(name = c("SE1", "SE2"), xmin = c(EL["SE1_Mhp", start], EL["SE2_Mhp", start]) / 1e3, xmax = c(EL["SE1_Mhp", end], EL["SE2_Mhp", end]) / 1e3)
e2r <- data.table(xmin = EL["E2_21kb", start] / 1e3, xmax = EL["E2_21kb", end] / 1e3)
els <- EL[c("E1", "E2_3kb", "E3", "E4")][, .(lab = c("E1", "E2", "E3", "E4"), x = (start + end) / 2 / 1e3)]
els[, hjust := c(0.5, 0.9, 1.0, 0.0)]                                              # 补充图条带：E3 右对齐 / E4 左对齐错开，不抬高
stopifnot(diff(els$x)[3] > 15, all(diff(els$x) > 0))
vp_x <- (VP0 + BIN / 2) / 1e3
X_BRK <- seq(WIN_KB[1], WIN_KB[2], by = 100)
sx <- function(mult = 0.005) scale_x_continuous(breaks = X_BRK, labels = function(x) format(x, big.mark = ","), expand = expansion(mult = mult))
strw <- function(s, pt, italic = FALSE) systemfonts::string_width(s, family = FONT, size = pt, italic = italic) * 25.4 / 72   # 字串宽（mm）

lo <- if (MAIN_ROWS == "oe") long[panel == "O/E"] else copy(long)
Y_MAX_OE <- ymax[panel == "O/E", ymax]                                             # 9.13（Mhp 峰，SE1 内 loop 锚点 bin）
Y_TOP  <- Y_MAX_OE * 1.32                                                          # 顶部留 32% 放两行注记
PANEL_H <- H_MAIN - 7.5 - (MAIN_ROWS == "both") * 1.5; if (MAIN_ROWS == "both") PANEL_H <- PANEL_H / 2
U_MM  <- Y_TOP / PANEL_H                                                           # 每 mm 对应的 O/E 数据单位
LH_MM <- 2.4                                                                       # 6 pt 行距（mm）
y1 <- Y_TOP - 0.5 * U_MM; y2 <- y1 - LH_MM * U_MM                                  # 第 1 / 2 行文字顶（vjust = 1）
KB_PER_MM <- (WIN_KB[2] - WIN_KB[1]) / (W_MAIN - 9.5 - RIGHT_MM)                   # 数据区宽 ≈ 116 − 左 9.5（轴名 + 刻度）− 右 5 mm 边距；≈ 5.9 kb/mm
SE1_LAB_X <- reg[name == "SE1", (xmin + xmax) / 2]
lab <- rbindlist(list(
  data.table(id = "SE1",   x = SE1_LAB_X, y = y1, hjust = 0.5, label = "SE1", parse = FALSE, col = "grey30", w = strw("SE1", FS_SMALL)),
  if (SHOW_RATIO) data.table(id = "ratio", x = SE1_LAB_X, y = y2, hjust = 0.5, label = RATIO_TXT, parse = FALSE, col = "grey30", w = strw(RATIO_TXT, FS_SMALL)),
  data.table(id = "SE2",   x = reg[name == "SE2", (xmin + xmax) / 2], y = y1, hjust = 0.5, label = "SE2", parse = FALSE, col = "grey30", w = strw("SE2", FS_SMALL))))
lab[, `:=`(x0 = x - hjust * w * KB_PER_MM, x1 = x + (1 - hjust) * w * KB_PER_MM, y_bot = y - 2.1 * U_MM)]   # 文字框（kb；6 pt 字高 ≈ 2.1 mm）
lab[, panel := factor("O/E", levels = levels(long$panel))]
# 自检 1：同一行的注记框互不重叠（≥ 1 mm 间隙）；自检 2：每个注记框下方的阶梯线不穿过文字（框内 bin 的 O/E 最大值 < 文字底）
for (yy in unique(lab$y)) { r <- lab[y == yy][order(x0)]; if (nrow(r) > 1) stopifnot(all(r$x0[-1] - r$x1[-nrow(r)] >= 1 * KB_PER_MM)) }
oe_max_in <- function(a, b) lo[panel == "O/E" & x >= a - 5 & x <= b, max(value, na.rm = TRUE)]
lab[, oe_under := mapply(oe_max_in, x0, x1)]
stopifnot(all(lab$oe_under < lab$y_bot - 0.3 * U_MM), all(lab$x0 >= WIN_KB[1]), all(lab$x1 <= WIN_KB[2]), all(lab$x1 < vp_x - 1 * KB_PER_MM | lab$x0 > vp_x + 1 * KB_PER_MM))   # 注记不压视点点线
cat(sprintf("主图注记：%s\n", paste(lab[, sprintf("%s [%.0f–%.0f kb] 下方 O/E 最大 %.2f < 文字底 %.2f", id, x0, x1, oe_under, y_bot)], collapse = "；")))
p_b <- ggplot() +
  geom_rect(data = reg, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf), fill = "grey93", colour = NA) +
  geom_vline(xintercept = vp_x, linetype = "22", linewidth = lw(0.5), colour = "black") +
  geom_step(data = lo, aes(x, value, colour = group), direction = "hv", linewidth = lw(0.6), na.rm = TRUE) +
  geom_text(data = lab[parse == FALSE], aes(x = x, y = y, label = label, hjust = hjust, colour = col), vjust = 1, size = pt2size(FS_SMALL), family = FONT, show.legend = FALSE) +
  (if (any(lab$parse)) geom_text(data = lab[parse == TRUE],  aes(x = x, y = y, label = label, hjust = hjust, colour = col), vjust = 1, size = pt2size(FS_SMALL), family = FONT, parse = TRUE, show.legend = FALSE) else NULL) +
  scale_colour_manual(values = c(grp_col, grey30 = "grey30", black = "black"), breaks = grp_levels, name = NULL) +
  sx(if (MAIN_ROWS == "oe") 0 else 0.005) +
  scale_y_continuous(breaks = function(l) { b <- scales::breaks_extended(4)(l); b[b <= Y_MAX_OE * 1.05] }, expand = c(0, 0)) +   # 顶部注记区不放刻度（"12" 会贴上缘被裁）
  labs(x = "chr13 position (kb)", y = if (MAIN_ROWS == "oe") "Distance-normalized contacts" else NULL) +
  theme_mhp() + theme(axis.text = element_text(colour = "black"), legend.position = "inside", legend.position.inside = c(0.955, 0.86),
                      axis.title.y = element_text(hjust = 0.3),
                      legend.key.size = unit(2.6, "mm"), legend.text = element_text(size = FS_SMALL), legend.background = element_blank(),
                      legend.key.spacing.y = unit(0.3, "mm"), plot.margin = margin(1.5, RIGHT_MM * 72.27 / 25.4, 1, 1))
if (MAIN_ROWS == "oe") {
  p_b <- p_b + coord_cartesian(xlim = WIN_KB, ylim = c(0, Y_TOP), clip = "off")
} else {                                                                           # 备选：两行（与补充图同），注记只在上行
  p_b <- p_b + facet_wrap(~ panel, ncol = 1, scales = "free_y", strip.position = "left",
                          labeller = as_labeller(c("O/E" = "'Obs. / exp.'", "CPM" = "atop('Contacts per', '10'^6*' cis pairs')"), label_parsed)) +
    coord_cartesian(xlim = WIN_KB, clip = "off") + scale_y_continuous(limits = c(0, NA), breaks = scales::breaks_extended(3), expand = expansion(mult = c(0, 0.32))) +
    theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 90, size = FS_BASE, lineheight = 0.9), panel.spacing.y = unit(1.5, "mm"))
}
sav(p_b, STEM_MAIN, PAN, W_MAIN, H_MAIN, panel_dir = PNL)                          # PNG 300 dpi 草稿 + rds（MHP_FINAL=1 时 SVG + TIFF）
svg_b <- file.path(PAN, paste0(STEM_MAIN, ".svg"))                                 # 另固定出 SVG（svglite 矢量；不出 quartz PDF）
ggsave(svg_b, p_b, width = W_MAIN, height = H_MAIN, units = "mm", device = .svg_dev)
stopifnot(file.exists(svg_b), file.info(svg_b)$size > 0, !any(grepl("<image", readLines(svg_b, warn = FALSE), fixed = TRUE)))   # SVG 内不得有 <image>
cat("✓", basename(svg_b), format(file.info(svg_b)$size, big.mark = ","), "字节（<image> 0 个）\n")
if (identical(Sys.getenv("MHP_FINAL"), "1")) {
  pdf_b <- file.path(PAN, paste0(STEM_MAIN, ".pdf")); ggsave(pdf_b, p_b, width = W_MAIN, height = H_MAIN, units = "mm", device = .pdf_dev)
  stopifnot(file.exists(pdf_b), length(grepRaw("/Subtype\\s*/Image", readBin(pdf_b, "raw", file.size(pdf_b)), all = TRUE)) == 0)
  cat("✓", basename(pdf_b), "（quartz PDF，/Subtype /Image 0 个）\n") }
png_b <- file.path(PAN, paste0(STEM_MAIN, ".png")); stopifnot(file.exists(png_b))
im_b <- png::readPNG(png_b); Hb <- nrow(im_b); Wb <- ncol(im_b)
blk_b <- im_b[, , 1] < 0.25 & im_b[, , 2] < 0.25 & im_b[, , 3] < 0.25
vcol <- which(colMeans(blk_b[round(Hb * 0.25):round(Hb * 0.85), ]) > 0.3 & seq_len(Wb) > Wb * 0.3)   # 点线列：右侧 70% 内黑色占比 > 30%（阶梯线为蓝 / 红，不计）
stopifnot(length(vcol) >= 1, diff(range(vcol)) <= 4)
VP_TXT_INK <- sum(blk_b[1:round(Hb * 0.2), round(Wb * 0.3):(min(vcol) - 4)])     # 点线左侧、面板顶部 20%：旧视点标签位置
stopifnot(VP_TXT_INK == 0); cat(sprintf("✓ 视点点线在 PNG 第 %d 列；其左侧顶部无黑色标签墨迹（视点文字标签已删）\n", min(vcol)))
g_b <- (im_b[, , 1] + im_b[, , 2] + im_b[, , 3]) / 3
EDGE_INK <- c(top = sum(g_b[1:2, ] < 0.9), bottom = sum(g_b[Hb - 0:1, ] < 0.9), left = sum(g_b[, 1:2] < 0.9), right = sum(g_b[, Wb - 0:1] < 0.9))
stopifnot(all(EDGE_INK == 0)); cat("✓ 面板 b PNG 四边最外 2 px 无墨迹（刻度标签未被裁切）\n")

if (MAKE_SUPP) {
  top <- data.table(panel = factor("O/E", levels = levels(long$panel)), y = ymax[panel == "O/E", ymax])
  p_main <- ggplot() +
    geom_rect(data = reg, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf), fill = "grey93", colour = NA) +
    geom_rect(data = e2r, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf), fill = BH_COL, alpha = 0.12, colour = NA) +
    geom_text(data = cbind(reg, top[rep(1, 2)]), aes(x = (xmin + xmax) / 2, y = y * 1.12, label = name), size = pt2size(FS_SMALL), family = FONT, colour = "grey30") +
    geom_vline(xintercept = vp_x, linetype = "22", linewidth = lw(0.5), colour = "black") +
    geom_text(data = top, aes(x = vp_x - 5, y = y * 1.12, label = "Viewpoint"), hjust = 1, size = pt2size(FS_SMALL), family = FONT, colour = "black") +   # 放虚线左侧，免与 SE2 标签相撞
    geom_step(data = long, aes(x, value, colour = group), direction = "hv", linewidth = lw(0.7), na.rm = TRUE) +
    scale_colour_manual(values = grp_col, name = NULL) +
    facet_wrap(~ panel, ncol = 1, scales = "free_y", strip.position = "left",
               labeller = as_labeller(c("O/E" = "'Obs. / exp.'", "CPM" = "atop('Contacts per', '10'^6*' cis pairs')"), label_parsed)) +   # 下行两行 plotmath；"chr13" 放图注（分母 = chr13 顺式接触对总数）
    sx() +
    scale_y_continuous(limits = c(0, NA), breaks = scales::breaks_extended(3), expand = expansion(mult = c(0, 0.28))) +
    coord_cartesian(xlim = WIN_KB, clip = "off") +
    labs(x = NULL, y = NULL) +
    theme_mhp() + theme(axis.text = element_text(colour = "black"), strip.placement = "outside", strip.text.y.left = element_text(angle = 90, size = FS_BASE, lineheight = 0.9),
                        axis.text.x = element_blank(), axis.ticks.x = element_blank(),                       # x 轴放到下方条带
                        panel.spacing.y = unit(1.5, "mm"), legend.position = "inside", legend.position.inside = c(0.93, 0.90),
                        legend.key.size = unit(2.6, "mm"), legend.text = element_text(size = FS_SMALL), legend.background = element_blank(),
                        legend.key.spacing.y = unit(0.3, "mm"), plot.margin = margin(2, 4, 0, 1))
  # E1–E4 条带：与主图共用 x 范围；刻度向上、标签在刻度下；灰 / 红底延续到条带，便于对位
  p_el <- ggplot() +
    geom_rect(data = reg, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf), fill = "grey93", colour = NA) +
    geom_rect(data = e2r, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf), fill = BH_COL, alpha = 0.12, colour = NA) +
    geom_segment(data = els, aes(x = x, xend = x, y = 0.55, yend = 1), linewidth = lw(0.8), colour = BH_COL) +
    geom_text(data = els, aes(x = x, y = 0.45, label = lab, hjust = hjust), vjust = 1, size = pt2size(FS_SMALL), family = FONT, colour = BH_COL) +
    sx() + scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    coord_cartesian(xlim = WIN_KB, clip = "off") +
    labs(x = "chr13 position (kb)", y = NULL) +
    theme_mhp() + theme(axis.text = element_text(colour = "black"), axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                        plot.margin = margin(0, 4, 1, 1))
  p_s <- p_main / p_el + plot_layout(heights = c(1, 0.12))                         # patchwork：面板列自动对齐，同一条 x 轴
  sav(p_s, STEM_SUPP, PAN, W_SUPP, H_SUPP, panel_dir = PNL)
  if (identical(Sys.getenv("MHP_FINAL"), "1")) {
    out <- file.path(PAN, paste0(STEM_SUPP, ".pdf")); ggsave(out, p_s, width = W_SUPP, height = H_SUPP, units = "mm", device = .pdf_dev)
    stopifnot(sum(grepl("/Subtype /Image", readLines(out, warn = FALSE, skipNul = TRUE), useBytes = TRUE)) == 0) }
}

# ---------------------------------------------------------------- ⑤ 数据表与关键数字
fwrite(W[, .(chr = "chr13", bin_start_1based = start + 1, bin_end = end, dist_to_viewpoint_bins = d, plotted = d >= D_MIN,
             raw_Con, raw_Mhp, bal_Con = signif(bal_Con, 5), bal_Mhp = signif(bal_Mhp, 5), oe_Con = signif(oe_Con, 4), oe_Mhp = signif(oe_Mhp, 4),
             cpm_Con = signif(rpm_Con, 4), cpm_Mhp = signif(rpm_Mhp, 4), in_SE1_Mhp, in_SE1_Con, in_SE2_Mhp, in_E1, in_E2_3kb, in_E3, in_E4, in_E2_21kb)],
       file.path(TAB, "07_4C_虚拟4C_作图数据.csv"))
key <- data.table(
  item = c("视点", "上行 O/E", "下行 Contacts per 10^6 cis pairs（补充图；原名 CPM）", "留空", "SE1（Mhp 区间，d ≥ 2）bin 数", "SE1 ICE 接触和 Mhp/Con", "SE1 中 Mhp ≥ Con 的 bin", "SE1 O/E 均值 Con → Mhp",
           "SE1 原始 CPM 和 Mhp/Con", "窗口内非 SE1/SE2 区 ICE 接触和 Mhp/Con（对照）", "O/E 峰位 Con / Mhp", "上游块 log2(Mhp/Con)（10–320 kb O/E 之和）",
           "chr13 全部视点百分位 / 接触密度相近视点百分位", "E2 3 kb bin（61,340,001–61,345,000）", "E4 loop 锚点 bin（61,385,001–61,390,000）", "E1 bin",
           "文库与统计", "两组猪号",
           "主图面板 b 尺寸 / 行数", "主图面板 b 行数选择理由", "主图面板 b 移除项", "主图面板 b 注记", "主图 SE1 比值口径",
           "主图视点标签",
           "图注建议句（主图面板 b）", "图注建议句（补充图两行版）"),
  value = c(sprintf("BHLHE40 TSS %s 所在 5 kb bin：chr13:%s-%s", format(TSS, big.mark = ","), format(VP0 + 1, big.mark = ","), format(VP0 + BIN, big.mark = ",")),
            "仅顺式 ICE 校正接触（与 cooler balance 逐值一致）÷ chr13 同距离期望值", "原始接触计数 / chr13 顺式总接触 × 1e6（未做 ICE；不是文库 CPM）；补充图 x 轴 kb（61,100 … 61,700）；E1–E4 标签在 x 轴上方独立条带",
            "视点 ± 1 bin（距离 < 10 kb）不画", se1$n_bins, sprintf("×%.3f", se1$ratio_bal), sprintf("%d / %d", se1$n_ge_bal, se1$n_bins),
            sprintf("%.2f → %.2f（%d/%d bin Mhp ≥ Con）", se1$mean_oe_Con, se1$mean_oe_Mhp, se1$n_ge_oe, se1$n_bins), sprintf("×%.3f", se1$ratio_rpm),
            sprintf("×%.3f（%d bin，%d 个 Mhp ≥ Con）", nse$ratio_bal, nse$n_bins, nse$n_ge_bal),
            PEAK_TXT,
            sprintf("%.3f（×%.2f）", vp$log2_Mhp_Con, 2^vp$log2_Mhp_Con), sprintf("第 %.1f / 第 %.1f 百分位（%s / %s 个视点）", pct_all, pct_sim, format(nrow(NB), big.mark = ","), format(nrow(sim), big.mark = ",")),
            sprintf("原始 %d → %d，O/E %.2f → %.2f（×%.2f）；★不是 loop 锚点", e2b$raw_Con, e2b$raw_Mhp, e2b$oe_Con, e2b$oe_Mhp, e2b$oe_Mhp / e2b$oe_Con),
            sprintf("原始 %d → %d，O/E %.2f → %.2f（×%.2f）", e4b$raw_Con, e4b$raw_Mhp, e4b$oe_Con, e4b$oe_Mhp, e4b$oe_Mhp / e4b$oe_Con),
            sprintf("原始 %d → %d，O/E ×%.2f（下降）", e1b$raw_Con, e1b$raw_Mhp, e1b$oe_Mhp / e1b$oe_Con),
            "Hi-C libraries pooled from two animals per group：每组 2 头猪合并为一个文库，每组 n = 1，无重复水平统计；图与图注只作描述，不写 P 值",
            "每组 2 头猪合并为一个文库",
            sprintf("%g × %g mm（行 1，面板 a 右侧）；%s", W_MAIN, H_MAIN, if (MAIN_ROWS == "oe") "只画 O/E 一行（MAIN_ROWS = \"oe\"）" else "两行（MAIN_ROWS = \"both\"）"),
            "44 mm 高画两行每行仅约 15 mm，6 pt 刻度 + 两行 plotmath 纵轴名放不下；O/E 已除去距离衰减，是唯一能把 176–320 kb 外的 SE1 与紧邻视点的 SE2 同尺度比较的量；原始接触行只体现距离衰减且与面板 c 的 Hi-C 三角热图重复 → 原始接触行保留在补充图两行版",
            "E1–E4 刻度条带与 E2 21 kb 红带（E1–E4 只在面板 c 引入）；保留 SE1 / SE2 灰底 + 标签、视点点线、Con / Mhp 图例",
            sprintf("顶部 6 pt 只标 SE1 / SE2 几何中心的 \"SE1\" / \"SE2\"（grey30）；纵轴名 \"Distance-normalized contacts\"（= obs./exp.）；x 轴无外扩（61,100–61,700 kb，与面板 c 同），右边距 %g mm（= 面板 c）；脚本自检注记框互不重叠、不压视点点线、注记下方无阶梯线穿过（各框下 O/E 最大 %.2f）",
                    RIGHT_MM, max(lab$oe_under)),
            sprintf("SE1（Mhp ROSE 区间）内 %d 个 5 kb bin（距视点 ≥ 2 bin）的仅顺式 ICE 校正接触之和 Mhp ÷ Con = %.3f（×%.2f）；描述性，每组 n = 1 个合并文库，无统计；图上不印，只作描述性比值；不是 O/E 均值之比（%.2f/%.2f = ×%.2f）也不是原始接触之比（×%.2f）",
                    se1$n_bins, se1$ratio_bal, se1$ratio_bal, se1$mean_oe_Mhp, se1$mean_oe_Con, se1$mean_oe_Mhp / se1$mean_oe_Con, se1$ratio_rpm),
            "视点只以点线表示、图注写明；PNG 自检：点线存在、其左侧顶部无黑色标签墨迹",
            sprintf("Virtual 4C from the BHLHE40 promoter viewpoint (Hi-C of lung tissue, 5-kb bins; viewpoint bin chr13:%s-%s, bins within 10 kb of the viewpoint omitted). y axis, distance-normalized contacts = observed/expected contacts (ICE-balanced contacts divided by the chr13 distance-matched expectation). Grey shading, SE1 and SE2 (Mhp ROSE intervals); dotted line, viewpoint at the BHLHE40 promoter. %s. Hi-C libraries pooled from two animals per group (one merged library per group, n = 1); values are descriptive and no statistical test was applied. Raw contacts per 10^6 chr13 cis contact pairs are shown in Supplementary Fig. X",
                    format(VP0 + 1, big.mark = ","), format(VP0 + BIN, big.mark = ","), PEAK_EN),
            sprintf("Virtual 4C from the BHLHE40 promoter viewpoint (5-kb bin chr13:%s-%s; bins within 10 kb of the viewpoint omitted). Top, ICE-balanced contacts divided by the chr13 distance-matched expectation; bottom, raw contacts per 10^6 chr13 cis contact pairs (not library CPM). Grey, SE1/SE2 (Mhp ROSE intervals); red, the 21,437-bp E2 interval; red ticks below the tracks, E1-E4. The bin containing the 3-kb E2 fragment rose from %d to %d raw contacts (observed/expected %.2f to %.2f), but no loop anchor falls within E2 in either group. Hi-C libraries pooled from two animals per group (one merged library per group); no replicate-level statistics",
                    format(VP0 + 1, big.mark = ","), format(VP0 + BIN, big.mark = ","), e2b$raw_Con, e2b$raw_Mhp, e2b$oe_Con, e2b$oe_Mhp)))
fwrite(key, file.path(TAB, "07_4C_虚拟4C_关键数字.csv"))
stopifnot(file.size(file.path(TAB, "07_4C_虚拟4C_关键数字.csv")) > 0)
print(key[!grepl("^图注建议句", item), .(item, value = substr(value, 1, 90))])
