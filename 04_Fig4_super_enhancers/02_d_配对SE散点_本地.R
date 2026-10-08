# ============================================================================
# 02_d_配对SE散点_本地.R — Fig. 4d：配对超级增强子的 H3K27ac 信号散点（Con 对 Mhp）
# 运行位置：本地；路径变量见 config/paths.R
# ============================================================================

source(file.path(Sys.getenv("MHP_CODE", "."), "config", "paths.R"))

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(ggrepel) })
source(file.path(MHP_CODE, "shared", "_palette_MHP.R"))

SE_DIR <- MHP_ROSE
FIG    <- file.path(MHP_RES, "6、SE/新Fig4_作图")
PAN    <- file.path(FIG, "Plots/面板"); TAB <- file.path(FIG, "Tables"); PNL <- file.path(FIG, "_panels")
for (d in c(PAN, TAB, PNL)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
BH_COL <- "#B2182B"                                                               # ★BHLHE40 高亮红：Fig 4 c / d / e 统一；待并入 _shared/_palette_MHP.R
W_MM <- 50; H_MM <- 56                                                            # 出图尺寸（mm）：整合面板建议 §2.1 的 d = 50 × 56

# ---------------------------------------------------------------- 读表与核对
if (!file.exists(file.path(TAB, "00b_SE配对_1680.csv"))) stop("缺少配对表，请先运行 00b_SE配对_本地重做.R")
P <- fread(file.path(TAB, "00b_SE配对_1680.csv"))
CO <- fread(file.path(SE_DIR, "DE se/SE_log2FC_all_results.csv"))                  # 公司 1,679 对：只核对
stopifnot(nrow(P) == 1680L, nrow(CO) == 1679L, anyDuplicated(P$EnhancerID_Con) == 0L, min(P$OverlapLength) >= 1000L,
          all(P$ActivityScore_Con > 0), all(P$ActivityScore_Mhp > 0), !is.unsorted(-P$AbsoluteScoreDifference), identical(P$rank_absdiff, seq_len(nrow(P))),
          sum(!P$in_company_table) == 1L, P[in_company_table == FALSE, Gene_Con_clean] == "ZFP36L1", P[in_company_table == FALSE, EnhancerID_Con] == "595_peak_287089_lociStitched")
K <- merge(CO[, .(EnhancerID_Con, EnhancerID_Mhp, ActivityScore_Con, ActivityScore_Mhp)], P[in_company_table == TRUE, .(EnhancerID_Con, EnhancerID_Mhp, ActivityScore_Con, ActivityScore_Mhp)],
           by = "EnhancerID_Con", suffixes = c("", ".new"))
stopifnot(nrow(K) == 1679L, identical(K$EnhancerID_Mhp, K$EnhancerID_Mhp.new), max(abs(K$ActivityScore_Con - K$ActivityScore_Con.new)) < 1e-6,
          max(abs(K$ActivityScore_Mhp - K$ActivityScore_Mhp.new)) < 1e-6, identical(P[in_company_table == TRUE, EnhancerID_Con], CO$EnhancerID_Con))   # 1,679 对逐行一致、顺序相同
P[, l2 := log2(ActivityScore_Mhp / ActivityScore_Con)]
stopifnot(max(abs(P$l2 - P$SE_log2FC_exact)) < 1e-9, max(abs(P$l2 - P$SE_log2FC)) < 1e-3)   # 与文件列一致（SE_log2FC 为 +1 伪计数版本）
P[, rank_gain_l2 := frank(-l2, ties.method = "first")]                            # 按倍数的增益名次
P[, rank_gain_diff := frank(ScoreDifference_Con_minus_Mhp, ties.method = "first")] # 按差值的增益名次（Con − Mhp 最负 = 增益最大）
P[, cls := fifelse(l2 >= 1, "Mhp", fifelse(l2 <= -1, "Con", "ns"))]
NC <- P[, .N, by = cls]; n_mhp <- NC[cls == "Mhp", N]; n_con <- NC[cls == "Con", N]; n_ns <- NC[cls == "ns", N]
stopifnot(n_mhp == 120L, n_con == 109L, n_ns == 1451L, n_mhp + n_con + n_ns == 1680L)   # ZFP36L1 对 log2FC 0.26 → 灰点；≥ 2 倍计数与 1,679 对时相同
rho <- cor(log10(P$ActivityScore_Con), log10(P$ActivityScore_Mhp), method = "spearman")
stopifnot(abs(rho - 0.6726) < 0.001)                                              # 1,679 对时 0.6721

# BHLHE40 两条 SE（SE1 = chr13:61,209,361–61,517,999，跨 ITPR1 3′ 段 ~190 kb 并延伸到 BHLHE40 TSS 上游 ~8 kb；SE2 = 61,531,355–61,583,348，只与 BHLHE40 3′ 末端外显子重叠、主体在下游）
B <- P[Gene_Con_clean == "BHLHE40" & Gene_Mhp_clean == "BHLHE40"][order(-ActivityScore_Mhp)]
stopifnot(nrow(B) == 2L, B$EnhancerID_Con[1] == "88_peak_78108_lociStitched", B$EnhancerID_Mhp[1] == "115_peak_71508_lociStitched",
          B$rank_absdiff[1] == 35L, B$rank_gain_diff[1] == 25L, B$rank_gain_l2[1] == 114L, B$rank_absdiff[2] == 378L,   # 1,679 对时 34 / 24 / 114 / 377：ZFP36L1 对（|差值| 109,789，增益）插到前面
          abs(B$l2[1] - 1.0397) < 1e-3, abs(B$ActivityScore_Mhp[1] / B$ActivityScore_Con[1] - 2.056) < 0.005,
          abs(B$l2[2] - 1.0796) < 1e-3, all(B$cls == "Mhp"))
# ZFP36L1 对（两组 ROSE rank 1；公司表漏掉的一对）：灰色斜体标签
Z <- P[Gene_Con_clean == "ZFP36L1"]; stopifnot(nrow(Z) == 1L, Z$cls == "ns", Z$in_company_table == FALSE, abs(Z$l2 - 0.2594) < 1e-3)
Z[, `:=`(kind = "ZFP36L1_rank1_pair", lab = 'italic("ZFP36L1")')]
B[, `:=`(kind = c("BHLHE40_SE1", "BHLHE40_SE2"), lab = c('bolditalic("BHLHE40")~"SE1"', 'bolditalic("BHLHE40")~"SE2"'))]

# 前 5 增益 SE（按倍数；标注 Mhp 侧最近基因）
T5 <- P[order(rank_gain_l2)][1:5][, `:=`(kind = "top5_gain_log2FC", lab = sprintf('italic("%s")', Gene_Mhp_clean))]
T5d <- P[order(rank_gain_diff)][1:5][, kind := "top5_gain_diff"]
stopifnot(identical(T5$Gene_Mhp_clean, c("TNS3", "ARID1B", "FGGY", "CALM1", "NAP1L5")), all(T5$cls == "Mhp"))
LAB <- rbind(B, T5, Z)
sel <- c("kind", "rank_gain_l2", "rank_gain_diff", "rank_absdiff", "EnhancerID_Con", "Gene_Con_clean", "ActivityScore_Con",
         "EnhancerID_Mhp", "Gene_Mhp_clean", "ActivityScore_Mhp", "OverlapLength", "l2")
fwrite(rbind(LAB[, ..sel], T5d[, ..sel])[, `:=`(l2 = round(l2, 3), ActivityScore_Con = round(ActivityScore_Con, 1), ActivityScore_Mhp = round(ActivityScore_Mhp, 1))],
       file.path(TAB, "02_d_标注点.csv"))
KN <- data.table(item = c("配对 SE 对数（≥1 kb 重叠；00b 本地重做）", "较公司表多出的对", "Mhp 高（log2FC ≥ 1）对数", "Con 高（log2FC ≤ −1）对数", "|log2FC| < 1 对数",
                          "Mhp 高（log2FC ≥ 0.58）对数", "Con 高（log2FC ≤ −0.58）对数", "Spearman ρ（log10 信号）",
                          "BHLHE40-SE1 Con 信号", "BHLHE40-SE1 Mhp 信号", "BHLHE40-SE1 log2FC", "BHLHE40-SE1 按 |Con−Mhp| 双向名次（非增益名次）", "BHLHE40-SE1 按倍数增益名次", "BHLHE40-SE1 按差值增益名次",
                          "BHLHE40-SE2 Con 信号", "BHLHE40-SE2 Mhp 信号", "BHLHE40-SE2 log2FC", "BHLHE40-SE2 按 |差值| 名次",
                          "ZFP36L1 对 Con 信号 / Mhp 信号 / log2FC / 按 |差值| 名次",
                          "配对中 Con 信号范围", "配对中 Mhp 信号范围", "配对 Mhp SE 去重数（76 个重复）", "配对 Con SE 去重数",
                          "正文名次句（建议）", "旧口径（公司 1,679 对）对照", "BHLHE40 高亮色（c/d/e 统一）"),
                 value = c(nrow(P), sprintf("%s ↔ %s（ZFP36L1，两组 ROSE rank 1，重叠 %s bp）", Z$EnhancerID_Con, Z$EnhancerID_Mhp, format(Z$OverlapLength, big.mark = ",")),
                           n_mhp, n_con, n_ns, sum(P$l2 >= 0.58), sum(P$l2 <= -0.58), round(rho, 3),
                           round(B$ActivityScore_Con[1], 1), round(B$ActivityScore_Mhp[1], 1), round(B$l2[1], 3), B$rank_absdiff[1], B$rank_gain_l2[1], B$rank_gain_diff[1],
                           round(B$ActivityScore_Con[2], 1), round(B$ActivityScore_Mhp[2], 1), round(B$l2[2], 3), B$rank_absdiff[2],
                           sprintf("%.1f / %.1f / %.3f / %d", Z$ActivityScore_Con, Z$ActivityScore_Mhp, Z$l2, Z$rank_absdiff),
                           paste(round(range(P$ActivityScore_Con)), collapse = "–"), paste(round(range(P$ActivityScore_Mhp)), collapse = "–"),
                           uniqueN(P$EnhancerID_Mhp), uniqueN(P$EnhancerID_Con),
                           sprintf("1,680 对配对 SE 中按信号差绝对值排第 %d（增益按差值第 %d、按倍数第 %d）", B$rank_absdiff[1], B$rank_gain_diff[1], B$rank_gain_l2[1]),
                           "1,679 对：|差值| 第 34、增益按差值第 24、按倍数第 114；ρ 0.672；≥ 2 倍上调 120 / 下调 109（不变）", BH_COL))
fwrite(KN, file.path(TAB, "02_d_配对SE关键数字.csv"))

# ---------------------------------------------------------------- 作图
COL <- c(grp_col, ns = "#C8C8C8")
P[, cls := factor(cls, levels = c("ns", "Con", "Mhp"))]; setorder(P, cls)     # 灰点先画，彩点在上
lim <- c(min(c(P$ActivityScore_Con, P$ActivityScore_Mhp)) * 0.45, max(c(P$ActivityScore_Con, P$ActivityScore_Mhp)) * 1.3)
LL <- log10(lim)
xlim <- c(lim[1] * 0.8, max(P$ActivityScore_Con) * 1.21); XLL <- log10(xlim)        # Con 最大现为 ZFP36L1 557,436（log10 5.75）→ 上限 5.83
stopifnot(XLL[2] > max(log10(P$ActivityScore_Con)) + 0.05, XLL[2] < LL[2])
mm_per_unit <- (W_MM - 9.3) / (XLL[2] - XLL[1])                                   # 面板每 log10 单位约 (W_MM − 9.3) / x 范围 mm（y 轴标题 / 刻度 / 边距 ≈ 9.3 mm）
# ZFP36L1 标签：右对齐置于点左侧（间隔 0.05 log10 单位），标签中心 = 点 − 0.05 − 半宽（7 字符斜体 6 pt ≈ 1.18 mm / 字符 + 0.6 mm 内边距）；白色描边
Z_HW <- (1.18 * 7 + 0.6) / 2 / mm_per_unit; Z_NX <- -(0.05 + Z_HW)
stopifnot(log10(Z$ActivityScore_Con) + Z_NX - Z_HW > 4.109 + 0.05,                # 标签左缘仍在点云区右侧远处（> 5.2）、不进左侧标签竖条
          log10(Z$ActivityScore_Mhp) - 0.10 > max(log10(P[Gene_Con_clean != "ZFP36L1", ActivityScore_Mhp])) + 0.05)   # 标签行（y = 点）之下 0.1 单位内无其他点（次高 KLF6 5.54）
# 面板内文字：左上 n 对 + ρ（顶行，ZFP36L1 对在右上、与其不同 x 区间）；右下两行 Mhp 高（红）/ Con 高（深蓝 smp_col Con_1）计数，位于所有点之下（y 下限已留出）
TXT <- data.table(x = c(xlim[1] * 1.06, xlim[2] / 1.06, xlim[2] / 1.06), y = c(lim[2] / 1.04, lim[1] * 1.045 * 10^0.10, lim[1] * 1.045),
                  hjust = c(0, 1, 1), vjust = c(1, 0, 0), colour = c("black", grp_col[["Mhp"]], smp_col[["Con_1"]]),
                  lab = c(sprintf('"n = %s pairs,"~rho~"= %.2f"', prettyNum(nrow(P), big.mark = ","), rho),
                          sprintf('"Higher in Mhp, n = %d"', n_mhp), sprintf('"Higher in Con, n = %d"', n_con)))
stopifnot(log10(min(P$ActivityScore_Mhp)) > log10(lim[1] * 1.045 * 10^0.10) + 0.05)   # 右下两行文字整体在最低点之下
stopifnot(log10(Z$ActivityScore_Mhp) + 0.07 < log10(lim[2] / 1.04) - 0.14 || log10(Z$ActivityScore_Con) + Z_NX - Z_HW > XLL[1] + 1.45)   # ZFP36L1 标签与顶行 n / ρ 文字（左起约 1.4 单位宽）不同行或不同 x 区间
# 前 5 增益 SE 标签：右对齐于 XR（log10），放在数据左侧空白竖条内（所有点 log10 x ≥ 4.109）；5 个槽位 y 紧贴这 5 个点（y 4.94–5.20）
XR <- 4.085; SL <- seq(4.95, 5.35, by = 0.10)
px <- log10(T5$ActivityScore_Con); py <- log10(T5$ActivityScore_Mhp)
T5_W <- (1.18 * max(nchar(T5$Gene_Mhp_clean)) + 0.6) / mm_per_unit                # 最长前 5 标签宽（6 字符斜体 6 pt ≈ 1.18 mm / 字符 + 0.6 mm 内边距）换算成 log10 单位
stopifnot(XR + 0.015 < min(log10(P$ActivityScore_Con)), XR - T5_W > XLL[1] + 0.05)   # 竖条内无点；最长标签左缘与 y 轴线间距 ≥ 0.05 单位（≈ 0.9 mm）
ori <- function(ax, ay, bx, by, cx, cy) sign((bx - ax) * (cy - ay) - (by - ay) * (cx - ax))
seg_cross <- function(a, b, c, d) ori(a[1], a[2], b[1], b[2], c[1], c[2]) != ori(a[1], a[2], b[1], b[2], d[1], d[2]) &&
                                   ori(c[1], c[2], d[1], d[2], a[1], a[2]) != ori(c[1], c[2], d[1], d[2], b[1], b[2])
perms <- as.matrix(expand.grid(rep(list(1:5), 5))); perms <- perms[apply(perms, 1, function(r) length(unique(r)) == 5), ]
score <- t(apply(perms, 1, function(pm) {
  e <- cbind(XR + 0.012, SL[pm]); s <- cbind(px, py); cr <- 0
  for (i in 1:4) for (j in (i + 1):5) cr <- cr + seg_cross(s[i, ], e[i, ], s[j, ], e[j, ])
  c(cr, sum(sqrt((e[, 1] - s[, 1])^2 + (e[, 2] - s[, 2])^2)))
}))
best <- perms[order(score[, 1], score[, 2])[1], ]; stopifnot(score[order(score[, 1], score[, 2])[1], 1] == 0)   # 必须存在不交叉方案
T5L <- copy(T5)[, sl := SL[best]]
seg_len <- sqrt((XR + 0.012 - px)^2 + (SL[best] - py)^2)
stopifnot(max(seg_len) < 0.5)                                                     # 引线均 < 0.5 log10 单位（≈ 11 mm）
# 引线是否擦过其他点（点到线段距离 < 0.02 log10 单位 ≈ 0.45 mm）：只回显供核对，不作硬性要求（点云内难以完全避开）
pt_seg_d <- function(qx, qy, a, b) { v <- b - a; t <- pmin(1, pmax(0, ((qx - a[1]) * v[1] + (qy - a[2]) * v[2]) / sum(v^2))); sqrt((qx - a[1] - t * v[1])^2 + (qy - a[2] - t * v[2])^2) }
n_graze <- sum(sapply(1:5, function(i) { d <- pt_seg_d(log10(P$ActivityScore_Con), log10(P$ActivityScore_Mhp), c(px[i], py[i]), c(XR + 0.012, SL[best][i]))
                                         sum(d < 0.02) - 1 }))
cat(sprintf("前 5 标签引线：槽位 %s；最长 %.2f log10 单位；被引线擦过的其他点 %d 个\n", paste(T5$Gene_Mhp_clean, "→", SL[best], collapse = ", "), max(seg_len), n_graze))
p <- ggplot(P, aes(ActivityScore_Con, ActivityScore_Mhp)) +
  geom_abline(slope = 1, intercept = 0, linewidth = lw(0.4), colour = "grey45") +
  geom_abline(slope = 1, intercept = c(log10(2), -log10(2)), linewidth = lw(0.3), linetype = "12", colour = "grey60") +   # 双 log10 轴上 2 倍 / 0.5 倍为平行线
  geom_point(aes(colour = cls), size = 0.55, stroke = 0, alpha = 0.85) +
  geom_point(data = B, shape = 21, size = 2.0, stroke = lw(0.6), colour = BH_COL, fill = NA) +       # BHLHE40 空心圈
  # 前 5 增益 SE：标签手动排成一列（log10 x = 4.085 右对齐，位于数据左侧空白竖条；y 4.95–5.35 紧贴点群；ggrepel 对 plotmath 盒子估计偏小会互相压叠），
  #   槽位分配穷举 5! 种，取引线不交叉且总长最短者；引线从点画到标签右缘
  geom_segment(data = T5L, aes(x = ActivityScore_Con, y = ActivityScore_Mhp, xend = 10^(XR + 0.012), yend = 10^sl), linewidth = lw(0.3), colour = "grey45") +
  geom_text(data = T5L, aes(10^XR, 10^sl, label = lab), parse = TRUE, hjust = 1, size = pt2size(FS_SMALL), family = FONT, colour = "black") +
  # BHLHE40：SE1 标签放到对角线右下方的稀疏区，SE2 标签右下（白色描边保证压在灰点上仍可读）
  geom_text_repel(data = B, aes(label = lab), parse = TRUE, size = pt2size(FS_SMALL), family = FONT, colour = BH_COL,
                  nudge_x = c(0.32, 0.30), nudge_y = c(-0.12, -0.11), direction = "both", bg.colour = "white", bg.r = 0.08,
                  segment.size = lw(0.3), segment.colour = "grey45", min.segment.length = 0, box.padding = 0.15, point.padding = 0.3,
                  max.overlaps = Inf, force = 1, seed = 7) +
  # ZFP36L1（两组 rank 1 的孤立灰点）：固定位（force = 0）、无引线、白色描边压过 2 倍点线
  geom_text_repel(data = Z, aes(label = lab), parse = TRUE, size = pt2size(FS_SMALL), family = FONT, colour = "grey30",
                  nudge_x = Z_NX, nudge_y = 0, force = 0, force_pull = 0, min.segment.length = Inf, bg.colour = "white", bg.r = 0.08,
                  box.padding = 0, point.padding = 0, max.overlaps = Inf, seed = 7) +
  geom_text(data = TXT, aes(x, y, label = lab, hjust = hjust, vjust = vjust), parse = TRUE, colour = TXT$colour, size = pt2size(FS_SMALL), family = FONT) +
  scale_colour_manual(values = COL, guide = "none") +
  scale_x_log10(breaks = 10^(4:5), minor_breaks = as.vector(outer(2:9, 10^(3:5))), labels = scales::label_log()) +
  scale_y_log10(breaks = 10^(4:5), minor_breaks = as.vector(outer(2:9, 10^(3:5))), labels = scales::label_log()) +
  guides(x = guide_axis(minor.ticks = TRUE), y = guide_axis(minor.ticks = TRUE)) +
  coord_fixed(xlim = xlim, ylim = lim, expand = FALSE, clip = "off") +                                  # coord 限范围：不丢点、不丢标签；x 上限裁掉右上空白
  labs(x = "Con SE H3K27ac signal", y = "Mhp SE H3K27ac signal") +
  theme_mhp() + theme(plot.margin = margin(2, 3, 2, 2), axis.minor.ticks.length = rel(0.5),
                      axis.minor.ticks.x.bottom = element_line(linewidth = lw(0.4)), axis.minor.ticks.y.left = element_line(linewidth = lw(0.4)))
# BHLHE40 SE1 标签右缘（中心 = 点 + nudge，半宽按 11 字符 × 0.62 em × 6 pt ≈ 12.3 mm 估）须在 x 上限之内；每 log10 单位 ≥ 18 mm（1,679 对版 20.9，现约 19.2）
stopifnot(log10(B$ActivityScore_Con[1]) + 0.32 + 12.3 / 2 / mm_per_unit < XLL[2], mm_per_unit > 18)
cat(sprintf("x 轴 log10 范围 %.3f–%.3f（y 轴 %.3f–%.3f）；约 %.1f mm / log10 单位；标签竖条 ≈ %.1f mm、最长前 5 标签 ≈ %.1f mm、与 y 轴线净空 ≈ %.1f mm\n",
            XLL[1], XLL[2], LL[1], LL[2], mm_per_unit, (4.109 - XLL[1]) * mm_per_unit, T5_W * mm_per_unit, (XR - T5_W - XLL[1]) * mm_per_unit))
sav(p, "Fig4d_配对SE散点", PAN, W_MM, H_MM, panel_dir = PNL)
if (identical(Sys.getenv("MHP_FINAL"), "1")) ggsave(file.path(PAN, "Fig4d_配对SE散点.pdf"), p, width = W_MM, height = H_MM, units = "mm", device = .pdf_dev)

cat(sprintf("n = %d 对；Mhp 高 %d / Con 高 %d / 其余 %d（|log2FC| ≥ 1）；Spearman ρ = %.3f\n", nrow(P), n_mhp, n_con, n_ns, rho))
print(LAB[, .(kind, Gene_Con_clean, Gene_Mhp_clean, Con = round(ActivityScore_Con), Mhp = round(ActivityScore_Mhp), l2 = round(l2, 2), rank_gain_l2, rank_gain_diff, rank_absdiff)])
