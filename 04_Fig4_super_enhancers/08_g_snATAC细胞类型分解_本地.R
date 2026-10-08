# ============================================================================
# 08_g_snATAC细胞类型分解_本地.R — Fig. 4g：按 snATAC-seq 伪 bulk 可及性把超级增强子指派到主导细胞类型并比较感染变化
# 运行位置：本地；路径变量见 config/paths.R
# ============================================================================

source(file.path(Sys.getenv("MHP_CODE", "."), "config", "paths.R"))

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(vipor) })
source(file.path(MHP_CODE, "shared", "_palette_MHP.R"))

ASSIGN <- Sys.getenv("FIG4G_ASSIGN", "balanced")
MAIN   <- "balanced"                                     # 主图方案；ASSIGN != MAIN 时输出改到 备选_指派方案/<ASSIGN>/
stopifnot(ASSIGN %in% c("balanced", "pooled"))
SCHEMES <- c("pooled", "Con_only", "Mhp_only", "balanced")   # 全部写入 Source Data 的四种方案
ATAC_TAB <- file.path(MHP_DIR$atac, "Tables/24_SE")
PROP_CSV <- file.path(MHP_DIR$atac, "Plots/10_CellProportion/CellProportion_byGroup.csv")
OUT <- file.path(MHP_RES, "6、SE/新Fig4_作图")
if (ASSIGN == MAIN) { PAN <- file.path(OUT, "Plots/面板"); TAB <- file.path(OUT, "Tables"); PNL <- file.path(OUT, "_panels")
} else              { PAN <- file.path(OUT, "备选_指派方案", ASSIGN, "Plots"); TAB <- file.path(OUT, "备选_指派方案", ASSIGN, "Tables"); PNL <- NULL }
for (d in c(PAN, TAB, PNL)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
SENS_CSV <- file.path(OUT, "Tables/15_g_敏感性_本地结果.csv")   # 15_g 独立复算（存在时交叉核对 balanced）
BH_COL <- "#B2182B"                                   # ★BHLHE40 高亮红：Fig 4 各面板统一；待并入 _shared/_palette_MHP.R
PC     <- 0.1                                         # 伪计数（CPM），记录口径（25_SE归因热图_本地.R）
EXCL   <- "B_cells"                                   # 不上图的类型（Mhp 组仅 3 个核）
BH_SE  <- "SE1977"                                    # BHLHE40 SE1 在合并 SE 表中的编号
cat(sprintf("主导类型指派方案 ASSIGN = %s（主图方案 = %s）；输出目录 %s\n", ASSIGN, MAIN, PAN))

# ---------------------------------------------------------------- ① 读表、排除 E2 子区间
a <- fread(file.path(ATAC_TAB, "SE_CPM_byCelltype_fixed.csv"))
g <- fread(file.path(ATAC_TAB, "SE_CPM_byCelltypeGroup_fixed.csv"))
stopifnot(nrow(a) == 2856L, nrow(g) == 2856L, identical(a$SE_ID, g$SE_ID), !anyNA(a), !anyNA(g),
          all(ct_atac_order %in% names(a)),
          all(c(paste0("Mhp_", ct_atac_order), paste0("Mock_", ct_atac_order)) %in% names(g)))
a <- a[SE_ID != "E2_BHLHE40"]; g <- g[SE_ID != "E2_BHLHE40"]
stopifnot(nrow(a) == 2855L, identical(a[, .N, by = src][order(src)]$N, c(1555L, 621L, 679L)))   # Both / Con_only / Mhp_only
N <- nrow(a); ii <- seq_len(N)

# ---------------------------------------------------------------- ② 主导类型：四种方案各算一遍（主图用 ASSIGN 方案）
M  <- as.matrix(a[, ..ct_atac_order])                                                       # 两组合并伪 bulk CPM
Mc <- as.matrix(g[, paste0("Mock_", ct_atac_order), with = FALSE]); colnames(Mc) <- ct_atac_order
Mh <- as.matrix(g[, paste0("Mhp_",  ct_atac_order), with = FALSE]); colnames(Mh) <- ct_atac_order
SCORE <- list(pooled = log2(M + 1), Con_only = log2(Mc + 1), Mhp_only = log2(Mh + 1), balanced = (log2(Mc + 1) + log2(Mh + 1)) / 2)
argmax_ct <- function(L) { top2 <- t(apply(L, 1, function(v) sort(v, decreasing = TRUE)[1:2])); stopifnot(all(top2[, 1] > top2[, 2]))   # 无并列
  list(j = max.col(L, ties.method = "first"), spec = top2[, 1] - top2[, 2]) }
DOM <- lapply(SCORE, argmax_ct)
CNT_DOC <- list(pooled   = c(412L, 358L, 302L, 419L, 117L, 427L, 420L, 400L), balanced = c(422L, 337L, 287L, 476L, 124L, 407L, 425L, 377L),
                Con_only = c(442L, 337L, 306L, 445L,  80L, 418L, 427L, 400L), Mhp_only = c(368L, 329L, 252L, 599L, 169L, 391L, 417L, 330L))
for (s in SCHEMES) stopifnot(identical(as.integer(table(factor(ct_atac_order[DOM[[s]]$j], levels = ct_atac_order))), CNT_DOC[[s]]))
j <- DOM[[ASSIGN]]$j; dom <- ct_atac_order[j]
a[, `:=`(assignment = ASSIGN, dominant = dom, specificity_log2 = DOM[[ASSIGN]]$spec, CPM_dominant_all = M[cbind(ii, j)])]

# ---------------------------------------------------------------- ③ 主导类型伪 bulk 内 log2(Mhp/Con)（四方案；主图列 = ASSIGN 方案）
LFC     <- log2((Mh + PC) / (Mc + PC))                                                     # 2,855 × 8：每个 SE 在每一类伪 bulk 内的 log2FC
LFC_pc1 <- log2((Mh + 1)  / (Mc + 1))                                                      # 【探针】口径，仅入表对照
for (s in SCHEMES) { js <- DOM[[s]]$j; a[, (paste0("dominant_", s)) := ct_atac_order[js]]; a[, (paste0("log2FC_", s)) := LFC[cbind(ii, js)]] }
a[, `:=`(CPM_dominant_Con = Mc[cbind(ii, j)], CPM_dominant_Mhp = Mh[cbind(ii, j)], log2FC_dominant = LFC[cbind(ii, j)], log2FC_dominant_pc1 = LFC_pc1[cbind(ii, j)])]
a[, rank_within_dominant := frank(-log2FC_dominant, ties.method = "min"), by = dominant]  # 1 = 类内增幅最大
a[, consensus_Con_Mhp := dominant_Con_only == dominant_Mhp_only]
a[, on_figure := dominant != EXCL]
a[, is_BHLHE40_SE1 := SE_ID == BH_SE]
stopifnot(identical(a$dominant, a[[paste0("dominant_", ASSIGN)]]), all(abs(a$log2FC_dominant - a[[paste0("log2FC_", ASSIGN)]]) < 1e-12))

# 各组细胞核数（图注 / Spearman 用；B 细胞 Mhp 3 / Mock 132）
prop <- fread(PROP_CSV)
stopifnot(setequal(prop$celltype, ct_atac_order), setequal(prop$group, c("Mock", "Mhp")))
prop[, group := to_grp(group)]
nuc <- dcast(prop, celltype ~ group, value.var = "count")
stopifnot(nuc[celltype == "B_cells", Mhp] == 3L, nuc[celltype == "B_cells", Con] == 132L, sum(nuc$Con) == 5044L, sum(nuc$Mhp) == 4566L)
nuc[, nuclei_ratio := Mhp / Con]

# 各方案 × 8 类汇总（长表）+ Spearman(7 类中位数, 核数比 Mhp/Con)
sum_scheme <- function(s) { js <- DOM[[s]]$j; d <- data.table(dominant = ct_atac_order[js], v = LFC[cbind(ii, js)], v1 = LFC_pc1[cbind(ii, js)])
  d[, .(n_SE = .N, median = median(v), mean = mean(v), q1 = quantile(v, .25), q3 = quantile(v, .75), min = min(v), max = max(v),
        pct_up = 100 * mean(v > 0), pct_down = 100 * mean(v < 0), median_pc1 = median(v1)), by = dominant][, assignment := s] }
sm_all <- rbindlist(lapply(SCHEMES, sum_scheme))
sm_all <- merge(sm_all, nuc[, .(dominant = celltype, nuclei_Con = Con, nuclei_Mhp = Mhp, nuclei_ratio)], by = "dominant")
sm_all[, `:=`(is_main = assignment == ASSIGN, on_figure = dominant != EXCL, label = lab_atac(dominant))]
sm_all[on_figure == TRUE, spearman_median_vs_nuclei_ratio := cor(median, nuclei_ratio, method = "spearman"), by = assignment]
GLOBAL_all <- sapply(SCHEMES, function(s) median(LFC[cbind(ii, DOM[[s]]$j)]))                # 各方案全局中位（含 B）

DOC <- c(Plasma_cells = 0.359, Neutrophils = 0.263, Macrophages = 0.031, B_cells = 0.018, Epithelial_cells = 0.012,
         T_NK_cells = -0.055, Endothelial_cells = -0.134, Fibroblasts = -0.187)
sm_all[, documented_median := ifelse(assignment == "pooled", DOC[dominant], NA_real_)]
sm_all[, diff_vs_documented := median - documented_median]
smp <- sm_all[assignment == "pooled"]
stopifnot(all(abs(smp$diff_vs_documented) < 0.0015), abs(GLOBAL_all[["pooled"]] - (-0.011)) < 0.0015,
          abs(smp[dominant == "Plasma_cells", median_pc1] - 0.28) < 0.006,                   # 【探针】+0.28 / +0.21 / +0.02 = 伪计数 1 口径
          abs(smp[dominant == "Neutrophils",  median_pc1] - 0.21) < 0.006,
          abs(smp[dominant == "Macrophages",  median_pc1] - 0.02) < 0.006,
          abs(smp[on_figure == TRUE, spearman_median_vs_nuclei_ratio][1] - 1) < 1e-9)      # pooled：中位数与核数比完全同序（15_g 的出发点）
# balanced 交叉核对：与 15_g_敏感性_本地结果.csv（独立复算，scheme 名含 "balanced"）逐类比 n 与中位数
sens_note <- "15_g 表不存在，未交叉核对（先跑 15_g_敏感性_本地.R 再重跑本脚本即可）"
if (file.exists(SENS_CSV)) {
  sens <- fread(SENS_CSV)[estimator == "pseudo-bulk CPM" & grepl("^balanced", scheme)]
  chk <- merge(sm_all[assignment == "balanced", .(dominant, n_SE, median)], sens[, .(dominant, n_15g = n_SE, med_15g = median)], by = "dominant")
  stopifnot(nrow(chk) == 8L, all(chk$n_SE == chk$n_15g), all(abs(chk$median - chk$med_15g) < 1e-3))
  sens_note <- sprintf("与 15_g_敏感性_本地结果.csv 独立复算逐类一致（8 类 n 全等、中位数差 < 1e-3；scheme = %s）", sens$scheme[1])
}
cat(sens_note, "\n")
sm <- sm_all[assignment == ASSIGN]                                                          # 主图方案汇总
GLOBAL <- GLOBAL_all[[ASSIGN]]
stopifnot(all(abs(a[on_figure == TRUE, log2FC_dominant]) < 2.05))                           # 上图 7 类范围 −1.79–1.97（pooled 与 balanced 同），x 轴 ±2.1 不需封顶（B 细胞 Mhp 多为 0，log2FC 可到 −7，不上图）

# BHLHE40 SE1（SE1977）核对：基因、主导类型（四方案均中性粒）、坐标（= Con ∪ Mhp ROSE SE1 并集）、CPM、log2FC、类内名次
bh <- a[SE_ID == BH_SE]
reg <- fread(file.path(ATAC_TAB, "SE_regions_ensembl.tsv"))
r1  <- reg[gene_clean == "BHLHE40" & chrom == "chr13" & start < 61.3e6]                    # 两组各自的 SE1（Mhp rank 48 / Con rank 304）
RANK_DOC <- c(pooled = 22L, balanced = 13L, Con_only = 5L, Mhp_only = 21L)
bh_rank <- sapply(SCHEMES, function(s) { v <- a[[paste0("log2FC_", s)]][a[[paste0("dominant_", s)]] == "Neutrophils"]; sum(v >= bh[[paste0("log2FC_", s)]]) })
stopifnot(nrow(bh) == 1L, bh$gene == "BHLHE40", all(unlist(bh[, paste0("dominant_", SCHEMES), with = FALSE]) == "Neutrophils"), bh$nPeak == 143L,
          nrow(r1) == 2L, setequal(r1$group, c("Mhp", "Con")), min(r1$start) == 61209361L, max(r1$end) == 61517999L,
          bh$width == max(r1$end) - min(r1$start) + 1L,                                     # fixed 表 width = end − start + 1 = 308,639（01_c / 06_g 按 STOP − START = 308,638，同一区间）
          abs(bh$Neutrophils - 5.5914) < 1e-3, abs(bh$Macrophages - 3.0308) < 1e-3,
          abs(bh$CPM_dominant_Con - 3.3328) < 1e-3, abs(bh$CPM_dominant_Mhp - 7.1187) < 1e-3, # 中性粒伪 bulk 3.33 → 7.12（+114%）
          abs(bh$log2FC_dominant - 1.072) < 0.002, identical(unname(bh_rank[SCHEMES]), unname(RANK_DOC[SCHEMES])), bh$rank_within_dominant == RANK_DOC[[ASSIGN]])
BH_COORD <- sprintf("chr13:%s-%s", format(min(r1$start), big.mark = ","), format(max(r1$end), big.mark = ","))
N_NEUT <- sapply(SCHEMES, function(s) sm_all[assignment == s & dominant == "Neutrophils", n_SE])
top2_of <- function(v) { o <- order(v, decreasing = TRUE)[1:2]; list(ct = names(v)[o], val = unname(v[o])) }
bh_i   <- which(a$SE_ID == BH_SE)
bh_sc  <- top2_of(setNames(SCORE[[ASSIGN]][bh_i, ], ct_atac_order))                       # ASSIGN 方案得分（balanced：中性粒 2.568 vs T/NK 2.021 → 特异性 0.547）
bh_cpm <- top2_of(setNames(M[bh_i, ], ct_atac_order))                                      # 两组合并伪 bulk CPM（中性粒 5.59 vs 浆 3.11；巨噬 3.03 第三）
TOP2_DOC <- c(pooled = "Plasma_cells", balanced = "T_NK_cells")
stopifnot(length(bh_i) == 1L, bh_sc$ct[1] == "Neutrophils", bh_sc$ct[2] == TOP2_DOC[[ASSIGN]], bh_cpm$ct[1] == "Neutrophils", bh_cpm$ct[2] == "Plasma_cells",
          abs(bh_sc$val[1] - bh_sc$val[2] - bh$specificity_log2) < 1e-9, abs(bh_cpm$val[2] - 3.1082) < 1e-3)

# ---------------------------------------------------------------- ④ 作图数据（7 类，按中位 log2FC 降序自上而下）
Y_LAB <- 0.31; Y_PT <- -0.15; SW <- 0.20; Y_MED <- 0.25                                    # 行名 / 点带中心 / 点带半宽 / 中位段半高（行高 = 1）
XL <- c(-2.1, 2.15)
ord <- sm[on_figure == TRUE][order(-median), dominant]
NR  <- length(ord)
TOP_DOC <- c(pooled = "Plasma_cells", balanced = "Neutrophils")                             # 首行：pooled 浆细胞 / balanced 中性粒；末行均为成纤维
stopifnot(NR == 7L, ord[1] == TOP_DOC[[ASSIGN]], ord[NR] == "Fibroblasts", "Neutrophils" %in% ord[1:2], "Plasma_cells" %in% ord[1:2])
rows <- data.table(dominant = ord, yc = NR:1, lab = lab_atac(ord), n_SE = sm[match(ord, dominant), n_SE], med = sm[match(ord, dominant), median])
d <- copy(a[on_figure == TRUE]); d[, yc := rows$yc[match(dominant, rows$dominant)]]                # 行心 y（rows 内的 dominant 与 d 的同名列作用域不同，故用 $ 取）
set.seed(1)
d[, off := vipor::offsetX(log2FC_dominant, width = SW, method = "quasirandom"), by = dominant]
stopifnot(all(abs(d$off) <= SW + 1e-8))
d[, y := yc + Y_PT + off]
bhp <- d[SE_ID == BH_SE]
X_BHLAB <- 1.42
MM_PER_UNIT <- (36 - 1.2) / diff(XL)
n_digits <- nchar(format(rows[dominant == "Neutrophils", n_SE], big.mark = ","))
n_left  <- XL[2] - 0.02 - n_digits * 0.556 * 6 * 25.4 / 72 / MM_PER_UNIT                   # n 列（中性粒行）左端（数据单位）
lab_left <- X_BHLAB - 7.0 * 6 * 25.4 / 72 / MM_PER_UNIT                                     # 标签左端（粗斜体约 7 em）
row_right <- XL[1] + 0.02 + 5.06 * 6 * 25.4 / 72 / MM_PER_UNIT                             # "Neutrophils" 右端（≈ 5.06 em）
stopifnot(X_BHLAB < 1.5, (n_left - X_BHLAB) * MM_PER_UNIT >= 2, (lab_left - row_right) * MM_PER_UNIT >= 2)
cat(sprintf("BHLHE40 SE1 标签右端 %.2f；与 n 列左端（≈ %.2f）间隔 ≈ %.1f mm；标签左端 ≈ %.2f 与行名右端（≈ %.2f）间隔 ≈ %.1f mm\n",
            X_BHLAB, n_left, (n_left - X_BHLAB) * MM_PER_UNIT, lab_left, row_right, (lab_left - row_right) * MM_PER_UNIT))
hdr <- data.table(x = XL + c(0.02, -0.02), y = NR + 0.85, lab = c("Dominant cell type", "n SEs"), hj = c(0, 1))

p <- ggplot() +
  # 每行零线（只穿过点带，不穿过行名）
  geom_segment(data = rows, aes(x = 0, xend = 0, y = yc + Y_PT - SW - 0.06, yend = yc + Y_PT + SW + 0.06),
               linetype = "22", linewidth = lw(0.4), colour = "grey55") +
  geom_point(data = d[SE_ID != BH_SE], aes(log2FC_dominant, y, colour = dominant), size = 0.45, stroke = 0, alpha = 0.7) +
  scale_colour_manual(values = ct_atac_col, guide = "none") +
  geom_segment(data = rows, aes(x = med, xend = med, y = yc + Y_PT - Y_MED, yend = yc + Y_PT + Y_MED),
               linewidth = lw(1.0), colour = "black", lineend = "butt") +
  # BHLHE40 SE1：红色实心点 + 引线 + 粗斜体标签（放在本行行名线上，中性粒名与 n 之间）
  geom_segment(data = bhp, aes(x = log2FC_dominant, xend = log2FC_dominant, y = y + 0.06, yend = yc + Y_LAB - 0.16),
               linewidth = lw(0.4), colour = BH_COL) +
  geom_point(data = bhp, aes(log2FC_dominant, y), shape = 21, fill = BH_COL, colour = "white", size = 1.4, stroke = 0.3) +
  geom_text(data = bhp, aes(X_BHLAB, yc + Y_LAB), label = 'bolditalic("BHLHE40")~"SE1"', parse = TRUE, hjust = 1,
            colour = BH_COL, size = pt2size(FS_SMALL), family = FONT) +
  # 行名（左）与 n（右）、表头
  geom_text(data = rows, aes(XL[1] + 0.02, yc + Y_LAB, label = lab), hjust = 0, size = pt2size(FS_SMALL), family = FONT, colour = "black") +
  geom_text(data = rows, aes(XL[2] - 0.02, yc + Y_LAB, label = format(n_SE, big.mark = ",")), hjust = 1, size = pt2size(FS_SMALL), family = FONT, colour = "grey30") +
  geom_text(data = hdr, aes(x, y, label = lab, hjust = hj), size = pt2size(FS_SMALL), family = FONT, colour = "grey30") +
  scale_x_continuous(limits = XL, breaks = -2:2, expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.45, NR + 1.05), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  labs(x = expression("snATAC " * log[2] * "FC (Mhp / Con)"), y = NULL) +                    # 单行（atop 双行会吃掉约 8 mm 高度；"pseudo-bulk / 主导类型 / 指派方案"写进表头与图注）
  theme_mhp() + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(),
                      plot.margin = margin(12, 2, 1, 1.5))
                                                              #   面板字母 g（8 pt 粗体，占 PNG 顶 0.6–2.6 mm）与左上角表头 "Dominant cell type" 相撞，表头顶须 ≥ 4.3 mm；
                                                              #   行距由约 6.5 降到约 6.0 mm，标签位置为数据坐标随之等比缩放

# ---------------------------------------------------------------- ⑤ 输出
STEM <- "Fig4g_SE感染变化_按细胞类型"
sav(p, STEM, PAN, 36, 56, panel_dir = PNL)
if (identical(Sys.getenv("MHP_FINAL"), "1")) ggsave(file.path(PAN, paste0(STEM, ".pdf")), p, width = 36, height = 56, units = "mm", device = .pdf_dev)
# 读回 PNG 自检：左侧 200 px（表头 / 行名区）从顶往下第一条深色行 = 表头顶，须 ≥ 4.3 mm，给拼版字母 g（8 pt 粗体，占顶 0.6–2.6 mm）留位
if (file.exists(file.path(PAN, paste0(STEM, ".png")))) {
  hdr_top_mm <- local({ g <- png::readPNG(file.path(PAN, paste0(STEM, ".png"))); gr <- (g[, , 1] + g[, , 2] + g[, , 3]) / 3
                        min(which(apply(gr[, 1:200], 1, min) < 0.6)) / 300 * 25.4 })
  cat(sprintf("PNG 顶部表头顶 = %.2f mm（须 ≥ 4.3）\n", hdr_top_mm)); stopifnot(hdr_top_mm >= 4.3)
}

# Source Data：2,855 行；主图列（dominant / log2FC_dominant …）= ASSIGN 方案；四方案的主导类型与 log2FC 各自成列
src_cols <- c("SE_ID", "gene", "src", "nPeak", "width", ct_atac_order)
alt_cols <- as.vector(rbind(paste0("dominant_", SCHEMES), paste0("log2FC_", SCHEMES)))
out <- a[, c(src_cols, "assignment", "dominant", "specificity_log2", "CPM_dominant_all", "CPM_dominant_Con", "CPM_dominant_Mhp",
             "log2FC_dominant", "log2FC_dominant_pc1", "rank_within_dominant", "on_figure", "is_BHLHE40_SE1", alt_cols, "consensus_Con_Mhp"), with = FALSE]
out[, dominant_label := lab_atac(dominant)]
rnd <- c("specificity_log2", "log2FC_dominant", "log2FC_dominant_pc1", paste0("log2FC_", SCHEMES))
out[, (rnd) := lapply(.SD, round, 4), .SDcols = rnd]
setcolorder(out, c(src_cols, "assignment", "dominant", "dominant_label"))
fwrite(out[order(dominant, rank_within_dominant)], file.path(TAB, "08_g_SE主导类型与感染log2FC.csv"))

sm_out <- copy(sm_all)[order(-is_main, match(assignment, SCHEMES), -median)]
sm_out[, order_on_figure := ifelse(is_main & on_figure, match(dominant, ord), NA_integer_)]
num_cols <- c("median", "mean", "q1", "q3", "min", "max", "pct_up", "pct_down", "median_pc1", "diff_vs_documented", "nuclei_ratio", "spearman_median_vs_nuclei_ratio")
sm_out[, (num_cols) := lapply(.SD, round, 4), .SDcols = num_cols]
fwrite(sm_out[, .(assignment, is_main, dominant, label, on_figure, order_on_figure, n_SE, nuclei_Con, nuclei_Mhp, nuclei_ratio, median_log2FC = median,
                  documented_median, diff_vs_documented, mean_log2FC = mean, q1, q3, min, max, pct_up, pct_down, median_log2FC_pseudocount1 = median_pc1,
                  spearman_median_vs_nuclei_ratio)],
       file.path(TAB, "08_g_各细胞类型汇总.csv"))

fmt3 <- function(x) sprintf("%+.3f", x)
cap_low <- function(x) ifelse(grepl("^T/NK", x), x, paste0(tolower(substr(x, 1, 1)), substring(x, 2)))   # 图注句内类型名小写首字母（T/NK 保留）
rho_of <- function(s) sm_all[assignment == s & on_figure == TRUE, spearman_median_vs_nuclei_ratio][1]
per_type <- function(s) { x <- sm_all[assignment == s & on_figure == TRUE][order(-median)]
  paste(sprintf("%s %s (n = %d, %.1f%% up)", x$label, fmt3(x$median), x$n_SE, x$pct_up), collapse = "; ") }
SCHEME_DEF <- c(pooled = "log2(pooled Con + Mhp CPM + 1)", Con_only = "log2(Con CPM + 1)", Mhp_only = "log2(Mhp CPM + 1)",
                balanced = "(log2(Con CPM + 1) + log2(Mhp CPM + 1)) / 2 (equal weight of the two groups)")
key <- data.table(
  item = c("SE 集合", "主图指派方案", "主导类型定义（各方案的得分；8 类中最高者，无并列）", "log2FC 定义（记录口径，四方案相同）", "伪计数（CPM）",
           "上图类型数 / SE 数", "不上图：B 细胞主导 SE 数（Mhp 核数 / Con 核数）",
           "全局中位 log2FC（2,855 个 SE，含 B）", sprintf("全局中位 log2FC（上图 %s 个，不含 B）", format(sum(sm[on_figure == TRUE, n_SE]), big.mark = ",")),
           "各类中位 log2FC（图序，自上而下；n；升高比例）",
           "各类中位 log2FC（Source Data，不上图）", "各类中位 log2FC（Source Data，不上图）", "各类中位 log2FC（Source Data，不上图）",
           "Spearman（7 类中位数 vs 核数比 Mhp/Con）",
           "记录值核对",
           "与 15_g 独立复算的交叉核对",
           "整合方案【探针】+0.28 / +0.21 / +0.02 的来源",
           "BHLHE40 SE1 编号 / 坐标（fixed 表 = 两组 ROSE SE1 并集）", "BHLHE40 SE1 主导类型 / 特异性（log2）",
           "BHLHE40 SE1 中性粒伪 bulk CPM Con → Mhp", "BHLHE40 SE1 log2FC（中性粒内名次；四方案）",
           "BHLHE40 SE1 巨噬伪 bulk CPM Con → Mhp（对照，不上图）",
           "浆细胞核数 Con / Mhp（图注注明 Con 组浆细胞核少）",
           "图注句（建议，英文）"),
  assignment = c("-", ASSIGN, "-", "-", "-", ASSIGN, ASSIGN, ASSIGN, ASSIGN, ASSIGN, setdiff(SCHEMES, ASSIGN), "全部方案", "pooled", "balanced", "pooled",
                 "-", "全部方案", ASSIGN, "全部方案", "-", ASSIGN, ASSIGN),
  value = c("2,855 non-redundant SE intervals = Con 2,307 ∪ Mhp 2,422 ROSE SEs reduced (Both 1,555 / Mhp_only 679 / Con_only 621); E2_BHLHE40 sub-interval excluded",
            sprintf("main figure = %s: dominant cell type = argmax over 8 types of %s; pooled / Con-only / Mhp-only assignments go to Source Data and the supplementary sketch (15_g). This run: ASSIGN = %s%s",
                    MAIN, SCHEME_DEF[[MAIN]], ASSIGN, if (ASSIGN == MAIN) "" else "（备选重现，输出在 备选_指派方案/ 下，不进主图）"),
            paste(sprintf("%s = %s", SCHEMES, SCHEME_DEF[SCHEMES]), collapse = "; "),
            "log2((CPM_Mhp + 0.1) / (CPM_Con + 0.1)) within the dominant cell type's pseudo-bulk", PC,
            sprintf("%d / %s", NR, format(sum(sm[on_figure == TRUE, n_SE]), big.mark = ",")),
            sprintf("%d (%d / %d)", sm[dominant == EXCL, n_SE], sm[dominant == EXCL, nuclei_Mhp], sm[dominant == EXCL, nuclei_Con]),
            fmt3(GLOBAL), fmt3(median(d$log2FC_dominant)),
            per_type(ASSIGN),
            per_type(setdiff(SCHEMES, ASSIGN)[1]), per_type(setdiff(SCHEMES, ASSIGN)[2]), per_type(setdiff(SCHEMES, ASSIGN)[3]),
            paste(sprintf("%s %.3f", SCHEMES, sapply(SCHEMES, rho_of)), collapse = "; "),
            paste(sprintf("%s %s vs %s", smp[order(-median), dominant], fmt3(smp[order(-median), median]), fmt3(smp[order(-median), documented_median])), collapse = "; "),
            sens_note,
            sprintf("pseudo-count 1 (log2(CPM + 1)) instead of 0.1, pooled assignment: plasma %s / neutrophil %s / macrophage %s",
                    fmt3(smp[dominant == "Plasma_cells", median_pc1]), fmt3(smp[dominant == "Neutrophils", median_pc1]), fmt3(smp[dominant == "Macrophages", median_pc1])),
            sprintf("%s / %s (width %s bp)", BH_SE, BH_COORD, format(bh$width, big.mark = ",")),
            sprintf("Neutrophils under all four schemes / %.3f (%s score %s %.3f vs %s %.3f; pooled CPM %s %.2f vs %s %.2f)", bh$specificity_log2, ASSIGN,
                    lab_atac(bh_sc$ct[1]), bh_sc$val[1], lab_atac(bh_sc$ct[2]), bh_sc$val[2], lab_atac(bh_cpm$ct[1]), bh_cpm$val[1], lab_atac(bh_cpm$ct[2]), bh_cpm$val[2]),
            sprintf("%.2f → %.2f (%+.0f%%)", bh$CPM_dominant_Con, bh$CPM_dominant_Mhp, 100 * (bh$CPM_dominant_Mhp / bh$CPM_dominant_Con - 1)),
            sprintf("%+.3f; rank %s", bh$log2FC_dominant, paste(sprintf("%s %d of %d", SCHEMES, bh_rank[SCHEMES], N_NEUT[SCHEMES]), collapse = " / ")),
            sprintf("%.2f → %.2f (%+.0f%%)", g[SE_ID == BH_SE, Mock_Macrophages], g[SE_ID == BH_SE, Mhp_Macrophages], 100 * (g[SE_ID == BH_SE, Mhp_Macrophages] / g[SE_ID == BH_SE, Mock_Macrophages] - 1)),
            sprintf("%d / %d", sm[dominant == "Plasma_cells", nuclei_Con], sm[dominant == "Plasma_cells", nuclei_Mhp]),
            sprintf(paste0("g, Infection-induced change in super-enhancer accessibility by dominant cell type (snATAC-seq). Each of the 2,855 non-redundant SE intervals (Con ∪ Mhp) ",
                           "was assigned to the cell type with the highest pseudo-bulk accessibility score among eight cell types%s; log2(Mhp/Con) was then computed within that cell type's pseudo-bulk ",
                           "(pseudo-count 0.1 CPM). Points, individual SEs; black bars, median; n, SEs per cell type (median: %s; genome-wide %s). B-cell-dominant SEs (n = %d) are omitted because the Mhp group ",
                           "contained only %d B-cell nuclei. The BHLHE40 SE1 (%s; neutrophil-dominant) is highlighted in red (log2FC %+.2f). Assignment on the pooled, Con-only or Mhp-only pseudo-bulk gives ",
                           "the medians listed in Source Data."),
                    if (ASSIGN == "balanced") " (mean of log2(CPM + 1) in the Con and Mhp pseudo-bulks, i.e. equal weight of the two groups)" else " (log2(CPM + 1) of the pooled Con + Mhp pseudo-bulk)",
                    paste(sprintf("%s %s", cap_low(rows$lab), fmt3(rows$med)), collapse = ", "), fmt3(GLOBAL),
                    sm[dominant == EXCL, n_SE], sm[dominant == EXCL, nuclei_Mhp], BH_COORD, bh$log2FC_dominant)))
stopifnot(nrow(key) == 24L)
fwrite(key, file.path(TAB, "08_g_关键数字.csv"))
stopifnot(file.size(file.path(TAB, "08_g_SE主导类型与感染log2FC.csv")) > 0, file.size(file.path(TAB, "08_g_各细胞类型汇总.csv")) > 0,
          file.size(file.path(TAB, "08_g_关键数字.csv")) > 0)

# 关键数字回显
print(sm_out[is_main == TRUE, .(dominant, on_figure, n_SE, nuclei_Con, nuclei_Mhp, median_log2FC = median, pct_up, median_log2FC_pseudocount1 = median_pc1)])
cat(sprintf("[%s] 全局中位 %+.4f（含 B）/ %+.4f（不含 B）；Spearman(中位数, 核数比) %s；BHLHE40 SE1 %s log2FC %+.3f，中性粒内第 %d/%d\n",
            ASSIGN, GLOBAL, median(d$log2FC_dominant), paste(sprintf("%s %.3f", SCHEMES, sapply(SCHEMES, rho_of)), collapse = " / "),
            BH_COORD, bh$log2FC_dominant, bh$rank_within_dominant, N_NEUT[[ASSIGN]]))
