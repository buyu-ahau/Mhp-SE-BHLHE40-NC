# ============================================================================
# 05_新Fig6_凋亡率_从FCS重算_本地.R — Fig. 6e：由原始 .fcs 文件按四象限门重算 Annexin V⁺ 凋亡率
# 运行位置：本地；路径变量见 config/paths.R
# ============================================================================

source(file.path(Sys.getenv("MHP_CODE", "."), "config", "paths.R"))

suppressPackageStartupMessages({ library(flowCore); library(readxl) })
source(file.path(MHP_CODE, "06_Fig6_E2_CRISPRi", "_fig6_common_本地.R"))

# ---------------------------------------------------------------- 每管四象限
TB <- list_tubes()
TB <- TB[, {
  r <- read_tube(file); q <- quad4(r$dt)
  .(n_total = r$n_total, n_P1 = r$n_P1, date = r$date, btim = r$btim, UL = q[["UL"]], UR = q[["UR"]], LL = q[["LL"]], LR = q[["LR"]])
}, by = .(file, raw_group, group, replicate)]
TB[, annexinV_pos := UR + LR]
stopifnot(nrow(TB) == 18, all(TB$n_P1 == 10000), abs(rowSums(TB[, .(UL, UR, LL, LR)]) - 100) < 1e-9)

# 核对 ①：现图 6 幅代表图的四象限（原图标注值）
IMG <- data.table(file = c("WT-1.fcs", "SE-KO-1.fcs", "SE-KO+OE-1.fcs", "WT+MHP-1.fcs", "SE-KO+MHP-1.fcs", "SE-KO+OE+MHP-3.fcs"),
                  UL_img = c(3.22, 1.81, 3.61, 1.82, 0.87, 2.81), UR_img = c(7.47, 12.39, 7.65, 11.83, 5.61, 9.49),
                  LL_img = c(79.36, 72.71, 78.92, 73.66, 64.39, 74.38), LR_img = c(9.95, 13.09, 9.82, 12.69, 29.13, 13.32))
TB <- merge(TB, IMG, by = "file", all.x = TRUE)
TB[, diff_img_sum := abs(UL - UL_img) + abs(UR - UR_img) + abs(LL - LL_img) + abs(LR - LR_img)]
stopifnot(max(TB$diff_img_sum, na.rm = TRUE) < 0.5)
a <- as.data.frame(suppressMessages(read_excel(APO_XLS, col_names = FALSE, .name_repair = "minimal")))
# 列位（探针实测）：1 象限名（Q1-UR 行 2、Q1-LR 行 5）；2 re1–re3；3–8 六组
stopifnot(identical(unname(unlist(a[1, 3:8])), c("WT", "KE-KO", "SE-KO+OE", "WT+MHP", "SE-KO+MHP", "SE-KO+OE+MHP")), a[2, 1] == "Q1-UR", a[5, 1] == "Q1-LR")
hand <- rbindlist(lapply(1:6, function(j) data.table(raw_group = names(RAW2KEY)[j], row_in_sheet = 1:3, UR_hand = 100 * as.numeric(a[2:4, j + 2]), LR_hand = 100 * as.numeric(a[5:7, j + 2]))))
hand[, annexinV_hand := UR_hand + LR_hand]
TB[, nearest_hand_row := sapply(seq_len(.N), function(i) { h <- hand[raw_group == TB$raw_group[i]]; h$row_in_sheet[which.min(abs(h$annexinV_hand - TB$annexinV_pos[i]))] })]
TB <- merge(TB, hand[, .(raw_group, nearest_hand_row = row_in_sheet, UR_hand, LR_hand, annexinV_hand)], by = c("raw_group", "nearest_hand_row"), all.x = TRUE)
TB[, diff_vs_hand := annexinV_pos - annexinV_hand]
TB <- TB[order(group, replicate)]

# ---------------------------------------------------------------- 统计与作图
d <- TB[, .(group, value = annexinV_pos)]
ST <- tukey_all(d, GRP6)[, `:=`(panel = "i", marker = "Annexin V+ (UR + LR, % of P1)")]
shown_key <- unlist(lapply(CMP6, function(a) c(paste(a, collapse = " "), paste(rev(a), collapse = " "))))
ST[, shown := paste(g1, g2) %in% shown_key][, p_adj_txt := fmtp(p_adj)]
pi <- pm_table(bar_open(d, pick_cmp(ST, CMP6, CMP6_H), "Annexin V<sup><span style='font-size:5pt'>+</span></sup> cells (%)", GRP6, COL6, COL6_PT), SIGN6) +
  theme(axis.title.y = ggtext::element_markdown(size = FS_BASE, family = FONT, colour = "black"))
sav6(pi, "Fig6i_凋亡率_AnnexinV")

# ---------------------------------------------------------------- 表
r2 <- function(v) round(v, 2)
fwrite(TB[, .(raw_group, group, replicate, file, date_in_FCS = date, btim, n_total, n_P1, UL = r2(UL), UR = r2(UR), LL = r2(LL), LR = r2(LR), annexinV_pos_plot = r2(annexinV_pos),
              UL_img, UR_img, LL_img, LR_img, diff_img_sum = r2(diff_img_sum), nearest_hand_row, UR_hand, LR_hand, annexinV_hand, diff_vs_hand = r2(diff_vs_hand))],
       file.path(TAB, "Fig6i_流式凋亡_每管四象限.csv"), bom = TRUE)
sm <- TB[, .(n = .N, n_unit = "tube (single acquisition)", mean_plot = r2(mean(annexinV_pos)), sd_plot = r2(sd(annexinV_pos)), mean_hand_table = r2(mean(annexinV_hand))), by = group]
sm[, oldfig := c(18.31, 26.43, 19.48, 24.89, 32.29, 20.78)[match(group, GRP6)]][, diff_vs_oldfig := r2(mean_plot - oldfig)]
fwrite(sm[order(match(group, GRP6))], file.path(TAB, "Fig6i_凋亡率_组均值_与现图核对.csv"), bom = TRUE)
fwrite(ST[, .(panel, marker, g1, g2, diff = round(diff, 4), lwr = round(lwr, 4), upr = round(upr, 4), p_adj = p_adj_txt, star, shown_in_figure = shown)],
       file.path(TAB, "Fig6i_凋亡率_统计_Tukey全部两两.csv"), bom = TRUE)
fwrite(cbind(panel = "i", attr(ST, "anova")), file.path(TAB, "Fig6i_凋亡率_ANOVA.csv"), bom = TRUE)
for (f in c("Fig6i_流式凋亡_每管四象限.csv", "Fig6i_凋亡率_组均值_与现图核对.csv", "Fig6i_凋亡率_统计_Tukey全部两两.csv", "Fig6i_凋亡率_ANOVA.csv")) chk_file(file.path(TAB, f))
print(TB[, .(group, replicate, UL = r2(UL), UR = r2(UR), LL = r2(LL), LR = r2(LR), annexinV = r2(annexinV_pos), diff_img_sum = r2(diff_img_sum), diff_vs_hand = r2(diff_vs_hand))])
print(sm[order(match(group, GRP6))]); print(ST[shown == TRUE, .(g1, g2, p_adj = p_adj_txt, star)])
cat(sprintf("代表图复现：四象限绝对差之和最大 %.2f 个百分点；与记录表最大差 %.2f\n", max(TB$diff_img_sum, na.rm = TRUE), max(abs(TB$diff_vs_hand))))
