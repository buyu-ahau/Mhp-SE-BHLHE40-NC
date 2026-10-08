# ============================================================================
# 04_新Fig6_CCK8_本地.R — Fig. 6d：六组 CCK-8 细胞活力
# 运行位置：本地；路径变量见 config/paths.R
# ============================================================================

source(file.path(Sys.getenv("MHP_CODE", "."), "config", "paths.R"))

suppressPackageStartupMessages({ library(readxl) })
source(file.path(MHP_CODE, "06_Fig6_E2_CRISPRi", "_fig6_common_本地.R"))
RAW_ORDER <- names(RAW2KEY)

x <- as.data.frame(suppressMessages(read_excel(RAW6, sheet = "CCK8", col_names = FALSE, .name_repair = "minimal")))
stopifnot(x[7, 2] == "ELx800", x[3, 2] == "板2", identical(unname(unlist(x[18, 4:9])), LETTERS[1:6]),
          identical(as.numeric(unlist(x[20, 3:14])), as.numeric(1:12)), identical(unname(unlist(x[21:28, 2])), LETTERS[1:8]))
od <- apply(x[21:28, 3:14], 2, as.numeric); dimnames(od) <- list(LETTERS[1:8], 1:12)   # 8 × 12 板面 OD
blank <- od[LETTERS[2:7], "1"]; b0 <- mean(blank)
D <- rbindlist(lapply(1:6, function(j) data.table(plate_col = j + 1, col_letter = LETTERS[j], raw_group = RAW_ORDER[j], plate_row = LETTERS[2:7],
                                                  OD_raw = od[LETTERS[2:7], as.character(j + 1)], blank_same_row = blank)))
D[, group := to_key6(raw_group)][, `:=`(OD_minus_blank = OD_raw - b0)]
chkE <- D[, .(m = mean(OD_raw)), by = group][order(group)]$m
stopifnot(nrow(D) == 36, !anyNA(D$OD_raw), max(abs(chkE - c(0.707, 0.205, 0.729, 0.541, 0.240, 1.280))) < 0.001)   # 现图 E（改图方案 §2）
cat(sprintf("空白孔（第 1 列 B–G）OD = %s，均值 %.4f；第 8 列 = %s；第 9–12 列均值 %.3f\n", paste(blank, collapse = "/"), b0,
            paste(od[LETTERS[2:7], "8"], collapse = "/"), mean(od[LETTERS[2:7], as.character(9:12)])))

# ---------------------------------------------------------------- 统计与作图（扣空白）
d <- D[, .(group, value = OD_minus_blank)]
ST <- tukey_all(d, GRP6)[, `:=`(panel = "h", marker = "CCK-8 OD450 − blank")]
ST_raw <- tukey_all(D[, .(group, value = OD_raw)], GRP6)[, `:=`(panel = "h_raw_OD", marker = "CCK-8 raw OD450 (not plotted)")]
shown_key <- unlist(lapply(CMP6, function(a) c(paste(a, collapse = " "), paste(rev(a), collapse = " "))))
STall <- rbind(ST, ST_raw)[, shown := panel == "h" & paste(g1, g2) %in% shown_key][, p_adj_txt := fmtp(p_adj)]
ph <- pm_table(bar_open(d, pick_cmp(ST, CMP6, CMP6_H), "CCK-8 (OD450 \u2212 blank)", GRP6, COL6, COL6_PT), SIGN6)
sav6(ph, "Fig6h_CCK8")

# ---------------------------------------------------------------- 表
fwrite(D[order(group, plate_row), .(raw_group, group, col_letter, plate_col, plate_row, OD_raw, blank_mean = round(b0, 4), OD_minus_blank_plot = round(OD_minus_blank, 4))],
       file.path(TAB, "Fig6h_CCK8_逐孔.csv"), bom = TRUE)
sm <- D[, .(n = .N, n_unit = "well (single plate)", mean_raw_OD = round(mean(OD_raw), 4), mean_plot = round(mean(OD_minus_blank), 4), sd_plot = round(sd(OD_minus_blank), 4)), by = group]
sm[, oldfig_raw_OD := c(0.707, 0.205, 0.729, 0.541, 0.240, 1.280)[match(group, GRP6)]][, diff_vs_oldfig := round(mean_raw_OD - oldfig_raw_OD, 4)]
sm[, ratio_to_Ctrl := round(mean_plot / mean_plot[group == "Ctrl"], 3)][, ratio_to_CtrlMhp := round(mean_plot / mean_plot[group == "Ctrl+Mhp"], 3)]
fwrite(sm[order(match(group, GRP6))], file.path(TAB, "Fig6h_CCK8_组均值_与现图核对.csv"), bom = TRUE)
fwrite(STall[, .(panel, marker, g1, g2, diff = round(diff, 4), lwr = round(lwr, 4), upr = round(upr, 4), p_adj = p_adj_txt, star, shown_in_figure = shown)],
       file.path(TAB, "Fig6h_CCK8_统计_Tukey全部两两.csv"), bom = TRUE)
fwrite(rbind(cbind(panel = "h", attr(ST, "anova")), cbind(panel = "h_raw_OD", attr(ST_raw, "anova"))), file.path(TAB, "Fig6h_CCK8_ANOVA.csv"), bom = TRUE)
for (f in c("Fig6h_CCK8_逐孔.csv", "Fig6h_CCK8_组均值_与现图核对.csv", "Fig6h_CCK8_统计_Tukey全部两两.csv", "Fig6h_CCK8_ANOVA.csv")) chk_file(file.path(TAB, f))
print(sm[order(match(group, GRP6))]); print(STall[shown == TRUE, .(g1, g2, p_adj = p_adj_txt, star)])
