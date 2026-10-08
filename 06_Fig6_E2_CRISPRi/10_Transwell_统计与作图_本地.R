# ============================================================================
# 10_Transwell_统计与作图_本地.R — Fig. 6f：六组 Transwell 染色面积的单因素方差分析（Tukey 多重比较）与柱图
# 运行位置：本地；路径变量见 config/paths.R
# ============================================================================

source(file.path(Sys.getenv("MHP_CODE", "."), "config", "paths.R"))

source(file.path(MHP_CODE, "06_Fig6_E2_CRISPRi", "_fig6_common_本地.R"))
QC <- file.path(OUT, "Plots/QC_Transwell"); dir.create(QC, recursive = TRUE, showWarnings = FALSE)
RES  <- fread(file.path(TAB, "Fig6_Transwell_逐图定量.csv"), encoding = "UTF-8")
SENS <- fread(file.path(TAB, "Fig6_Transwell_敏感性网格_逐图.csv"), encoding = "UTF-8")
stopifnot(nrow(RES) == 18, all(RES$group %in% GRP6), RES[, .N, by = group][, all(N == 3)], nrow(SENS) == 18 * 11)
RES[, group := factor(group, levels = GRP6)]; SENS[, group := factor(group, levels = GRP6)]

# ---------------------------------------------------------------- 主口径：ANOVA + Tukey + 柱图
d <- RES[, .(group, value = pct_area)]
ST <- tukey_all(d, GRP6)[, `:=`(panel = "Transwell", marker = "crystal violet-stained area (% of field), pores removed")]
shown_key <- unlist(lapply(CMP6, function(a) c(paste(a, collapse = " "), paste(rev(a), collapse = " "))))
ST[, shown := paste(g1, g2) %in% shown_key][, p_adj_txt := fmtp(p_adj)]
p_main <- pm_table(bar_open(d, pick_cmp(ST, CMP6, CMP6_H), "Migrated cells (% area)", GRP6, COL6, COL6_PT), SIGN6)
sav6(p_main, "Fig6_Transwell_迁移面积")
# 次要指标：深染核对象数 / mm²（相触核不拆分，只用于核对方向）
dn <- RES[, .(group, value = nuclear_obj_per_mm2)]
STn <- tukey_all(dn, GRP6)[, `:=`(panel = "Transwell_nuclear_objects", marker = "dark nuclear objects per mm2 (touching nuclei not split)")][, shown := FALSE][, p_adj_txt := fmtp(p_adj)]
p_nuc <- pm_table(bar_open(dn, pick_cmp(STn, CMP6, CMP6_H), expression(paste("Nuclear objects (", mm^-2, ")")), GRP6, COL6, COL6_PT), SIGN6)
sav6(p_nuc, "Fig6_Transwell_核对象数")

# ---------------------------------------------------------------- 组均值与比值
sm <- RES[, .(n = .N, n_unit = "independent well (1 field/well)", mean_pct_area = mean(pct_area), sd_pct_area = sd(pct_area),
              values = paste(round(pct_area, 2), collapse = " / "), mean_nuclear_obj_per_mm2 = mean(nuclear_obj_per_mm2), sd_nuclear_obj_per_mm2 = sd(nuclear_obj_per_mm2),
              mean_imagej_like = mean(pct_imagej_like)), by = group][order(group)]
m_of <- function(g) sm[group == g, mean_pct_area]
sm[, ratio_to_Ctrl := mean_pct_area / m_of("Ctrl")][, ratio_to_CtrlMhp := mean_pct_area / m_of("Ctrl+Mhp")]
sm[, note := "ratios depend on quantification variant (see Fig6_Transwell_倍数范围.csv); do not cite 2-decimal fold changes in text/legend"]
fwrite(sm[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)], file.path(TAB, "Fig6_Transwell_组均值.csv"), bom = TRUE)
fwrite(rbind(ST, STn)[, .(panel, marker, g1, g2, diff = round(diff, 3), lwr = round(lwr, 3), upr = round(upr, 3), p_adj = p_adj_txt, star, shown_in_figure = shown)],
       file.path(TAB, "Fig6_Transwell_统计_Tukey全部两两.csv"), bom = TRUE)
fwrite(rbind(cbind(panel = "Transwell", attr(ST, "anova")), cbind(panel = "Transwell_nuclear_objects", attr(STn, "anova"))), file.path(TAB, "Fig6_Transwell_ANOVA.csv"), bom = TRUE)

# ---------------------------------------------------------------- 敏感性：11 个变体 × （4 对关键比较 + 1 对未作图比较）+ 去碎片核对
CMPX <- c(CMP6, list(c("CRISPRi", "CRISPRi+OE"))); CMPX_H <- c(CMP6_H, 1); CMPX_LAB <- c(CMP6_LAB, "OE rescue, uninfected (not shown in figure)")
VAR <- unique(SENS[, .(variant, type, k_s, pore_max)])
FRAG_V <- "main_minus_objects_lt200px (check, not one of the 11 variants)"
sk_one <- function(dv, vrow) { st <- tukey_all(dv, GRP6); an <- attr(st, "anova")
  cbind(vrow, pick_cmp(st, CMPX, CMPX_H)[, .(comparison = CMPX_LAB, g1, g2, diff = round(diff, 3), p_adj, star, sig = p_adj < 0.05)], F_anova = round(an$F, 2), p_anova = an$p_anova) }
SK <- rbindlist(c(lapply(seq_len(nrow(VAR)), function(i) sk_one(SENS[variant == VAR$variant[i], .(group, value = pct_area)], VAR[i])),
                  list(sk_one(RES[, .(group, value = pct_area_minus_lt200px)], data.table(variant = FRAG_V, type = "fragments_removed_check", k_s = 0.75, pore_max = 400L)))))
SK[, is_main := variant == sprintf("k_s=%.2f|pore_max=%d", 0.75, 400)]
MAIN <- SK[is_main == TRUE, .(comparison, main_sign = sign(diff), main_sig = sig)]
FRG  <- SK[type == "fragments_removed_check", .(comparison, fragments_removed_p = fmtp(p_adj), fragments_removed_star = star)]
VER <- merge(SK[type != "fragments_removed_check"], MAIN, by = "comparison")[, .(n_variants = .N, n_same_direction = sum(sign(diff) == main_sign), n_significant = sum(sig), n_same_significance = sum(sig == main_sig),
                                             p_range = sprintf("%s – %s", fmtp(min(p_adj)), fmtp(max(p_adj))),
                                             p_range_significant_variants = if (any(sig)) sprintf("%s – %s (%d/%d significant)", fmtp(min(p_adj[sig])), fmtp(max(p_adj[sig])), sum(sig), .N) else sprintf("none significant (0/%d)", .N),
                                             p_nonsignificant_variants = if (any(!sig)) paste(sprintf("%s: %s", variant[!sig], fmtp(p_adj[!sig])), collapse = "; ") else "",
                                             main_diff = diff[is_main], main_p = fmtp(p_adj[is_main]), main_star = star[is_main]), by = .(comparison, g1, g2)]
VER[, verdict := fifelse(n_same_direction == n_variants & n_same_significance == n_variants, "stable: same direction and same significance in all 11 variants",
                         fifelse(n_same_direction == n_variants, "direction stable; significance changes in some variants", "direction changes in some variants"))]
VER <- merge(VER, FRG, by = "comparison")[match(CMPX_LAB, comparison)]
fwrite(SK[, .(variant, type, k_s, pore_max, is_main, comparison, shown_in_figure = comparison %in% CMP6_LAB, g1, g2, diff, p_adj = fmtp(p_adj), star, significant = sig, F_anova, p_anova = fmtp(p_anova))],
       file.path(TAB, "Fig6_Transwell_敏感性_关键比较.csv"), bom = TRUE)
fwrite(VER[, c(.SD, list(shown_in_figure = comparison %in% CMP6_LAB))], file.path(TAB, "Fig6_Transwell_敏感性_判定.csv"), bom = TRUE)

# ---------------------------------------------------------------- 倍数范围：11 个变体 + 去碎片核对（正文 / 图注只引用范围）
RAT <- list(`CRISPRi / Ctrl` = c("CRISPRi", "Ctrl"), `Ctrl+Mhp / Ctrl` = c("Ctrl+Mhp", "Ctrl"), `CRISPRi+Mhp / Ctrl+Mhp` = c("CRISPRi+Mhp", "Ctrl+Mhp"),
            `CRISPRi+OE+Mhp / CRISPRi+Mhp` = c("CRISPRi+OE+Mhp", "CRISPRi+Mhp"), `CRISPRi+OE / CRISPRi` = c("CRISPRi+OE", "CRISPRi"))
gm <- rbind(SENS[, .(mean = mean(pct_area)), by = .(variant, type, group)], RES[, .(variant = FRAG_V, type = "fragments_removed_check", mean = mean(pct_area_minus_lt200px)), by = group])
FOLD <- gm[, c(list(is_main = variant[1] == "k_s=0.75|pore_max=400"), lapply(RAT, function(r) mean[group == r[1]] / mean[group == r[2]])), by = .(variant, type)]
FOLD_RANGE <- data.table(ratio = names(RAT), main = sapply(names(RAT), function(k) FOLD[is_main == TRUE][[k]]),
                         min_11_variants = sapply(names(RAT), function(k) min(FOLD[type != "fragments_removed_check"][[k]])), max_11_variants = sapply(names(RAT), function(k) max(FOLD[type != "fragments_removed_check"][[k]])),
                         fragments_removed = sapply(names(RAT), function(k) FOLD[type == "fragments_removed_check"][[k]]))
FOLD_RANGE[, `:=`(min_all = pmin(min_11_variants, fragments_removed), max_all = pmax(max_11_variants, fragments_removed))]
FOLD_RANGE[, `:=`(min_pore_removed_9 = sapply(ratio, function(k) min(FOLD[type == "grid"][[k]])), max_pore_removed_9 = sapply(ratio, function(k) max(FOLD[type == "grid"][[k]])),
                  min_no_pore_removal_2 = sapply(ratio, function(k) min(FOLD[type %in% c("no_pore", "baseline")][[k]])), max_no_pore_removal_2 = sapply(ratio, function(k) max(FOLD[type %in% c("no_pore", "baseline")][[k]])))]
FOLD_RANGE[, independent_verification := c("0.35-0.42 (OD deconvolution / Lab a*, 10 thresholds; independent check on first-version data)", "", "",
                                           "2.6-5.1 (OD deconvolution / Lab a*, 10 thresholds; independent check on first-version data)", "")]
FOLD_RANGE[, text_wording := sprintf("%.2f-%.2f overall (pore-removed grid %.2f-%.2f; no pore removal %.2f-%.2f; fragment-removed %.2f); report group mean ± s.d. or a range, never the main value %.2f alone",
                                      min_all, max_all, min_pore_removed_9, max_pore_removed_9, min_no_pore_removal_2, max_no_pore_removal_2, fragments_removed, main)]
fwrite(rbind(FOLD[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)][, row := "per_variant"], fill = TRUE,
             FOLD_RANGE[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)][, row := "range_summary"]), file.path(TAB, "Fig6_Transwell_倍数范围.csv"), bom = TRUE)
# QC 图：各变体的组均值 ± SD
SM <- SENS[, .(mean = mean(pct_area), sd = sd(pct_area)), by = .(variant, type, group)]
SM[, variant_short := fifelse(type == "baseline", "ImageJ-like (gray<=162)", fifelse(type == "no_pore", "colour rule, no pore removal", sub("k_s=", "k_s ", sub("\\|pore_max=", ", pore<=", variant))))]
SM[, is_main := variant == "k_s=0.75|pore_max=400"]
pq <- ggplot(SM, aes(group, mean, colour = variant_short, group = variant_short)) + geom_line(aes(linewidth = is_main)) + geom_point(size = 1) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.15, linewidth = lw(0.4)) + scale_linewidth_manual(values = c(`FALSE` = lw(0.5), `TRUE` = lw(1.4)), guide = "none") +
  scale_colour_manual(values = c(scales::hue_pal()(9), "black", "grey50"), name = NULL) + labs(x = NULL, y = "Stained area (% of field), mean ± s.d. (n = 3 wells)", title = "Sensitivity grid: group means by variant (thick = main)") +
  theme_mhp() + theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "right", legend.key.size = unit(2.5, "mm"))
ggsave(file.path(QC, "QC_敏感性_组均值.png"), pq, width = 150, height = 80, units = "mm", dpi = 300, bg = "white", device = ragg::agg_png)

# ---------------------------------------------------------------- 关键数字 + 图注建议句
k4 <- ST[shown == TRUE]; k4 <- k4[match(sapply(CMP6, function(a) paste(sort(a), collapse = "|")), paste(pmin(g1, g2), pmax(g1, g2), sep = "|"))]
kn <- STn[(g1 == "Ctrl+Mhp" & g2 == "CRISPRi+Mhp") | (g1 == "CRISPRi+Mhp" & g2 == "Ctrl+Mhp")]; stopifnot(nrow(kn) == 1)
KN <- cbind(pick_cmp(STn, CMP6, CMP6_H)[, .(g1, g2, p_nuc = p_adj, star_nuc = star)], pick_cmp(ST, CMP6, CMP6_H)[, .(p_area = p_adj, star_area = star)])
KN[, differs := star_nuc != star_area]
kn_txt  <- KN[, paste(sprintf("%s vs %s P = %s (%s)", g1, g2, fmtp(p_nuc), star_nuc), collapse = "; ")]
kn_diff <- KN[differs == TRUE, paste(sprintf("%s vs %s (nuclear P = %s, %s; area P = %s, %s)", g1, g2, fmtp(p_nuc), star_nuc, fmtp(p_area), star_area), collapse = "; ")]
kn_ns   <- KN[differs == TRUE & star_nuc == "ns", paste(sprintf("%s vs %s, P = %s", g1, g2, fmtp(p_nuc)), collapse = "; ")]
kx <- VER[comparison == CMPX_LAB[5]]; k1 <- VER[comparison == CMP6_LAB[4]]
fr <- function(k) FOLD_RANGE[ratio == k, sprintf("%.2f–%.2f", min_all, max_all)]
q18 <- RES[qc_flag == "flagged"]; stopifnot(nrow(q18) == 1)
key <- data.table(item = c(paste0("mean_sd_", GRP6), paste0("P_", CMP6_LAB), "ANOVA",
                           "fold_range_CRISPRi_vs_Ctrl", "fold_range_CtrlMhp_vs_Ctrl", "fold_range_CRISPRiMhp_vs_CtrlMhp", "fold_range_OEMhp_vs_CRISPRiMhp", "fold_wording_rule",
                           "sensitivity_verdict", "sensitivity_Ctrl_vs_CRISPRi_P", "not_shown_CRISPRi_vs_CRISPRiOE", "nuclear_objects_panel_note", "QC_flagged_image",
                           "n_per_group", "quantification_unit", "legend_sentence_EN", "legend_sentence_nuclear_objects_EN"),
                  value = c(sm[, sprintf("%.1f ± %.1f %% (%s)", mean_pct_area, sd_pct_area, values)], k4[, sprintf("P = %s (%s), diff %.1f", p_adj_txt, star, diff)],
                            attr(ST, "anova")[, sprintf("F(%d, %d) = %.2f, P = %s", df1, df2, F, fmtp(p_anova))],
                            sprintf("%s (11 variants + fragment-removed check; independent colour paths 0.35–0.42) — about one-quarter to one-half of Ctrl depending on method (0.31–0.32 with pore removal)", fr("CRISPRi / Ctrl")),
                            sprintf("%s (Ctrl vs Ctrl+Mhp ns)", fr("Ctrl+Mhp / Ctrl")), sprintf("%s", fr("CRISPRi+Mhp / Ctrl+Mhp")),
                            sprintf("%s (independent colour paths 2.6–5.1)", fr("CRISPRi+OE+Mhp / CRISPRi+Mhp")),
                            "Text and legend: give group mean ± s.d. or the ranges above; do not cite 2-decimal fold changes (sub-200 px pore-ring fragments are ~1/5 of the value in sparse groups but ~2% in CRISPRi+OE+Mhp, so fold changes depend on the method)",
                            paste(VER[, sprintf("%s: %s", comparison, verdict)], collapse = "; "),
                            sprintf("Ctrl vs CRISPRi: %s; non-significant: %s; fragment-removed check P = %s", k1$p_range_significant_variants, k1$p_nonsignificant_variants, k1$fragments_removed_p),
                            sprintf("CRISPRi vs CRISPRi+OE (uninfected rescue, not plotted): main P = %s (%s); %s; fragment-removed check P = %s → significance is method-dependent. Do NOT write that BHLHE40 rescue works only under Mhp or that there is no rescue in uninfected cells.",
                                    kx$main_p, kx$main_star, kx$p_range_significant_variants, kx$fragments_removed_p),
                            sprintf("Nuclear-object count (secondary; touching nuclei not split, dense groups under-estimated): %s. Differs from the area metric: %s. If the panel goes to Supplementary, state these pairs as ns in the legend.",
                                    kn_txt, kn_diff),
                            q18[, sprintf("%s #%d (%s): %s Main value %.2f%% after shadow exclusion.", group, replicate, image_file, qc_note, pct_area)],
                            "3 independent wells per group, one 10x field (1.79 x 1.00 mm) per well (independence of wells relies on the author's confirmation; chamber / plate layout pending, checklist Fig6-8)",
                            "crystal violet-stained area as % of the field; membrane pores excluded by colour/size (rule failed on one Ctrl image with dye-filled pores, flagged); illumination-shadow tiles excluded from numerator and denominator (post hoc, one image); one parameter set fixed on blinded images",
                            paste0("Transwell migration of 3D4/21 cells. Migrated cells on the lower membrane surface were stained with crystal violet and quantified as the stained area ",
                                   "(% of one 10x field per well; membrane pores excluded by an automated colour/size rule with one parameter set fixed on blinded images; ",
                                   "one control image with dye-filled pores and a corner illumination shadow is flagged in the Source Data). ",
                                   "Data are mean ± s.d.; n = 3 independent wells per group; one-way ANOVA with Tukey's HSD test; *P < 0.05, **P < 0.01, ***P < 0.001, ns not significant. Scale bar, 100 µm."),
                            sprintf("Number of darkly stained nuclear objects per mm² (touching nuclei were not separated, so dense groups are under-estimated); not significant: %s.", kn_ns)))
fwrite(key, file.path(TAB, "Fig6_Transwell_关键数字与图注.csv"), bom = TRUE)
for (f in file.path(TAB, c("Fig6_Transwell_组均值.csv", "Fig6_Transwell_统计_Tukey全部两两.csv", "Fig6_Transwell_ANOVA.csv", "Fig6_Transwell_敏感性_关键比较.csv", "Fig6_Transwell_敏感性_判定.csv",
                           "Fig6_Transwell_倍数范围.csv", "Fig6_Transwell_关键数字与图注.csv"))) chk_file(f)
chk_file(file.path(QC, "QC_敏感性_组均值.png"))
print(sm[, .(group, n, mean = round(mean_pct_area, 2), sd = round(sd_pct_area, 2), values, ratio_to_Ctrl = round(ratio_to_Ctrl, 2), ratio_to_CtrlMhp = round(ratio_to_CtrlMhp, 2))])
print(attr(ST, "anova")); print(ST[shown == TRUE, .(g1, g2, diff = round(diff, 2), p_adj = p_adj_txt, star)])
print(VER[, .(comparison, main_diff, main_p, main_star, n_same_direction, n_significant, p_range_significant_variants, fragments_removed_p, verdict)])
print(FOLD_RANGE[, .(ratio, main = round(main, 3), min_11 = round(min_11_variants, 3), max_11 = round(max_11_variants, 3), frag_removed = round(fragments_removed, 3))])
print(KN[, .(g1, g2, p_nuc = fmtp(p_nuc), star_nuc, p_area = fmtp(p_area), star_area, differs)])
print(FOLD_RANGE[, .(ratio, text_wording)])
