# ============================================================================
# 05_e_E2_P2G8_CPM柱图_本地.R — Fig. 5f：E2 内与 BHLHE40 表达显著关联（peak-to-gene）的 8 个 peak 的伪 bulk 可及性（文件名中的面板字母为旧编号）
# 运行位置：本地；路径变量见 config/paths.R
# ============================================================================

source(file.path(Sys.getenv("MHP_CODE", "."), "config", "paths.R"))

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
source(file.path(MHP_CODE, "shared", "_palette_MHP.R"))   # theme_mhp、sav、lw、pt2size、FONT、grp_col、ct_atac_order、lab_atac、to_grp、MHP_RES

CSV  <- file.path(MHP_RES, "9.scatac 补充分析/Tables/24_SE/Fig5_byCelltypeGroup_fixed.csv")
P2G  <- file.path(MHP_RES, "9.scatac 补充分析/Tables/24_SE/P2G_links.csv")
OUT <- file.path(MHP_RES, "7、HI-C结果/新Fig5_作图")
PAN <- file.path(OUT, "Plots/面板"); TAB <- file.path(OUT, "Tables"); PNL <- file.path(OUT, "_panels")
for (d in c(PAN, TAB, PNL)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
VAL <- "P2G8"; STEM <- "Fig5f_E2_P2G8_CPM"; W_MM <- 44; H_MM <- 44
XLAB <- "8 *BHLHE40*-linked<br>E2 peaks (CPM)"
XLAB_TXT <- "8 BHLHE40-linked / E2 peaks (CPM)"                                     # 关键数字表用的纯文本版
suppressPackageStartupMessages(library(ggtext))                                    # ★rds 含 element_markdown，拼版脚本读 rds 前也需 library(ggtext)
CT_KEEP <- setdiff(ct_atac_order, "B_cells")                                        # B 细胞剔除
ANNOT <- c("Macrophages", "Neutrophils")                                            # 只在髓系两类标百分比
E2_REGION <- c(61329246L, 61350682L)

# ---------------------------------------------------------------- ① 读表 + 核对
d0 <- fread(CSV)
stopifnot(nrow(d0) == 16, identical(names(d0), c("group_ct", "E2_20peak", "P2G8", "Core61335772", "BHLHE40")))
d0[, `:=`(grp_raw = sub("_.*$", "", group_ct), ct = sub("^(Mock|Mhp)_", "", group_ct))]
stopifnot(all(d0$grp_raw %in% c("Mock", "Mhp")), all(d0$ct %in% ct_atac_order), d0[, .N, by = ct][, all(N == 2)])
d0[, group := to_grp(grp_raw)]
d <- d0[ct %in% CT_KEEP][, .(ct, group, value = get(VAL))]
d[, ct := factor(ct, levels = rev(CT_KEEP))]                                        # 横向：第一类在最上
stopifnot(nrow(d) == 14, !anyNA(d$value), all(d$value >= 0))
w <- dcast(d, ct ~ group, value.var = "value")[, `:=`(pct = (Mhp - Con) / Con * 100, fold = Mhp / Con)]
w[, lab := ifelse(ct %in% ANNOT, sprintf("%+.0f%%", pct), NA_character_)]
stopifnot(abs(w[ct == "Neutrophils", Con] - 16.04) < 0.01, abs(w[ct == "Neutrophils", Mhp] - 32.85) < 0.01,   # 与方案 §2.2 / 图注草稿一致
          abs(w[ct == "Macrophages", Con] - 6.14) < 0.01, abs(w[ct == "Macrophages", Mhp] - 9.49) < 0.01)

# 这 8 个 peak 的清单（E2 区间内、gene == BHLHE40）
p2g <- fread(P2G)
stopifnot(all(c("peak_chr", "peak_start", "peak_end", "gene", "gene_tss", "cor", "FDR") %in% names(p2g)))
pk <- p2g[peak_chr == "13" & gene == "BHLHE40" & peak_start >= E2_REGION[1] & peak_end <= E2_REGION[2]][order(peak_start)]
stopifnot(nrow(pk) == 8, all(pk$gene_tss == 61526116), min(pk$cor) > 0.45, max(pk$FDR) <= 1e-4)
cat(sprintf("E2 内 P2G-BHLHE40 peak %d 个：r %.3f–%.3f，FDR ≤ %.3g\n", nrow(pk), min(pk$cor), max(pk$cor), max(pk$FDR)))

# ---------------------------------------------------------------- ② 作图（横向成对柱；Con 在上、Mhp 在下）
xmax <- max(d$value)
lab <- w[!is.na(lab), .(ct, x = pmax(Con, Mhp) + xmax * 0.03, lab)]
p <- ggplot(d, aes(x = value, y = ct, fill = group)) +
  geom_col(position = position_dodge2(width = 0.72, reverse = TRUE, padding = 0.08), width = 0.72) +   # reverse：Con 在上
  geom_text(data = lab, aes(x = x, y = ct, label = lab), inherit.aes = FALSE, hjust = 0, vjust = 0.5, size = pt2size(FS_SMALL), family = FONT, colour = "black") +
  scale_fill_manual(values = grp_col, name = NULL) +
  scale_y_discrete(labels = lab_atac) +
  scale_x_continuous(limits = c(0, NA), breaks = scales::breaks_extended(4), expand = expansion(mult = c(0, 0.36))) +   # 右侧留白放百分比标签
  labs(x = XLAB, y = NULL) +
  theme_mhp() + theme(axis.text = element_text(colour = "black"), axis.line.y = element_blank(), axis.ticks.y = element_blank(),
                      axis.title.x = ggtext::element_markdown(size = FS_BASE, lineheight = 1.05, halign = 0.5, margin = margin(t = 1)),   # 两行轴标题（ggtext）
                      legend.position = "inside", legend.position.inside = c(0.80, 0.13), legend.key.size = unit(2.6, "mm"),
                      legend.text = element_text(size = FS_SMALL), legend.key.spacing.y = unit(0.4, "mm"), legend.background = element_blank(),
                      plot.margin = margin(2, 2, 0.5, 1))
b <- ggplot_build(p)$layout$panel_params[[1]]
PW_EST <- W_MM - 15 - 2                                                            # 面板宽估计（mm）：总宽 − 类型名列（"Endothelial cells" 6 pt ≈ 14 mm + 刻度）− 边距
lab[, w_data := nchar(lab) * 0.55 * FS_SMALL * 25.4 / 72 * b$x.range[2] / PW_EST]  # 标签宽（数据单位；字宽按 0.55 em 估）
stopifnot(b$x.range[1] == 0, all(b$x.range[2] > lab$x + lab$w_data))              # x 从 0 起，百分比标签整段在面板内

# ---------------------------------------------------------------- ③ 输出
sav(p, STEM, PAN, W_MM, H_MM, panel_dir = PNL)
if (identical(Sys.getenv("MHP_FINAL"), "1")) {
  out <- file.path(PAN, paste0(STEM, ".pdf")); ggsave(out, p, width = W_MM, height = H_MM, units = "mm", device = .pdf_dev)
  stopifnot(sum(grepl("/Subtype /Image", readLines(out, warn = FALSE, skipNul = TRUE), useBytes = TRUE)) == 0) }
fwrite(w[order(-as.integer(ct)), .(celltype = as.character(ct), display = lab_atac(as.character(ct)), metric = VAL, Con = signif(Con, 5), Mhp = signif(Mhp, 5),
                                   fold = round(fold, 3), pct_change = round(pct, 1), annotated = !is.na(lab))],
       file.path(TAB, "05_e_E2_P2G8_作图数据.csv"))
fwrite(pk[, .(peak_chr, peak_start, peak_end, gene, gene_tss, cor = signif(cor, 4), FDR = signif(FDR, 3),
              offset_from_E2_start_kb = round((peak_start - E2_REGION[1]) / 1000, 2))], file.path(TAB, "05_e_E2_P2G8_peak清单.csv"))
key <- data.table(
  item = c("横轴量", "横轴名（图上）", "8 个 peak 的坐标（chr13）", "P2G 相关 r 范围 / FDR 上限", "Con → Mhp（中性粒）", "Con → Mhp（巨噬）", "其余 5 类变化",
           "剔除", "n 的性质（无统计）", "与面板 e 的口径差异", "图注建议句"),
  value = c("E2 内 8 个 P2G-BHLHE40 peak 的伪 bulk 可及性均值（CPM；ArchR getGroupSE PeakMatrix, divideN = TRUE，按组 CPM）",
            XLAB_TXT,
            paste(pk[, sprintf("%s-%s", format(peak_start, big.mark = ","), format(peak_end, big.mark = ","))], collapse = "; "),
            sprintf("r %.3f–%.3f；FDR %.3g–%.3g（图注草稿写「FDR ≤ 1.55e-86」实为最小值，上限应写 %.2g）", min(pk$cor), max(pk$cor), min(pk$FDR), max(pk$FDR), max(pk$FDR)),
            w[ct == "Neutrophils", sprintf("%.2f → %.2f（%+.0f%%）", Con, Mhp, pct)], w[ct == "Macrophages", sprintf("%.2f → %.2f（%+.0f%%）", Con, Mhp, pct)],
            paste(w[!ct %in% ANNOT][order(-as.integer(ct)), sprintf("%s %+.0f%%", lab_atac(as.character(ct)), pct)], collapse = "; "),
            sprintf("B cells（Mhp 组仅 3 核）"), "每组 = 3 只动物合并伪 bulk，每类型 × 组一个值；不做检验，只作描述",
            "本面板 = 8 个 P2G peak 的伪 bulk CPM（中性粒 +105%、巨噬 +55%）；面板 e = 20 peak 插补值（+95% / +35%），两种口径不同",
            sprintf("Pseudo-bulk accessibility (CPM) of the eight E2 peaks whose accessibility co-varies with BHLHE40 expression (peak-to-gene links, r = %.3f-%.3f, FDR <= %.2g), by cell type and condition (this 8-peak CPM is a different quantity from the imputed 20-peak score in e); B cells excluded (n = 3 nuclei in Mhp); percentages = (Mhp - Con)/Con for myeloid populations (neutrophils %.2f -> %.2f, %+.0f%%; macrophages %.2f -> %.2f, %+.0f%%); all other cell types < 5 CPM",
                    min(pk$cor), max(pk$cor), max(pk$FDR), w[ct == "Neutrophils", Con], w[ct == "Neutrophils", Mhp], w[ct == "Neutrophils", pct],
                    w[ct == "Macrophages", Con], w[ct == "Macrophages", Mhp], w[ct == "Macrophages", pct])))
fwrite(key, file.path(TAB, "05_e_E2_P2G8_关键数字.csv"))
stopifnot(file.size(file.path(TAB, "05_e_E2_P2G8_关键数字.csv")) > 0)
print(w[order(-as.integer(ct)), .(ct, Con = round(Con, 2), Mhp = round(Mhp, 2), pct = round(pct))])
