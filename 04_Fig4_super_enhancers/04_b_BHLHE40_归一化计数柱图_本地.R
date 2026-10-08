# ============================================================================
# 04_b_BHLHE40_归一化计数柱图_本地.R — Fig. 4b：肺组织 bulk RNA-seq 的 DESeq2 分析与 BHLHE40 归一化计数
# 运行位置：本地；路径变量见 config/paths.R
# ============================================================================

source(file.path(Sys.getenv("MHP_CODE", "."), "config", "paths.R"))

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
source(file.path(MHP_CODE, "shared", "_palette_MHP.R"))   # theme_mhp、sav、lw、pt2size、FONT、grp_col、smp_col、MHP_RES

RC_IN  <- file.path(MHP_RES, "1、RNA-seq结果/2、DEGs/gene.read_count.xls")
DEG_IN <- file.path(MHP_RES, "1、RNA-seq结果/2、DEGs/Con-vs-Mhp.DEseq2.before-filter.xls")
OUT <- file.path(MHP_RES, "6、SE/新Fig4_作图")
PAN <- file.path(OUT, "Plots/面板"); TAB <- file.path(OUT, "Tables"); PNL <- file.path(OUT, "_panels")
for (d in c(PAN, TAB, PNL)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
SMP <- c("Con1", "Con2", "Con3", "Mhp1", "Mhp2", "Mhp3")
GENE <- "BHLHE40"
# 散点描边用同色系深一档（取自调色板样本色 Con_1 / Mhp_1），柱边框用 grp_col（Fig 3 的 COL / COL_PT 做法）
PT_COL <- c("Con" = smp_col[["Con_1"]], "Mhp" = smp_col[["Mhp_1"]])
SF_REF <- c(Con1 = 1.0326, Con2 = 1.0491, Con3 = 1.1228, Mhp1 = 0.9542, Mhp2 = 0.8349, Mhp3 = 1.0621)

# ---------------------------------------------------------------- 读入原始计数 → DESeq2 中位数比值 size factor
rc <- fread(RC_IN, sep = "\t", quote = "", select = c("Gene_id", SMP))
stopifnot(nrow(rc) == 27577, !anyDuplicated(rc$Gene_id), GENE %in% rc$Gene_id, all(SMP %in% names(rc)))
cnt <- as.matrix(rc[, ..SMP]); rownames(cnt) <- rc$Gene_id
stopifnot(all(cnt >= 0), all(cnt == round(cnt)))                                # 原始计数为非负整数
size_factor <- function(m) {                                                    # DESeq2 默认 type = "ratio"：几何均值全为有限值的基因参与
  lg <- log(m); gm <- rowMeans(lg); ok <- is.finite(gm)
  apply(m[ok, , drop = FALSE], 2, function(x) exp(median(log(x) - gm[ok])))
}
sf <- size_factor(cnt)
stopifnot(all(abs(sf[SMP] - SF_REF[SMP]) < 5e-4))
if (requireNamespace("DESeq2", quietly = TRUE)) {                               # 交叉核对：与 DESeq2 官方实现一致（可选）
  dds <- suppressMessages(DESeq2::DESeqDataSetFromMatrix(cnt, data.frame(g = factor(substr(SMP, 1, 3))), ~ g))
  dds <- DESeq2::estimateSizeFactors(dds)
  stopifnot(max(abs(DESeq2::sizeFactors(dds)[SMP] - sf[SMP])) < 1e-6)
}
d <- data.table(gene = GENE, sample = SMP, group = factor(substr(SMP, 1, 3), levels = grp_levels),
                raw_count = cnt[GENE, SMP], size_factor = sf[SMP])
d[, norm_count := raw_count / size_factor]
stopifnot(nrow(d) == 6, d[, .N, by = group][, all(N == 3)], !anyNA(d$norm_count))
sm <- d[, .(n = .N, mean = mean(norm_count), sd = sd(norm_count)), by = group]
fold <- sm[group == "Mhp", mean] / sm[group == "Con", mean]

# ---------------------------------------------------------------- DESeq2 基因级 FDR → 星号；柱比与 log2FC 同源核对
dg <- fread(DEG_IN, sep = "\t", quote = "", select = c("GeneID", "Gene_name", "logFC", "Pvalue", "FDR", "Up/Down-Regulation"))
stopifnot(nrow(dg) == 23006, sum(dg$GeneID == GENE) == 1)                 # 基因级全表（非 68,921 行转录本级）
de <- dg[GeneID == GENE]
stopifnot(abs(de$logFC - 1.665) < 0.001, de$FDR < 1e-70, de$`Up/Down-Regulation` == "Up")
stopifnot(abs(fold / 2^de$logFC - 1) < 0.01)                                     # 3.17 = 2^1.665：柱比与星号所引统计量同一口径
star <- function(p) fcase(p < 0.001, "***", p < 0.01, "**", p < 0.05, "*", default = "ns")
st <- star(de$FDR)
welch <- t.test(d[group == "Mhp", norm_count], d[group == "Con", norm_count])   # 仅供参考，不上图

# ---------------------------------------------------------------- 作图（Fig 3 bar4 样式：三边空心柱 + 向上 SD + 空心圆点 + 横线星号）
ymax <- max(c(sm$mean + sm$sd, d$norm_count))
box <- sm[, .(x = as.integer(group) + c(-0.34, -0.34, 0.34, 0.34), y = c(0, mean, mean, 0)), by = group]
br <- data.table(x1 = 1, x2 = 2, h = ymax * 1.10, lab = st)
p <- ggplot(d, aes(group, norm_count)) +
  geom_path(data = box, aes(x = x, y = y, colour = group, group = group), linewidth = lw(1.0), lineend = "butt", linejoin = "mitre") +
  geom_errorbar(data = sm, aes(y = mean, ymin = mean, ymax = mean + sd, colour = group), width = 0.28, linewidth = lw(0.8)) +
  scale_colour_manual(values = grp_col, guide = "none") +
  ggnewscale::new_scale_colour() +
  geom_point(aes(colour = group), position = position_jitter(width = 0.12, height = 0, seed = 1),
             size = 1.5, stroke = 0.55, shape = 21, fill = "white") +          # 空心圆：描边即点本身，故保留 stroke
  scale_colour_manual(values = PT_COL, guide = "none") +
  geom_segment(data = br, aes(x = x1 + 0.08, xend = x2 - 0.08, y = h, yend = h), inherit.aes = FALSE, linewidth = lw(0.6)) +
  geom_text(data = br, aes(x = (x1 + x2) / 2, y = h, label = lab), inherit.aes = FALSE, vjust = 0.25,
            size = pt2size(FS_BASE + 1), family = FONT) +
  scale_x_discrete(expand = expansion(add = 0.55)) +                            # 柱两侧留白（与 Fig 3 / Fig 2h 一致）
  scale_y_continuous(limits = c(0, NA), breaks = scales::breaks_extended(4), expand = expansion(mult = c(0, 0.22))) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = expression(italic("BHLHE40") ~ "(normalized counts)")) +
  theme_mhp() + theme(axis.text = element_text(colour = "black"), plot.margin = margin(3, 2, 1, 1))
b <- ggplot_build(p)$layout$panel_params[[1]]
stopifnot(b$x.range[1] < 0.66 - 0.1, b$x.range[2] > 2.34 + 0.1, b$y.range[1] == 0)   # x 两侧确有留白、y 从 0 起

# ---------------------------------------------------------------- 输出
sav(p, "Fig4b_BHLHE40_归一化计数柱图", PAN, 22, 50, panel_dir = PNL)
if (identical(Sys.getenv("MHP_FINAL"), "1")) ggsave(file.path(PAN, "Fig4b_BHLHE40_归一化计数柱图.pdf"), p, width = 22, height = 50, units = "mm", device = .pdf_dev)
fwrite(d[, .(gene, sample, group, raw_count, size_factor = round(size_factor, 4), norm_count = round(norm_count, 2))],
       file.path(TAB, "04_b_BHLHE40_归一化计数_作图数据.csv"))
key <- data.table(item = c("纵轴量", "size factor 方法", "size factor Con1–3 / Mhp1–3",
                           "Con 归一化计数均值", "Con SD", "Mhp 归一化计数均值", "Mhp SD", "倍数 Mhp/Con（归一化计数）", "2^DESeq2 log2FC",
                           "n 每组（头猪）", "DESeq2 log2FC（基因级）", "DESeq2 P（Wald）", "DESeq2 FDR（BH）", "星号（按 DESeq2 FDR）",
                           "归一化计数 Welch t P（参考，不上图）", "read count 表行数", "DESeq2 表行数", "read count 表路径", "DESeq2 表路径"),
                  value = c("DESeq2 归一化计数 = raw count / size factor", "DESeq2 median-of-ratios（全部 27,577 基因；与 DESeq2::estimateSizeFactors 一致）",
                            paste(round(sf[SMP], 4), collapse = "; "),
                            round(sm[group == "Con", mean], 2), round(sm[group == "Con", sd], 2), round(sm[group == "Mhp", mean], 2),
                            round(sm[group == "Mhp", sd], 2), round(fold, 3), round(2^de$logFC, 3), 3, round(de$logFC, 4), signif(de$Pvalue, 3), signif(de$FDR, 3),
                            st, signif(welch$p.value, 3), nrow(rc), nrow(dg), RC_IN, DEG_IN))
fwrite(key, file.path(TAB, "04_b_BHLHE40_归一化计数_关键数字.csv"))
stopifnot(file.size(file.path(TAB, "04_b_BHLHE40_归一化计数_关键数字.csv")) > 0, st == "***")
print(sm); cat(sprintf("BHLHE40: 归一化计数 %.1f → %.1f（×%.2f；2^log2FC = %.2f）；DESeq2 log2FC %.3f，FDR %.2e → %s；Welch P = %.2e\n",
                       sm[group == "Con", mean], sm[group == "Mhp", mean], fold, 2^de$logFC, de$logFC, de$FDR, st, welch$p.value))
