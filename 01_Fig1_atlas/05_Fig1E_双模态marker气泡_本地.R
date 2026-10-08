# ============================================================================
# 05_Fig1E_双模态marker气泡_本地.R — Fig. 1e、1f：scRNA-seq（11 类）与 snATAC-seq（8 类）marker 基因气泡图
# 运行位置：本地；路径变量见 config/paths.R
# ============================================================================

source(file.path(Sys.getenv("MHP_CODE", "."), "config", "paths.R"))

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
source(file.path(MHP_CODE, "shared", "_palette_MHP.R"))

TB   <- file.path(MHP_DIR$fig1, "Tables")
OUT  <- file.path(MHP_DIR$fig1, "Plots/02_双模态marker")
PANL <- file.path(MHP_DIR$fig1, "_panels")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ---------- 1. 读表并算口径 ----------
cnt <- fread(file.path(TB, "30_Fig1v2_RNA导出/RNA_pb_counts_byCelltype.csv.gz"))
pct <- fread(file.path(TB, "30_Fig1v2_RNA导出/RNA_pct_byCelltype.csv.gz"))
gs  <- fread(file.path(TB, "31_Fig1v2_ATAC导出/ATAC_GeneScore_mean_byCelltype.csv.gz"))
stopifnot(identical(cnt$symbol, pct$symbol), all(ct_rna_order %in% names(cnt)), all(ct_atac_order %in% names(gs)))

TRD_ID <- "ENSSSCG00000056250"; TRD_SYM <- "TRDC"
i_trd  <- which(cnt$ensembl_id == TRD_ID)
stopifnot(length(i_trd) == 1, !TRD_SYM %in% cnt$symbol, TRD_ID %in% gs$ensembl_id)   # RNA 唯一、无重名、ATAC 有对应行
cnt[i_trd, symbol := TRD_SYM]; pct[i_trd, symbol := TRD_SYM]
FIXED <- list("NK cells" = c(TRD_SYM, "SOX13"))

M   <- as.matrix(cnt[, ct_rna_order, with = FALSE]); rownames(M) <- cnt$symbol
CPM <- sweep(M, 2, colSums(M), "/") * 1e6                                 # 伪 bulk CPM（11 类）
PC  <- as.matrix(pct[, ct_rna_order, with = FALSE]); rownames(PC) <- cnt$symbol
GS  <- as.matrix(gs[, ct_atac_order, with = FALSE]); rownames(GS) <- gs$ensembl_id   # GeneScore 类均值（8 类），行名 = ensembl

# RNA 基因 → ATAC 行：按 ensembl 匹配；"a;b" 一对多拆开取 GeneScore 总和最高者
atac_row <- rep(NA_character_, nrow(cnt)); multi <- grepl(";", cnt$ensembl_id)
atac_row[!multi] <- ifelse(cnt$ensembl_id[!multi] %in% rownames(GS), cnt$ensembl_id[!multi], NA)
for (i in which(multi)) { ids <- intersect(strsplit(cnt$ensembl_id[i], ";")[[1]], rownames(GS))
  if (length(ids)) atac_row[i] <- ids[which.max(rowSums(GS[ids, , drop = FALSE]))] }
cat("RNA 基因", nrow(cnt), "| 匹配到 ATAC 行:", sum(!is.na(atac_row)), "\n")

# 每基因统计量：RNA 最高类型、领先幅度 gap（最高 − 次高，按 e/max 缩放后）、最高类型内表达比例；ATAC 对应类排名与 log2 fold
zr  <- CPM / pmax(apply(CPM, 1, max), 1e-9)
top <- colnames(zr)[apply(zr, 1, which.max)]
gap <- apply(zr, 1, function(x) { v <- sort(x, decreasing = TRUE); v[1] - v[2] })
pct_top <- PC[cbind(seq_len(nrow(PC)), match(top, colnames(PC)))]
GA  <- matrix(NA_real_, nrow(cnt), ncol(GS), dimnames = list(cnt$symbol, colnames(GS)))
hit <- !is.na(atac_row); GA[hit, ] <- GS[atac_row[hit], , drop = FALSE]
fold <- log2((GA + .01) / (rowMeans(GA) + .01))                          # 同 08：相对 8 类均值的 log2 fold
es_a <- pmax(pmin(fold / 2, 1), 0)
atac_type <- ct_map$atac[match(top, ct_map$rna)]                           # RNA 最高类型对应的 ATAC 类型（DC 为 NA）
rk <- t(apply(-GA, 1, rank)); rk[is.na(GA)] <- NA
a_rank <- rk[cbind(seq_len(nrow(rk)), match(atac_type, colnames(rk)))]
a_fold <- fold[cbind(seq_len(nrow(fold)), match(atac_type, colnames(fold)))]
a_es   <- es_a[cbind(seq_len(nrow(es_a)), match(atac_type, colnames(es_a)))]

st <- data.table(symbol = cnt$symbol, ensembl_id = cnt$ensembl_id, atac_ensembl = atac_row, top = top,
                 gap = gap, pct = pct_top, cpm = apply(CPM, 1, max), atac_rank = a_rank, atac_fold = a_fold,
                 atac_es = ifelse(is.na(a_es), 0, a_es))

# ---------- 2. 候选池（自动 + 知识约束）----------
CLASSIC <- c("CD163","C1QA","MRC1","CD68","MARCO","CSF3R","CXCL8","S100A8","FLT3","HLA-DRA","FCER1A","CD4","CD40LG","IL7R",
             "CD8A","CD8B","GZMA","NCR1","KLRB1","CD79A","MS4A1","CD19","MZB1","EPCAM","SFTPC","KRT8","NKX2-1",
             "CDH5","PECAM1","CLDN5","KDR","COL1A1","COL3A1","DCN","BGN")
EXT <- c("LYZ","C1QB","C1QC","MSR1","CD14","APOE",                                          # 巨噬
         "TREM1","G0S2","S100A12","SELL","IL1R2","MMP9","LTF","CXCR2",                       # 中性粒
         "IRF8","CD86","CD74","CST3","CLEC9A","XCR1",                                        # DC
         "ICOS","CD5","CD28","CD2","CD3E","CD3D","TCF7","LEF1","CTLA4",                      # CD4 T
         "KLRK1","CCL5","GNLY","PRF1","GZMB","GZMK","CD3G",                                  # CD8 T / 细胞毒
         "KLRD1","KLRF1","NCAM1","EOMES","TBX21","XCL1","KLRC1",                             # NK
         "CD79B","PAX5","BANK1","CD22",                                                      # B
         "XBP1","PRDM1","IRF4","TNFRSF17","SDC1","CD38",                                     # 浆
         "SFTPB","SFTPA1","KRT18","NAPSA","CDH1","SCGB1A1","FOXJ1",                          # 上皮
         "VWF","TEK","FLT1","ESAM","EMCN","CD34",                                            # 内皮
         "PDGFRA","LUM","COL1A2","FBLN1","PDGFRB")                                           # 成纤维
BLACK <- c("NKG7", "JCHAIN")                                                                 # 已知不宜作 marker
GAP_MIN <- 0.35; PCT_MIN <- 0.40; N_FILL <- 30

st[, cn := fifelse(symbol %in% CLASSIC, 1, fifelse(symbol %in% EXT, 0.6, 0))]
st[, fit := gap >= GAP_MIN & pct >= PCT_MIN & (is.na(atac_type) | (!is.na(atac_rank) & atac_rank <= 2))]
usable <- st[!is.na(atac_ensembl) & !symbol %in% BLACK & !grepl("^(ENSSSCG|LOC)", symbol) & cpm >= 50]
# 三层候选：经典 + 扩展（凡 RNA 有最高类型者都进池）+ 每类型按打分补选 N_FILL 个（无合格经典候选的类型靠这层）
usable[, bs := 0.35 * cn + 0.45 * gap + 0.35 * pct + 0.10 * atac_es + 1e-4 * gap]     # 末项：平局按 gap 高者
pool <- rbind(usable[cn > 0], usable[cn == 0][order(-fit, -bs)][, head(.SD, N_FILL), by = top])
pool <- unique(rbind(pool, usable[symbol %in% unlist(FIXED)]), by = "symbol")
stopifnot(all(unlist(FIXED) %in% pool$symbol))
cat("候选池:", nrow(pool), "基因 | 经典", sum(pool$cn == 1), "| 扩展", sum(pool$cn == 0.6), "| 全表补选", sum(pool$cn == 0), "\n")

# ---------- 3. 全局最优分配：22 个槽位（11 类 × 2）× 候选基因，匈牙利算法 ----------
SLOT <- rep(ct_rna_order, each = 2)
SC <- matrix(0, nrow(pool), length(SLOT), dimnames = list(pool$symbol, paste0(SLOT, "#", 1:2)))
for (j in seq_along(SLOT)) { r <- pool$top == SLOT[j]                     # 只能分配给自己 RNA 最高的类型
  SC[r, j] <- ifelse(pool$fit[r], 2.0, 0.5) + pool$bs[r]                 # 合格档 2.0 / 不合格档 0.5（同 08）
  if (SLOT[j] %in% names(FIXED)) {
    fx <- pool$symbol %in% FIXED[[SLOT[j]]]
    stopifnot(sum(fx) == 2, all(pool$top[fx] == SLOT[j]))                  # 指定基因必须恰为 2 个，且 RNA 最高类型就是该类
    SC[!fx, j] <- 0; SC[fx, j] <- 100 + pool$bs[fx] } }
stopifnot(all(colSums(SC > 0) >= 1))
a <- clue::solve_LSAP(t(SC), maximum = TRUE)                                # 22 槽位 × 候选基因
sel <- pool[as.integer(a)][, celltype := SLOT]
sel <- sel[order(match(celltype, ct_rna_order), -bs)]
stopifnot(nrow(sel) == 22, !anyDuplicated(sel$symbol), all(sel$top == sel$celltype))

# 入选理由（中文，写入清单）
src <- c("1" = "经典候选", "0.6" = "扩展文献marker", "0" = "全表补选")
sel[, source := fifelse(symbol %in% unlist(FIXED), "人工指定", src[as.character(cn)])]
sel[, 入选理由 := sprintf("%s；RNA 最高=%s，gap=%.2f，pct=%.2f%s%s", source, celltype, gap, pct,
                          ifelse(is.na(atac_type[match(symbol, st$symbol)]), "；DC 在 snATAC 无对应类，不设 ATAC 门槛",
                                 sprintf("；ATAC %s 排名 %d（log2 fold %.2f）", ct_map$atac[match(celltype, ct_map$rna)], atac_rank, atac_fold)),
                          fifelse(fit, "；合格", fifelse(symbol %in% unlist(FIXED), "；不合格（ATAC 门槛未达，人工指定保留）", "；不合格（无合格候选时的最优补选）")))]
sel[celltype == "NK cells", 入选理由 := paste0(入选理由,
    "；注：数据值 NK cells 簇经 marker 核实为 γδ T 细胞（图上显示 γδ T cells，数据值不变）——",
    "TRD 链 ENSSSCG00000056250（猪 Ensembl 无 symbol，图上显示 TRDC）在该簇 68% 细胞表达（CD8 簇 11%、CD4 簇 4%），",
    "CD3E+ 81%，而 NCR1+ 仅 1%、KLRB1+ 4%，并特异表达 SOX13/BLK/RHEX；ATAC T_NK_cells 的 TRD 链 GeneScore 为 8 类最高（log2 fold 0.48）。",
    "本行两个 marker 人工指定为 TRDC + SOX13，不参与匈牙利分配")]
n_fit <- sum(sel$fit); cat("合格 marker:", n_fit, "/ 22\n")
print(sel[, .(celltype, symbol, gap = round(gap, 2), pct = round(pct, 2), cpm = round(cpm), atac_rank, atac_fold = round(atac_fold, 2), src = source, fit)], nrows = 30)
out_csv <- sel[, .(celltype, symbol, ensembl_id, atac_ensembl, rna_gap = round(gap, 3), rna_pct = round(pct, 3), rna_cpm = round(cpm, 1),
                   atac_rank, atac_log2fold = round(atac_fold, 3), score = round(bs, 4), qualified = fit, source, 入选理由)]
fwrite(out_csv, file.path(TB, "Fig1E_marker清单.csv"), bom = TRUE)         # BOM 便于 Excel 正确识别中文
cat("✓ Fig1E_marker清单.csv", nrow(out_csv), "行\n")

# ---------- 4. 作图数据：行 = 22 基因（按类型分组），列 = RNA 11 类 / ATAC 8 类 ----------
G   <- sel$symbol
LIN <- c(rep("Myeloid", 3), rep("Lymphoid", 5), rep("Structural", 3)); names(LIN) <- ct_rna_order   # 髓系 / 淋巴 / 结构
build <- function(MODE) {
  if (MODE == "RNA") { E <- CPM[G, ct_rna_order]; P <- PC[G, ct_rna_order]; lab <- lab_rna(ct_rna_order) }
  else               { E <- GS[sel$atac_ensembl, ct_atac_order]; rownames(E) <- G; P <- E; lab <- lab_atac(ct_atac_order) }
  dt <- CJ(gene = factor(G, levels = rev(G)), ct = factor(lab, levels = lab), sorted = FALSE)
  dt[, `:=`(e = E[cbind(match(as.character(gene), rownames(E)), as.integer(ct))],
            p = P[cbind(match(as.character(gene), rownames(P)), as.integer(ct))])]
  if (MODE == "ATAC") dt[, es := { f <- log2((e + .01) / (mean(e) + .01)); pmax(pmin(f / 2, 1), 0) }, by = gene]   # 绝不按行缩放
  else                dt[, es := e / max(e), by = gene]
  dt }
dR <- build("RNA"); dA <- build("ATAC")
stopifnot(!anyNA(dR$e), !anyNA(dA$e))

# ---------- 5. 版式 ----------
PR <- c("#EEF3F8", "#9EC5E0", "#2D6FA8", "#12365C")                       # RNA 色带（蓝，沿用 08；本脚本自定义）
PA <- c("#FDF2EA", "#F3C39C", "#D2703C", "#7E2C10")                       # ATAC 色带（橙棕，沿用 08；本脚本自定义）
FS_LEG <- 5; SEP_COL <- "grey75"; PT_MAX <- 2.3
nG <- length(G)
gene_brk <- cumsum(table(factor(sel$celltype, levels = ct_rna_order)))    # 各类型基因数累计（每类 2 → 2,4,…,22）
lin_row_brk <- gene_brk[c("DC cells", "Plasma cells")]                    # 髓系|淋巴、淋巴|结构 的行分隔
base <- function(legpos) theme_mhp(base = FS_BASE) + theme(
  axis.line = element_blank(), axis.title = element_blank(),
  axis.ticks = element_line(linewidth = lw(0.4), colour = "black"), axis.ticks.length = unit(0.5, "mm"),
  axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = FS_SMALL, colour = "black"),
  panel.border = element_rect(fill = NA, colour = "black", linewidth = lw(0.5)),
  plot.title = element_text(size = FS_BASE, face = "bold", hjust = 0.5, margin = margin(0, 0, 1.5, 0, "mm")),
  legend.position = legpos, legend.key.size = unit(2.4, "mm"), legend.text = element_text(size = FS_LEG),
  legend.title = element_text(size = FS_SMALL), legend.margin = margin(0, 0, 0, 0), legend.spacing.y = unit(2, "mm"),   # 审图：1 → 2 mm，色带下标 "0" 与下一图例标题不再贴近
  legend.box.spacing = unit(1.5, "mm"), legend.ticks.length = unit(0.4, "mm"),
  plot.margin = margin(1, 0.5, 0.5, 0.5, "mm"))
# 顶部类型色条 + 三谱系分隔线（行：基因分组；列：类型分组）
deco <- function(cols, col_brk) list(
  annotate("rect", xmin = seq_along(cols) - 0.42, xmax = seq_along(cols) + 0.42, ymin = nG + 0.9, ymax = nG + 1.35, fill = unname(cols), colour = NA),
  geom_hline(yintercept = nG - lin_row_brk + 0.5, linewidth = lw(0.4), colour = SEP_COL),
  geom_vline(xintercept = col_brk + 0.5, linewidth = lw(0.4), colour = SEP_COL),
  coord_cartesian(clip = "off"))
yexp <- expansion(add = c(0.6, 1.6)); xexp <- expansion(add = 0.6)
cbar <- function() guide_colourbar(order = 1, theme = theme(legend.key.width = unit(1.8, "mm"), legend.key.height = unit(9, "mm")))

p1 <- ggplot(dR, aes(ct, gene)) + deco(ct_rna_col, cumsum(table(LIN[ct_rna_order])[c("Myeloid", "Lymphoid")])) +
  geom_point(aes(size = p, colour = es), stroke = 0) +
  scale_size_continuous(range = c(0.05, PT_MAX), limits = c(0, 1), breaks = c(.25, .5, .75), labels = c("25", "50", "75"), name = "% cells\nexpressing") +
  scale_colour_gradientn(colours = PR, limits = c(0, 1), breaks = c(0, .5, 1), labels = c("0", "0.5", "1"), name = "Scaled\nexpression") +
  scale_y_discrete(position = "right", expand = yexp) + scale_x_discrete(expand = xexp) +
  labs(title = "scRNA-seq") + base("left") +
  guides(colour = cbar(), size = guide_legend(order = 2, override.aes = list(colour = PR[3]))) +
  theme(axis.text.y.right = element_text(size = FS_SMALL, face = "italic", colour = "black", hjust = 0.5, margin = margin(0, 0, 0, 1, "mm")),
        axis.ticks.y = element_blank())

p2 <- ggplot(dA, aes(ct, gene)) + deco(ct_atac_col, cumsum(table(LIN[ct_map$rna[match(ct_atac_order, ct_map$atac)]])[c("Myeloid", "Lymphoid")])) +
  geom_point(aes(size = es, colour = es), stroke = 0) +
  scale_size_continuous(range = c(0.05, PT_MAX), limits = c(0, 1), guide = "none") +
  scale_colour_gradientn(colours = PA, limits = c(0, 1), breaks = c(0, .5, 1), labels = c("0", "1", "2+"), name = "Gene score\n(log2 FC)") +
  scale_y_discrete(expand = yexp) + scale_x_discrete(expand = xexp) +
  labs(title = "snATAC-seq") + base("right") + guides(colour = cbar()) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

p <- p1 + p2 + plot_layout(widths = c(length(ct_rna_order), length(ct_atac_order)))
sav(p, "Fig1E_双模态marker气泡", OUT, w_mm = 118, h_mm = 74, final = FALSE, panel_dir = PANL)

# ---------- 6. 自检数字 ----------
diagR <- sum(sel$celltype == colnames(zr)[apply(zr[G, , drop = FALSE], 1, which.max)])
mapped <- !is.na(ct_map$atac[match(sel$celltype, ct_map$rna)])
diagA <- sum(sel$atac_rank[mapped] == 1, na.rm = TRUE)
cat(sprintf("RNA 对角(基因在自己类型最高) %d/22 | ATAC 对应类排名第一 %d/%d、前二 %d/%d | 合格 %d/22\n",
            diagR, diagA, sum(mapped), sum(sel$atac_rank[mapped] <= 2, na.rm = TRUE), sum(mapped), n_fit))
