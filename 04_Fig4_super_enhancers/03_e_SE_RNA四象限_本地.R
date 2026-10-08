# ============================================================================
# 03_e_SE_RNA四象限_本地.R — Fig. 4e：超级增强子 H3K27ac 变化与对应基因 mRNA 变化的四象限图
# 运行位置：本地；路径变量见 config/paths.R
# ============================================================================

source(file.path(Sys.getenv("MHP_CODE", "."), "config", "paths.R"))

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(ggrepel) })
source(file.path(MHP_CODE, "shared", "_palette_MHP.R"))   # theme_mhp、sav、lw、pt2size、FONT、grp_col、MHP_RES

SE_DIR <- MHP_ROSE
RNA_F  <- file.path(MHP_RES, "1、RNA-seq结果/2、DEGs/Con-vs-Mhp.DEseq2.before-filter.xls")
FIG    <- file.path(MHP_RES, "6、SE/新Fig4_作图"); PAN <- file.path(FIG, "Plots/面板"); TAB <- file.path(FIG, "Tables"); PNL <- file.path(FIG, "_panels")
for (d in c(PAN, TAB, PNL)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
W_MM <- 50; H_MM <- 56                                                # 出图尺寸（mm）：整合面板建议 §2.1 的 e = 50 × 56
THR <- 0.58; LAB_THR <- 1; FDR_CUT <- 0.05; SEED <- 20260927
BH_COL <- "#B2182B"        # ★BHLHE40 高亮红：Fig 4 c / d / e 三面板统一（Brewer RdBu 深红，与面板 f 热图色标同源）；待并入 _shared/_palette_MHP.R
LABEL_MODE <- Sys.getenv("FIG4E_LABEL", "author")
stopifnot(LABEL_MODE %in% c("author", "deg"))
LAB_GENES <- c("BHLHE40", "TNIP3", "IL1R2", "IL1R1", "DUSP4")

# ---------------------------------------------------------------- 1. 读入与核验
PAIR_F <- file.path(TAB, "00b_SE配对_1680.csv")
if (!file.exists(PAIR_F)) stop("缺少配对表，请先运行 00b_SE配对_本地重做.R：", PAIR_F)
SE <- fread(PAIR_F)
stopifnot(nrow(SE) == 1680, !anyDuplicated(SE$EnhancerID_Con),
          all(c("EnhancerID_Con", "EnhancerID_Mhp", "ActivityScore_Con", "ActivityScore_Mhp", "SE_log2FC", "Gene_Con_clean", "Gene_Mhp_clean", "in_company_table") %in% names(SE)),
          max(abs(SE$SE_log2FC - log2((SE$ActivityScore_Mhp + 1) / (SE$ActivityScore_Con + 1)))) < 1e-6,
          sum(SE$in_company_table == FALSE) == 1, SE[in_company_table == FALSE, Gene_Con_clean] == "ZFP36L1")
CO <- fread(file.path(SE_DIR, "DE se/SE_log2FC_all_results.csv"))                          # 公司 1,679 对：只核对（配对与 SE_log2FC 逐行一致）
KC <- merge(CO[, .(EnhancerID_Con, EnhancerID_Mhp, SE_log2FC)], SE[in_company_table == TRUE, .(EnhancerID_Con, EnhancerID_Mhp, SE_log2FC)], by = "EnhancerID_Con", suffixes = c("", ".new"))
stopifnot(nrow(CO) == 1679, nrow(KC) == 1679, identical(KC$EnhancerID_Mhp, KC$EnhancerID_Mhp.new), max(abs(KC$SE_log2FC - KC$SE_log2FC.new)) < 1e-9)
G <- fread(RNA_F, sep = "\t", quote = "", select = 1:13, colClasses = list(character = c(1, 12, 13)))
setnames(G, "Up/Down-Regulation", "RNA_dir")
# 基因级全表核验：行数、GeneID 无转录本号、Up / Down 计数、BHLHE40 与 DEG 表一致（logFC 1.665、FDR 7.52e-76）
stopifnot(nrow(G) == 23006, !anyDuplicated(G$GeneID), !any(grepl("^(XM_|XR_|NM_|NR_)[0-9]", G$GeneID)),
          sum(G$RNA_dir == "Up") == 1230, sum(G$RNA_dir == "Down") == 909, !anyNA(G$logFC), !anyNA(G$FDR))
bh <- G[GeneID == "BHLHE40"]
stopifnot(nrow(bh) == 1, abs(bh$logFC - 1.664776) < 1e-5, abs(bh$FDR / 7.52e-76 - 1) < 0.01, bh$RNA_dir == "Up")

# ---------------------------------------------------------------- 2. SE ↔ 基因配对（Con 侧最近基因优先，退而 Mhp 侧）
SE[, con_hit := Gene_Con_clean %in% G$GeneID][, mhp_hit := Gene_Mhp_clean %in% G$GeneID]
SE[, gene := fifelse(con_hit, Gene_Con_clean, fifelse(mhp_hit, Gene_Mhp_clean, NA_character_))]
SE[, gene_side := fifelse(con_hit, "Con", fifelse(mhp_hit, "Mhp", "none"))]
n_drop <- sum(is.na(SE$gene))
M <- merge(SE[!is.na(gene)], G[, .(gene = GeneID, RNA_log2FC = logFC, RNA_P = Pvalue, RNA_FDR = FDR, RNA_dir)], by = "gene")
stopifnot(nrow(M) + n_drop == 1680, nrow(M) == 1642, uniqueN(M$gene) == 1546, sum(M$gene_side == "Mhp") == 7, n_drop == 38)   # 1,679 对时 1,641 / 1,545
M[, quadrant := fcase(SE_log2FC >=  THR & RNA_log2FC >=  THR, "Q1",
                      SE_log2FC >=  THR & RNA_log2FC <= -THR, "Q2",
                      SE_log2FC <= -THR & RNA_log2FC <= -THR, "Q3",
                      SE_log2FC <= -THR & RNA_log2FC >=  THR, "Q4", default = "NS")]
M[, is_DEG := RNA_FDR < FDR_CUT & abs(RNA_log2FC) >= 1]            # 手稿 DEG 口径（公司 known DEG 表：|log2FC| ≥ 1、FDR < 0.05）

# 与公司 141 行配对 DEG 及其 0.58 阈值 Q1 清单核对：本表为其超集，基因与 RNA log2FC 逐行一致
CD <- fread(file.path(SE_DIR, "SE-RNA/RNA_seq_vs_SE_complete_data.csv"))
K  <- merge(CD[, .(EnhancerID_Con, EnhancerID_Mhp, Gene_for_merge, RNA_logFC)],
            M[, .(EnhancerID_Con, EnhancerID_Mhp, gene, RNA_log2FC, quadrant)], by = c("EnhancerID_Con", "EnhancerID_Mhp"))
stopifnot(nrow(K) == 141, all(K$gene == K$Gene_for_merge), max(abs(K$RNA_logFC - K$RNA_log2FC)) < 1e-6)
q1_company <- readLines(file.path(SE_DIR, "SE-RNA/Q1_Both_up_threshold_058.txt"), warn = FALSE); q1_company <- unique(q1_company[nzchar(q1_company)])
stopifnot(length(q1_company) == 21, all(q1_company %in% M[quadrant == "Q1", gene]))
# 公司 141 行中的 Q1（23 行 / 21 基因）在本表中仍为 Q1
stopifnot(sum(K$quadrant == "Q1") == 23)

# ---------------------------------------------------------------- 3. 象限计数、相关
cnt     <- M[, .N, by = quadrant][order(quadrant)]
cnt_fdr <- M[RNA_FDR < FDR_CUT, .N, by = quadrant][order(quadrant)]
cnt_deg <- M[is_DEG == TRUE, .N, by = quadrant][order(quadrant)]
n_of <- function(tb, q) { v <- tb[quadrant == q, N]; if (length(v) == 0) 0L else v }
stopifnot(n_of(cnt, "Q1") == 56, n_of(cnt, "Q2") == 21, n_of(cnt, "Q3") == 90, n_of(cnt, "Q4") == 19, n_of(cnt, "NS") == 1456)   # ZFP36L1 对落在阈值带内：NS 1455 → 1456，四象限不变
stopifnot(M[gene == "ZFP36L1", quadrant] == "NS", abs(M[gene == "ZFP36L1", SE_log2FC] - 0.2594) < 1e-3)
ct  <- cor.test(M$SE_log2FC, M$RNA_log2FC, method = "spearman", exact = FALSE)
rho <- unname(ct$estimate); rho_p <- ct$p.value
stopifnot(rho > 0.3, rho < 0.45, rho_p < 1e-20)
bhr <- M[gene == "BHLHE40"]
stopifnot(nrow(bhr) == 2, all(bhr$quadrant == "Q1"), all(abs(sort(bhr$SE_log2FC) - c(1.039693, 1.079533)) < 1e-5))

# ---------------------------------------------------------------- 4. 标签集
M[, is_bh := gene == "BHLHE40"]
cand <- M[quadrant == "Q1" & RNA_FDR < FDR_CUT & SE_log2FC >= LAB_THR & RNA_log2FC >= LAB_THR][order(-SE_log2FC)][!duplicated(gene)]   # 同一基因多条 SE（BHLHE40、DUSP4）只取 SE_log2FC 最大的一条
stopifnot(nrow(cand) == 13, all(LAB_GENES %in% cand$gene))
lab_pool <- if (LABEL_MODE == "author") cand[gene %in% LAB_GENES] else M[quadrant == "Q1" & is_DEG == TRUE][order(-SE_log2FC)][!duplicated(gene)]
M[, lab := fifelse(paste(gene, EnhancerID_Con) %in% paste(lab_pool$gene, lab_pool$EnhancerID_Con), gene, "")]
stopifnot("BHLHE40" %in% M$lab, sum(M$lab == "BHLHE40") == 1)
if (LABEL_MODE == "author") stopifnot(setequal(M[nzchar(lab), lab], LAB_GENES), sum(nzchar(M$lab)) == length(LAB_GENES), all(M[nzchar(lab), quadrant] == "Q1"),
                                      abs(M[lab == "DUSP4", SE_log2FC] - 2.4134) < 1e-3)   # DUSP4 两条 SE（2.41 / 1.50，mRNA 同为 1.45）只标 SE log2FC 较大的一条
n_lab <- sum(nzchar(M$lab))
# 标签表：BHLHE40 两条 SE 合并为一个锚点（两点均值；两点 x 差 0.04，图上 0.3 mm），其余每基因一行；plotmath 斜体、BHLHE40 粗斜体
LABS <- rbind(M[nzchar(lab) & !is_bh, .(gene = lab, x = SE_log2FC, y = RNA_log2FC)],
              M[is_bh == TRUE, .(gene = "BHLHE40", x = mean(SE_log2FC), y = mean(RNA_log2FC))])
LABS[, `:=`(expr  = fifelse(gene == "BHLHE40", 'bolditalic("BHLHE40")', sprintf('italic("%s")', gene)),
            col   = fifelse(gene == "BHLHE40", BH_COL, "black"),
            psize = fifelse(gene == "BHLHE40", 2.4, 1.0),                    # 引线止于圈缘 / 红点缘
            nch   = fifelse(gene == "BHLHE40", 8L, nchar(gene)))]                # 粗体略宽，按 8 字符估计
POS <- list(TNIP3 = c(2.25, 3.80), IL1R2 = c(2.35, 3.30), BHLHE40 = c(-1.50, 2.90), IL1R1 = c(1.95, 2.25), DUSP4 = c(2.69, 1.90))
BH_NOTE <- c(-1.50, 2.55)                                            # "(SE1/SE2)" 小注位置：紧贴 BHLHE40 标签下方（该处无 Q1 点）
XL <- c(-3.1, 3.3); YL <- c(-4.3, 4.3)                               # 坐标范围（作图段沿用）
MM_PER_UNIT   <- (W_MM - 8.5) / diff(XL)                             # 面板宽 ≈ 出图宽 − y 轴标题 / 刻度 / 边距 8.5 mm → ≈ 6.5 mm / log2 单位
MM_PER_UNIT_Y <- (H_MM - 9.5) / diff(YL)                             # 面板高 ≈ 出图高 − x 轴标题 / 刻度 / 边距 9.5 mm → ≈ 5.4 mm / log2 单位（仅用于引线几何自检）
half_w <- function(nch) (1.18 * nch + 0.6) / MM_PER_UNIT / 2         # 6 pt Arial 斜体 ≈ 1.18 mm / 字符（0.56 em）+ 0.6 mm 内边距，换算成数据单位的半宽
ori <- function(ax, ay, bx, by, cx, cy) sign((bx - ax) * (cy - ay) - (by - ay) * (cx - ax))
seg_cross <- function(a, b, c, d) ori(a[1], a[2], b[1], b[2], c[1], c[2]) != ori(a[1], a[2], b[1], b[2], d[1], d[2]) &&
                                   ori(c[1], c[2], d[1], d[2], a[1], a[2]) != ori(c[1], c[2], d[1], d[2], b[1], b[2])
seg_hits_box <- function(a, b, bx) {                                   # 线段 a–b 是否碰到矩形 bx = c(x1, y1, x2, y2)（端点在框内或与任一边相交）
  inb <- function(p) p[1] >= bx[1] && p[1] <= bx[3] && p[2] >= bx[2] && p[2] <= bx[4]
  cn <- list(c(bx[1], bx[2]), c(bx[3], bx[2]), c(bx[3], bx[4]), c(bx[1], bx[4]))
  inb(a) || inb(b) || any(vapply(1:4, function(k) seg_cross(a, b, cn[[k]], cn[[k %% 4 + 1]]), logical(1)))
}
TARGET <- POS
LABS[, `:=`(nx = 0, ny = 0)]
for (g in names(TARGET)) { i <- match(g, LABS$gene); if (!is.na(i)) { LABS$nx[i] <- TARGET[[g]][1] - LABS$x[i]; LABS$ny[i] <- TARGET[[g]][2] - LABS$y[i] } }
if (LABEL_MODE == "author") {
  stopifnot(setequal(names(POS), LABS$gene))
  # ① 预置位自检：文字框（半高 0.14）不跨越四条 ±THR 虚线、不压住任何 Q1 红点（本基因锚点除外）、不出面板、右缘 ≤ 3.2、互不重叠
  LB <- LABS[, .(gene, cx = x + nx, cy = y + ny, hw = half_w(nch), hh = 0.14)]
  LB <- rbind(LB, data.table(gene = "BHLHE40_note", cx = BH_NOTE[1], cy = BH_NOTE[2], hw = half_w(9L), hh = 0.14))   # 小注同样参加自检
  q1 <- M[quadrant == "Q1", .(gene, x = SE_log2FC, y = RNA_log2FC)]
  for (k in seq_len(nrow(LB))) {
    inside <- q1[abs(x - LB$cx[k]) < LB$hw[k] + 0.05 & abs(y - LB$cy[k]) < LB$hh[k] + 0.06 & gene != sub("_note$", "", LB$gene[k])]
    if (nrow(inside)) stop(sprintf("标签 %s 压住 Q1 点：%s", LB$gene[k], paste(inside$gene, collapse = ", ")))
    for (t in c(-THR, THR)) if ((LB$cy[k] - LB$hh[k] < t & LB$cy[k] + LB$hh[k] > t) || (LB$cx[k] - LB$hw[k] < t & LB$cx[k] + LB$hw[k] > t))
      stop(sprintf("标签 %s 的文字框跨越 %.2f 阈值虚线", LB$gene[k], t))
    stopifnot(LB$cx[k] - LB$hw[k] > -3.05, LB$cx[k] + LB$hw[k] <= 3.2, LB$cy[k] + LB$hh[k] < 4.2)   # 右缘 ≤ 3.2（距面板右缘 ≥ 0.1 单位 ≈ 0.65 mm）
  }
  for (i in seq_len(nrow(LB) - 1)) for (j in (i + 1):nrow(LB))
    if (!(abs(LB$cx[i] - LB$cx[j]) >= LB$hw[i] + LB$hw[j] || abs(LB$cy[i] - LB$cy[j]) >= LB$hh[i] + LB$hh[j] + 0.04)) stop(sprintf("标签 %s 与 %s 重叠", LB$gene[i], LB$gene[j]))
  n_grey_under <- sum(sapply(seq_len(nrow(LB)), function(k) M[quadrant != "Q1" & abs(SE_log2FC - LB$cx[k]) < LB$hw[k] & abs(RNA_log2FC - LB$cy[k]) < LB$hh[k], .N]))
  mmx <- function(v) v * MM_PER_UNIT; mmy <- function(v) v * MM_PER_UNIT_Y
  LB[, `:=`(bx1 = mmx(cx - hw), by1 = mmy(cy - hh), bx2 = mmx(cx + hw), by2 = mmy(cy + hh))]
  LD <- LABS[, .(gene, px = mmx(x), py = mmy(y), prad = psize / 2)]
  LD <- cbind(LD, t(sapply(seq_len(nrow(LD)), function(i) { b <- LB[gene == LD$gene[i]]
    ggrepel:::select_line_connection(c(LD$px[i], LD$py[i]), c(b$bx1, b$by1, b$bx2, b$by2)) })))
  setnames(LD, c("V1", "V2"), c("sx", "sy"))
  LD[, len_mm := sqrt((sx - px)^2 + (sy - py)^2) - prad]
  n_cross <- 0L; n_box <- 0L
  for (i in seq_len(nrow(LD) - 1)) for (j in (i + 1):nrow(LD))
    if (seg_cross(c(LD$px[i], LD$py[i]), c(LD$sx[i], LD$sy[i]), c(LD$px[j], LD$py[j]), c(LD$sx[j], LD$sy[j]))) { n_cross <- n_cross + 1L; cat("引线交叉：", LD$gene[i], "×", LD$gene[j], "\n") }
  for (i in seq_len(nrow(LD))) for (k in which(LB$gene != LD$gene[i] & !(LD$gene[i] == "BHLHE40" & LB$gene == "BHLHE40_note")))
    if (seg_hits_box(c(LD$px[i], LD$py[i]), c(LD$sx[i], LD$sy[i]), unlist(LB[k, .(bx1, by1, bx2, by2)]))) { n_box <- n_box + 1L; cat("引线穿过标签框：", LD$gene[i], "→", LB$gene[k], "\n") }
  stopifnot(n_cross == 0, n_box == 0, all(LD$len_mm >= 0.8))
  # BHLHE40 引线另查不穿过自己的小注（小注在标签正下方，引线从标签右缘出发）
  bn <- LB[gene == "BHLHE40_note"]; bl <- LD[gene == "BHLHE40"]
  stopifnot(!seg_hits_box(c(bl$px, bl$py), c(bl$sx, bl$sy), c(bn$bx1, bn$by1, bn$bx2, bn$by2)))
  cat(sprintf("标签预置位自检通过：%d 个标签框 + 1 个小注不跨四条 ±0.58 阈值线、不压 Q1 红点、互不重叠；被标签覆盖的灰点 %d 个（可接受）\n", nrow(LABS), n_grey_under))
  cat(sprintf("引线自检通过：%d 条引线两两不交叉、不穿过其他标签框；可见长度 %s mm\n", nrow(LD), paste(sprintf("%s %.1f", LD$gene, LD$len_mm), collapse = " / ")))
  print(LB[, .(gene, cx = round(cx, 2), cy = round(cy, 2))])
  LEADER_CHECK <- sprintf("%d 条引线两两不交叉、不穿过其他标签框 / 小注；标签框不跨四条 ±0.58 阈值虚线、不压 Q1 红点、互不重叠；引线可见长度 %s mm（ggrepel select_line_connection 规则、毫米坐标重建）",
                          nrow(LD), paste(sprintf("%s %.1f", LD$gene, LD$len_mm), collapse = " / "))
  LABEL_POS <- paste(LB[, sprintf("%s (%.2f, %.2f)", gene, cx, cy)], collapse = "; ")
}

# ---------------------------------------------------------------- 5. 作图（50 × 56 mm；XL / YL 在第 4 节已定义）
stopifnot(all(M$SE_log2FC > XL[1] & M$SE_log2FC < XL[2]), all(M$RNA_log2FC > YL[1] & M$RNA_log2FC < YL[2]))
corner <- function(q, x, y, hj, vj, col = "grey30") annotate("text", x = x, y = y, hjust = hj, vjust = vj, label = sprintf("italic(n) == %d", n_of(cnt, q)),
                                                             parse = TRUE, size = pt2size(FS_SMALL), family = FONT, colour = col)
bh_pt <- LABS[gene == "BHLHE40"]
p <- ggplot(M, aes(SE_log2FC, RNA_log2FC)) +
  geom_hline(yintercept = c(-THR, THR), linetype = "22", linewidth = lw(0.4), colour = "grey55") +
  geom_vline(xintercept = c(-THR, THR), linetype = "22", linewidth = lw(0.4), colour = "grey55") +
  geom_point(data = M[quadrant != "Q1"], colour = "#C6C6C6", size = 0.65, stroke = 0) +
  geom_point(data = M[quadrant == "Q1"], colour = grp_col[["Mhp"]], size = 1.0, stroke = 0) +
  geom_point(data = bh_pt, aes(x, y), shape = 21, fill = NA, colour = BH_COL, size = 2.4, stroke = lw(0.6)) +   # 一个圈套住两条 SE
  annotate("text", x = BH_NOTE[1], y = BH_NOTE[2], label = "(SE1/SE2)", colour = BH_COL, size = pt2size(FS_SMALL), family = FONT) +
  geom_text_repel(data = LABS, aes(x, y, label = expr, colour = I(col), point.size = psize), parse = TRUE,
                  size = pt2size(FS_SMALL), family = FONT, min.segment.length = 0, segment.size = lw(0.3), segment.colour = "grey45",
                  nudge_x = LABS$nx, nudge_y = LABS$ny,
                  # author 模式：5 个标签全部预置位，force = 0 使标签固定在预置位（互不重叠 / 引线不交叉由第 4 节自检保证），ggrepel 只负责从标签框边缘画引线到点缘
                  box.padding = 0.15, point.padding = 0.02, force = if (LABEL_MODE == "author") 0 else 1, force_pull = if (LABEL_MODE == "author") 0 else 0.3,
                  max.overlaps = Inf, xlim = c(XL[1] + 0.05, 3.25), ylim = c(YL[1] + 0.05, 4.2), max.time = 5, max.iter = 1e5, seed = SEED, show.legend = FALSE) +
  corner("Q4", -Inf, Inf, -0.1, 1.4) + corner("Q3", -Inf, -Inf, -0.1, -0.5) +                        # 左侧两象限：外角
  corner("Q1", THR + 0.12, YL[2] - 0.05, 0, 1, grp_col[["Mhp"]]) + corner("Q2", THR + 0.12, YL[1] + 0.05, 0, 0) +   # 右侧两象限：阈值竖线右侧内角
  annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.6, parse = TRUE, size = pt2size(FS_SMALL), family = FONT, colour = "grey30",
           label = sprintf("rho == %.2f", rho)) +                                                     # ρ 单独在右下角，不与计数叠行
  scale_x_continuous(limits = XL, breaks = -3:3, expand = c(0, 0)) +
  scale_y_continuous(limits = YL, breaks = seq(-4, 4, 2), expand = c(0, 0)) +
  labs(x = expression("SE H3K27ac "*log[2]*"FC (Mhp / Con)"), y = expression("mRNA "*log[2]*"FC (Mhp / Con)")) +
  theme_mhp() + theme(legend.position = "none", plot.margin = margin(1.5, 2, 1, 1, "mm"))
stem <- if (LABEL_MODE == "author") "Fig4e_SE_RNA四象限" else "Fig4e_SE_RNA四象限_全标签"
out <- sav(p, stem, PAN, W_MM, H_MM, panel_dir = if (LABEL_MODE == "author") PNL else NULL)
if (identical(Sys.getenv("MHP_FINAL"), "1") && LABEL_MODE == "author") ggsave(file.path(PAN, paste0(stem, ".pdf")), p, width = W_MM, height = H_MM, units = "mm", device = .pdf_dev)
stopifnot(file.size(out) > 20000)

# ---------------------------------------------------------------- 6. 输出表
setcolorder(M, c("gene", "gene_side", "quadrant", "is_DEG", "SE_log2FC", "RNA_log2FC", "RNA_FDR", "RNA_P", "RNA_dir",
                 "EnhancerID_Con", "EnhancerID_Mhp", "ActivityScore_Con", "ActivityScore_Mhp", "OverlapLength", "Gene_Con_clean", "Gene_Mhp_clean", "lab"))
OUTM <- M[order(quadrant, -RNA_log2FC), .(gene, gene_side, quadrant, is_DEG, SE_log2FC = round(SE_log2FC, 4), RNA_log2FC = round(RNA_log2FC, 4),
                                          RNA_FDR = signif(RNA_FDR, 3), RNA_P = signif(RNA_P, 3), RNA_dir, EnhancerID_Con, EnhancerID_Mhp,
                                          ActivityScore_Con = round(ActivityScore_Con, 2), ActivityScore_Mhp = round(ActivityScore_Mhp, 2),
                                          OverlapLength, Gene_Con_clean, Gene_Mhp_clean, labelled = nzchar(lab))]
if (LABEL_MODE == "author") fwrite(OUTM, file.path(TAB, "03_e_四象限_全部配对.csv"))
key <- rbindlist(list(
  data.table(item = "配对 SE 总数（≥1 kb 重叠；00b 本地重做，含公司表漏掉的 ZFP36L1 对）", value = nrow(SE)),
  data.table(item = "ZFP36L1 对（公司表漏对）在本图的位置", value = sprintf("SE_log2FC %.3f / mRNA log2FC %.3f → 阈值带内灰点（不计入象限）", M[gene == "ZFP36L1", SE_log2FC], M[gene == "ZFP36L1", RNA_log2FC])),
  data.table(item = "旧口径（公司 1,679 对）对照", value = "1,641 点 / 1,545 基因；Q1 56 / Q2 21 / Q3 90 / Q4 19 / 阈值带内 1,455；ρ 0.367（P 1.75e-53）"),
  data.table(item = "两侧最近基因都不在 RNA 表而弃用的配对", value = n_drop),
  data.table(item = "画出的点数（配对 SE）", value = nrow(M)),
  data.table(item = "其中 Con 侧最近基因命中 / 改用 Mhp 侧", value = sprintf("%d / %d", sum(M$gene_side == "Con"), sum(M$gene_side == "Mhp"))),
  data.table(item = "去重基因数", value = uniqueN(M$gene)),
  data.table(item = "象限阈值 |log2FC|", value = THR),
  data.table(item = "Q1 SE↑mRNA↑ 配对数 / 基因数", value = sprintf("%d / %d", n_of(cnt, "Q1"), uniqueN(M[quadrant == "Q1", gene]))),
  data.table(item = "Q2 SE↑mRNA↓ 配对数 / 基因数", value = sprintf("%d / %d", n_of(cnt, "Q2"), uniqueN(M[quadrant == "Q2", gene]))),
  data.table(item = "Q3 SE↓mRNA↓ 配对数 / 基因数", value = sprintf("%d / %d", n_of(cnt, "Q3"), uniqueN(M[quadrant == "Q3", gene]))),
  data.table(item = "Q4 SE↓mRNA↑ 配对数 / 基因数", value = sprintf("%d / %d", n_of(cnt, "Q4"), uniqueN(M[quadrant == "Q4", gene]))),
  data.table(item = "阈值带内（不计象限）配对数", value = n_of(cnt, "NS")),
  data.table(item = "Q1–Q4 配对数（RNA FDR < 0.05 过滤后）", value = paste(sprintf("%s %d", c("Q1", "Q2", "Q3", "Q4"), sapply(c("Q1", "Q2", "Q3", "Q4"), n_of, tb = cnt_fdr)), collapse = "; ")),
  data.table(item = "Q1–Q4 配对数（手稿 DEG 口径 |log2FC|≥1 & FDR<0.05）", value = paste(sprintf("%s %d", c("Q1", "Q2", "Q3", "Q4"), sapply(c("Q1", "Q2", "Q3", "Q4"), n_of, tb = cnt_deg)), collapse = "; ")),
  data.table(item = "Q1 基因数（RNA FDR < 0.05）", value = uniqueN(M[quadrant == "Q1" & RNA_FDR < FDR_CUT, gene])),
  data.table(item = "Q1 基因数（手稿 DEG 口径，= 补充热图 05 的行数）", value = uniqueN(M[quadrant == "Q1" & is_DEG == TRUE, gene])),
  data.table(item = "Q1 中 SE 与 mRNA log2FC 均 ≥ 1 且 FDR < 0.05 的基因（标签候选池）", value = paste(cand[order(-RNA_log2FC), gene], collapse = ";")),
  data.table(item = "图上标注的基因（按 mRNA log2FC 降序）", value = paste(lab_pool[order(-RNA_log2FC), gene], collapse = ";")),
  data.table(item = "图上标注的基因数", value = n_lab),
  data.table(item = "标注基因 SE log2FC / mRNA log2FC（DUSP4 取 SE log2FC 较大的一条；BHLHE40 为两条 SE）",
             value = paste(M[nzchar(lab) | is_bh][order(-RNA_log2FC, -SE_log2FC), sprintf("%s %.2f / %.2f", gene, SE_log2FC, RNA_log2FC)], collapse = "; ")),
  data.table(item = "标签中心位置（数据坐标）", value = if (LABEL_MODE == "author") LABEL_POS else "deg 模式：ggrepel 自动"),
  data.table(item = "标签 / 引线自检", value = if (LABEL_MODE == "author") LEADER_CHECK else "deg 模式未做"),
  data.table(item = "Spearman rho（全部点）", value = round(rho, 4)),
  data.table(item = "Spearman P", value = signif(rho_p, 3)),
  data.table(item = "BHLHE40 SE1 (88_peak_78108) SE_log2FC", value = round(bhr[EnhancerID_Con == "88_peak_78108_lociStitched", SE_log2FC], 4)),
  data.table(item = "BHLHE40 SE2 (19_peak_78161) SE_log2FC", value = round(bhr[EnhancerID_Con == "19_peak_78161_lociStitched", SE_log2FC], 4)),
  data.table(item = "BHLHE40 mRNA log2FC / FDR", value = sprintf("%.4f / %.3g", bh$logFC, bh$FDR)),
  data.table(item = "Q1 中 RNA FDR ≥ 0.05 仍着红并计入 n 的点", value = paste(M[quadrant == "Q1" & RNA_FDR >= FDR_CUT][order(-RNA_log2FC), sprintf("%s(FDR %.2f)", gene, RNA_FDR)], collapse = ";")),
  data.table(item = "Mhp SE 重复配对：出现 ≥ 2 次的 Mhp SE 数 / 重复配对数", value = sprintf("%d / %d", M[, .N, by = EnhancerID_Mhp][N >= 2, .N], nrow(M) - uniqueN(M$EnhancerID_Mhp))),
  data.table(item = "图注句（计数口径）", value = sprintf("n = pairs counted by |log2FC| >= %.2f alone (no FDR filter); with RNA FDR < 0.05, Q1 %d / Q2 %d / Q3 %d / Q4 %d", THR,
                                                     n_of(cnt_fdr, "Q1"), n_of(cnt_fdr, "Q2"), n_of(cnt_fdr, "Q3"), n_of(cnt_fdr, "Q4"))),
  data.table(item = "图注句（点的定义）", value = sprintf("Each point is one Con-Mhp SE pair (%d pairs; Con SE unique, %d Mhp SEs paired with >= 2 Con SEs); BHLHE40 circle encloses both of its SEs (SE1 / SE2)", nrow(M), M[, .N, by = EnhancerID_Mhp][N >= 2, .N])),
  data.table(item = "BHLHE40 高亮色（c/d/e 统一）", value = BH_COL),
  data.table(item = "公司 141 配对 DEG 在本表中逐行一致", value = sprintf("%d / 141", nrow(K))),
  data.table(item = "公司 0.58 阈值 Q1 清单（21 基因）全部落在本图 Q1", value = all(q1_company %in% M[quadrant == "Q1", gene])),
  data.table(item = "标签模式", value = LABEL_MODE)))
if (LABEL_MODE == "author") fwrite(key, file.path(TAB, "03_e_关键数字.csv"))
cat("\n===== 面板 e 关键数字 =====\n"); print(key, nrows = 100)
cat("\nQ1（DEG 口径）基因按 mRNA log2FC 降序:\n")
print(M[quadrant == "Q1" & is_DEG == TRUE][order(-RNA_log2FC)][!duplicated(gene), .(gene, SE_log2FC = round(SE_log2FC, 2), RNA_log2FC = round(RNA_log2FC, 2), RNA_FDR = signif(RNA_FDR, 2), labelled = nzchar(lab))], nrows = 100)
if (LABEL_MODE == "author") stopifnot(file.size(file.path(TAB, "03_e_四象限_全部配对.csv")) > 0, file.size(file.path(TAB, "03_e_关键数字.csv")) > 0)
