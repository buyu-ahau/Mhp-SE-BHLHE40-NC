# ============================================================================
# 00b_SE配对_本地重做.R — Fig. 4d、4e 的输入：Con 与 Mhp 两组 ROSE 超级增强子按区间重叠配对
# 运行位置：本地；路径变量见 config/paths.R
# ============================================================================

source(file.path(Sys.getenv("MHP_CODE", "."), "config", "paths.R"))

suppressPackageStartupMessages({ library(data.table) })
source(file.path(MHP_CODE, "shared", "_palette_MHP.R"))   # MHP_RES、grp_levels

SE_DIR <- MHP_ROSE
TAB    <- file.path(MHP_RES, "6、SE/新Fig4_作图/Tables"); dir.create(TAB, recursive = TRUE, showWarnings = FALSE)
MIN_OV <- 1000L                                                                # 公司规则：重叠 ≥ 1,000 bp
DIR_THR <- 0.5                                                                 # 公司 Change_Direction 阈值（由公司表反推，见头注）
COMPANY_COLS <- c("EnhancerID_Con", "CLOSEST_GENE_Con", "ActivityScore_Con", "EnhancerID_Mhp", "CLOSEST_GENE_Mhp", "ActivityScore_Mhp",
                  "ScoreDifference_Con_minus_Mhp", "AbsoluteScoreDifference", "OverlapLength", "SE_log2FC", "SE_log2FC_exact", "Change_Direction",
                  "Gene_Con_clean", "Gene_Mhp_clean")

# ---------------------------------------------------------------- ① 读 ROSE SE 表 + 最近基因（与 01_c 同法）
ROSE_COLS <- c("REGION_ID", "CHROM", "START", "STOP", "NUM_LOCI", "CONSTITUENT_SIZE", "signal", "enhancerRank", "isSuper")
E2G_COLS  <- c(ROSE_COLS, "OVERLAP_GENES", "PROXIMAL_GENES", "CLOSEST_GENE", "enhancerRank_2", "isSuper_2")
read_rose <- function(g) {
  f <- file.path(SE_DIR, paste0(g, "_SuperEnhancers.table.txt")); hd <- readLines(f, n = 6)
  stopifnot(grepl(paste0("^#", g, " Enhancers"), hd[1]), grepl("^REGION_ID\t", hd[6]))       # 第 6 行才是列名：第 7 行（rank 1）是数据，不能被当作表头
  se <- fread(f, skip = "REGION_ID", sep = "\t"); setnames(se, ROSE_COLS)
  e2g <- fread(file.path(SE_DIR, paste0(g, "_SuperEnhancers_ENHANCER_TO_GENE.txt")), skip = 1, header = FALSE, sep = "\t", quote = "",
               colClasses = list(character = c(1, 2, 10, 11, 12)))
  stopifnot(ncol(e2g) == 14); setnames(e2g, E2G_COLS)
  stopifnot(setequal(se$REGION_ID, e2g$REGION_ID), !anyDuplicated(se$REGION_ID), all(se$isSuper == 1L), all(se$STOP > se$START),
            identical(se$enhancerRank, seq_len(nrow(se))), !is.unsorted(-se$signal))                          # rank = 1..n 且按信号降序
  se <- merge(se, e2g[, .(REGION_ID, CLOSEST_GENE)], by = "REGION_ID")[order(enhancerRank)]
  stopifnot(!anyNA(se$CLOSEST_GENE), all(nzchar(se$CLOSEST_GENE)))
  se[, group := g][]
}
SE <- lapply(setNames(grp_levels, grp_levels), read_rose)
stopifnot(nrow(SE$Con) == 2307L, nrow(SE$Mhp) == 2422L,
          SE$Con[enhancerRank == 1, REGION_ID] == "595_peak_287089_lociStitched", SE$Mhp[enhancerRank == 1, REGION_ID] == "568_peak_269553_lociStitched")
n_scaf <- sapply(SE, function(x) sum(grepl("^NW_", x$CHROM)))
cat(sprintf("ROSE SE：Con %d（scaffold %d）/ Mhp %d（scaffold %d）；两组 rank 1 均为 ZFP36L1 位点 chr7:91.60–92.81 Mb\n",
            nrow(SE$Con), n_scaf[["Con"]], nrow(SE$Mhp), n_scaf[["Mhp"]]))

# ---------------------------------------------------------------- ② 配对函数：每个 Con SE → 重叠最长的 Mhp SE（≥ MIN_OV）
pair_se <- function(con, mhp, min_ov = MIN_OV) {
  cand <- merge(con[, .(EnhancerID_Con = REGION_ID, CHROM, Con_START = START, Con_STOP = STOP, Con_signal = signal, Con_rank = enhancerRank, Con_gene = CLOSEST_GENE)],
                mhp[, .(EnhancerID_Mhp = REGION_ID, CHROM, Mhp_START = START, Mhp_STOP = STOP, Mhp_signal = signal, Mhp_rank = enhancerRank, Mhp_gene = CLOSEST_GENE)],
                by = "CHROM", allow.cartesian = TRUE)                                                           # 同染色体（含 NW_ scaffold 字面匹配）全部组合
  cand[, OverlapLength := pmin(Con_STOP, Mhp_STOP) - pmax(Con_START, Mhp_START)]                              # BED 半开区间：min(STOP) − max(START)
  cand <- cand[OverlapLength >= min_ov]
  cand[, n_Mhp_candidates := .N, by = EnhancerID_Con]
  cand[, tie_longest := sum(OverlapLength == max(OverlapLength)) > 1L, by = EnhancerID_Con]                    # 最长重叠并列（→ 取信号更高者）
  setorder(cand, EnhancerID_Con, -OverlapLength, -Mhp_signal, Mhp_rank)
  cand[, .SD[1L], by = EnhancerID_Con]
}
recompute <- function(p) {                                                                                       # 公司 14 列复算
  p[, `:=`(CLOSEST_GENE_Con = Con_gene, ActivityScore_Con = Con_signal, CLOSEST_GENE_Mhp = Mhp_gene, ActivityScore_Mhp = Mhp_signal)]
  p[, `:=`(ScoreDifference_Con_minus_Mhp = ActivityScore_Con - ActivityScore_Mhp)]
  p[, `:=`(AbsoluteScoreDifference = abs(ScoreDifference_Con_minus_Mhp),
           SE_log2FC = log2((ActivityScore_Mhp + 1) / (ActivityScore_Con + 1)), SE_log2FC_exact = log2(ActivityScore_Mhp / ActivityScore_Con))]
  p[, Change_Direction := fifelse(SE_log2FC >= DIR_THR, "Mhp_High", fifelse(SE_log2FC <= -DIR_THR, "Con_High", "No_Significant_Change"))]
  p[, `:=`(Gene_Con_clean = sub("^gene-", "", CLOSEST_GENE_Con), Gene_Mhp_clean = sub("^gene-", "", CLOSEST_GENE_Mhp))]
  setorder(p, -AbsoluteScoreDifference, EnhancerID_Con)                                                         # 与公司表同序：|差值| 降序
  p[, rank_absdiff := .I][]
}
P <- recompute(pair_se(SE$Con, SE$Mhp))
stopifnot(nrow(P) == 1680L, !anyDuplicated(P$EnhancerID_Con), min(P$OverlapLength) >= MIN_OV, all(P$ActivityScore_Con > 0), all(P$ActivityScore_Mhp > 0),
          !any(P$tie_longest))                                                                                   # 无最长重叠并列

# ---------------------------------------------------------------- ③ 与公司表逐行核对：1,679 对完全一致，多出的 1 对 = ZFP36L1
CO <- fread(file.path(SE_DIR, "DE se/SE_log2FC_all_results.csv"))
L  <- readLines(file.path(SE_DIR, "overlap_score_diff_results.txt"), encoding = "UTF-8")
hdr <- grep("EnhancerID_Con", L)[1]; body <- L[hdr:length(L)]; body <- body[nzchar(trimws(body))]
OV <- as.data.table(read.table(text = body, header = TRUE, stringsAsFactors = FALSE))
stopifnot(nrow(CO) == 1679L, nrow(OV) == 1679L, identical(names(CO), COMPANY_COLS), identical(CO$EnhancerID_Con, OV$EnhancerID_Con),
          identical(CO$EnhancerID_Mhp, OV$EnhancerID_Mhp), max(abs(CO$OverlapLength - OV$OverlapLength)) == 0)
P[, in_company_table := EnhancerID_Con %in% CO$EnhancerID_Con]
extra <- P[in_company_table == FALSE]
stopifnot(nrow(extra) == 1L, extra$EnhancerID_Con == "595_peak_287089_lociStitched", extra$EnhancerID_Mhp == "568_peak_269553_lociStitched",
          extra$Gene_Con_clean == "ZFP36L1", extra$Gene_Mhp_clean == "ZFP36L1", extra$OverlapLength == 1201343L, extra$Con_rank == 1L, extra$Mhp_rank == 1L)
K <- merge(CO, P[in_company_table == TRUE, c(COMPANY_COLS, "rank_absdiff"), with = FALSE], by = "EnhancerID_Con", suffixes = c("", ".new"))
stopifnot(nrow(K) == 1679L, identical(K$EnhancerID_Mhp, K$EnhancerID_Mhp.new), identical(K$CLOSEST_GENE_Con, K$CLOSEST_GENE_Con.new),
          identical(K$CLOSEST_GENE_Mhp, K$CLOSEST_GENE_Mhp.new), identical(K$Change_Direction, K$Change_Direction.new),
          identical(K$Gene_Con_clean, K$Gene_Con_clean.new), identical(K$Gene_Mhp_clean, K$Gene_Mhp_clean.new), identical(K$OverlapLength, K$OverlapLength.new))
num_cols <- c("ActivityScore_Con", "ActivityScore_Mhp", "ScoreDifference_Con_minus_Mhp", "AbsoluteScoreDifference", "SE_log2FC", "SE_log2FC_exact")
max_dev <- sapply(num_cols, function(v) max(abs(K[[v]] - K[[paste0(v, ".new")]])))
stopifnot(all(max_dev < 1e-6))                                                                                    # 数值列逐对一致（公司 csv 保留 4 位小数以上）
# 行序：去掉 ZFP36L1 对后，本表顺序须与公司表完全相同（公司表按 |差值| 降序）
stopifnot(identical(P[in_company_table == TRUE, EnhancerID_Con], CO$EnhancerID_Con))
# 漏对原因验证：两表各去掉第 1 行数据（rank 1）再配对 → 与公司表完全一致？
P_drop1 <- recompute(pair_se(SE$Con[enhancerRank > 1], SE$Mhp[enhancerRank > 1]))
drop1_identical <- nrow(P_drop1) == 1679L && identical(P_drop1$EnhancerID_Con, CO$EnhancerID_Con) && identical(P_drop1$EnhancerID_Mhp, CO$EnhancerID_Mhp)
cat(sprintf("本地重做 %d 对；其中 %d 对与公司表逐行逐值一致（数值列最大偏差 %.1e）；多出 1 对 = %s ↔ %s（%s，重叠 %s bp，Con %s / Mhp %s，log2FC %.3f）\n",
            nrow(P), nrow(K), max(max_dev), extra$EnhancerID_Con, extra$EnhancerID_Mhp, extra$Gene_Con_clean, format(extra$OverlapLength, big.mark = ","),
            format(round(extra$ActivityScore_Con), big.mark = ","), format(round(extra$ActivityScore_Mhp), big.mark = ","), extra$SE_log2FC))
cat(sprintf("漏对原因验证：两张 ROSE 表各去掉 rank 1 行再配对 → 与公司 1,679 对完全一致 = %s（支持『读表时把第 1 行数据当作列名』的解释）\n", drop1_identical))

# ---------------------------------------------------------------- ④ 其他统计（供关键数字表 / 图注）
n_dup_mhp <- sum(duplicated(P$EnhancerID_Mhp)); n_scaf_pair <- sum(grepl("^NW_", P$CHROM))
n_multi <- sum(P$n_Mhp_candidates > 1L)
b1 <- P[EnhancerID_Con == "88_peak_78108_lociStitched"]; b2 <- P[EnhancerID_Con == "19_peak_78161_lociStitched"]
stopifnot(nrow(b1) == 1L, nrow(b2) == 1L, b1$Gene_Con_clean == "BHLHE40", b2$Gene_Con_clean == "BHLHE40")
rho <- cor(log10(P$ActivityScore_Con), log10(P$ActivityScore_Mhp), method = "spearman")
rho_old <- cor(log10(CO$ActivityScore_Con), log10(CO$ActivityScore_Mhp), method = "spearman")
cat(sprintf("Mhp SE 重复配对 %d 个；scaffold 配对 %d 对；有 ≥ 2 个候选 Mhp SE 的 Con SE %d 个（均取最长重叠，无并列）；Spearman ρ %.4f（公司 1,679 对 %.4f）\n",
            n_dup_mhp, n_scaf_pair, n_multi, rho, rho_old))
cat(sprintf("BHLHE40 SE1 按 |差值| 名次 %d（公司表 %d）、SE2 %d（公司表 %d）\n", b1$rank_absdiff, match("88_peak_78108_lociStitched", CO$EnhancerID_Con),
            b2$rank_absdiff, match("19_peak_78161_lociStitched", CO$EnhancerID_Con)))

# ---------------------------------------------------------------- ⑤ 输出
OUT <- P[, c(COMPANY_COLS, "CHROM", "Con_START", "Con_STOP", "Con_rank", "Mhp_START", "Mhp_STOP", "Mhp_rank", "n_Mhp_candidates", "tie_longest",
             "rank_absdiff", "in_company_table"), with = FALSE]
fwrite(OUT, file.path(TAB, "00b_SE配对_1680.csv"))
key <- data.table(item = c("配对规则", "重叠长度约定", "并列处理", "scaffold 处理", "本地配对对数", "与公司表逐行逐值一致的对数", "数值列最大偏差",
                           "多出的 1 对", "多出 1 对的信号 Con / Mhp / log2FC（+1 伪计数）/ Change_Direction", "漏对原因验证（两表各去掉 rank 1 行再配对 = 公司表）",
                           "Mhp SE 重复配对数（一个 Mhp SE 对多个 Con SE）", "去重 Mhp SE / Con SE", "scaffold（NW_*）配对数", "有 ≥ 2 个候选 Mhp SE 的 Con SE 数", "最长重叠并列数",
                           "Change_Direction 阈值（由公司表反推）", "Change_Direction 计数 Mhp_High / Con_High / No_Significant_Change",
                           "Spearman ρ（log10 信号，1,680 对）", "Spearman ρ（公司 1,679 对）", "BHLHE40 SE1 按 |差值| 名次（1,680 对 / 公司 1,679 对）",
                           "BHLHE40 SE2 按 |差值| 名次（1,680 对 / 公司 1,679 对）", "行序"),
                  value = c(sprintf("每个 Con SE → 同 CHROM 上重叠 ≥ %d bp 的 Mhp SE 中取重叠最长者（公司说明同）", MIN_OV), "min(STOP) − max(START)（BED 半开区间；与公司 OverlapLength 逐对一致）",
                            "并列时取 Mhp 信号更高者（enhancerRank 更小）；本数据无并列", "NW_* 按 CHROM 字面匹配、参与配对（公司表含 22 对）", nrow(P), nrow(K), sprintf("%.1e", max(max_dev)),
                            sprintf("%s ↔ %s（%s，重叠 %s bp；两组 ROSE rank 1）", extra$EnhancerID_Con, extra$EnhancerID_Mhp, extra$Gene_Con_clean, format(extra$OverlapLength, big.mark = ",")),
                            sprintf("%.1f / %.1f / %.4f / %s", extra$ActivityScore_Con, extra$ActivityScore_Mhp, extra$SE_log2FC, extra$Change_Direction), drop1_identical,
                            n_dup_mhp, sprintf("%d / %d", uniqueN(P$EnhancerID_Mhp), uniqueN(P$EnhancerID_Con)), n_scaf_pair, n_multi, sum(P$tie_longest),
                            DIR_THR, sprintf("%d / %d / %d", sum(P$Change_Direction == "Mhp_High"), sum(P$Change_Direction == "Con_High"), sum(P$Change_Direction == "No_Significant_Change")),
                            round(rho, 4), round(rho_old, 4), sprintf("%d / %d", b1$rank_absdiff, match("88_peak_78108_lociStitched", CO$EnhancerID_Con)),
                            sprintf("%d / %d", b2$rank_absdiff, match("19_peak_78161_lociStitched", CO$EnhancerID_Con)), "AbsoluteScoreDifference 降序（= 公司表序；rank_absdiff = 行号）"))
fwrite(key, file.path(TAB, "00b_SE配对_关键数字.csv"))
stopifnot(file.size(file.path(TAB, "00b_SE配对_1680.csv")) > 1e5, file.size(file.path(TAB, "00b_SE配对_关键数字.csv")) > 0)
cat("\n===== 00b 关键数字 =====\n"); print(key, nrows = 100)
cat(sprintf("\n✓ Tables/00b_SE配对_1680.csv（%d 行 × %d 列）、00b_SE配对_关键数字.csv\n", nrow(OUT), ncol(OUT)))
