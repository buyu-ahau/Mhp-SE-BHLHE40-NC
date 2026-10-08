# ============================================================================
# 01_a_loop分类_本地.R — Fig. 5a：Hi-C loop 按两端锚点所在元件分类
# 运行位置：本地；路径变量见 config/paths.R
# ============================================================================

source(file.path(Sys.getenv("MHP_CODE", "."), "config", "paths.R"))

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(GenomicRanges); library(GenomicFeatures); library(AnnotationDbi)
})
source(file.path(MHP_CODE, "shared", "_palette_MHP.R"))   # theme_mhp、sav、lw、pt2size、FONT、grp_levels、MHP_RES、.pdf_dev

SE_DIR   <- MHP_ROSE
TE_DIR   <- MHP_TE
HIC_DIR  <- file.path(MHP_RES, "7、HI-C结果")
LOOP_DIR <- file.path(MHP_HIC, "single_result/05.loop")
TXDB_IN  <- file.path(HIC_DIR, "Sus_scrofa.Sscrofa11.1.114.chr_txdb.sqlite")
OLD_IN   <- file.path(HIC_DIR, "HI-C联合SE/Fig4C_Loop_Categories_Data_Fixed.txt")
OUT <- file.path(HIC_DIR, "新Fig5_作图")
PAN <- file.path(OUT, "Plots/面板"); TAB <- file.path(OUT, "Tables"); PNL <- file.path(OUT, "_panels")
for (d in c(PAN, TAB, PNL)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
STEM <- "Fig5a_loop两端元件分类"

LOOP_PREFIX <- c(Con = "_Con", Mhp = "Mhp")          # bedpe 文件名前缀（公司命名：Con 组带下划线）
N_LOOP_REF  <- c(Con = 22441L, Mhp = 19203L)         # loop_stat.txt / 探针记录值
N_SE_REF    <- c(Con = 2307L,  Mhp = 2422L)          # Fig 4c 图上 n = 2,307 / 2,422
N_TE_REF    <- c(Con = 101894L, Mhp = 102376L)
PROM_FLANK  <- 2500L                                 # 启动子 = TSS ± 2.5 kb
MAIN_CHR    <- paste0("chr", 1:19)                   # Hi-C 有 loop 鉴定的主染色体（chr19 = X）
LEN_CHR18   <- 55982971; LEN_X <- 125939595
ELEM_RANK   <- c("SE", "TE", "P", "Other")           # 类别名内两端的书写顺序（SE–TE 而非 TE–SE）
CAT_LEVELS  <- c("SE-P", "SE-SE", "SE-TE", "SE-Other", "TE-P", "TE-TE", "TE-Other", "P-P", "P-Other", "Other-Other")
SCHEMES     <- list(main = c("SE", "P", "TE"), altA = c("SE", "TE", "P"), altB = c("P", "SE", "TE"))   # 锚点标签优先级（高 → 低）
LAB_MIN_PCT <- 5                                     # 段内标百分比的最小段高（%）
SE_PCT_SCHEME <- "main"                              # 图上 "有 ≥ 1 启动子 loop 的 SE 占比" 的口径：main = 严格（SE–P 类 loop，36.9 / 37.3%）；
stopifnot(SE_PCT_SCHEME %in% c("main", "relaxed"))

# ---------------------------------------------------------------- 读入工具
bed2gr <- function(chr, s0, e1) GRanges(chr, IRanges(as.integer(s0) + 1L, as.integer(e1)))   # 0-based 半开 → 1-based 闭
read_se <- function(g) {                                                                       # ROSE 超级增强子表
  f <- file.path(SE_DIR, paste0(g, "_SuperEnhancers.table.txt"))
  hdr <- readLines(f, n = 6)
  cutoff <- as.numeric(sub("^#Using cutoff of ([0-9.]+) .*$", "\\1", hdr[4]))
  d <- fread(f, skip = "REGION_ID", sep = "\t")
  setnames(d, 7, "signal")
  stopifnot(is.finite(cutoff), ncol(d) == 9, nrow(d) == N_SE_REF[[g]], !anyDuplicated(d$REGION_ID), all(d$isSuper == 1),
            all(d$STOP > d$START), all(d$signal > cutoff), all(grepl("^(chr[0-9]+|NW_)", d$CHROM)))
  bed <- fread(file.path(SE_DIR, paste0(g, "_SuperEnhancers.bed")), header = FALSE)
  stopifnot(nrow(bed) == nrow(d), identical(bed$V4, d$REGION_ID), identical(bed$V2, d$START), identical(bed$V3, d$STOP))   # bed 与表同坐标
  gr <- bed2gr(d$CHROM, d$START, d$STOP); mcols(gr)$id <- d$REGION_ID; mcols(gr)$signal <- d$signal; mcols(gr)$rank <- d$enhancerRank
  attr(gr, "cutoff") <- cutoff
  gr }
read_te <- function(g, n_se, cutoff) {                                                         # ROSE 普通增强子表（isSuper = 0）
  f <- file.path(TE_DIR, paste0(g, "_AllEnhancers_ENHANCER_TO_GENE.txt"))
  d <- fread(f, skip = 1, header = FALSE, sep = "\t", fill = TRUE, quote = "",
             col.names = c("REGION_ID", "CHROM", "START", "STOP", "NUM_LOCI", "CONSTITUENT_SIZE", "signal", "enhancerRank", "isSuper",
                           "OVERLAP_GENES", "PROXIMAL_GENES", "CLOSEST_GENE", "enhancerRank2", "isSuper2"))
  stopifnot(ncol(d) == 14, nrow(d) == N_TE_REF[[g]], !anyDuplicated(d$REGION_ID), all(d$isSuper == 0), all(d$isSuper2 == 0),
            all(d$STOP > d$START), min(d$enhancerRank) == n_se + 1L, max(d$enhancerRank) == n_se + nrow(d),
            all(d$signal <= cutoff + 1e-6), all(grepl("^(chr[0-9]+|NW_)", d$CHROM)))
  gr <- bed2gr(d$CHROM, d$START, d$STOP); mcols(gr)$id <- d$REGION_ID
  gr }
read_loops <- function(g) {                                                                    # 公司 combine.bedpe
  f <- file.path(LOOP_DIR, paste0(LOOP_PREFIX[[g]], ".combine.bedpe"))
  d <- fread(f, header = FALSE, sep = "\t", col.names = c("chr1", "s1", "e1", "chr2", "s2", "e2"))
  stopifnot(ncol(d) == 6, nrow(d) == N_LOOP_REF[[g]], !anyDuplicated(d), all(d$chr1 == d$chr2), all(d$chr1 %in% MAIN_CHR),
            all((d$e1 - d$s1) %in% c(5000L, 10000L, 20000L)), all((d$e2 - d$s2) %in% c(5000L, 10000L, 20000L)), all(d$s1 < d$s2))
  d[, `:=`(group = g, loop_id = sprintf("%s_loop%05d", g, .I))]
  d }

# ---------------------------------------------------------------- 启动子（Ensembl 114 转录本 TSS ± 2.5 kb → chr 命名）
txdb <- loadDb(TXDB_IN)
pr <- promoters(txdb, upstream = PROM_FLANK, downstream = PROM_FLANK, columns = c("tx_name", "gene_id"))
stopifnot(length(pr) == 58674, all(width(pr) == 2 * PROM_FLANK))
sn <- as.character(seqnames(pr)); keep <- sn %in% c(as.character(1:18), "X")
pr <- pr[keep]; sn <- sn[keep]
prom <- GRanges(paste0("chr", ifelse(sn == "X", "19", sn)), IRanges(pmax(start(pr), 1L), end(pr)))
mcols(prom)$tx_name <- mcols(pr)$tx_name; mcols(prom)$gene_id <- vapply(mcols(pr)$gene_id, `[`, character(1), 1)
stopifnot(length(prom) == 58674 - 182 - 37, all(seqnames(prom) %in% MAIN_CHR))           # 去掉 Y（182）与 MT（37）转录本
stopifnot(max(end(prom[seqnames(prom) == "chr19"])) <= LEN_X, max(end(prom[seqnames(prom) == "chr19"])) > LEN_CHR18)   # chr19 = X
cat(sprintf("启动子窗口：%d 条转录本（TSS ± %d bp），合并后覆盖 %d 段 / %.1f Mb\n",
            length(prom), PROM_FLANK, length(reduce(prom)), sum(width(reduce(prom))) / 1e6))

# ---------------------------------------------------------------- 逐组：读 SE / TE / loop → 锚点标记 → 分类
ov <- function(q, s) suppressWarnings(overlapsAny(q, s, minoverlap = 1L))
fo <- function(q, s) suppressWarnings(findOverlaps(q, s, minoverlap = 1L))
same_elem <- function(a1, a2, elem) {              # 两端锚点是否落在同一个元件内（≥ 1 个共同元件，如 SE 内部 loop）
  h1 <- fo(a1, elem); h2 <- fo(a2, elem)
  k1 <- paste(queryHits(h1), subjectHits(h1)); k2 <- paste(queryHits(h2), subjectHits(h2))
  out <- logical(length(a1)); out[as.integer(sub(" .*", "", intersect(k1, k2)))] <- TRUE
  out }
lab_anchor <- function(fl, order) {                  # fl = list(SE=, TE=, P=) 逻辑向量；order = 优先级高 → 低；未命中 = Other
  out <- rep("Other", length(fl$SE))
  for (k in rev(order)) out[fl[[k]]] <- k             # 低优先级先写、高优先级后覆盖
  out }
pair_cat <- function(a, b) {                         # 两端标签 → 无序类别名（按 ELEM_RANK 顺序书写）
  ra <- match(a, ELEM_RANK); rb <- match(b, ELEM_RANK)
  paste(ifelse(ra <= rb, a, b), ifelse(ra <= rb, b, a), sep = "-") }

loops_all <- list(); se_tab_all <- list(); elem_n <- list()
for (g in grp_levels) {
  se <- read_se(g); te <- read_te(g, length(se), attr(se, "cutoff")); lp <- read_loops(g)
  stopifnot(length(fo(se, te)) == 0)                                              # 同一次 stitching：SE 与 TE 不重叠
  se_main <- se[as.character(seqnames(se)) %in% MAIN_CHR]
  stopifnot(max(end(se[seqnames(se) == "chr19"])) <= LEN_X, max(end(se[seqnames(se) == "chr19"])) > LEN_CHR18)
  a1 <- bed2gr(lp$chr1, lp$s1, lp$e1); a2 <- bed2gr(lp$chr2, lp$s2, lp$e2)
  stopifnot(all(start(a2) > end(a1)))                                                       # 两锚点不重叠、锚点 1 在上游
  fl1 <- list(SE = ov(a1, se), TE = ov(a1, te), P = ov(a1, prom))
  fl2 <- list(SE = ov(a2, se), TE = ov(a2, te), P = ov(a2, prom))
  # SE 与 TE 本身不重叠，但 10 / 20 kb 的锚点可以同时跨到相邻的 SE 与 TE（stitching 间隔 > 12.5 kb），由优先级规则消歧
  cat(sprintf("  同时命中 SE 与 TE 的锚点：%d / %d；同时命中 SE 与启动子：%d；TE 与启动子：%d\n", sum(fl1$SE & fl1$TE) + sum(fl2$SE & fl2$TE), 2L * nrow(lp),
              sum(fl1$SE & fl1$P) + sum(fl2$SE & fl2$P), sum(fl1$TE & fl1$P) + sum(fl2$TE & fl2$P)))
  lp[, `:=`(a1_in_SE = fl1$SE, a1_in_TE = fl1$TE, a1_in_P = fl1$P, a2_in_SE = fl2$SE, a2_in_TE = fl2$TE, a2_in_P = fl2$P,
            same_SE = same_elem(a1, a2, se), same_TE = same_elem(a1, a2, te))]                # SE 中位宽约 100 kb，两端同在一个 SE 内的 loop 不少
  for (s in names(SCHEMES)) {
    l1 <- lab_anchor(fl1, SCHEMES[[s]]); l2 <- lab_anchor(fl2, SCHEMES[[s]])
    lp[, (paste0("a1_label_", s)) := l1]; lp[, (paste0("a2_label_", s)) := l2]; lp[, (paste0("category_", s)) := pair_cat(l1, l2)]
  }
  stopifnot(all(lp$category_main %in% CAT_LEVELS), all(lp$category_main[lp$same_SE] == "SE-SE"), all(lp$category_main[lp$same_TE] %in% c("TE-TE", "TE-P")))
  cat(sprintf("  SE–SE 类 %d 条，其中两端在同一 SE 内 %d 条（%.1f%%）；TE–TE 类中两端同一 TE 内 %d 条\n", sum(lp$category_main == "SE-SE"), sum(lp$same_SE),
              100 * sum(lp$same_SE) / sum(lp$category_main == "SE-SE"), sum(lp$same_TE & lp$category_main == "TE-TE")))
  # ---- 每个 SE 是否有 ≥ 1 启动子 loop（主口径：另一端标签 = P，即在启动子内且不在任何 SE 内；放宽：另一端只要在启动子内）
  se_p_main  <- lp$category_main == "SE-P"
  hit_main   <- c(queryHits(fo(se_main, a1[se_p_main & lp$a1_label_main == "SE"])),
                  queryHits(fo(se_main, a2[se_p_main & lp$a2_label_main == "SE"])))
  relax      <- (fl1$SE & fl2$P) | (fl2$SE & fl1$P)
  hit_relax  <- c(queryHits(fo(se_main, a1[relax & fl1$SE])), queryHits(fo(se_main, a2[relax & fl2$SE])))
  hit_any    <- c(queryHits(fo(se_main, a1)), queryHits(fo(se_main, a2)))
  se_tab_all[[g]] <- data.table(group = g, SE_id = mcols(se_main)$id, chr = as.character(seqnames(se_main)),
                                start_1based = start(se_main), end = end(se_main), width_bp = width(se_main),
                                signal = mcols(se_main)$signal, rank = mcols(se_main)$rank,
                                n_loops_any = tabulate(hit_any, length(se_main)),
                                n_promoter_loops_main = tabulate(hit_main, length(se_main)),
                                n_promoter_loops_relaxed = tabulate(hit_relax, length(se_main)))
  elem_n[[g]] <- data.table(group = g, n_SE_all = length(se), n_SE_main_chr = length(se_main), n_SE_scaffold = length(se) - length(se_main),
                            median_SE_width_kb = median(width(se)) / 1e3, n_TE = length(te), median_TE_width_kb = median(width(te)) / 1e3,
                            n_loops = nrow(lp), ROSE_cutoff = attr(se, "cutoff"))
  loops_all[[g]] <- lp
  cat(sprintf("%s：SE %d（主染色体 %d）、TE %d、loop %d；锚点在 SE / TE / 启动子内的比例 %.1f%% / %.1f%% / %.1f%%\n", g, length(se), length(se_main),
              length(te), nrow(lp), 100 * mean(c(fl1$SE, fl2$SE)), 100 * mean(c(fl1$TE, fl2$TE)), 100 * mean(c(fl1$P, fl2$P))))
}
loops <- rbindlist(loops_all); se_tab <- rbindlist(se_tab_all); elem_n <- rbindlist(elem_n)
stopifnot(nrow(loops) == sum(N_LOOP_REF), nrow(se_tab) == sum(elem_n$n_SE_main_chr))

# ---------------------------------------------------------------- 汇总：组 × 口径 × 类别
summ <- rbindlist(lapply(names(SCHEMES), function(s) {
  x <- loops[, .N, by = .(group, category = get(paste0("category_", s)))]
  x <- merge(CJ(group = grp_levels, category = CAT_LEVELS), x, by = c("group", "category"), all.x = TRUE)
  x[is.na(N), N := 0L]
  x[, `:=`(scheme = s, priority = paste(c(SCHEMES[[s]], "Other"), collapse = " > "), pct = 100 * N / sum(N)), by = group]
  x }))
summ[, category := factor(category, levels = CAT_LEVELS)]
summ[, group := factor(group, levels = grp_levels)]
setorder(summ, scheme, group, category)
stopifnot(summ[, .(s = sum(N)), by = .(scheme, group)][, all(s == N_LOOP_REF[as.character(group)])])
main <- summ[scheme == "main"]

old <- fread(OLD_IN)
old[, category := fcase(Interaction_Type == "Super Enhancer - Promoter", "SE-P", Interaction_Type == "Super Enhancer - Super Enhancer", "SE-SE",
                        Interaction_Type == "Super Enhancer - Typical Enhancer", "SE-TE", Interaction_Type == "Super Enhancer - Other", "SE-Other",
                        Interaction_Type == "Promoter - Typical Enhancer", "TE-P", Interaction_Type == "Typical Enhancer - Typical Enhancer", "TE-TE",
                        Interaction_Type == "Typical Enhancer - Other", "TE-Other", Interaction_Type == "Promoter - Promoter", "P-P",
                        Interaction_Type == "Promoter - Other", "P-Other", Interaction_Type == "Other - Other", "Other-Other")]
stopifnot(!anyNA(old$category), sum(old$Count) == N_LOOP_REF[["Mhp"]])
cmp <- merge(main[group == "Mhp", .(category, N_new = N, pct_new = round(pct, 2))],
             old[, .(category, N_old = Count, pct_old = round(Percent, 2))], by = "category")
cmp[, category := factor(category, levels = CAT_LEVELS)]; setorder(cmp, category)
cat("\nMhp 主口径 vs 公司旧版 Fig4C_Fixed（旧版启动子注释来源不明，仅作量级对照）：\n"); print(cmp)

# ---------------------------------------------------------------- 关键数字
se_frac <- se_tab[, .(n_SE_main_chr = .N, n_SE_with_promoter_loop = sum(n_promoter_loops_main > 0),
                      n_SE_with_promoter_loop_relaxed = sum(n_promoter_loops_relaxed > 0), n_SE_with_any_loop = sum(n_loops_any > 0)), by = group]
se_frac <- merge(se_frac, elem_n[, .(group, n_SE_all)], by = "group")
se_frac[, `:=`(pct_SE_with_promoter_loop = 100 * n_SE_with_promoter_loop / n_SE_main_chr,
               pct_SE_with_promoter_loop_relaxed = 100 * n_SE_with_promoter_loop_relaxed / n_SE_main_chr,
               pct_SE_with_promoter_loop_allSE_denom = 100 * n_SE_with_promoter_loop / n_SE_all,
               pct_SE_with_any_loop = 100 * n_SE_with_any_loop / n_SE_main_chr)]
se_frac[, group := factor(group, levels = grp_levels)]; setorder(se_frac, group)
sePct <- setNames(se_frac$pct_SE_with_promoter_loop, as.character(se_frac$group))
sePct_relax <- setNames(se_frac$pct_SE_with_promoter_loop_relaxed, as.character(se_frac$group))
cat("\nSE 中有 ≥ 1 启动子 loop 的占比（主口径 / 放宽口径）：\n"); print(se_frac)
# ---- 图上 / 图注采用的口径（SE_PCT_SCHEME 开关）：两种口径的数字都进关键数字表，图上只写一种并把定义写清
sePct_fig <- if (SE_PCT_SCHEME == "main") sePct else sePct_relax
seN_fig   <- if (SE_PCT_SCHEME == "main") se_frac$n_SE_with_promoter_loop else se_frac$n_SE_with_promoter_loop_relaxed
names(seN_fig) <- as.character(se_frac$group)
fig_lab_line <- if (SE_PCT_SCHEME == "main") "SE–P loop" else "promoter loop"        # 图上第二行：主口径写类别名 SE–P，与图例一致
cap_def <- if (SE_PCT_SCHEME == "main") {
  "had at least one SE–P class loop (that is, the other anchor lies in a promoter and outside any SE; with the relaxed definition, any promoter-overlapping partner anchor, the fractions are %.1f%% and %.1f%%)"
} else {
  "had at least one loop whose other anchor overlaps a promoter (irrespective of whether that anchor also lies in an SE or TE; restricted to SE–P class loops the fractions are %.1f%% and %.1f%%)"
}
cap_alt <- if (SE_PCT_SCHEME == "main") sePct_relax else sePct
stopifnot(abs(sePct[["Con"]] - 36.9) < 0.1, abs(sePct[["Mhp"]] - 37.3) < 0.1, abs(sePct_relax[["Con"]] - 54.7) < 0.1, abs(sePct_relax[["Mhp"]] - 56.9) < 0.1)

kv <- function(item, Con, Mhp, note = "") data.table(item = item, Con = as.character(Con), Mhp = as.character(Mhp), note = note)   # 关键数字表统一为字符
key <- rbindlist(list(
  kv(item = "loop 总数（条顶数字）", Con = N_LOOP_REF[["Con"]], Mhp = N_LOOP_REF[["Mhp"]], note = "公司 combine.bedpe 行数 = loop_stat.txt；每组合并文库 n = 1"),
  kv(item = "SE 数（全部 / 主染色体 chr1–chr19）", Con = sprintf("%d / %d", elem_n[group == "Con", n_SE_all], elem_n[group == "Con", n_SE_main_chr]),
             Mhp = sprintf("%d / %d", elem_n[group == "Mhp", n_SE_all], elem_n[group == "Mhp", n_SE_main_chr]), note = "scaffold 上的 SE 无 loop 鉴定，占比分母用主染色体 SE"),
  kv(item = "TE 数", Con = elem_n[group == "Con", n_TE], Mhp = elem_n[group == "Mhp", n_TE], note = "AllEnhancers 表 isSuper = 0 全部行；与 SE 同一次 stitching（不重叠）"),
  kv(item = "SE / TE 中位宽度 (kb)", Con = sprintf("%.1f / %.1f", elem_n[group == "Con", median_SE_width_kb], elem_n[group == "Con", median_TE_width_kb]),
             Mhp = sprintf("%.1f / %.1f", elem_n[group == "Mhp", median_SE_width_kb], elem_n[group == "Mhp", median_TE_width_kb]), note = "长度效应：不做每元件 loop 数的 SE vs TE 比较"),
  kv(item = "启动子窗口", Con = length(prom), Mhp = length(prom), note = sprintf("Ensembl 114 转录本 TSS ± %d bp（chr1–chr19），两组共用", PROM_FLANK)),
  rbindlist(lapply(CAT_LEVELS, function(k) kv(item = paste0("主口径 ", k, "：n（%）"),
    Con = sprintf("%d (%.1f%%)", main[group == "Con" & category == k, N], main[group == "Con" & category == k, pct]),
    Mhp = sprintf("%d (%.1f%%)", main[group == "Mhp" & category == k, N], main[group == "Mhp" & category == k, pct]),
    note = "锚点标签优先级 SE > P > TE > Other；类别 = 两端标签无序对"))),
  rbindlist(lapply(c("altA", "altB"), function(s) rbindlist(lapply(c("SE-P", "SE-SE", "SE-TE", "SE-Other"), function(k) kv(
    item = sprintf("敏感性口径 %s（%s）%s：n（%%）", s, summ[scheme == s, priority][1], k),
    Con = sprintf("%d (%.1f%%)", summ[scheme == s & group == "Con" & category == k, N], summ[scheme == s & group == "Con" & category == k, pct]),
    Mhp = sprintf("%d (%.1f%%)", summ[scheme == s & group == "Mhp" & category == k, N], summ[scheme == s & group == "Mhp" & category == k, pct]),
    note = "只列 SE 相关四类；全表见 01_a_loop分类_汇总.csv"))))),
  kv(item = "SE 相关 loop 合计（主口径 SE–P + SE–SE + SE–TE + SE–Other）：n（%）",
             Con = sprintf("%d (%.1f%%)", main[group == "Con" & grepl("^SE-", category), sum(N)], main[group == "Con" & grepl("^SE-", category), sum(pct)]),
             Mhp = sprintf("%d (%.1f%%)", main[group == "Mhp" & grepl("^SE-", category), sum(N)], main[group == "Mhp" & grepl("^SE-", category), sum(pct)]), note = ""),
  kv(item = "SE–SE 中两端落在同一 SE 内（SE 内部 loop）：n（占 SE–SE %）",
     Con = sprintf("%d (%.1f%%)", loops[group == "Con", sum(same_SE)], 100 * loops[group == "Con", sum(same_SE) / sum(category_main == "SE-SE")]),
     Mhp = sprintf("%d (%.1f%%)", loops[group == "Mhp", sum(same_SE)], 100 * loops[group == "Mhp", sum(same_SE) / sum(category_main == "SE-SE")]),
     note = "图中仍计入 SE–SE；公司旧版 Mhp SE–SE = 1,082 与本表扣除内部 loop 后的数目接近，推测旧版排除了 SE 内部 loop"),
  kv(item = "有 ≥ 1 启动子 loop 的 SE 占比（主口径）", Con = sprintf("%.1f%% (%d/%d)", sePct[["Con"]], se_frac[group == "Con", n_SE_with_promoter_loop], se_frac[group == "Con", n_SE_main_chr]),
             Mhp = sprintf("%.1f%% (%d/%d)", sePct[["Mhp"]], se_frac[group == "Mhp", n_SE_with_promoter_loop], se_frac[group == "Mhp", n_SE_main_chr]),
             note = "主口径 = SE–P 类 loop（另一端在 SE > P 优先级下标为 P，即在启动子内且不落在任何 SE 内）的 SE 端锚点与该 SE 重叠 / 主染色体 SE 数"),
  kv(item = "同上，放宽口径（另一端在启动子内即可，不论是否同时在 SE / TE 内）", Con = sprintf("%.1f%% (%d/%d)", sePct_relax[["Con"]], se_frac[group == "Con", n_SE_with_promoter_loop_relaxed], se_frac[group == "Con", n_SE_main_chr]),
             Mhp = sprintf("%.1f%% (%d/%d)", sePct_relax[["Mhp"]], se_frac[group == "Mhp", n_SE_with_promoter_loop_relaxed], se_frac[group == "Mhp", n_SE_main_chr]), note = ""),
  kv(item = "图上采用的口径（SE_PCT_SCHEME）", Con = sprintf("%.1f%%", sePct_fig[["Con"]]), Mhp = sprintf("%.1f%%", sePct_fig[["Mhp"]]),
             note = sprintf("%s；图上第二行写 \"%s\"", SE_PCT_SCHEME, fig_lab_line)),
  kv(item = "同上，分母改为全部 SE（含 scaffold）", Con = sprintf("%.1f%%", se_frac[group == "Con", pct_SE_with_promoter_loop_allSE_denom]),
             Mhp = sprintf("%.1f%%", se_frac[group == "Mhp", pct_SE_with_promoter_loop_allSE_denom]), note = ""),
  kv(item = "有 ≥ 1 条任意 loop 的 SE 占比", Con = sprintf("%.1f%%", se_frac[group == "Con", pct_SE_with_any_loop]), Mhp = sprintf("%.1f%%", se_frac[group == "Mhp", pct_SE_with_any_loop]), note = "仅供图注 / 正文参考"),
  kv(item = "公司旧版 Fig4C_Fixed（Mhp）SE–P：n（%）", Con = "—", Mhp = sprintf("%d (%.1f%%)", old[category == "SE-P", Count], old[category == "SE-P", Percent]),
             note = "旧版无脚本、启动子注释来源不明；本脚本主口径由其类别合并规律反推为 SE > P > TE"),
  kv(item = "图注建议句", Con = "", Mhp = "",
             note = sprintf(paste0("Genome-wide loops called from the merged Hi-C library of each group (Con, n = %s; Mhp, n = %s) classified by the elements at their two anchors ",
                                   "(SE and TE, each group's own ROSE calls; P, promoter = Ensembl TSS ± 2.5 kb; each anchor takes a single label with priority SE > P > TE, ",
                                   "so an anchor inside an SE is labelled SE even when it also overlaps a promoter). ",
                                   "SE–P loops: Con %.1f%%, Mhp %.1f%%. The SE–SE class includes loops with both anchors inside the same SE (Con %s, Mhp %s; SE median width ~100 kb). ",
                                   "%.1f%% (Con) and %.1f%% (Mhp) of SEs on chr1–X ", cap_def, ". One merged library per group; descriptive only, no statistical test."),
                            format(N_LOOP_REF[["Con"]], big.mark = ","), format(N_LOOP_REF[["Mhp"]], big.mark = ","),
                            main[group == "Con" & category == "SE-P", pct], main[group == "Mhp" & category == "SE-P", pct],
                            format(loops[group == "Con", sum(same_SE)], big.mark = ","), format(loops[group == "Mhp", sum(same_SE)], big.mark = ","),
                            sePct_fig[["Con"]], sePct_fig[["Mhp"]], cap_alt[["Con"]], cap_alt[["Mhp"]]))
))
stopifnot(identical(names(key), c("item", "Con", "Mhp", "note")), !anyNA(key$Con))

# ---------------------------------------------------------------- 写表
fwrite(loops[, c("group", "loop_id", "chr1", "s1", "e1", "chr2", "s2", "e2",
                 "a1_in_SE", "a1_in_TE", "a1_in_P", "a2_in_SE", "a2_in_TE", "a2_in_P",
                 "same_SE", "same_TE", "a1_label_main", "a2_label_main", "category_main", "category_altA", "category_altB"), with = FALSE],
       file.path(TAB, "01_a_loop分类_逐loop.csv"))
fwrite(summ[, .(scheme, priority, group, category, N, pct = round(pct, 3))], file.path(TAB, "01_a_loop分类_汇总.csv"))
fwrite(se_tab, file.path(TAB, "01_a_SE启动子loop_逐SE.csv"))
fwrite(key, file.path(TAB, "01_a_关键数字.csv"))
for (f in c("01_a_loop分类_逐loop.csv", "01_a_loop分类_汇总.csv", "01_a_SE启动子loop_逐SE.csv", "01_a_关键数字.csv")) stopifnot(file.size(file.path(TAB, f)) > 200)

# ---------------------------------------------------------------- 作图（主口径）
cat_lab <- setNames(gsub("-", "–", CAT_LEVELS), CAT_LEVELS)                 # 图例用 en-dash
cat_col <- c("SE-P" = "#99000D", "SE-SE" = "#CB181D", "SE-TE" = "#FB6A4A", "SE-Other" = "#FCBBA1",   # SE 红系（深 → 浅）
             "TE-P" = "#D95F02", "TE-TE" = "#FDAE6B", "TE-Other" = "#FEE6CE",                          # TE 橙系
             "P-P" = "#2166AC", "P-Other" = "#92C5DE",                                                 # P 蓝系
             "Other-Other" = "#C8C8C8")                                                                # 两端皆无元件：灰
stopifnot(identical(names(cat_col), CAT_LEVELS))
lum <- function(hex) { m <- grDevices::col2rgb(hex) / 255; unname(0.2126 * m[1, ] + 0.7152 * m[2, ] + 0.0722 * m[3, ]) }
txt_col <- ifelse(lum(cat_col) < 0.5, "white", "black"); names(txt_col) <- CAT_LEVELS
pd <- copy(main); setorder(pd, group, category)
pd[, ymid := cumsum(pct) - pct / 2, by = group]
pd[, `:=`(lab = ifelse(pct >= LAB_MIN_PCT, sprintf("%.1f", pct), ""), tcol = txt_col[as.character(category)])]
tot <- pd[, .(y = 100, lab = paste0("n = ", format(sum(N), big.mark = ","))), by = group]
se_txt <- sprintf("SEs with \u2265 1\n%s\nCon %.1f%%\nMhp %.1f%%", fig_lab_line, sePct_fig[["Con"]], sePct_fig[["Mhp"]])   # \u4e3b\u53e3\u5f84\u5199 "SE\u2013P loop"\uff08= \u56fe\u4f8b\u7c7b\u522b\u540d\uff09\uff0c\u653e\u5bbd\u53e3\u5f84\u624d\u5199 "promoter loop"
X_SETXT <- 3.42                                                                    # 右下角文字的右缘（数据坐标；面板 x 范围 0.55–2.45，图例在其右侧）

p <- ggplot(pd, aes(x = group, y = pct, fill = category)) +
  geom_col(width = 0.7, colour = "white", linewidth = lw(0.3), position = position_stack(reverse = TRUE)) +   # 第一类 SE–P 在最底
  geom_text(aes(y = ymid, label = lab, colour = tcol), size = pt2size(5), family = FONT, show.legend = FALSE) +
  geom_text(data = tot, aes(x = group, y = y, label = lab), inherit.aes = FALSE, vjust = -0.45, size = pt2size(FS_SMALL), family = FONT) +
  annotation_custom(grid::textGrob(se_txt, x = 0, y = 0, hjust = 1, vjust = 0, gp = grid::gpar(fontsize = FS_SMALL, fontfamily = FONT, lineheight = 0.95)),
                    xmin = X_SETXT, xmax = X_SETXT, ymin = 0, ymax = 0) +          # 画在面板外、图例下方（clip = "off"；annotation_custom 不参与坐标训练）
  scale_colour_identity() +
  scale_fill_manual(values = cat_col, breaks = rev(CAT_LEVELS), labels = cat_lab[rev(CAT_LEVELS)], name = NULL) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25), expand = expansion(mult = c(0, 0.10))) +
  scale_x_discrete(expand = expansion(add = 0.45)) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = "Loops (%)") +
  theme_mhp() +
  theme(legend.position = "right", legend.justification = c(0, 1), legend.key.size = unit(2.4, "mm"), legend.key.spacing.y = unit(0.35, "mm"),
        legend.margin = margin(0, 0, 0, 1), legend.box.margin = margin(0, 0, 0, 0), legend.text = element_text(size = FS_SMALL, margin = margin(l = 1)),
        axis.text.x = element_text(size = FS_BASE), plot.margin = margin(2, 1, 1, 1))
# 自检：ggplot 实际堆叠的段中点 = 手算 ymid（防止标签错段）
lb <- as.data.table(layer_data(p, 1))[, .(x, fill, ymid_built = (ymin + ymax) / 2)]
chk <- merge(pd[, .(x = as.integer(group), fill = unname(cat_col[as.character(category)]), ymid)], lb, by = c("x", "fill"))
stopifnot(nrow(chk) == nrow(pd), max(abs(chk$ymid - chk$ymid_built)) < 1e-6)

sav(p, STEM, PAN, 60, 44, panel_dir = PNL)
if (identical(Sys.getenv("MHP_FINAL"), "1")) {                                     # 无栅格面板：另存 quartz PDF（供 12_ 核对与打包）
  ggsave(file.path(PAN, paste0(STEM, ".pdf")), p, width = 60, height = 44, units = "mm", device = .pdf_dev)
  stopifnot(file.exists(file.path(PAN, paste0(STEM, ".pdf"))))
}
cat(sprintf("\n图上数字：Con n = %s，Mhp n = %s；SE–P %.1f%% / %.1f%%；SE 有 ≥ 1 %s（口径 %s）%.1f%% / %.1f%%（另一口径 %.1f%% / %.1f%%）\n",
            format(N_LOOP_REF[["Con"]], big.mark = ","), format(N_LOOP_REF[["Mhp"]], big.mark = ","),
            main[group == "Con" & category == "SE-P", pct], main[group == "Mhp" & category == "SE-P", pct], fig_lab_line, SE_PCT_SCHEME,
            sePct_fig[["Con"]], sePct_fig[["Mhp"]], cap_alt[["Con"]], cap_alt[["Mhp"]]))
