# ============================================================================
# 07_infection_DA.R — Fig. 2c、2f 的输入：snATAC-seq 各细胞类型 Mhp 与 Con 的差异分析（基因活性分数、peak、chromVAR 偏离值）
# 运行位置：服务器（conda 环境 Archr）；路径变量见 config/paths.R
# ============================================================================

source(file.path(Sys.getenv("MHP_CODE", "."), "config", "paths.R"))

# ============= D0：配置 =============
TABLE_DIR <- file.path(PROJ_DIR, "outputs", "07_infection_DA", "Tables")
SE_DIR    <- file.path(PROJ_DIR, "outputs", "07_infection_DA", "SE_rds")
PLOT_DIR7 <- file.path(PROJ_DIR, "outputs", "07_infection_DA", "Plots")
for (dd in c(TABLE_DIR, SE_DIR, PLOT_DIR7)) dir.create(dd, recursive = TRUE, showWarnings = FALSE)

if (!"GroupCellType" %in% colnames(proj@cellColData)) {
  proj$GroupCellType <- paste(proj$group, proj$celltype, sep = "_")
  saveArchRProject(proj, OUT_DIR, load = FALSE)
}
gc_counts <- table(proj$GroupCellType); print(gc_counts)

peak_set <- getPeakSet(proj)
peak_anno_df <- data.frame(
  peak = paste0(seqnames(peak_set), ":", start(peak_set), "-", end(peak_set)),
  seqnames = as.character(seqnames(peak_set)), start = start(peak_set), end = end(peak_set),
  peakType = peak_set$peakType, nearestGene = peak_set$nearestGene, distToTSS = peak_set$distToTSS,
  stringsAsFactors = FALSE)
MIN_NCELLS <- 30
get_col <- function(se, name) as.numeric(as.matrix(assay(se, name))[, 1])

# ============= D1：Phase 1 — 三矩阵 × 8 类 =============
for (mi in list(c("PeakMatrix","DAPeak"), c("MotifMatrix","DAMotif"), c("GeneScoreMatrix","DAGeneScore"))) {
  mtx <- mi[1]; lab <- mi[2]
  cat("\n======== ", mtx, " ========\n")
  for (ct in ct_order) {
    g1 <- paste0("Mhp_", ct); g2 <- paste0("Mock_", ct)
    n1 <- as.integer(gc_counts[g1]); n2 <- as.integer(gc_counts[g2])
    n1 <- ifelse(is.na(n1), 0, n1); n2 <- ifelse(is.na(n2), 0, n2)
    csv_fn <- file.path(TABLE_DIR, paste0(ct, "_MhpvsMock_", lab, ".csv"))
    rds_fn <- file.path(SE_DIR,    paste0(ct, "_MhpvsMock_", lab, ".rds"))
    if (n1 < MIN_NCELLS || n2 < MIN_NCELLS) { cat(sprintf("  跳过 %-18s Mhp=%d Mock=%d\n", ct, n1, n2)); next }
    if (file.exists(csv_fn) && file.exists(rds_fn)) { cat(sprintf("  已完成 %-18s\n", ct)); next }
    cat(sprintf("  跑 %-18s Mhp=%d Mock=%d ... ", ct, n1, n2)); t0 <- Sys.time()
    se <- tryCatch(getMarkerFeatures(
      ArchRProj = proj, useMatrix = mtx, groupBy = "GroupCellType",
      useGroups = g1, bgdGroups = g2,
      useSeqnames = if (mtx == "MotifMatrix") "z" else NULL,
      bias = c("TSSEnrichment", "log10(nFrags)"), testMethod = "wilcoxon"),
      error = function(e) { cat("失败:", conditionMessage(e), "\n"); NULL })
    if (is.null(se)) next
    saveRDS(se, rds_fn)
    rd <- as.data.frame(rowData(se))
    if (mtx == "PeakMatrix") {
      stopifnot(nrow(rd) == nrow(peak_anno_df))
      df <- cbind(peak_anno_df, data.frame(
        Mean_Mhp = get_col(se, "Mean"), Mean_Mock = get_col(se, "MeanBGD"),
        Log2FC = get_col(se, "Log2FC"), MeanDiff = get_col(se, "MeanDiff"),
        FDR = get_col(se, "FDR"), Pval = get_col(se, "Pval"), AUC = get_col(se, "AUC")))
    } else if (mtx == "MotifMatrix") {
      df <- data.frame(motif = rd$name,
        Mean_Mhp = get_col(se, "Mean"), Mean_Mock = get_col(se, "MeanBGD"), MeanDiff = get_col(se, "MeanDiff"),
        FDR = get_col(se, "FDR"), Pval = get_col(se, "Pval"), AUC = get_col(se, "AUC"))
    } else {
      df <- data.frame(gene = rd$name, seqnames = as.character(rd$seqnames), start = rd$start, end = rd$end,
        Mean_Mhp = get_col(se, "Mean"), Mean_Mock = get_col(se, "MeanBGD"),
        Log2FC = get_col(se, "Log2FC"), MeanDiff = get_col(se, "MeanDiff"),
        FDR = get_col(se, "FDR"), Pval = get_col(se, "Pval"), AUC = get_col(se, "AUC"))
    }
    write.csv(df, csv_fn, row.names = FALSE)
    cat(sprintf("✓ %d 行  %.1f min\n", nrow(df), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
}
cat("\n✅ Phase 1 完成 |", format(Sys.time()), "\n")


# ============= D1p：并行版 Phase 1（照 GGJ/HC 07 的 mclapply 写法；与 D1 二选一，都会跳过已完成）=============
# getMarkerFeatures 1-vs-1 是单核，addArchRThreads 无效 → 单核 + mclapply 跨 combo 并行，21 combo ≈ 2 min
# 中断串行版 D1 是安全的（getMarkerFeatures 只读 Arrow）。补丁在 ArchR 命名空间，fork 子进程继承。
library(parallel); addArchRThreads(threads = 1); Sys.setenv(HDF5_USE_FILE_LOCKING = "FALSE")
tasks <- list()
for (mi in list(c("PeakMatrix","DAPeak"), c("MotifMatrix","DAMotif"), c("GeneScoreMatrix","DAGeneScore"))) {
  for (ct in ct_order) {
    g1 <- paste0("Mhp_", ct); g2 <- paste0("Mock_", ct)
    n1 <- as.integer(gc_counts[g1]); n2 <- as.integer(gc_counts[g2]); n1 <- ifelse(is.na(n1), 0, n1); n2 <- ifelse(is.na(n2), 0, n2)
    csv_fn <- file.path(TABLE_DIR, paste0(ct, "_MhpvsMock_", mi[2], ".csv")); rds_fn <- file.path(SE_DIR, paste0(ct, "_MhpvsMock_", mi[2], ".rds"))
    if (n1 < MIN_NCELLS || n2 < MIN_NCELLS) next
    if (file.exists(csv_fn) && file.exists(rds_fn)) next
    tasks[[length(tasks) + 1]] <- list(mtx = mi[1], g1 = g1, g2 = g2, ct = ct, csv_fn = csv_fn, rds_fn = rds_fn)
  }
}
cat("待跑 combos:", length(tasks), "\n")
res <- mclapply(tasks, function(tk) {
  se <- tryCatch(getMarkerFeatures(proj, useMatrix = tk$mtx, groupBy = "GroupCellType", useGroups = tk$g1, bgdGroups = tk$g2,
    useSeqnames = if (tk$mtx == "MotifMatrix") "z" else NULL, bias = c("TSSEnrichment", "log10(nFrags)"), testMethod = "wilcoxon"),
    error = function(e) conditionMessage(e))
  if (is.character(se)) return(paste("FAIL", tk$mtx, tk$ct, "|", se))
  saveRDS(se, tk$rds_fn)
  gv <- function(nm) as.numeric(as.matrix(assay(se, nm))[, 1]); rd <- as.data.frame(rowData(se))
  if (tk$mtx == "PeakMatrix") {
    df <- cbind(peak_anno_df, data.frame(Mean_Mhp = gv("Mean"), Mean_Mock = gv("MeanBGD"), Log2FC = gv("Log2FC"), MeanDiff = gv("MeanDiff"), FDR = gv("FDR"), Pval = gv("Pval"), AUC = gv("AUC")))
  } else if (tk$mtx == "MotifMatrix") {
    df <- data.frame(motif = rd$name, Mean_Mhp = gv("Mean"), Mean_Mock = gv("MeanBGD"), MeanDiff = gv("MeanDiff"), FDR = gv("FDR"), Pval = gv("Pval"), AUC = gv("AUC"))
  } else {
    df <- data.frame(gene = rd$name, seqnames = as.character(rd$seqnames), start = rd$start, end = rd$end, Mean_Mhp = gv("Mean"), Mean_Mock = gv("MeanBGD"), Log2FC = gv("Log2FC"), MeanDiff = gv("MeanDiff"), FDR = gv("FDR"), Pval = gv("Pval"), AUC = gv("AUC"))
  }
  write.csv(df, tk$csv_fn, row.names = FALSE); paste("OK", tk$mtx, tk$ct, nrow(df))
}, mc.cores = min(21, length(tasks)))
cat(paste(unlist(res), collapse = "\n"), "\n"); cat("失败:", sum(grepl("^FAIL", unlist(res))), "/", length(tasks), "\n")
addArchRThreads(threads = 62)

# ============= D2：汇总 =============
summ <- do.call(rbind, lapply(ct_order, function(ct) {
  f <- file.path(TABLE_DIR, paste0(ct, "_MhpvsMock_DAPeak.csv")); if (!file.exists(f)) return(NULL)
  d <- read.csv(f)
  data.frame(celltype = ct, n_tested = nrow(d),
    up_loose   = sum(d$FDR < 0.05 & d$Log2FC >  0.5, na.rm = TRUE), down_loose = sum(d$FDR < 0.05 & d$Log2FC < -0.5, na.rm = TRUE),
    up_strict  = sum(d$FDR < 0.01 & d$Log2FC >  1,   na.rm = TRUE), down_strict = sum(d$FDR < 0.01 & d$Log2FC < -1,  na.rm = TRUE))
}))
print(summ); write.csv(summ, file.path(TABLE_DIR, "_DAR_summary.csv"), row.names = FALSE)

dm <- read.csv(file.path(TABLE_DIR, "Macrophages_MhpvsMock_DAMotif.csv"))
dms <- dm[dm$FDR < 0.05, ]; dms <- dms[order(-dms$MeanDiff), ]
cat("巨噬 DA motif 升:\n"); print(head(dms[, c("motif","MeanDiff","FDR")], 10)); cat("降:\n"); print(tail(dms[, c("motif","MeanDiff","FDR")], 10))

bh <- do.call(rbind, lapply(ct_order, function(ct) {
  f <- file.path(TABLE_DIR, paste0(ct, "_MhpvsMock_DAGeneScore.csv")); if (!file.exists(f)) return(NULL)
  d <- read.csv(f); r <- d[d$gene == "BHLHE40", ][1, ]
  data.frame(celltype = ct, Mean_Mock = r$Mean_Mock, Mean_Mhp = r$Mean_Mhp, Log2FC = r$Log2FC, FDR = r$FDR)
}))
print(bh); write.csv(bh, file.path(TABLE_DIR, "_BHLHE40_GeneScore_MhpvsMock.csv"), row.names = FALSE)

dp <- read.csv(file.path(TABLE_DIR, "Macrophages_MhpvsMock_DAPeak.csv"))
e2p <- dp[dp$seqnames == "13" & dp$end >= 61329246 & dp$start <= 61350682, c("peak","peakType","Log2FC","FDR","Mean_Mock","Mean_Mhp")]
print(e2p); write.csv(e2p, file.path(TABLE_DIR, "_E2_peaks_Macrophages_DA.csv"), row.names = FALSE)


# ============= S1：BHLHE40 GeneScore 按样本 n=3v3 检验（细胞水平 Wilcoxon 是伪重复）=============
if (!exists("gs")) { gsm <- getMatrixFromProject(proj, useMatrix = "GeneScoreMatrix"); gs <- assay(gsm); rownames(gs) <- rowData(gsm)$name; rm(gsm); gc(verbose = FALSE) }
bv <- data.frame(v = as.numeric(gs["BHLHE40", ]), celltype = proj$celltype, group = proj$group, Sample = proj$Sample)
bs <- bv %>% group_by(celltype, Sample, group) %>% summarise(v = mean(v), .groups = "drop")
bt <- bs %>% group_by(celltype) %>% summarise(mock = mean(v[group == "Mock"]), mhp = mean(v[group == "Mhp"]), log2FC = log2(mhp / mock),
  t_p = tryCatch(t.test(v[group == "Mhp"], v[group == "Mock"])$p.value, error = function(e) NA_real_), .groups = "drop") %>% arrange(t_p)
print(bt); write.csv(bt, file.path(TABLE_DIR, "_BHLHE40_GeneScore_replicate_ttest.csv"), row.names = FALSE)

# ============= S2：巨噬 DA motif 关注家族排名 =============
dm <- read.csv(file.path(TABLE_DIR, "Macrophages_MhpvsMock_DAMotif.csv")); dm <- dm[order(-dm$MeanDiff), ]
dm$rank_up <- seq_len(nrow(dm)); dm$rank_down <- rev(seq_len(nrow(dm)))
print(dm[grepl("^(IRF|STAT|SPI1|SPIB|SPIC|NFKB|REL|CEBP|BHLHE40|KLF2|KLF4)", dm$motif), c("motif","MeanDiff","FDR","rank_up","rank_down")])

# ============= Phase 2：Fig3 D/E/F（统一 pal_div 两端色 #B2182B / #2166AC）=============
# D：DAR 数量双向柱图
summ <- read.csv(file.path(TABLE_DIR, "_DAR_summary.csv"))
sd <- bind_rows(summ %>% transmute(celltype, dir = "Up in Mhp", n = up_loose), summ %>% transmute(celltype, dir = "Down in Mhp", n = -down_loose)) %>%
  mutate(celltype = factor(celltype, levels = rev(ct_order[ct_order %in% summ$celltype])), dir = factor(dir, levels = c("Up in Mhp", "Down in Mhp")))
p_d <- ggplot(sd, aes(n, celltype, fill = dir)) + geom_col(width = 0.7, color = "black", linewidth = 0.3) +
  geom_text(aes(label = abs(n), hjust = ifelse(n >= 0, -0.15, 1.15)), size = 3.3) + geom_vline(xintercept = 0, linewidth = 0.4) +
  scale_fill_manual(values = c("Up in Mhp" = "#B2182B", "Down in Mhp" = "#2166AC")) + scale_x_continuous(labels = abs, expand = expansion(mult = 0.2)) +
  labs(x = "Number of DARs (FDR < 0.05, |log2FC| > 0.5)", y = NULL, fill = NULL) + theme_classic(base_size = 12) + theme(axis.text = element_text(color = "black"), legend.position = "top")
ggsave(file.path(PLOT_DIR7, "MHP_ATAC_DAR_counts.pdf"), p_d, width = 6, height = 4.2); ggsave(file.path(PLOT_DIR7, "MHP_ATAC_DAR_counts.png"), p_d, width = 6, height = 4.2, dpi = 300, bg = "white")

# E：巨噬 GeneScore 火山图，标 BHLHE40（bioplot 火山规范：padj=0 用最小非零替换、X 对称、上下调计数）
dg <- read.csv(file.path(TABLE_DIR, "Macrophages_MhpvsMock_DAGeneScore.csv"))
dg <- dg[!is.na(dg$FDR) & !grepl("^ENSSSCG|^U[0-9]|^SNOR|^ssc-mir|^5S_|^7SK|^Y_RNA|^Metazoa|RNase", dg$gene), ]
mn <- min(dg$FDR[dg$FDR > 0]); dg$FDR[dg$FDR == 0] <- mn; dg$nlp <- -log10(dg$FDR)
dg$change <- ifelse(dg$FDR < 0.05 & dg$Log2FC > 0.5, "UP", ifelse(dg$FDR < 0.05 & dg$Log2FC < -0.5, "DOWN", "NS"))
lab <- bind_rows(dg %>% filter(change == "UP") %>% arrange(FDR, -Log2FC) %>% head(12), dg %>% filter(change == "DOWN") %>% arrange(FDR, Log2FC) %>% head(8), dg %>% filter(gene == "BHLHE40")) %>% distinct(gene, .keep_all = TRUE)
xr <- max(abs(dg$Log2FC)) * 1.05; has_rep <- requireNamespace("ggrepel", quietly = TRUE)
p_e <- ggplot(dg, aes(Log2FC, nlp)) + geom_point(data = filter(dg, change == "NS"), color = "grey85", size = 1.1, alpha = 0.5) +
  geom_point(data = filter(dg, change != "NS"), aes(color = change), size = 1.5) + scale_color_manual(values = c(UP = "#B2182B", DOWN = "#2166AC"), guide = "none") +
  geom_point(data = filter(dg, gene == "BHLHE40"), shape = 21, size = 4.5, fill = "#B2182B", color = "black", stroke = 0.8) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey50", linewidth = 0.4) + geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50", linewidth = 0.4) +
  annotate("text", x = -xr * 0.85, y = max(dg$nlp) * 0.97, label = paste(sum(dg$change == "DOWN"), "down"), color = "#2166AC", fontface = "bold", size = 4) +
  annotate("text", x = xr * 0.85, y = max(dg$nlp) * 0.97, label = paste(sum(dg$change == "UP"), "up"), color = "#B2182B", fontface = "bold", size = 4) +
  labs(x = "log2(Mhp / Mock)  gene score", y = "-log10(FDR)", title = "Macrophages") + coord_cartesian(xlim = c(-xr, xr)) +
  theme_classic(base_size = 12) + theme(axis.text = element_text(color = "black"), plot.title = element_text(hjust = 0.5, face = "bold"))
p_e <- if (has_rep) p_e + ggrepel::geom_text_repel(data = lab, aes(label = gene), size = 3.2, fontface = ifelse(lab$gene == "BHLHE40", "bold.italic", "italic"), max.overlaps = 30, box.padding = 0.4, segment.size = 0.3) else p_e + geom_text(data = lab, aes(label = gene), size = 3, fontface = "italic", vjust = -0.8)
ggsave(file.path(PLOT_DIR7, "MHP_ATAC_Macrophage_GeneScore_volcano.pdf"), p_e, width = 6, height = 5.5); ggsave(file.path(PLOT_DIR7, "MHP_ATAC_Macrophage_GeneScore_volcano.png"), p_e, width = 6, height = 5.5, dpi = 300, bg = "white")

# F：巨噬差异 motif 双向柱图（top15 升 / top15 降）
dmf <- dm %>% filter(FDR < 0.05)
mf <- bind_rows(dmf %>% arrange(-MeanDiff) %>% head(15) %>% mutate(dir = "Up"), dmf %>% arrange(MeanDiff) %>% head(15) %>% mutate(dir = "Down")) %>%
  mutate(name = sub("_[0-9]+$", "", motif), name = factor(name, levels = rev(unique(name))))
p_f <- ggplot(mf, aes(MeanDiff, name, fill = dir)) + geom_col(width = 0.72, color = "black", linewidth = 0.25) + geom_vline(xintercept = 0, linewidth = 0.4) +
  scale_fill_manual(values = c(Up = "#B2182B", Down = "#2166AC"), guide = "none") +
  labs(x = "ChromVAR deviation difference (Mhp - Mock)", y = NULL, title = "Macrophages: differential motif activity") +
  theme_classic(base_size = 11) + theme(axis.text = element_text(color = "black"), axis.text.y = element_text(face = "italic"), plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave(file.path(PLOT_DIR7, "MHP_ATAC_Macrophage_DAmotif_bar.pdf"), p_f, width = 5.5, height = 7); ggsave(file.path(PLOT_DIR7, "MHP_ATAC_Macrophage_DAmotif_bar.png"), p_f, width = 5.5, height = 7, dpi = 300, bg = "white")

