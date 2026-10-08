# ============================================================================
# 05_SCENIC_B_感染差异regulon_本地.R — Fig. 2b：各细胞类型 Mhp 与 Con 的差异调控子活性（细胞水平双侧 Wilcoxon 检验，Benjamini–Hochberg 校正）
# 运行位置：本地；路径变量见 config/paths.R
# ============================================================================

source(file.path(Sys.getenv("MHP_CODE", "."), "config", "paths.R"))

suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(ggrepel) })

BASE <- file.path(MHP_RES, "8、scRNA-seq/5.转录因子分析_v2_2026-09")
TAB <- file.path(BASE, "Tables"); PLT <- file.path(BASE, "Plots")
MIN_CELLS <- 20; LFC_CUT <- 0.25; FDR_CUT <- 0.05; EPS <- 1e-3      # EPS: log2FC 的伪计数（AUC 量级）
DE_CTS <- c("Macrophages", "Neutrophils", "DC cells", "CD4 cells", "CD8 cells", "NK cells",
            "B cells", "Plasma cells", "Epithelial cells")
KEY_TF <- c("BHLHE40", "JUN", "FOS", "ATF3", "JUNB", "NFKB1", "RELA", "REL", "CEBPB", "SPI1", "STAT3", "IRF4")

o <- readRDS(file.path(TAB, "A_objects.rds")); A <- o$A; meta <- o$meta
grp <- meta$Groups; smp <- meta$orig.ident; ct <- meta$celltype_new
stopifnot(all(grp %in% c("Day0", "Day28")))

# ---------------------------------------------------------------- 对数尺度的 Wilcoxon 与 BH
# 返回自然对数 p（双侧），与 wilcox.test(x, y, exact = FALSE, correct = TRUE) 同法
logp_wilcox <- function(x, y) {
  n1 <- length(x); n0 <- length(y); N <- n1 + n0
  r <- rank(c(x, y)); U <- sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2
  tt <- table(r); tie <- sum(tt^3 - tt)
  sig <- sqrt(n1 * n0 / 12 * ((N + 1) - tie / (N * (N - 1))))
  d <- U - n1 * n0 / 2; z <- (d - 0.5 * sign(d)) / sig
  if (sig == 0) return(0)
  min(0, log(2) + pnorm(-abs(z), log.p = TRUE))
}
# BH 在对数尺度：q_(i) = min_{j>=i} p_(j)·m/j
logBH <- function(lp) {
  m <- length(lp); o <- order(lp, decreasing = TRUE)
  q <- lp[o] + log(m) - log(m - seq_len(m) + 1)
  q <- cummin(q); out <- numeric(m); out[o] <- pmin(q, 0); out
}

# ---------------------------------------------------------------- 逐细胞类型检验
run_ct <- function(k) {
  i0 <- which(ct == k & grp == "Day0"); i1 <- which(ct == k & grp == "Day28")
  m0 <- colMeans(A[i0, , drop = FALSE]); m1 <- colMeans(A[i1, , drop = FALSE])
  lp <- vapply(seq_len(ncol(A)), function(j) logp_wilcox(A[i1, j], A[i0, j]), numeric(1))
  # 样本均值（细胞数不足的样本置 NA）
  sm <- sapply(sort(unique(smp[ct == k])), function(s) { ii <- which(ct == k & smp == s)
    if (length(ii) < MIN_CELLS) rep(NA_real_, ncol(A)) else colMeans(A[ii, , drop = FALSE]) })
  s0 <- sm[, grepl("^Day0_", colnames(sm)), drop = FALSE]; s1 <- sm[, grepl("^Day28_", colnames(sm)), drop = FALSE]
  n0s <- sum(!is.na(s0[1, ])); n1s <- sum(!is.na(s1[1, ]))
  up <- apply(s1, 1, min, na.rm = TRUE) > apply(s0, 1, max, na.rm = TRUE)
  dn <- apply(s1, 1, max, na.rm = TRUE) < apply(s0, 1, min, na.rm = TRUE)
  cons <- if (n0s >= 2 && n1s >= 2) ifelse(up, "up", ifelse(dn, "down", "inconsistent")) else NA_character_
  data.frame(celltype = k, regulon = colnames(A), TF = sub("\\(\\+\\)$", "", colnames(A)),
             n_Day0 = length(i0), n_Day28 = length(i1), samples_Day0 = n0s, samples_Day28 = n1s,
             meanAUC_Day0 = m0, meanAUC_Day28 = m1, dAUC = m1 - m0, log2FC = log2((m1 + EPS) / (m0 + EPS)),
             logp = lp, sample_direction = cons, stringsAsFactors = FALSE)
}
res <- do.call(rbind, lapply(DE_CTS, run_ct)) %>%
  group_by(celltype) %>% mutate(logFDR = logBH(logp)) %>% ungroup() %>%
  mutate(p = exp(logp), FDR = exp(logFDR), mlog10FDR = -logFDR / log(10)) %>%
  mutate(sig_cell = FDR < FDR_CUT & abs(log2FC) > LFC_CUT,
         robust = sig_cell & !is.na(sample_direction) & sample_direction == ifelse(log2FC > 0, "up", "down"),
         class = case_when(robust & log2FC > 0 ~ "Robust up", robust & log2FC < 0 ~ "Robust down",
                           sig_cell ~ "Cell-level only", TRUE ~ "NS")) %>%
  group_by(celltype) %>% mutate(rank_log2FC = rank(-log2FC, ties.method = "min")) %>% ungroup()

# 与 wilcox.test 比对（巨噬，未下溢部分）
chk_i0 <- which(ct == "Macrophages" & grp == "Day0"); chk_i1 <- which(ct == "Macrophages" & grp == "Day28")
chk <- res %>% filter(celltype == "Macrophages") %>% mutate(j = match(regulon, colnames(A)))
p_ref <- vapply(chk$j, function(j) suppressWarnings(wilcox.test(A[chk_i1, j], A[chk_i0, j], exact = FALSE)$p.value), numeric(1))
ok <- !is.na(p_ref) & p_ref > 1e-300; cat("wilcox.test 给出 NaN 的 regulon:", sum(is.na(p_ref)), "\n")
cat(sprintf("与 wilcox.test 比对（巨噬，%d 个未下溢）: log10 p 最大绝对差 = %.2e\n", sum(ok), max(abs(log10(p_ref[ok]) - chk$logp[ok] / log(10)))))

sumt <- res %>% group_by(celltype) %>%
  summarise(cells = paste0(n_Day0[1], "/", n_Day28[1]), samples = paste0(samples_Day0[1], "/", samples_Day28[1]),
            sig_cell = sum(sig_cell), robust_up = sum(class == "Robust up"), robust_down = sum(class == "Robust down"), .groups = "drop") %>%
  mutate(celltype = factor(celltype, levels = DE_CTS)) %>% arrange(celltype)
write.csv(res, file.path(TAB, "B_diff_regulon_allCelltypes.csv"), row.names = FALSE)
write.csv(sumt, file.path(TAB, "B_summary.csv"), row.names = FALSE)
saveRDS(res, file.path(TAB, "B_objects.rds"))

mac <- res %>% filter(celltype == "Macrophages")

# ---------------------------------------------------------------- 核对数字
cat("\n===== 各细胞类型差异 regulon 数 =====\n"); print(as.data.frame(sumt), row.names = FALSE)
cat("\n===== 巨噬 稳健上调 前 15（按 log2FC）=====\n")
print(as.data.frame(mac %>% filter(class == "Robust up") %>% arrange(-log2FC) %>% head(15) %>%
      transmute(TF, log2FC = round(log2FC, 2), mlog10FDR = round(mlog10FDR, 1))), row.names = FALSE)
cat("\n===== 巨噬 稳健下调 前 10（按 log2FC）=====\n")
print(as.data.frame(mac %>% filter(class == "Robust down") %>% arrange(log2FC) %>% head(10) %>%
      transmute(TF, log2FC = round(log2FC, 2), mlog10FDR = round(mlog10FDR, 1))), row.names = FALSE)
cat("\n===== 关键 TF（巨噬 / 中性粒）=====\n")
print(as.data.frame(res %>% filter(TF %in% KEY_TF, celltype %in% c("Macrophages", "Neutrophils")) %>%
      transmute(celltype, TF, log2FC = round(log2FC, 2), mlog10FDR = round(mlog10FDR, 1), sample_direction, class) %>%
      arrange(celltype, -log2FC)), row.names = FALSE)
for (k in c("Macrophages", "Neutrophils")) {
  d <- res %>% filter(celltype == k); b <- d %>% filter(TF == "BHLHE40")
  cat(sprintf("\nBHLHE40 @ %s: log2FC %.2f, -log10FDR %.1f, %s；按 log2FC 在 409 个中排第 %d，在 %d 个稳健上调中排第 %d\n",
      k, b$log2FC, b$mlog10FDR, b$class, b$rank_log2FC, sum(d$class == "Robust up"),
      sum(d$class == "Robust up" & d$log2FC > b$log2FC) + 1))
}
