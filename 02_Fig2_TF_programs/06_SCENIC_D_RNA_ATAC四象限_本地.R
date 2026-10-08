# ============================================================================
# 06_SCENIC_D_RNA_ATAC四象限_本地.R — Fig. 2e：调控子活性变化与 chromVAR 基序可及性变化的四象限图
# 运行位置：本地；路径变量见 config/paths.R
# ============================================================================

source(file.path(Sys.getenv("MHP_CODE", "."), "config", "paths.R"))

suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(ggrepel) })
source(file.path(MHP_CODE, "shared", "_palette_MHP.R"))   # sav()、.pdf_dev

BASE <- file.path(MHP_RES, "8、scRNA-seq/5.转录因子分析_v2_2026-09")
TAB <- file.path(BASE, "Tables"); PLT <- file.path(BASE, "Plots")
CV  <- file.path(MHP_RES, "9.scatac 补充分析/Tables/14_DA_GeneLists")
ATAC_FDR_CUT <- 0.05                            # 常量不能与列名 ATAC_FDR 同名，否则 mutate 内优先解析为列、比较恒为假
CTS <- c("Macrophages", "Neutrophils")          # scRNA 与 snATAC 同名，可直接对应
LAB_SEED <- as.integer(Sys.getenv("LAB_SEED", "42"))   # ggrepel 随机种子：标签布局随种子变化，择优后写死

fam_of <- function(tf) case_when(
  tf %in% c("NFKB1", "NFKB2", "RELA", "RELB", "REL") ~ "NF-κB",
  tf %in% c("FOS", "FOSB", "FOSL1", "FOSL2", "JUN", "JUNB", "JUND", "ATF3", "BATF", "BATF3", "ATF2", "ATF7", "JDP2") ~ "AP-1",
  grepl("^IRF[0-9]$", tf) ~ "IRF",
  grepl("^STAT[0-9]", tf) ~ "STAT",
  TRUE ~ "Other")
fam_col <- c("AP-1" = "#D55E00", "NF-κB" = "#CC79A7", "IRF" = "#0072B2", "STAT" = "#009E73")   # Okabe-Ito，色盲友好

B <- readRDS(file.path(TAB, "B_objects.rds"))

make_one <- function(k) {
  cv <- read.csv(file.path(CV, paste0(k, "_MhpvsMock_DAMotif.csv")), stringsAsFactors = FALSE) %>%
    mutate(TF = sub("_[0-9]+$", "", motif))
  stopifnot(!anyDuplicated(cv$TF))
  m <- inner_join(B %>% filter(celltype == k) %>% select(TF, log2FC, class, mlog10FDR),
                  cv %>% select(TF, motif, MeanDiff, ATAC_FDR = FDR), by = "TF") %>%
    mutate(family = fam_of(TF),
           rna_sig = class %in% c("Robust up", "Robust down"),
           both = rna_sig & ATAC_FDR < ATAC_FDR_CUT,
           concordant = sign(log2FC) == sign(MeanDiff),
           celltype = k)
  rho <- cor(m$log2FC, m$MeanDiff, method = "spearman")
  famsum <- m %>% filter(family != "Other") %>% group_by(family) %>%
    summarise(n = n(), both_sig = sum(both), both_concordant = sum(both & concordant),
              RNA_median = round(median(log2FC), 2), ATAC_median = round(median(MeanDiff), 2), .groups = "drop") %>%
    mutate(celltype = k)
  allsum <- m %>% summarise(family = "All TFs", n = n(), both_sig = sum(both), both_concordant = sum(both & concordant),
                            RNA_median = round(median(log2FC), 2), ATAC_median = round(median(MeanDiff), 2)) %>% mutate(celltype = k)

  hl <- m %>% filter(family != "Other") %>% mutate(lab = TF)
  bh <- m %>% filter(TF == "BHLHE40") %>% mutate(lab = TF)
  # 标签层包含全部点（未标注点标签为空），彩色点按画出的大小设 point.size → 标签不压彩色点；灰点不避让
  labd <- m %>% select(log2FC, MeanDiff, TF, family) %>%
    left_join(hl %>% select(TF, lab), by = "TF") %>%
    mutate(lab = ifelse(TF == "BHLHE40", TF, ifelse(is.na(lab), "", lab)),
           psize = case_when(TF == "BHLHE40" ~ 3.8, family != "Other" ~ 2.6, TRUE ~ 0),   # 灰点不避让（避让会把标签推远、引线过长）
           lab_col = ifelse(TF == "BHLHE40", "key", family), face = ifelse(TF == "BHLHE40", "bold.italic", "italic"))
  # 自动避让解不开的拥挤处手动推开（换 4 个种子布局相同）：巨噬 ATF2 / ATF3 / IRF6，中性粒 IRF1 / IRF3 / STAT5A；其余引线交叉留待 Illustrator
  NUDGE <- list(Macrophages = list(ATF2 = c(0.10, -0.45), ATF3 = c(0.35, 0.30), IRF6 = c(0.45, -0.25)),
                Neutrophils = list(IRF3 = c(-0.35, -0.35), IRF1 = c(0.40, -0.20), STAT5A = c(0.30, 0.25), BHLHE40 = c(0.20, 0.45)))
  labd$nx <- 0; labd$ny <- 0
  for (tf in names(NUDGE[[k]])) { i <- labd$TF == tf; labd$nx[i] <- NUDGE[[k]][[tf]][1]; labd$ny[i] <- NUDGE[[k]][[tf]][2] }
  p <- ggplot(m, aes(log2FC, MeanDiff)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60") +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey60") +
    geom_point(data = m %>% filter(family == "Other"), colour = "#D0D0D0", size = 1.1, stroke = 0) +
    geom_point(data = hl, colour = "white", size = 3.1, stroke = 0) +     # 白底衬，点重叠时仍分得开
    geom_point(data = hl, aes(colour = family), size = 2.4, stroke = 0, shape = 16) +
    geom_point(data = bh, shape = 21, fill = NA, colour = "#A11122", size = 3.8, stroke = 0.9) +
    geom_text_repel(data = labd, aes(label = lab, colour = lab_col, fontface = I(face), point.size = psize),
                    size = 2.6, min.segment.length = 0, segment.size = 0.25, segment.colour = "grey50",
                    nudge_x = labd$nx, nudge_y = labd$ny,
                    box.padding = 0.35, point.padding = 0.2, force = 3, force_pull = 1, max.overlaps = Inf,
                    max.time = 5, max.iter = 1e5, seed = LAB_SEED, show.legend = FALSE) +
    scale_colour_manual(values = c(fam_col, key = "#A11122"), breaks = names(fam_col), name = NULL) +
    guides(colour = guide_legend(override.aes = list(shape = 16, size = 3))) +
    labs(x = expression("Regulon activity, "*log[2]*" FC (Mhp / Con)"),
         y = expression("Motif accessibility, "*Delta*" chromVAR z (Mhp − Con)")) +
    theme_bw() +
    theme(panel.grid = element_blank(), axis.title = element_text(size = 11), axis.text = element_text(size = 10, colour = "black"),
          legend.position = "top", legend.box = "vertical", legend.justification = "left", legend.spacing.y = unit(0, "mm"),
          legend.text = element_text(size = 9), legend.key.size = unit(3.5, "mm"), legend.margin = margin(0, 0, 0, 0))
  fn <- file.path(PLT, paste0("D_quadrant_", k))
  # 两段式出图（palette 的 sav）：默认只出 PNG 草稿；MHP_FINAL=1 出 SVG + TIFF，并另出 PDF（quartz 支持 κ/Δ/−；本图无栅格）
  sav(p, basename(fn), PLT, 127, 118)
  if (identical(Sys.getenv("MHP_FINAL"), "1")) {
    ggsave(paste0(fn, ".pdf"), p, width = 127, height = 118, units = "mm", device = .pdf_dev)
    stopifnot(sum(grepl("/Subtype /Image", readLines(paste0(fn, ".pdf"), warn = FALSE, skipNul = TRUE), useBytes = TRUE)) == 0)
  }
  cat("✓", basename(fn), sprintf("  ρ = %.3f\n", rho))
  list(m = m, fam = bind_rows(famsum, allsum), rho = rho)
}

out <- lapply(CTS, make_one)
paired <- bind_rows(lapply(out, `[[`, "m"))
fams <- bind_rows(lapply(out, `[[`, "fam"))
write.csv(paired, file.path(TAB, "D_RNA_ATAC_TF_paired.csv"), row.names = FALSE)
write.csv(fams, file.path(TAB, "D_family_summary.csv"), row.names = FALSE)

cat("\n===== 图注用数字（两层都显著 / 其中方向一致）=====\n")
print(as.data.frame(fams %>% select(celltype, family, n, both_sig, both_concordant, RNA_median, ATAC_median)), row.names = FALSE)
cat("\n===== 四个炎症家族中两层都显著的成员 =====\n")
print(as.data.frame(paired %>% filter(family != "Other", both) %>%
      transmute(celltype, family, TF, RNA_log2FC = round(log2FC, 2), ATAC_dZ = round(MeanDiff, 2), concordant) %>%
      arrange(celltype, family, -RNA_log2FC)), row.names = FALSE)
cat("\n===== BHLHE40 =====\n")
print(as.data.frame(paired %>% filter(TF == "BHLHE40") %>% transmute(celltype, RNA_log2FC = round(log2FC, 2), class,
      ATAC_dZ = round(MeanDiff, 3), ATAC_FDR = signif(ATAC_FDR, 2))), row.names = FALSE)
