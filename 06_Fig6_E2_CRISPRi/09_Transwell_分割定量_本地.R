# ============================================================================
# 09_Transwell_分割定量_本地.R — Fig. 6f 的输入：Transwell 结晶紫染色面积的盲化图像分割定量（盲化拷贝与对应表见 新Fig6_作图/Transwell盲化存档_2026-09-30/）
# 运行位置：本地；路径变量见 config/paths.R
# ============================================================================

source(file.path(Sys.getenv("MHP_CODE", "."), "config", "paths.R"))

source(file.path(MHP_CODE, "06_Fig6_E2_CRISPRi", "_fig6_common_本地.R"))
source(file.path(MHP_CODE, "06_Fig6_E2_CRISPRi", "_transwell_seg_本地.R"))
suppressPackageStartupMessages(library(magick))
SC <- file.path(OUT, "Transwell盲化存档_2026-09-30")
BL <- file.path(SC, "blind"); MK <- file.path(SC, "masks"); OV <- file.path(SC, "ov_blind"); QC <- file.path(OUT, "Plots/QC_Transwell")
for (d in c(MK, OV, QC)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
files <- sort(list.files(BL, "\\.tif$", full.names = TRUE)); stopifnot(length(files) == 18)
P <- SEG_PAR
GRID <- CJ(k_s = P$k_s_steps, pore_max = P$pore_max_steps)[, variant := sprintf("k_s=%.2f|pore_max=%d", k_s, pore_max)][, type := "grid"]
GRID <- rbind(GRID, data.table(k_s = P$k_s, pore_max = NA_integer_, variant = "colour_rule_no_pore_removal", type = "no_pore"),
              data.table(k_s = NA_real_, pore_max = NA_integer_, variant = "ImageJ_like_gray_le_162_no_pore_removal", type = "baseline"))
GRID[, is_main := variant == sprintf("k_s=%.2f|pore_max=%d", P$k_s, P$pore_max)]; stopifnot(sum(GRID$is_main) == 1)

# ---------------------------------------------------------------- 逐盲图：主口径 + 网格
t_all <- Sys.time()
RES <- list(); SENS <- list()
for (f in files) { id <- sub("\\.tif$", "", basename(f)); t0 <- Sys.time()
  a <- read_rgb(f); pre <- seg_prepare(a, P)
  rows <- list()
  for (i in seq_len(nrow(GRID))) { g <- GRID[i]
    if (g$type == "baseline") { pct <- pre$pct_imagej; r <- NULL }
    else { r <- seg_run(pre, P, k_s = g$k_s, pore_max = if (g$type == "no_pore") P$pore_max else g$pore_max, remove_pores = g$type == "grid"); pct <- r$pct }
    rows[[i]] <- data.table(blind_id = id, variant = g$variant, type = g$type, k_s = g$k_s, pore_max = g$pore_max, pct_area = pct,
                            pct_candidate_S = if (is.null(r)) NA_real_ else r$pct_S, pct_pore_removed = if (is.null(r)) NA_real_ else r$pct_pore)
    if (g$is_main) { main <- r
      n_nuc <- count_nuclei(pre, r$pore, P)
      saveRDS(list(mask = r$mask, pore = r$pore, excl = pre$excl), file.path(MK, paste0(id, ".rds")), compress = "xz")
      overlay_png(a, r$mask, r$pore, file.path(OV, paste0(id, ".png")), excl = pre$excl)
      RES[[id]] <- data.table(blind_id = id, pct_area = r$pct, pct_shadow_excluded = 100 * mean(pre$excl), valid_px = pre$n_valid,
                              B_tile90_smoothed_min = min(pre$shadow_tile), B_tile90_smoothed_max = max(pre$shadow_tile),
                              pct_candidate_S = r$pct_S, pct_pore_removed = r$pct_pore, n_pore_obj = pre$n_pore_obj,
                              n_dark_obj_gt_pore_max = pre$n_big_dark, area_dark_obj_gt_pore_max_px = pre$area_big_dark, n_nuclear_obj = n_nuc,
                              R_bg = pre$bg[["R"]], G_bg = pre$bg[["G"]], B_bg = pre$bg[["B"]], pct_imagej_like = pre$pct_imagej) } }
  SENS[[id]] <- rbindlist(rows)
  cat(sprintf("%s  主口径 %.2f%%  基线 %.2f%%  核 %d  孔对象 %d  (%.0f s)\n", id, RES[[id]]$pct_area, pre$pct_imagej, RES[[id]]$n_nuclear_obj, pre$n_pore_obj, as.numeric(Sys.time() - t0, units = "secs"))) }
RES <- rbindlist(RES); SENS <- rbindlist(SENS)
cat("18 张 × 11 变体完成，用时", round(as.numeric(Sys.time() - t_all, units = "mins"), 1), "min\n")

# ---------------------------------------------------------------- 解盲
KEY <- fread(file.path(SC, "blind_key.csv"), encoding = "UTF-8"); stopifnot(nrow(KEY) == 18, setequal(KEY$blind_id, RES$blind_id))
KEY[, group := factor(group, levels = GRP6)][, image_file := sub(".*transell/", "", path)]
RES  <- merge(KEY[, .(blind_id, raw_group, group, replicate, image_file)], RES, by = "blind_id")[order(group, replicate)]
SENS <- merge(KEY[, .(blind_id, raw_group, group, replicate, image_file)], SENS, by = "blind_id")[order(group, replicate, match(variant, GRID$variant))]
# 与 lead 探针（probe_pctarea.R，8-bit 灰度 ≤ 162 整幅）核对：WT 23.6/22.3/16.2、SE-KO 10.6/13.4/9.6、SE-KO+OE 18.6/19.9/26.2、WT+MHP 33.3/24.5/28.2、SE-KO+MHP 11.7/13.9/14.0、SE-KO+OE+MHP 36.4/41.4/43.8
LEAD <- data.table(raw_group = rep(names(RAW2KEY), each = 3), replicate = rep(1:3, 6),
                   lead_probe = c(23.6, 22.3, 16.2, 10.6, 13.4, 9.6, 18.6, 19.9, 26.2, 33.3, 24.5, 28.2, 11.7, 13.9, 14.0, 36.4, 41.4, 43.8))
RES <- merge(RES, LEAD, by = c("raw_group", "replicate"))[order(group, replicate)]
stopifnot(max(abs(RES$pct_imagej_like - RES$lead_probe)) < 0.1)
cat("解盲：基线（灰度 ≤ 162）与 lead 探针 18 个值全部一致（|差| < 0.1 个百分点）\n")
RES[, `:=`(um_per_px = UM_PER_PX, field_area_mm2 = round(1920 * 1080 * UM_PER_PX^2 / 1e6, 4), valid_area_mm2 = round(valid_px * UM_PER_PX^2 / 1e6, 4))]
RES[, nuclear_obj_per_mm2 := round(n_nuclear_obj / valid_area_mm2, 1)]

OBJ <- rbindlist(lapply(RES$blind_id, function(id) { mk <- readRDS(file.path(MK, paste0(id, ".rds"))); M <- mk$mask; L <- label(M); sz <- comp_size(L); n <- sum(!mk$excl)
  data.table(blind_id = id, n_mask_obj = length(sz), n_mask_obj_lt200px = sum(sz < 200), pct_mask_obj_lt200px = 100 * sum(sz[sz < 200]) / n,
             pct_mask_obj_200_1000px = 100 * sum(sz[sz >= 200 & sz < 1000]) / n, pct_mask_obj_ge1000px = 100 * sum(sz[sz >= 1000]) / n) }))
RES <- merge(RES, OBJ, by = "blind_id")[order(group, replicate)]
RES[, `:=`(frag_share_of_value = pct_mask_obj_lt200px / pct_area, pct_area_minus_lt200px = pct_area - pct_mask_obj_lt200px)]
cat(sprintf("掩膜中 < 200 px 碎片占视野：全部 18 图 %.2f–%.2f%%（中位 %.2f%%）；各组均值 %s\n", min(RES$pct_mask_obj_lt200px), max(RES$pct_mask_obj_lt200px), median(RES$pct_mask_obj_lt200px),
            RES[, .(m = round(mean(pct_mask_obj_lt200px), 2)), by = group][, paste(group, m, collapse = " / ")]))
cat("碎片占各图数值的比例（%）：", RES[, .(s = paste(round(100 * frag_share_of_value), collapse = "/")), by = group][, paste(group, s, collapse = "；")], "\n")

# ---------------------------------------------------------------- QC 标记（解盲后、目检 + 客观指标；只标注，不改主口径以外的数值）
# 孔规则失效：孔对象数 < 全部 18 图中位数的一半（img_18 为 98，其余 17 张 244–507）；暗影：触发⑥；折痕 / 照明波动：目检（img_18 斜向折痕，片内照明 ±12%，其余 ±3–5%）
RES[, qc_pore_rule_failed := n_pore_obj < 0.5 * median(n_pore_obj)][, qc_shadow_excluded := pct_shadow_excluded > 0]
NOTE_18 <- paste("Pores filled with indigo dye (RGB ~55/0/155, B above 0.60 x B_bg) -> pore rule largely failed (98 pore objects vs 244-507 in the other 17 images),",
                 "many pores counted as stain -> value still over-estimated; bottom-right illumination shadow (3.3% of field) excluded from numerator and denominator",
                 "(post hoc rule; whole-image value without exclusion 22.39%); diagonal membrane crease; within-image illumination +/-12% (others +/-3-5%).",
                 "Not used as representative image.")
RES[, qc_flag := fifelse(qc_pore_rule_failed | qc_shadow_excluded, "flagged", "pass")]
RES[, qc_note := fifelse(image_file == "重复/WT/3.tif", NOTE_18, "")]
stopifnot(RES[qc_flag == "flagged", all(image_file == "重复/WT/3.tif")], RES[image_file == "重复/WT/3.tif", qc_pore_rule_failed & qc_shadow_excluded])
cat("QC 标记：", RES[qc_flag == "flagged", sprintf("%s #%d（%s）孔对象 %d、暗影排除 %.2f%% 视野", group, replicate, image_file, n_pore_obj, pct_shadow_excluded)], "\n")

# ---------------------------------------------------------------- QC 叠加图改名（组_重复）+ 18 图拼版（按组排）
lab_of <- function(r) sprintf("%s  #%d  (%s)  %.1f%%", r$group, r$replicate, r$blind_id, r$pct_area)
for (i in seq_len(nrow(RES))) { r <- RES[i]; to <- file.path(QC, sprintf("QC_叠加_%s_%d.png", gsub("\\+", "_", r$group), r$replicate))
  file.copy(file.path(OV, paste0(r$blind_id, ".png")), to, overwrite = TRUE) }
tiles <- lapply(seq_len(nrow(RES)), function(i) { r <- RES[i]
  image_annotate(image_border(image_scale(image_read(file.path(OV, paste0(r$blind_id, ".png"))), "640x360"), "white", "0x22"), lab_of(r), size = 18, font = "Arial", color = "black", gravity = "south") })
rows <- lapply(seq(1, 18, 3), function(i) image_append(do.call(c, tiles[i:(i + 2)])))
sheet <- image_append(do.call(c, rows), stack = TRUE)
image_write(sheet, file.path(QC, "QC_叠加_拼版18图.png"), format = "png")

# ---------------------------------------------------------------- 表
PARTAB <- data.table(parameter = c("background", "stain_rule", "k_s", "pore_rule", "k_p", "pore_area_px", "pore_dilation", "opening", "min_object_px", "fill_holes_le_px",
                                   "nuclear_rule", "baseline", "field", "um_per_px", "tuning", "shadow_exclusion"),
                     value = c("per-channel histogram mode of pixels > 100 (R_bg, G_bg, B_bg)", "G <= k_s x G_bg", P$k_s, "connected component of B <= k_p x B_bg with area in pore_area_px",
                               P$k_p, sprintf("%d-%d", P$pore_min, P$pore_max), sprintf("%d px (%dx%d square) before subtraction", P$dil_px, 2 * P$dil_px + 1, 2 * P$dil_px + 1), "3x3 square",
                               P$a_min, P$h_max, sprintf("G <= %.2f x G_bg, minus dilated pores, opening 3x3, >= %d px; touching nuclei not split", P$k_n, P$nuc_min),
                               sprintf("ImageJ-like: unweighted mean(R,G,B) rounded to 8-bit <= %d, whole image, no pore removal", P$imagej_thr), "whole 1920 x 1080 image", UM_PER_PX,
                               "one parameter set fixed on 18 blinded images (random IDs, seed 20260928) from pooled histograms + overlays; unblinded after all values were computed",
                               sprintf("Illumination-shadow exclusion: %d-px tiles, B-channel %.0fth percentile / B_bg, 3x3 tile median; tiles < %.2f excluded from numerator and denominator (other 17 images min 0.979; only img_18 triggered); applied to all non-baseline variants",
                                       P$shadow_tile, 100 * P$shadow_q, P$k_shadow)),
                     sensitivity_steps = c("", "", paste(P$k_s_steps, collapse = " / "), "", "", paste(P$pore_max_steps, collapse = " / "), "", "", "", "", "", "", "", "", "", ""))
fwrite(RES[, .(raw_group, group, replicate, image_file, blind_id, pct_area, qc_flag, qc_pore_rule_failed, qc_shadow_excluded, pct_shadow_excluded, qc_note,
               n_nuclear_obj, nuclear_obj_per_mm2, pct_candidate_S, pct_pore_removed, n_pore_obj,
               n_dark_obj_gt_pore_max, area_dark_obj_gt_pore_max_px, n_mask_obj, n_mask_obj_lt200px, pct_mask_obj_lt200px, frag_share_of_value, pct_area_minus_lt200px,
               pct_mask_obj_200_1000px, pct_mask_obj_ge1000px, B_tile90_smoothed_min, B_tile90_smoothed_max,
               R_bg, G_bg, B_bg, pct_imagej_like, lead_probe, um_per_px, field_area_mm2, valid_area_mm2)],
       file.path(TAB, "Fig6_Transwell_逐图定量.csv"), bom = TRUE)
fwrite(SENS[, .(raw_group, group, replicate, image_file, blind_id, variant, type, k_s, pore_max, pct_area, pct_candidate_S, pct_pore_removed)], file.path(TAB, "Fig6_Transwell_敏感性网格_逐图.csv"), bom = TRUE)
fwrite(PARTAB, file.path(TAB, "Fig6_Transwell_分割参数.csv"), bom = TRUE)
fwrite(KEY[order(group, replicate), .(blind_id, raw_group, group, replicate, image_file)], file.path(TAB, "Fig6_Transwell_盲化对应表.csv"), bom = TRUE)
for (f in c(file.path(TAB, c("Fig6_Transwell_逐图定量.csv", "Fig6_Transwell_敏感性网格_逐图.csv", "Fig6_Transwell_分割参数.csv", "Fig6_Transwell_盲化对应表.csv")),
            file.path(QC, "QC_叠加_拼版18图.png"), file.path(QC, sprintf("QC_叠加_%s_%d.png", gsub("\\+", "_", RES$group), RES$replicate)))) chk_file(f)
print(RES[, .(group, replicate, blind_id, pct_area = round(pct_area, 2), excl = round(pct_shadow_excluded, 2), frag = round(pct_mask_obj_lt200px, 2), frag_share = round(frag_share_of_value, 3),
               n_nuclear_obj, n_pore_obj, G_bg, B_bg, imagej = round(pct_imagej_like, 1), qc_flag)])
print(RES[, .(n = .N, mean = round(mean(pct_area), 2), sd = round(sd(pct_area), 2), imagej_mean = round(mean(pct_imagej_like), 2)), by = group])
