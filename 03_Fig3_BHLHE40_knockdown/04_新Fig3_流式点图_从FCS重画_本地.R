# ============================================================================
# 04_新Fig3_流式点图_从FCS重画_本地.R — Fig. 3h、3i：由原始 .fcs 文件重画 Annexin V-FITC / PI 与 EdU 代表性流式点图
# 运行位置：本地；路径变量见 config/paths.R
# ============================================================================

source(file.path(Sys.getenv("MHP_CODE", "."), "config", "paths.R"))

suppressPackageStartupMessages({ library(CytoML); library(flowWorkspace); library(flowCore); library(data.table); library(ggplot2); library(patchwork) })
source(file.path(MHP_CODE, "shared", "_palette_MHP.R"))   # theme_mhp、sav、FS_*、FONT、.pdf_dev

D   <- file.path(MHP_RES, "10.验证结果/BHLHE40补充验证/20260618肺炎项目数据汇总")
OUT <- file.path(MHP_RES, "10.验证结果/BHLHE40补充验证/新Fig3_作图")
PAN <- file.path(OUT, "Plots/面板")
GRP <- c("si-NC", "si-NC+Mhp", "si-BHLHE40", "si-BHLHE40+Mhp")
PSEUDO <- c("#2B3B9C", "#2E86C8", "#35B6B0", "#8CCB4A", "#F2D33A", "#F28C28", "#D7301F")   # 仿 FlowJo 伪彩：蓝 → 绿 → 黄 → 红

# ---------------------------------------------------------------- 工具
link_import <- function(sub, wsp, grp) {                     # 以 $FIL 原名建软链接后导入工作区
  fcs <- list.files(file.path(D, sub, "FCS"), "\\.fcs$", recursive = TRUE, full.names = TRUE)
  lk <- file.path(tempdir(), paste0("fcs_", gsub("[^A-Za-z0-9]", "", sub))); unlink(lk, recursive = TRUE); dir.create(lk)
  for (f in fcs) file.symlink(f, file.path(lk, read.FCSheader(f)[[1]][["$FIL"]]))
  suppressMessages(flowjo_to_gatingset(open_flowjo_xml(file.path(D, sub, "FCS", wsp)), name = grp, path = lk))
}
raw_single <- function(gs, fil, chans) {                     # 取某管 Single Cells 事件的原始荧光值（反变换 CytoML 的变换）
  sn <- grep(paste0("^", gsub("([.+])", "\\\\\\1", fil), "_"), sampleNames(gs), value = TRUE); stopifnot(length(sn) == 1)
  ex <- exprs(gh_pop_get_data(gs[[sn]], "Single Cells")); inv <- gh_get_transformations(gs[[sn]], inverse = TRUE)
  out <- sapply(chans, function(ch) if (ch %in% names(inv)) inv[[ch]](ex[, ch]) else ex[, ch]); as.data.table(out)
}
dens_order <- function(x, y) {                               # 局部密度（2D 核密度）→ 颜色，低密度先画
  k <- MASS::kde2d(x, y, n = 150, lims = c(range(x), range(y)))
  ix <- findInterval(x, k$x, all.inside = TRUE); iy <- findInterval(y, k$y, all.inside = TRUE)
  d <- k$z[cbind(ix, iy)]; list(d = d, o = order(d))
}
bx_ticks <- function(bx) list(breaks = bx(c(0, 1e3, 1e4, 1e5)), labels = c(expression(0), expression(10^3), expression(10^4), expression(10^5)))
theme_flow <- function() theme_mhp() + theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = lw(0.6)), axis.line = element_blank(),
  plot.title = element_text(size = FS_BASE, face = "plain", hjust = 0.5, margin = margin(b = 2)), axis.text = element_text(colour = "black", size = FS_SMALL),
  plot.margin = margin(2, 3, 2, 2))

# ---------------------------------------------------------------- 凋亡（Annexin V-FITC / PI）
gsA <- link_import("3.2 细胞流式凋亡", "17-Jun-2026  48 .wsp", "SW48")
bxF <- flowjo_biexp(channelRange = 256, maxValue = 262144, pos = 4.65, neg = 0, widthBasis = -39.8107170553)
bxP <- flowjo_biexp(channelRange = 256, maxValue = 262144, pos = 4.418539922, neg = 0, widthBasis = -31.6227766017)
QX <- bxF(1073.1439194519692); QY <- bxP(751.098551652953)                # 四象限分界（工作区 RectangleGate，原始值 → 显示坐标）
REP_A <- c("si-NC" = "Specimen_002_Tube_017.fcs", "si-NC+Mhp" = "Specimen_002_Tube_015.fcs",      # 各组重复 1（与导出表第 1 行一致）
           "si-BHLHE40" = "Specimen_002_Tube_008.fcs", "si-BHLHE40+Mhp" = "Specimen_002_Tube_004.fcs")
PCT_A <- data.table(group = GRP, Q2 = c("0.90", "17.4", "6.24", "34.0"), Q3 = c("0.98", "2.96", "1.44", "3.55"))  # FlowJo 导出表（数据统计.xls）重复 1，按原样字符串显示
tkF <- bx_ticks(bxF); tkP <- bx_ticks(bxP)
mk_apop <- function(g) {
  r <- raw_single(gsA, REP_A[[g]], c("FITC-A", "PE-Texas Red-A"))
  dt <- data.table(x = bxF(r[["FITC-A"]]), y = bxP(r[["PE-Texas Red-A"]]))
  dd <- dens_order(dt$x, dt$y); dt <- dt[dd$o][, d := dd$d[dd$o]]
  pc <- PCT_A[group == g]
  ggplot(dt, aes(x, y, colour = d)) +
    geom_point(size = 0.25, stroke = 0, shape = 16) +
    geom_vline(xintercept = QX, linewidth = lw(0.6)) + geom_hline(yintercept = QY, linewidth = lw(0.6)) +
    annotate("text", x = 250, y = 250, label = pc$Q2, hjust = 1, vjust = 1, size = pt2size(FS_SMALL), family = FONT) +
    annotate("text", x = 250, y = 6, label = pc$Q3, hjust = 1, vjust = 0, size = pt2size(FS_SMALL), family = FONT) +
    scale_colour_gradientn(colours = PSEUDO, guide = "none") +
    scale_x_continuous(limits = c(0, 256), breaks = tkF$breaks, labels = tkF$labels, expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 256), breaks = tkP$breaks, labels = tkP$labels, expand = c(0, 0)) +
    labs(title = g, x = "Annexin V-FITC", y = "PI") + coord_fixed() + theme_flow()
}
PA <- lapply(GRP, mk_apop)

# ---------------------------------------------------------------- EdU（FSC-A / EdU）
gsE <- link_import("3.4 流式EDU", "EDU.wsp", "EDU")
bxA <- flowjo_biexp(channelRange = 256, maxValue = 262144, pos = 4.418539922, neg = 0, widthBasis = -100)
F1 <- c(15.565851112, 139.7875499389); F2 <- c(171.434148888, 174.2124500611); DS <- 165.4262341595998   # 椭圆门（显示坐标 0–256）
ell <- local({ ctr <- (F1 + F2) / 2; a <- DS / 2; cc <- sqrt(sum((F2 - F1)^2)) / 2; b <- sqrt(a^2 - cc^2); th <- atan2(F2[2] - F1[2], F2[1] - F1[1])
  t <- seq(0, 2 * pi, length.out = 240); data.table(x = ctr[1] + a * cos(t) * cos(th) - b * sin(t) * sin(th), y = ctr[2] + a * cos(t) * sin(th) + b * sin(t) * cos(th)) })
ell[, x := x / 256 * 262144]                                 # FSC-A 为线性：显示坐标 → 原始值（作图横轴用原始 FSC-A）
REP_E <- c("si-NC" = "Specimen_001_Tube_3_018.fcs", "si-NC+Mhp" = "Specimen_001_Tube_3_011.fcs",
           "si-BHLHE40" = "Specimen_001_Tube_2-1_002.fcs", "si-BHLHE40+Mhp" = "Specimen_001_Tube_3_015.fcs")   # 按导出表重复 1
PCT_E <- c("si-NC" = "82.1", "si-NC+Mhp" = "53.8", "si-BHLHE40" = "67.8", "si-BHLHE40+Mhp" = "12.8")          # EDU 数据统计.xls 重复 1
tkA <- bx_ticks(bxA)
mk_edu <- function(g) {
  r <- raw_single(gsE, REP_E[[g]], c("FSC-A", "AmCyan-A"))
  dt <- data.table(x = r[["FSC-A"]], y = bxA(r[["AmCyan-A"]]))
  dd <- dens_order(dt$x / 1024, dt$y); dt <- dt[dd$o][, d := dd$d[dd$o]]
  ggplot(dt, aes(x, y, colour = d)) +
    geom_point(size = 0.25, stroke = 0, shape = 16) +
    geom_path(data = ell, aes(x, y), inherit.aes = FALSE, linewidth = lw(0.7)) +
    annotate("text", x = mean(range(ell$x)), y = max(ell$y) + 8, label = PCT_E[[g]], vjust = 0, size = pt2size(FS_SMALL), family = FONT) +
    scale_colour_gradientn(colours = PSEUDO, guide = "none") +
    scale_x_continuous(limits = c(0, 262144), breaks = c(0, 1e5, 2e5), labels = c("0", "100K", "200K"), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 256), breaks = tkA$breaks, labels = tkA$labels, expand = c(0, 0)) +
    labs(title = g, x = "FSC-A", y = "EdU") + coord_fixed(ratio = 262144 / 256) + theme_flow()
}
PE <- lapply(GRP, mk_edu)

# ---------------------------------------------------------------- 每管统计（R 复现 vs FlowJo 导出表）→ CSV
tube_map <- function(sub) {                                  # 本地文件 ↔ $FIL ↔ 组（本地文件夹名）
  fcs <- list.files(file.path(D, sub, "FCS"), "\\.fcs$", recursive = TRUE, full.names = TRUE)
  data.table(group = basename(dirname(fcs)), local_file = file.path(basename(dirname(fcs)), basename(fcs)),
             FIL = vapply(fcs, function(f) read.FCSheader(f)[[1]][["$FIL"]], ""))
}
read_fj <- function(f) {                                     # FlowJo 导出表：组名只写在每组第 1 行 → 向下填充；"0.90 %" → 0.90
  x <- as.data.table(readxl::read_excel(f, col_names = FALSE, .name_repair = "minimal")); setnames(x, paste0("V", seq_len(ncol(x))))
  x <- x[!is.na(V2)]; x[, V1 := V1[1], by = cumsum(!is.na(V1))]; x
}
num_pct <- function(s) as.numeric(sub("\\s*%$", "", s))

# 凋亡：四象限（Q1 Annexin V−PI+、Q2 Annexin V+PI+、Q3 Annexin V+PI−、Q4 双阴），以 Single Cells 为分母；导出表无文件名 → 按本地 1/2/3.fcs 对应重复 1/2/3
fjA <- read_fj(file.path(D, "3.2 细胞流式凋亡/数据统计.xls"))[, .(group = V1, replicate = as.integer(V2), Q2_FlowJo = num_pct(V3), Q3_FlowJo = num_pct(V4))]
TA <- tube_map("3.2 细胞流式凋亡")[, replicate := as.integer(sub("\\.fcs$", "", basename(local_file)))]
TA <- TA[, {
  r <- raw_single(gsA, FIL, c("FITC-A", "PE-Texas Red-A")); x <- bxF(r[["FITC-A"]]); y <- bxP(r[["PE-Texas Red-A"]])
  .(n_single_cells = length(x), Q1_R = 100 * mean(x < QX & y >= QY), Q2_R = 100 * mean(x >= QX & y >= QY),
    Q3_R = 100 * mean(x >= QX & y < QY), Q4_R = 100 * mean(x < QX & y < QY))
}, by = .(group, replicate, local_file, FIL)]
TA <- merge(TA, fjA, by = c("group", "replicate"), all = TRUE)
TA[, `:=`(apoptosis_R = Q2_R + Q3_R, apoptosis_FlowJo = Q2_FlowJo + Q3_FlowJo, representative = FIL %in% REP_A)]
TA[, diff_apoptosis := apoptosis_R - apoptosis_FlowJo]
stopifnot(nrow(TA) == 12, !anyNA(TA), max(abs(TA$diff_apoptosis)) < 0.2)

# EdU：椭圆门（显示坐标，两焦点距离和 ≤ DS），以 Single Cells 为分母；导出表带文件名 → 按 $FIL 对应；另核对改正后的工作区组名
fjE <- read_fj(file.path(D, "3.4 流式EDU/EDU 数据统计.xls"))[, .(group_FlowJo_table = V1, replicate = as.integer(V2), FIL = V3, EdU_FlowJo = num_pct(V4))]
wsE <- open_flowjo_xml(file.path(D, "3.4 流式EDU/FCS/EDU.wsp"))
wgE <- merge(as.data.table(fj_ws_get_sample_groups(wsE))[groupName %in% GRP], as.data.table(fj_ws_get_samples(wsE)), by = "sampleID")[, .(FIL = name, group_wsp = groupName)]
TE <- tube_map("3.4 流式EDU")[, {
  r <- raw_single(gsE, FIL, c("FSC-A", "AmCyan-A")); x <- r[["FSC-A"]] / 1024; y <- bxA(r[["AmCyan-A"]])   # FSC-A 线性：原始值 / 1024 = 显示坐标
  .(n_single_cells = length(x), EdU_R = 100 * mean(sqrt((x - F1[1])^2 + (y - F1[2])^2) + sqrt((x - F2[1])^2 + (y - F2[2])^2) <= DS))
}, by = .(group, local_file, FIL)]
TE <- merge(merge(TE, fjE, by = "FIL", all = TRUE), wgE, by = "FIL", all.x = TRUE)
TE[, `:=`(diff_EdU = EdU_R - EdU_FlowJo, representative = FIL %in% REP_E)]
stopifnot(nrow(TE) == 12, !anyNA(TE), TE$group == TE$group_FlowJo_table, TE$group == TE$group_wsp, max(abs(TE$diff_EdU)) < 0.5)

TAB <- file.path(OUT, "Tables"); FIG_TAB <- file.path(MHP_RES, "Figure/Fig3_v2_BHLHE40功能_2026-09/数据表")
ord <- function(dt) dt[order(match(group, GRP), replicate)]
csvA <- file.path(TAB, "Fig3E_流式凋亡_每管统计.csv"); csvE <- file.path(TAB, "Fig3H_流式EdU_每管统计.csv")
fwrite(ord(TA)[, .(group, replicate, local_file, FIL, representative, n_single_cells,
                   Q1_R = round(Q1_R, 2), Q2_R = round(Q2_R, 2), Q3_R = round(Q3_R, 2), Q4_R = round(Q4_R, 2),
                   Q2_FlowJo, Q3_FlowJo, apoptosis_R = round(apoptosis_R, 2), apoptosis_FlowJo, diff_apoptosis = round(diff_apoptosis, 2))], csvA, bom = TRUE)
fwrite(ord(TE)[, .(group, replicate, local_file, FIL, group_wsp, representative, n_single_cells,
                   EdU_R = round(EdU_R, 2), EdU_FlowJo, diff_EdU = round(diff_EdU, 2))], csvE, bom = TRUE)
for (f in c(csvA, csvE)) { stopifnot(file.size(f) > 0, file.copy(f, file.path(FIG_TAB, basename(f)), overwrite = TRUE)); cat("✓", basename(f), file.size(f), "字节\n") }
cat(sprintf("凋亡 Q2+Q3：R 与 FlowJo 最大差 %.2f 个百分点；EdU：最大差 %.2f 个百分点\n", max(abs(TA$diff_apoptosis)), max(abs(TE$diff_EdU))))

# ---------------------------------------------------------------- 输出（4 组一行，只第一张保留纵轴标题）
row4 <- function(L) { for (i in 2:4) L[[i]] <- L[[i]] + labs(y = NULL); wrap_plots(L, nrow = 1) }
W <- 4 * 30; H <- 36
sav(row4(PA), "Fig3E_流式点图_AnnexinV-PI", PAN, W, H)
sav(row4(PE), "Fig3H_流式点图_EdU", PAN, W, H)
if (identical(Sys.getenv("MHP_FINAL"), "1")) for (x in list(list(row4(PA), "Fig3E_流式点图_AnnexinV-PI"), list(row4(PE), "Fig3H_流式点图_EdU"))) {
  out <- file.path(PAN, paste0(x[[2]], ".pdf")); ggsave(out, x[[1]], width = W, height = H, units = "mm", device = .pdf_dev)
  stopifnot(sum(grepl("/Subtype /Image", readLines(out, warn = FALSE, skipNul = TRUE), useBytes = TRUE)) == 0)
}
cat("✓ 流式点图完成\n")
