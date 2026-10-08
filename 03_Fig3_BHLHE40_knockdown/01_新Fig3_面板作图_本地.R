# ============================================================================
# 01_新Fig3_面板作图_本地.R — Fig. 3 定量面板（qPCR、免疫印迹灰度比值、CCK-8、ELISA、Annexin V⁺ 与 EdU⁺ 比例），由 Prism 文件中的原始重复值作图与统计
# 运行位置：本地；路径变量见 config/paths.R
# ============================================================================

source(file.path(Sys.getenv("MHP_CODE", "."), "config", "paths.R"))

suppressPackageStartupMessages({ library(xml2); library(data.table); library(readxl); library(ggplot2); library(patchwork) })
source(file.path(MHP_CODE, "shared", "_palette_MHP.R"))   # theme_mhp、sav、lw、pt2size、FONT、.pdf_dev

D   <- file.path(MHP_RES, "10.验证结果/BHLHE40补充验证/20260618肺炎项目数据汇总")
OUT <- file.path(MHP_RES, "10.验证结果/BHLHE40补充验证/新Fig3_作图")
PAN <- file.path(OUT, "Plots/面板"); dir.create(PAN, recursive = TRUE, showWarnings = FALSE); dir.create(file.path(OUT, "Tables"), showWarnings = FALSE)
GRP <- c("si-NC", "si-NC+Mhp", "si-BHLHE40", "si-BHLHE40+Mhp")
COL    <- c("si-NC" = "#707070", "si-NC+Mhp" = "#C2507B", "si-BHLHE40" = "#E9A06B", "si-BHLHE40+Mhp" = "#7B3F8C")
COL_PT <- c("si-NC" = "#1A1A1A", "si-NC+Mhp" = "#9C2A58", "si-BHLHE40" = "#C9763D", "si-BHLHE40+Mhp" = "#55256A")
CMP <- list(c("si-NC", "si-NC+Mhp"), c("si-NC+Mhp", "si-BHLHE40+Mhp"), c("si-NC", "si-BHLHE40"))

# ---------------------------------------------------------------- 读 pzfx
pz_tables <- function(f) { x <- read_xml(f); xml_ns_strip(x); xml_find_all(x, "//Table") }
one_way <- function(t) {                                   # OneWay / Column 表：每组一列，重复为行
  rbindlist(lapply(xml_find_all(t, "./YColumn"), function(yc) {
    v <- suppressWarnings(as.numeric(xml_text(xml_find_all(yc, "./Subcolumn/d")))); v <- v[!is.na(v)]
    data.table(group = trimws(xml_text(xml_find_first(yc, "./Title"))), rep = seq_along(v), value = v) }))
}
two_way <- function(t) {                                   # TwoWay 表：行 = 时间点，子列 = 重复
  rows <- xml_text(xml_find_all(t, "./RowTitlesColumn/Subcolumn/d"))
  rbindlist(lapply(xml_find_all(t, "./YColumn"), function(yc) {
    subs <- xml_find_all(yc, "./Subcolumn")
    rbindlist(lapply(seq_along(subs), function(r) { v <- suppressWarnings(as.numeric(xml_text(xml_find_all(subs[[r]], "./d"))))
      data.table(group = trimws(xml_text(xml_find_first(yc, "./Title"))), time = rows[seq_along(v)], rep = r, value = v) })) }))
}
tb <- function(sub) pz_tables(file.path(D, sub))

read_wb_ma <- function() {
  nm <- c("roi", "Area", "Mean", "Min", "Max", "IntDen", "RawIntDen")
  x <- strsplit(iconv(readLines(file.path(D, "1.2 WB敲低验证/灰度值数据.xls"), warn = FALSE), "latin1", "UTF-8"), "\t", fixed = TRUE)
  stopifnot(length(x) >= 25, identical(trimws(x[[1]][2:7]), nm[2:7]))
  a <- setnames(rbindlist(lapply(x[2:25], function(r) as.data.table(as.list(as.numeric(r[1:7]))))), nm)
  stopifnot(identical(a$roi, as.numeric(1:24)), all(a$Area == 468), !anyNA(a$Mean))
  a[, `:=`(set = "1.2", target = ifelse(roi <= 12, "BHLHE40 protein", "β-actin"))]
  b <- as.data.table(read_excel(file.path(D, "3.1 WB检测凋亡/灰度值数据.xlsx"), sheet = "Results", col_names = FALSE, .name_repair = "minimal"))
  stopifnot(identical(unname(unlist(b[1, 2:7])), nm[2:7]))
  b <- setnames(b[2:49, 1:7], nm)[, lapply(.SD, as.numeric)]                 # 先命名再转数值（readxl minimal 列名为空，.SD 不能用）
  stopifnot(identical(b$roi, as.numeric(1:48)), all(b$Area == 341), !anyNA(b$Mean))
  b[, `:=`(set = "3.1", target = c("Bax", "Bcl-2", "Cleaved caspase-3", "β-actin")[(roi - 1) %/% 12 + 1])]
  o <- rbind(a, b)[, `:=`(rep = ((roi - 1) %% 12) %/% 4 + 1, lane = (roi - 1) %% 4 + 1, MA = Mean * Area)]
  act <- o[target == "β-actin", .(set, rep, lane, actin_area = Area, actin_mean = Mean, actin_ma = MA, actin_rec = IntDen)]
  w <- merge(o[target != "β-actin"], act, by = c("set", "rep", "lane"))
  w[, `:=`(group = GRP[lane], value = MA / actin_ma, rec = IntDen / actin_rec)]
  stopifnot(nrow(w) == 48, !anyNA(w$value))
  # ROI 对应核对：记录值之比 = pzfx（1.2 Table 0；3.1 Data 1/2/3 = Bax / Bcl-2 / Cleaved caspase-3）
  pz <- rbind(one_way(tb("1.2 WB敲低验证/灰度值数据分析.pzfx")[[1]])[, marker := "BHLHE40 protein"],
              rbindlist(lapply(1:3, function(i) one_way(tb("3.1 WB检测凋亡/灰度值数据分析.pzfx")[[i]])[, marker := c("Bax", "Bcl-2", "Cleaved caspase-3")[i]])))
  ck <- merge(pz, w[, .(marker = target, group, rep, rec)], by = c("marker", "group", "rep"))
  stopifnot(nrow(pz) == 48, nrow(ck) == 48, all(abs(ck$value - ck$rec) < 1e-6))
  # 行序与 pzfx 一致（marker → 组 → 重复）：保证 set.seed(1) 下 geom_jitter 的抖动位置与原面板逐点相同
  w <- w[pz[, .(target = marker, group, rep)], on = .(target, group, rep)]
  w[, .(marker = target, set, group, rep, lane, target_area = Area, target_mean = Mean, target_ma = MA, actin_area, actin_mean, actin_ma, value)]
}
WB <- read_wb_ma()
dat <- list()
dat$A <- one_way(tb("1.1 qPCR敲低验证/GraphPad.pzfx")[[1]])[, `:=`(panel = "A", marker = "BHLHE40 mRNA")]
dat$B <- WB[marker == "BHLHE40 protein", .(group, rep, value)][, `:=`(panel = "B", marker = "BHLHE40 protein")]
q <- tb("2.1 qPCR检测炎症因子/GraphPad.pzfx"); qn <- c("IL-1β", "IL-6", "IL-8", "TNF-α")
qa <- rbindlist(lapply(1:4, function(i) two_way(q[[i]])[, marker := qn[i]]))[!is.na(value)]
dat$C <- qa[time == "24h"][, `:=`(panel = "C", time = NULL)]
e <- tb("2.2 ELISA检测炎症因子/GraphPad.pzfx"); en <- c("TNF-α", "IL-6", "IL-1β", "IL-8")   # ELISA 表序与 qPCR 不同
dat$D <- rbindlist(lapply(1:4, function(i) one_way(e[[i]])[, marker := en[i]]))[, panel := "D"]
dat$E <- one_way(tb("3.2 细胞流式凋亡/数据分析.pzfx")[[1]])[, `:=`(panel = "E", marker = "Apoptosis")]
wn <- c("Bax", "Bcl-2", "Cleaved caspase-3")
dat$F <- WB[marker %in% wn, .(group, rep, value, marker)][, panel := "F"]
dat$H <- one_way(tb("3.4 流式EDU/数据分析.pzfx")[[1]])[, `:=`(panel = "H", marker = "EdU")]
X <- rbindlist(dat, use.names = TRUE)
stopifnot(all(X$group %in% GRP))
X[, group := factor(group, levels = GRP)]
# 表 ↔ 指标核对（si-NC+Mhp 组均值，来自放置评估报告）
chk <- function(p, m, v) stopifnot(abs(X[panel == p & marker == m & group == "si-NC+Mhp", mean(value)] - v) < 0.02 * v)
chk("C", "IL-1β", 2.63); chk("C", "IL-6", 2.68); chk("C", "IL-8", 2.35); chk("C", "TNF-α", 2.15)
chk("D", "IL-1β", 165.55); chk("D", "IL-6", 26.07); chk("D", "IL-8", 45.41); chk("D", "TNF-α", 173.77)
chk("E", "Apoptosis", 17.43); chk("H", "EdU", 56.00)
chk("B", "BHLHE40 protein", 0.909); chk("F", "Bax", 1.252); chk("F", "Bcl-2", 0.825); chk("F", "Cleaved caspase-3", 1.091)

# CCK-8：原始 OD 块（列 2 = blank，列 3–6 = 四组）减 blank 均值
ck <- as.data.frame(suppressMessages(read_excel(file.path(D, "3.3 CCK8检测细胞增殖/CCK8数据统计.xlsx"), sheet = 1, col_names = FALSE)))
G <- rbindlist(lapply(c("0h", "24h", "48h", "72h"), function(tp) {
  r0 <- which(ck[[1]] == tp); stopifnot(length(r0) == 1, identical(unname(unlist(ck[r0, 3:6])), GRP))
  blk <- as.data.frame(ck[(r0 + 1):(r0 + 5), 2:6]); m <- apply(blk, 2, as.numeric); b <- mean(m[, 1])
  data.table(time = as.numeric(sub("h", "", tp)), group = rep(GRP, each = 5), rep = rep(1:5, 4), value = as.vector(m[, 2:5]) - b)
}))
G[, group := factor(group, levels = GRP)]
stopifnot(abs(G[time == 72 & group == "si-NC", mean(value)] - 3.467) < 0.01)

# ---------------------------------------------------------------- 统计：单因素方差分析 + Tukey
tukey <- function(d) {
  d <- copy(d)[, g := factor(paste0("G", as.integer(group)))]
  tk <- TukeyHSD(aov(value ~ g, data = d))$g
  rbindlist(lapply(CMP, function(cp) { a <- paste0("G", match(cp, GRP)); nm <- c(paste(a[2], a[1], sep = "-"), paste(a[1], a[2], sep = "-"))
    k <- intersect(nm, rownames(tk)); data.table(g1 = cp[1], g2 = cp[2], p = tk[k, "p adj"]) }))
}
star <- function(p) fcase(p < 0.001, "***", p < 0.01, "**", p < 0.05, "*", default = "ns")
fmtp <- function(p) ifelse(p < 1e-4, "<0.0001", sprintf("%.4f", p))   # 图上写精确 P 值（参考图样式）
ST <- X[, tukey(.SD), by = .(panel, marker)][, star := star(p)]
ST72 <- G[time == 72, tukey(.SD)][, `:=`(panel = "G", marker = "CCK-8 72 h", star = star(p))]
fwrite(rbind(ST, ST72, use.names = TRUE)[, p := signif(p, 3)], file.path(OUT, "Tables/新Fig3_统计.csv"))
fwrite(rbind(X, copy(G)[, `:=`(panel = "G", marker = paste0("CCK-8 ", time, " h"))][, time := NULL], use.names = TRUE, fill = TRUE),   # copy：避免 := 原地删掉 G 的 time 列
       file.path(OUT, "Tables/新Fig3_作图数据.csv"))

# ---------------------------------------------------------------- 作图
pm_grid <- function(p, ymax) {                             # 横轴下方 Mhp / si-BHLHE40 两行 +/− 表（画在面板外，靠 coord 的 clip = "off"）
  gy <- -c(0.10, 0.20) * ymax
  lab <- data.table(x = rep(1:4, 2), y = rep(gy, each = 4), t = c("−", "+", "−", "+", "−", "−", "+", "+"))
  p + geom_text(data = lab, aes(x, y, label = t), inherit.aes = FALSE, size = pt2size(FS_BASE), family = FONT) +
    annotate("text", x = 0.45, y = gy, label = c("Mhp", "si-BHLHE40"), hjust = 1, size = pt2size(FS_SMALL), family = FONT) +
    coord_cartesian(ylim = c(0, ymax * 1.55), expand = FALSE, clip = "off") +   # 只在显示层面限范围：scale 设 limits 会把负 y 的 +/− 文字删掉
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.margin = margin(4, 4, 18, 26))
}
brk <- function(st, ymax) {
  br <- copy(st)[, `:=`(x1 = match(g1, GRP), x2 = match(g2, GRP))][, h := ymax * c(1.10, 1.26, 1.42)[seq_len(.N)]][, lab := star(p)]
  list(geom_segment(data = br, aes(x = x1 + 0.08, xend = x2 - 0.08, y = h, yend = h), inherit.aes = FALSE, linewidth = lw(0.6)),
       geom_text(data = br[lab != "ns"], aes(x = (x1 + x2) / 2, y = h, label = lab), inherit.aes = FALSE, vjust = 0.25, size = pt2size(FS_BASE + 1), family = FONT),
       geom_text(data = br[lab == "ns"], aes(x = (x1 + x2) / 2, y = h, label = lab), inherit.aes = FALSE, vjust = -0.35, size = pt2size(FS_SMALL), family = FONT))
}
theme_f3 <- function(title) list(labs(title = title, x = NULL), theme_mhp(),
  theme(plot.title = element_text(size = FS_BASE, face = "plain", hjust = 0.5), axis.text = element_text(colour = "black")))

# ① mRNA / WB 定量 / 流式比例：空心柱（三边折线，底边由 x 轴承担）+ 实心散点 + 向上误差棒（参考图 GraphPad 样式）
bar4 <- function(d, st, title, ylab) {
  sm <- d[, .(mean = mean(value), sd = sd(value)), by = group]
  ymax <- max(c(sm$mean + sm$sd, d$value))
  box <- sm[, .(x = as.integer(group) + c(-0.32, -0.32, 0.32, 0.32), y = c(0, mean, mean, 0)), by = group]
  p <- ggplot(d, aes(group, value)) +
    geom_path(data = box, aes(x = x, y = y, colour = group, group = group), linewidth = lw(1.0), lineend = "butt", linejoin = "mitre") +
    geom_errorbar(data = sm, aes(y = mean, ymin = mean, ymax = mean + sd, colour = group), width = 0.26, linewidth = lw(0.8)) +
    scale_colour_manual(values = COL, guide = "none") +
    ggnewscale::new_scale_colour() +
    geom_jitter(aes(colour = group), width = 0.12, height = 0, size = 1.5, stroke = 0.55, shape = 21, fill = "white") +
    scale_colour_manual(values = COL_PT, guide = "none") +
    brk(st, ymax) + scale_x_discrete() + scale_y_continuous(breaks = scales::breaks_extended(5)) +
    theme_f3(title) + labs(y = ylab)
  pm_grid(p, ymax)
}
# ② ELISA 分泌量：小提琴 + 空心圆点 + 中位数短线（参考图 o–s）
vio4 <- function(d, st, title, ylab) {
  ymax <- max(d[, max(value) + 3 * bw.nrd0(value), by = group]$V1)   # 小提琴（trim = FALSE）尾部约伸出 3 个带宽，括号放在其上方
  p <- ggplot(d, aes(group, value)) +
    geom_violin(aes(colour = group), fill = NA, trim = FALSE, scale = "width", width = 0.75, linewidth = lw(0.8)) +
    stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", width = 0.32, linewidth = lw(0.5), colour = "black") +
    geom_jitter(aes(colour = group), shape = 21, fill = "white", stroke = 0.5, size = 1.4, width = 0.08, height = 0) +
    scale_colour_manual(values = COL, guide = "none") +
    brk(st, ymax) + scale_x_discrete() + scale_y_continuous(breaks = scales::breaks_extended(5)) +
    theme_f3(title) + labs(y = ylab)
  pm_grid(p, ymax)
}
set.seed(1)
P <- list()
P$A <- bar4(X[panel == "A"], ST[panel == "A"], NULL, expression(paste("Relative ", italic("BHLHE40"), " mRNA")))
P$B <- bar4(X[panel == "B"], ST[panel == "B"], NULL, "BHLHE40 / β-actin")
for (m in qn) P[[paste0("C_", m)]] <- bar4(X[panel == "C" & marker == m], ST[panel == "C" & marker == m], NULL, paste("Relative", m, "mRNA"))
for (m in qn) P[[paste0("D_", m)]] <- vio4(X[panel == "D" & marker == m], ST[panel == "D" & marker == m], NULL, paste(m, "(pg/mL)"))   # ELISA → 小提琴
P$E <- bar4(X[panel == "E"], ST[panel == "E"], NULL, expression("Annexin V"^"+"*" cells (%)"))   # plotmath 上标：Unicode 上标 + 在 quartz PDF 可能缺字形
WLAB <- c("Bax" = "Bax", "Bcl-2" = "Bcl-2", "Cleaved caspase-3" = "C-caspase-3")
for (m in wn) P[[paste0("F_", m)]] <- bar4(X[panel == "F" & marker == m], ST[panel == "F" & marker == m], NULL, paste(WLAB[[m]], "/ β-actin"))
P$H <- bar4(X[panel == "H"], ST[panel == "H"], NULL, expression("EdU"^"+"*" cells (%)"))

# G：CCK-8 生长曲线（均值 ± SD，n = 5 孔）；72 h 处标两组比较
gs <- G[, .(mean = mean(value), sd = sd(value)), by = .(time, group)]
b72 <- ST72[1:2][, `:=`(y1 = gs[time == 72][match(g1, group), mean], y2 = gs[time == 72][match(g2, group), mean], x = c(76, 80))][, lab := star(p)]   # 两组比较共用 si-NC+Mhp 端点，错开 x；标签先算好（rds 在别的脚本里画时无 fmtp）
# ③ 时间曲线（参考图 w/x）：各重复为淡色大圆点，均值为白底描边点并连线
P$G <- ggplot(gs, aes(time, mean, colour = group)) +
  geom_point(data = G, aes(time, value, colour = group), size = 2.2, alpha = 0.30, stroke = 0, shape = 16, inherit.aes = FALSE) +
  geom_line(linewidth = lw(0.9)) + geom_point(size = 1.7, stroke = 0.7, shape = 21, fill = "white") +
  geom_segment(data = b72, aes(x = x, xend = x, y = y1, yend = y2), inherit.aes = FALSE, linewidth = lw(0.5)) +
  geom_segment(data = rbind(b72[, .(x, y = y1)], b72[, .(x, y = y2)]), aes(x = x, xend = x - 1.2, y = y, yend = y), inherit.aes = FALSE, linewidth = lw(0.5)) +   # 括号端点短横线
  geom_text(data = b72, aes(x = x + 1.2, y = (y1 + y2) / 2, label = lab), inherit.aes = FALSE, hjust = 0, vjust = 0.75,   # 横排星号，放竖线右侧
            size = pt2size(FS_BASE + 1), family = FONT) +
  scale_colour_manual(values = COL, name = NULL, guide = "none") +                     # 图例单独成面板（原图例压住 si-NC 曲线）
  scale_x_continuous(breaks = c(0, 24, 48, 72), expand = expansion(add = c(3, 16))) +
  labs(x = "Time (h)", y = expression(paste("CCK-8 (", OD[450], ")"))) +
  theme_mhp() + theme(axis.text = element_text(colour = "black"))

# 组别图例（四组通用，拼版时放一次）
LEG <- cowplot::get_legend(ggplot(gs, aes(time, mean, colour = group)) + geom_line(linewidth = lw(0.9)) +
  geom_point(size = 1.7, stroke = 0.7, shape = 21, fill = "white") + scale_colour_manual(values = COL, name = NULL) +
  theme_mhp() + theme(legend.key.size = unit(3, "mm"), legend.text = element_text(size = FS_BASE)))
P$LEG <- cowplot::ggdraw(LEG)

qa[, `:=`(group = factor(group, levels = GRP), h = as.numeric(sub("h", "", time)), marker = factor(marker, levels = qn))]
qs <- qa[, .(mean = mean(value)), by = .(marker, h, group)]
qst <- qa[, tukey(.SD[, .(group, value)]), by = .(marker, h)][g1 == "si-NC+Mhp" & g2 == "si-BHLHE40+Mhp"][, star := star(p)]
qst <- merge(qst, qa[group %in% c("si-NC+Mhp", "si-BHLHE40+Mhp"), .(ytop = max(value)), by = .(marker, h)], by = c("marker", "h"))
qst[, y := ytop * 1.07]
fwrite(qst[, .(marker, time_h = h, cmp = "si-NC+Mhp vs si-BHLHE40+Mhp", p = signif(p, 3), star)], file.path(OUT, "Tables/补充_四因子mRNA时间序列_统计.csv"))
PS <- ggplot(qs, aes(h, mean, colour = group)) +
  geom_point(data = qa, aes(h, value, colour = group), size = 2.0, alpha = 0.30, stroke = 0, shape = 16, inherit.aes = FALSE) +
  geom_line(linewidth = lw(0.9)) + geom_point(size = 1.6, stroke = 0.7, shape = 21, fill = "white") +
  geom_text(data = qst, aes(h, y, label = star), inherit.aes = FALSE, vjust = ifelse(qst$star == "ns", -0.2, 0.3),
            size = pt2size(ifelse(qst$star == "ns", FS_SMALL, FS_BASE + 1)), family = FONT) +
  facet_wrap(~ marker, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = COL, name = NULL, breaks = GRP) +
  scale_x_continuous(breaks = c(6, 12, 24), expand = expansion(add = c(2.5, 2.5))) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
  labs(x = "Time after infection (h)", y = "Relative mRNA") +
  theme_mhp() + theme(strip.background = element_blank(), strip.text = element_text(size = FS_BASE), legend.position = "top",
                      legend.key.size = unit(3, "mm"), axis.text = element_text(colour = "black"))

W1 <- 34; H1 <- 44
NM <- c(A = "Fig3A_BHLHE40_mRNA_敲低效率", B = "Fig3B_BHLHE40蛋白_WB定量",
        setNames(paste0("Fig3C_", c("IL-1β", "IL-6", "IL-8", "TNF-α"), "_mRNA_24h"), paste0("C_", qn)),
        setNames(paste0("Fig3D_", c("IL-1β", "IL-6", "IL-8", "TNF-α"), "_ELISA"), paste0("D_", qn)),
        E = "Fig3E_凋亡率_AnnexinV",
        setNames(paste0("Fig3F_", c("Bax", "Bcl-2", "CleavedCaspase-3"), "_WB定量"), paste0("F_", wn)),
        G = "Fig3G_CCK8生长曲线", H = "Fig3H_EdU阳性率", LEG = "Fig3_组别图例")
stopifnot(setequal(names(NM), names(P)))
PW <- function(k) switch(k, G = 58, LEG = 36, W1); PH <- function(k) if (k == "LEG") 20 else H1
ONLY <- trimws(strsplit(Sys.getenv("MHP_PANELS"), ",", fixed = TRUE)[[1]]); ONLY <- ONLY[nzchar(ONLY)]
stopifnot(all(ONLY %in% names(P)))
WRITE <- if (length(ONLY)) ONLY else names(P)
if (length(ONLY)) cat("MHP_PANELS：只写", length(WRITE), "个面板 →", paste(NM[WRITE], collapse = "、"), "\n")
TMP <- file.path(tempdir(), "未写面板")                                    # 未列出的面板渲染到此（R 退出即清理），只为保持 RNG 取种顺序（见头注 ★）
DST <- function(k) if (k %in% WRITE) PAN else TMP
for (k in names(P)) sav(P[[k]], NM[[k]], DST(k), PW(k), PH(k))
dir.create(file.path(OUT, "_panels"), showWarnings = FALSE)                  # 供 02 拼版脚本直接摆放
saveRDS(list(P = P, W = sapply(names(P), PW), H = sapply(names(P), PH), NM = NM), file.path(OUT, "_panels/新Fig3_面板.rds"))
sav(PS, "补充_四因子mRNA_6-12-24h", if (length(ONLY)) TMP else PAN, 150, 55)   # MHP_PANELS 时补充图也只渲染到临时目录（保持渲染顺序）
if (identical(Sys.getenv("MHP_FINAL"), "1")) {
  pdf1 <- function(p, stem, w, h, dir = PAN) { out <- file.path(dir, paste0(stem, ".pdf")); ggsave(out, p, width = w, height = h, units = "mm", device = .pdf_dev)
    stopifnot(sum(grepl("/Subtype /Image", readLines(out, warn = FALSE, skipNul = TRUE), useBytes = TRUE)) == 0) }
  for (k in names(P)) pdf1(P[[k]], NM[[k]], PW(k), PH(k), DST(k))
  pdf1(PS, "补充_四因子mRNA_6-12-24h", 150, 55, if (length(ONLY)) TMP else PAN)
  cat("✓ PDF", length(WRITE) + !length(ONLY), "个（栅格数均为 0）\n")
}
# MHP_PANELS 时核对：临时目录中未写面板的 PNG / SVG / TIFF 与 Plots/面板/ 现有文件逐字节相同 → 渲染顺序与全量运行一致
#   （PDF 含时间戳不比；补充图无抖动层、不取 RNG，不参与此核对）
if (length(ONLY)) {
  ext <- if (identical(Sys.getenv("MHP_FINAL"), "1")) c("svg", "tiff") else "png"
  stems <- unname(NM[setdiff(names(P), WRITE)])
  cmp <- CJ(stem = stems, ext = ext)[, same := unname(tools::md5sum(file.path(TMP, paste0(stem, ".", ext)))) ==
                                              unname(tools::md5sum(file.path(PAN, paste0(stem, ".", ext))))]
  cat(sprintf("未写面板与现有文件逐字节相同：%d / %d（%s）\n", sum(cmp$same, na.rm = TRUE), nrow(cmp), paste(ext, collapse = "、")))
  if (!all(cmp$same %in% TRUE)) print(cmp[same %in% FALSE | is.na(same)])
  stopifnot(all(cmp$same %in% TRUE))
}
# 拼版预览改由 02_新Fig3_预览拼版_本地.R 完成（含 WB / 流式原图）
print(ST[, .(panel, marker, cmp = paste(g1, "vs", g2), p = signif(p, 3), star)])
print(ST72[, .(cmp = paste(g1, "vs", g2), p = signif(p, 3), star)])
