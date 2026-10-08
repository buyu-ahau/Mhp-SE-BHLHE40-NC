# ============================================================================
# _fig6_common_本地.R — Fig. 6 各脚本共用的数据路径、分组与作图定义（由其他脚本 source，不单独运行）
# 运行位置：本地；路径变量见 config/paths.R
# ============================================================================

source(file.path(Sys.getenv("MHP_CODE", "."), "config", "paths.R"))

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
source(file.path(MHP_CODE, "shared", "_palette_MHP.R"))   # theme_mhp、sav、lw、pt2size、FONT、FS_*

# ---------- ① 路径 ----------
R5      <- file.path(MHP_RES, "10.验证结果/5.敲低过表达细胞实验")
RAW6    <- file.path(R5, "补充数据/原始数据.xlsx")                                   # 六组：PCR / PCR1–3 / WB / CCK8 sheet
FCS_DIR <- file.path(R5, "补充数据/凋亡流式/上机数据/Exp_20251222_1")                # 18 个样本管 .fcs（CytoFLEX）
APO_XLS <- file.path(R5, "20251228 布宇老师细胞实验结果/凋亡流式/凋亡.xlsx")
X8      <- file.path(MHP_RES, "10.验证结果/2.crisper/PCR&WB/实验结果 X8385 布宇")   # sgRNA 筛选
OUT <- file.path(R5, "新Fig6_作图"); PAN <- file.path(OUT, "Plots/面板"); TAB <- file.path(OUT, "Tables"); PNL <- file.path(OUT, "_panels")
for (d in c(PAN, TAB, PNL)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
chk_file <- function(f) { stopifnot(file.exists(f), file.size(f) > 0); cat("✓", basename(f), file.size(f), "字节\n"); invisible(f) }   # 中文文件名落盘核验

# ---------- ② 六组编码 ----------
GRP6 <- c("Ctrl", "CRISPRi", "CRISPRi+OE", "Ctrl+Mhp", "CRISPRi+Mhp", "CRISPRi+OE+Mhp")
RAW2KEY <- c("WT" = "Ctrl", "SE-KO" = "CRISPRi", "SE-KO+OE" = "CRISPRi+OE",
             "WT+MHP" = "Ctrl+Mhp", "SE-KO+MHP" = "CRISPRi+Mhp", "SE-KO+OE+MHP" = "CRISPRi+OE+Mhp")
to_key6 <- function(x) {                                    # 原始组名 → 图中键（未知组名直接报错，避免静默错配）
  k <- unname(RAW2KEY[as.character(x)])
  if (anyNA(k)) stop("未知组名: ", paste(unique(x[is.na(k)]), collapse = ", "))
  factor(k, levels = GRP6) }
SIGN6 <- matrix(c("−", "−", "−", "+", "+", "+",          # Mhp
                  "−", "+", "+", "−", "+", "+",               # E2 CRISPRi
                  "−", "−", "+", "−", "−", "+"),    # BHLHE40 OE
                nrow = 3, byrow = TRUE, dimnames = list(c("Mhp", "E2 CRISPRi", "BHLHE40 OE"), GRP6))

# ---------- ③ 配色 ----------
COL6    <- c("Ctrl" = "#666666", "CRISPRi" = "#999999", "CRISPRi+OE" = "#CCCCCC",
             "Ctrl+Mhp" = "#A30001", "CRISPRi+Mhp" = "#E64B35", "CRISPRi+OE+Mhp" = "#FE8264")
COL6_PT <- c("Ctrl" = "#1A1A1A", "CRISPRi" = "#5A5A5A", "CRISPRi+OE" = "#8C8C8C",
             "Ctrl+Mhp" = "#5E0000", "CRISPRi+Mhp" = "#A8321C", "CRISPRi+OE+Mhp" = "#C9633E")   # 散点描边：比柱边框深一档

# ---------- ④ 统计 ----------
star <- function(p) fcase(p < 0.001, "***", p < 0.01, "**", p < 0.05, "*", default = "ns")
fmtp <- function(p) ifelse(p < 1e-4, "<0.0001", formatC(p, digits = 3, format = "g"))   # Source Data 口径：< 1e-4 一律写 "<0.0001"
# 单因素方差分析 + Tukey：返回全部两两比较（g1 在 lev 中序号 < g2；diff = mean(g2) − mean(g1)）
tukey_all <- function(d, lev) {
  d <- copy(d)[, g := factor(paste0("G", match(as.character(group), lev)), levels = paste0("G", seq_along(lev)))]
  stopifnot(!anyNA(d$g))
  fit <- aov(value ~ g, data = d); tk <- TukeyHSD(fit)$g; av <- summary(fit)[[1]]
  out <- as.data.table(tk, keep.rownames = "pair")
  out[, c("ia", "ib") := lapply(tstrsplit(sub("^G", "", pair), "-G"), as.integer)]   # 行名 "G5-G2"：ia = 5, ib = 2，diff = mean(G5) − mean(G2)
  out[, `:=`(i1 = pmin(ia, ib), i2 = pmax(ia, ib))][, sgn := ifelse(ia > ib, 1, -1)]
  out <- out[, .(g1 = lev[i1], g2 = lev[i2], diff = sgn * diff, lwr = ifelse(sgn > 0, lwr, -upr), upr = ifelse(sgn > 0, upr, -lwr), p_adj = `p adj`)]
  out[, star := star(p_adj)]
  attr(out, "anova") <- data.table(F = av[["F value"]][1], df1 = av[["Df"]][1], df2 = av[["Df"]][2], p_anova = av[["Pr(>F)"]][1])
  out[order(match(g1, lev), match(g2, lev))]
}
# 从全部两两比较中挑出图上要标的几对（cmp 为 list(c(a, b), …)），并附括号层级 h_level
pick_cmp <- function(st, cmp, h_level) {
  rbindlist(lapply(seq_along(cmp), function(i) { a <- cmp[[i]]
    r <- st[(g1 == a[1] & g2 == a[2]) | (g1 == a[2] & g2 == a[1])]; stopifnot(nrow(r) == 1)
    r[, .(g1 = a[1], g2 = a[2], diff, p_adj, star, h_level = h_level[i], shown = TRUE)] }))
}
CMP6      <- list(c("Ctrl", "Ctrl+Mhp"), c("Ctrl+Mhp", "CRISPRi+Mhp"), c("CRISPRi+Mhp", "CRISPRi+OE+Mhp"), c("Ctrl", "CRISPRi"))
CMP6_H    <- c(3, 2, 1, 1)                                  # 括号层级：1–4 最高、4–5 中、5–6 与 1–2 最低（互不重叠）
CMP6_LAB  <- c("Mhp effect in control", "CRISPRi vs control under Mhp", "OE rescue under Mhp", "CRISPRi vs control, uninfected")

# ---------- ⑤ 作图 ----------
# 显著性括号：横线 + 上方星号（ns 用小字）；层级 1 / 2 / 3 → ymax × 1.10 / 1.26 / 1.42
brk <- function(br, ymax, lev) {
  br <- copy(br)[, `:=`(x1 = pmin(match(g1, lev), match(g2, lev)), x2 = pmax(match(g1, lev), match(g2, lev)))][, h := ymax * (0.94 + 0.16 * h_level)][, lab := star(p_adj)]
  list(geom_segment(data = br, aes(x = x1 + 0.08, xend = x2 - 0.08, y = h, yend = h), inherit.aes = FALSE, linewidth = lw(0.6)),
       geom_text(data = br[lab != "ns"], aes(x = (x1 + x2) / 2, y = h, label = lab), inherit.aes = FALSE, vjust = 0.25, size = pt2size(FS_BASE + 1), family = FONT),
       geom_text(data = br[lab == "ns"], aes(x = (x1 + x2) / 2, y = h, label = lab), inherit.aes = FALSE, vjust = -0.35, size = pt2size(FS_SMALL), family = FONT))
}
# 三边开口柱 + 空心圆点 + 均值 + SD（向上）；返回 list(p, ymax)，再交给 pm_table（+/− 表）或 xlab_text（文字横轴）收尾
bar_open <- function(d, br, ylab, lev, col, col_pt, seed = 1) {
  d <- copy(d)[, group := factor(as.character(group), levels = lev)]
  sm <- d[, .(mean = mean(value), sd = sd(value)), by = group]
  ymax <- max(c(sm$mean + sm$sd, d$value))
  box <- sm[, .(x = as.integer(group) + c(-0.32, -0.32, 0.32, 0.32), y = c(0, mean, mean, 0)), by = group]
  set.seed(seed)
  p <- ggplot(d, aes(group, value)) +
    geom_path(data = box, aes(x = x, y = y, colour = group, group = group), linewidth = lw(1.0), lineend = "butt", linejoin = "mitre") +
    geom_errorbar(data = sm, aes(y = mean, ymin = mean, ymax = mean + sd, colour = group), width = 0.26, linewidth = lw(0.8)) +
    scale_colour_manual(values = col, guide = "none") +
    ggnewscale::new_scale_colour() +
    geom_jitter(aes(colour = group), width = 0.12, height = 0, size = 1.4, stroke = 0.55, shape = 21, fill = "white") +
    scale_colour_manual(values = col_pt, guide = "none") +
    brk(br, ymax, lev) + scale_x_discrete(drop = FALSE) + scale_y_continuous(breaks = scales::breaks_extended(5)) +
    labs(x = NULL, y = ylab) + theme_mhp() + theme(axis.text = element_text(colour = "black"))
  list(p = p, ymax = ymax)
}
# 横轴下方 +/− 表（行数 = nrow(sign)），画在面板外（coord clip = "off"）；只在显示层面限 y 范围
pm_table <- function(b, sign, top = 1.55) {
  nr <- nrow(sign); nc <- ncol(sign); ymax <- b$ymax
  gy <- -(seq_len(nr) * 0.10) * ymax
  lab <- data.table(x = rep(seq_len(nc), nr), y = rep(gy, each = nc), t = as.vector(t(sign)))
  b$p + geom_text(data = lab, aes(x, y, label = t), inherit.aes = FALSE, size = pt2size(FS_BASE), family = FONT) +
    annotate("text", x = 0.45, y = gy, label = rownames(sign), hjust = 1, size = pt2size(FS_SMALL), family = FONT) +
    coord_cartesian(ylim = c(0, ymax * top), expand = FALSE, clip = "off") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.margin = margin(4, 4, 8 * nr + 2, 30))
}
# 文字横轴（sgRNA 筛选五组：Mock / NC / sg1 / sg2 / sg3）
xlab_text <- function(b, labels, top = 1.55) {
  b$p + scale_x_discrete(drop = FALSE, labels = labels) +
    coord_cartesian(ylim = c(0, b$ymax * top), expand = FALSE, clip = "off") +
    theme(axis.text.x = element_text(size = FS_SMALL, colour = "black"), plot.margin = margin(4, 4, 4, 4))
}
W6 <- 42; H6 <- 46                                          # 六组柱面板尺寸（mm）；新 Fig 3 四组为 34 × 44，六组加宽
# 存 PNG 草稿（sav 默认）+ rds（供拼版脚本只摆放不重画）
sav6 <- function(p, stem, w = W6, h = H6) { out <- sav(p, stem, PAN, w, h); saveRDS(p, file.path(PNL, paste0(stem, ".rds"))); for (f in out) chk_file(f); invisible(out) }

# ---------- ⑥ CytoFLEX FCS：门与象限（来源：00b 脚本从 Exp_20251222_1.xit 解析）----------
# P1 前门：FSC-A / SSC-A 多边形 16 顶点（原始值）。.xit 中另有一版第 10 顶点为 (2637820, 441220)，两版对 18 管 P1 计数均 = 10000（采集停止条件）
P1_POLY <- matrix(c(136886.6, 1696537.6,  904218.6, 4162340.8, 1415773.2, 5596990.0, 3120955.2, 7166137.5,
                    5337692.0, 8017960.5, 6218702.5, 6986806.5, 5906086.0, 5910819.5, 4797717.5, 3534681.8,
                    3120955.2, 1875868.8, 2664695.8,  940824.4, 1842068.6,  620550.8, 1159995.9,  127390.1,
                     676861.0,   -7108.3,  -62051.3, -231272.2, -118890.7,  217055.7,  -62051.3,  710216.3),
                  ncol = 2, byrow = TRUE, dimnames = list(NULL, c("FSC-A", "SSC-A")))
# 补偿：FCS 头 $SPILLOVER 的 FL1-A / FL2-A 子块（FITC → PE 14.32%，PE → FITC 10%）；补偿后值 = 原始值 %*% solve(SPILL)
SPILL6 <- matrix(c(1, 0.1432, 0.1, 1), 2, byrow = TRUE, dimnames = list(c("FL1-A", "FL2-A"), c("FL1-A", "FL2-A")))
QX <- 438840.8; QY <- 35817.0                               # 四象限分界（补偿后 FITC-A / PE-A 原始值）
TUBE_RE <- "^(WT|SE-KO|SE-KO\\+OE)(\\+MHP)?-([123])\\.fcs$"
list_tubes <- function() {                                  # 18 个样本管：文件名 ↔ 原始组名 ↔ 图中键 ↔ 重复号（AV / PI / BLANK / Tube1 / FXRKO Model 不计）
  f <- list.files(FCS_DIR, "\\.fcs$"); f <- f[grepl(TUBE_RE, f)]
  t <- data.table(file = f, raw_group = sub("-[123]\\.fcs$", "", f), replicate = as.integer(sub(".*-([123])\\.fcs$", "\\1", f)))
  t[, group := to_key6(raw_group)]; stopifnot(nrow(t) == 18); t[order(group, replicate)] }
read_tube <- function(file) {                               # 读一管：P1 内事件的补偿后 FITC-A / PE-A
  ff <- flowCore::read.FCS(file.path(FCS_DIR, file), transformation = FALSE, truncate_max_range = FALSE)
  sp <- flowCore::keyword(ff)[["$SPILLOVER"]]; j <- match(c("FL1-A", "FL2-A"), colnames(sp))
  stopifnot(max(abs(sp[j, j] - SPILL6)) < 1e-6)             # 每管自带补偿矩阵须与固定矩阵一致
  e <- flowCore::exprs(ff)
  inP1 <- sp::point.in.polygon(e[, "FSC-A"], e[, "SSC-A"], P1_POLY[, 1], P1_POLY[, 2]) > 0
  x <- e[inP1, c("FL1-A", "FL2-A"), drop = FALSE] %*% solve(SPILL6)
  list(n_total = nrow(e), n_P1 = sum(inP1), dt = data.table(FITC = x[, 1], PE = x[, 2]),
       date = flowCore::keyword(ff)[["$DATE"]], btim = flowCore::keyword(ff)[["$BTIM"]])
}
quad4 <- function(dt) c(UL = 100 * mean(dt$FITC <  QX & dt$PE >= QY), UR = 100 * mean(dt$FITC >= QX & dt$PE >= QY),
                        LL = 100 * mean(dt$FITC <  QX & dt$PE <  QY), LR = 100 * mean(dt$FITC >= QX & dt$PE <  QY))

stopifnot(length(GRP6) == 6, identical(colnames(SIGN6), GRP6), setequal(names(COL6), GRP6), setequal(names(COL6_PT), GRP6),
          all(unlist(CMP6) %in% GRP6), nrow(P1_POLY) == 16, all(RAW2KEY %in% GRP6))
