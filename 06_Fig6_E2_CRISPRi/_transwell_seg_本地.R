# ============================================================================
# _transwell_seg_本地.R — Transwell 结晶紫图像分割函数与参数（由 09 号脚本 source，不单独运行）
# 运行位置：本地；路径变量见 config/paths.R
# ============================================================================

source(file.path(Sys.getenv("MHP_CODE", "."), "config", "paths.R"))

suppressPackageStartupMessages({ library(tiff); library(png); library(data.table); library(spatstat.geom) })

# ---------- 参数（主口径 + 敏感性步长）----------
SEG_PAR <- list(k_s = 0.75, k_s_steps = c(0.70, 0.75, 0.80),           # 染色阈值：G ≤ k_s × G_bg
                k_p = 0.60,                                            # 孔阈值：B ≤ k_p × B_bg
                pore_min = 30, pore_max = 400, pore_max_steps = c(250, 400, 600),   # 孔对象面积范围（px）
                dil_px = 2,                                            # 孔向外膨胀（5 × 5）
                a_min = 80,                                            # 最小保留对象（px；≈ 69 µm²，单个细胞 ≥ 200 px）
                h_max = 600,                                           # 填孔上限（px；> pore_max + 膨胀环）
                k_n = 0.35, nuc_min = 30,                              # 核对象：G ≤ 0.35 × G_bg，≥ 30 px
                imagej_thr = 162,                                      # 基线：8-bit 灰度 ≤ 162
                shadow_tile = 60, shadow_q = 0.90, k_shadow = 0.95)    # 暗影排除：60 px 块 B 通道 90 分位 / B_bg，3 × 3 平滑后 < 0.95
UM_PER_PX <- 0.9301587301587303                                        # OME-XML PhysicalSizeX（DP74，10 ×，Zoom 0.63）；08 脚本用 200 µm 标尺核过

# ---------- 基本运算 ----------
peak <- function(x, lo = 100) { h <- tabulate(round(x[x > lo]) + 1L, 256); which.max(h) - 1L }   # 直方图峰（只看 > lo 的像素）
sq_off <- function(r) { g <- as.matrix(expand.grid(dr = -r:r, dc = -r:r)); g[!(g[, 1] == 0 & g[, 2] == 0), , drop = FALSE] }
OFF3 <- sq_off(1)
dil <- function(M, offs = OFF3) {                          # 膨胀（方形结构元），逻辑矩阵
  n <- nrow(M); m <- ncol(M); out <- M
  for (i in seq_len(nrow(offs))) { dr <- offs[i, 1]; dc <- offs[i, 2]
    r1 <- max(1, 1 + dr); r2 <- min(n, n + dr); c1 <- max(1, 1 + dc); c2 <- min(m, m + dc)
    out[r1:r2, c1:c2] <- out[r1:r2, c1:c2] | M[(r1 - dr):(r2 - dr), (c1 - dc):(c2 - dc)] }
  out }
ero   <- function(M, offs = OFF3) !dil(!M, offs)
open3 <- function(M) dil(ero(M))
label <- function(M) {                                     # 8 连通域标记；0 = 背景
  L <- as.integer(as.matrix(connected(im(M), background = FALSE, method = "C"))); L[is.na(L)] <- 0L; matrix(L, nrow(M), ncol(M)) }
comp_size <- function(L) { s <- tabulate(L[L > 0]); if (!length(s)) 0L else s }
rm_small <- function(M, amin) { L <- label(M); sz <- comp_size(L); M & (L > 0) & (sz[pmax(L, 1)] >= amin) }
fill_holes <- function(M, hmax) {                          # 填充不与图像边界相连、面积 ≤ hmax 的背景连通域
  L <- label(!M); sz <- comp_size(L); n <- nrow(L); m <- ncol(L)
  bord <- unique(c(L[1, ], L[n, ], L[, 1], L[, m])); ok <- sz <= hmax; ok[bord[bord > 0]] <- FALSE
  M | ((L > 0) & ok[pmax(L, 1)]) }

# ---------- 分割 ----------
read_rgb <- function(f) { a <- suppressWarnings(readTIFF(f)); stopifnot(length(dim(a)) == 3, dim(a)[3] >= 3, max(a) <= 1); a[, , 1:3] }
shadow_map <- function(B, bgB, P = SEG_PAR) {             # ⑥ 暗影：块 B 通道 90 分位 / B_bg → 3 × 3 块中位数平滑 → < k_shadow 的块
  t <- P$shadow_tile; nr <- nrow(B) %/% t; nc <- ncol(B) %/% t; stopifnot(nr * t == nrow(B), nc * t == ncol(B))
  q <- matrix(NA_real_, nr, nc)
  for (i in seq_len(nr)) for (j in seq_len(nc)) q[i, j] <- quantile(B[((i - 1) * t + 1):(i * t), ((j - 1) * t + 1):(j * t)], P$shadow_q, names = FALSE) / bgB
  s <- q; for (i in seq_len(nr)) for (j in seq_len(nc)) s[i, j] <- median(q[max(1, i - 1):min(nr, i + 1), max(1, j - 1):min(nc, j + 1)])
  list(tile = s, excl = kronecker(s < P$k_shadow, matrix(TRUE, t, t)) > 0) }
seg_prepare <- function(a, P = SEG_PAR) {                  # 与 k_s / pore_max 无关的部分：背景峰、B 暗连通域、暗影排除区、基线
  R <- a[, , 1] * 255; G <- a[, , 2] * 255; B <- a[, , 3] * 255
  bg <- c(R = peak(R), G = peak(G), B = peak(B))
  L <- label(B <= P$k_p * bg[["B"]]); sz <- comp_size(L)
  sh <- shadow_map(B, bg[["B"]], P)
  list(G = G, bg = bg, L = L, sz = sz, npx = length(G), excl = sh$excl, n_valid = sum(!sh$excl), shadow_tile = sh$tile,
       pct_imagej = 100 * mean(floor((R + G + B) / 3 + 0.5) <= P$imagej_thr),
       n_pore_obj = sum(sz >= P$pore_min & sz <= P$pore_max), n_big_dark = sum(sz > P$pore_max), area_big_dark = sum(sz[sz > P$pore_max]))
}
seg_run <- function(pre, P = SEG_PAR, k_s = P$k_s, pore_max = P$pore_max, remove_pores = TRUE) {
  S <- pre$G <= k_s * pre$bg[["G"]]
  if (remove_pores) { pore <- (pre$L > 0) & (pre$sz[pmax(pre$L, 1)] >= P$pore_min) & (pre$sz[pmax(pre$L, 1)] <= pore_max); M <- S & !dil(pore, sq_off(P$dil_px)) }
  else { pore <- matrix(FALSE, nrow(S), ncol(S)); M <- S }
  M <- open3(M); M <- rm_small(M, P$a_min); M <- fill_holes(M, P$h_max)
  M <- M & !pre$excl; nv <- pre$n_valid                    # 暗影排除区：分子分母同时剔除（无暗影的图 excl 全 FALSE，结果与未加⑥时完全相同）
  list(mask = M, pore = pore & !pre$excl, S = S & !pre$excl, pct = 100 * sum(M) / nv, pct_S = 100 * sum(S & !pre$excl) / nv, pct_pore = 100 * sum(pore & !pre$excl) / nv)
}
count_nuclei <- function(pre, pore, P = SEG_PAR) {          # 深染核对象数（相触的核不拆分，只作次要指标；暗影排除区不计）
  N <- (pre$G <= P$k_n * pre$bg[["G"]]) & !dil(pore, sq_off(P$dil_px)) & !pre$excl; N <- open3(N); N <- rm_small(N, P$nuc_min); max(label(N)) }

# ---------- QC 叠加图：绿 = 掩膜轮廓（2 px），青 = 剔除的孔（半透明填充），橙 = 暗影排除区（半透明），缩至 1/2 ----------
overlay_png <- function(a, M, pore, file, scale = 2, excl = NULL) {
  o <- a; edge <- M & !ero(ero(M))
  o[, , 1][pore] <- o[, , 1][pore] * 0.3; o[, , 2][pore] <- 0.5 + 0.5 * o[, , 2][pore]; o[, , 3][pore] <- 0.5 + 0.5 * o[, , 3][pore]
  if (!is.null(excl) && any(excl)) { o[, , 1][excl] <- 0.5 + 0.5 * o[, , 1][excl]; o[, , 2][excl] <- 0.25 + 0.5 * o[, , 2][excl]; o[, , 3][excl] <- 0.5 * o[, , 3][excl] }
  o[, , 1][edge] <- 0; o[, , 2][edge] <- 0.9; o[, , 3][edge] <- 0
  if (scale > 1) { n <- dim(o)[1] %/% scale; m <- dim(o)[2] %/% scale; o <- o[seq(1, n * scale, scale), seq(1, m * scale, scale), ] }
  writePNG(o, file); invisible(file) }
