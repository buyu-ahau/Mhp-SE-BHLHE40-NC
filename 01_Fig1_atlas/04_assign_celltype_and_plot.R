# ============================================================================
# 04_assign_celltype_and_plot.R — snATAC-seq：细胞类型注释与 UMAP / t-SNE 嵌入（Fig. 1d 的数据来源）
# 运行位置：服务器（conda 环境 Archr）；路径变量见 config/paths.R
# ============================================================================

source(file.path(Sys.getenv("MHP_CODE", "."), "config", "paths.R"))

suppressPackageStartupMessages({
  library(ArchR); library(BSgenome.Sscrofa.Ensembl.v11.1)
  library(ggplot2); library(patchwork); library(dplyr)
})
addArchRChrPrefix(chrPrefix = FALSE); addArchRThreads(threads = 62); set.seed(42)

PROJ_DIR  <- MHP_ARCHR
OUT_DIR   <- file.path(PROJ_DIR, "ArchR_filtered")
PLOT_DIR2 <- file.path(PROJ_DIR, "outputs", "04_FinalAnnotation")
dir.create(PLOT_DIR2, recursive = TRUE, showWarnings = FALSE)

proj <- loadArchRProject(OUT_DIR, showLogo = FALSE)
cat("cells:", nCells(proj), "\n"); print(table(proj$celltype, proj$group))


# ============= 段1：合并（只跑一次；已跑过则跳过）=============
if (!"celltype_fine" %in% colnames(getCellColData(proj))) {
  proj$celltype_fine <- proj$celltype
  proj$celltype <- ifelse(proj$celltype == "Macrophages_Inf", "Macrophages", proj$celltype)
  saveArchRProject(proj, OUT_DIR, load = FALSE)
  cat("✓ 已合并并保存\n")
} else cat("celltype_fine 已存在，跳过合并\n")
print(table(proj$celltype, proj$group))

# 9 类配色：髓系绿 / 淋巴紫 / 上皮蓝 / 内皮 / 间质棕
ct_col <- c("Macrophages"="#f69886","Neutrophils"="#d984b7","T_NK_cells"="#ef3f33","Plasma_cells"="#895da7",
            "B_cells"="#26a870","Epithelial_cells"="#d4c1ee","AT2_cells"="#d4c1ee","Endothelial_cells"="#a8d490",
            "Fibroblasts"="#fdc182")
ct_order <- names(ct_col)           # 变量名沿用，实际 9 个
colors_fine   <- c(ct_col, "Macrophages_Inf" = "#E76F51")
ct_fine_order <- c("Macrophages","Macrophages_Inf","Neutrophils","T_NK_cells","Plasma_cells",
                   "B_cells","Epithelial_cells","AT2_cells","Endothelial_cells","Fibroblasts")


cat("embeddings:", paste(names(proj@embeddings), collapse = ","), "\n")
umap_df <- data.frame(
  UMAP_1   = proj@embeddings$UMAP$df[, 1],
  UMAP_2   = proj@embeddings$UMAP$df[, 2],
  celltype = factor(proj$celltype, levels = ct_order),
  fine     = factor(proj$celltype_fine, levels = ct_fine_order),
  group    = factor(proj$group, levels = c("Mock", "Mhp")))
x_min <- min(umap_df$UMAP_1); x_max <- max(umap_df$UMAP_1)
y_min <- min(umap_df$UMAP_2); y_max <- max(umap_df$UMAP_2)
x_range <- x_max - x_min; y_range <- y_max - y_min
x_lower <- x_min - x_range*0.08; x_upper <- x_max + x_range*0.05
y_lower <- y_min - y_range*0.08; y_upper <- y_max + y_range*0.05
arrow_len_x <- x_range*0.15; arrow_len_y <- y_range*0.15
origin_x <- x_min - x_range*0.02; origin_y <- y_min - y_range*0.02

theme_no_axes <- theme_void() + theme(
  plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
  legend.text = element_text(size = 11), legend.title = element_blank(),
  legend.key.size = unit(5, "mm"), aspect.ratio = 1, plot.margin = margin(10, 10, 10, 10))
add_L_axis <- function(p) p +
  annotate("segment", x = origin_x, y = origin_y, xend = origin_x + arrow_len_x, yend = origin_y,
           arrow = arrow(length = unit(0.2, "cm"), type = "closed"), linewidth = 1) +
  annotate("segment", x = origin_x, y = origin_y, xend = origin_x, yend = origin_y + arrow_len_y,
           arrow = arrow(length = unit(0.2, "cm"), type = "closed"), linewidth = 1) +
  annotate("text", x = origin_x + arrow_len_x/2, y = origin_y - y_range*0.04,
           label = "UMAP1", size = 4.5, fontface = "bold") +
  annotate("text", x = origin_x - x_range*0.04, y = origin_y + arrow_len_y/2,
           label = "UMAP2", angle = 90, size = 4.5, fontface = "bold")

umap_base <- function(df, col_vec) ggplot(df, aes(UMAP_1, UMAP_2, color = celltype)) +
  geom_point(size = 0.3, stroke = 0) + scale_color_manual(values = col_vec) +
  coord_fixed(1, xlim = c(x_lower, x_upper), ylim = c(y_lower, y_upper), clip = "off") + theme_no_axes


# ============= 段3：主图（9 类）+ 高亮图（那 137 个在哪）=============
p_single <- add_L_axis(umap_base(umap_df, ct_col) + theme(legend.position = "right") +
  guides(color = guide_legend(override.aes = list(size = 4))))
ggsave(file.path(PLOT_DIR2, "MHP_ATAC_UMAP_CellType.pdf"), p_single, width = 8, height = 7)
ggsave(file.path(PLOT_DIR2, "MHP_ATAC_UMAP_CellType.png"), p_single, width = 8, height = 7, dpi = 300)

hl <- umap_df
hl$hl <- ifelse(proj$celltype_fine == "Macrophages_Inf", "Macrophages_Inf",
         ifelse(proj$celltype == "Macrophages", "Macrophages", "Other"))
hl$hl <- factor(hl$hl, levels = c("Other", "Macrophages", "Macrophages_Inf"))
hl <- hl[order(hl$hl), ]                 # 橘色最后画，压在最上层
p_hl <- ggplot(hl, aes(UMAP_1, UMAP_2, color = hl)) +
  geom_point(size = 0.35, stroke = 0) +
  scale_color_manual(values = c(Other = "grey88", Macrophages = "#f69886", Macrophages_Inf = "#E76F51")) +
  coord_fixed(1, xlim = c(x_lower, x_upper), ylim = c(y_lower, y_upper), clip = "off") +
  ggtitle("Macrophages_Inf (n = 137)") + theme_no_axes + theme(legend.position = "right") +
  guides(color = guide_legend(override.aes = list(size = 4)))
p_hl <- add_L_axis(p_hl)
ggsave(file.path(PLOT_DIR2, "MHP_ATAC_UMAP_MacInf_highlight.pdf"), p_hl, width = 8, height = 7)
ggsave(file.path(PLOT_DIR2, "MHP_ATAC_UMAP_MacInf_highlight.png"), p_hl, width = 8, height = 7, dpi = 300)


# ============= 段4：按 Mock / Mhp 分面 =============
p_all <- add_L_axis(umap_base(umap_df, ct_col) + ggtitle("All") + theme(legend.position = "right") +
  guides(color = guide_legend(override.aes = list(size = 3.5))))
plot_grp <- function(g) umap_base(umap_df[umap_df$group == g, ], ct_col) +
  ggtitle(sprintf("%s (n = %d)", g, sum(umap_df$group == g))) + guides(color = "none")
p_split <- p_all + plot_grp("Mock") + plot_grp("Mhp") +
  plot_layout(ncol = 3, guides = "collect") & theme(legend.position = "right")
ggsave(file.path(PLOT_DIR2, "MHP_ATAC_UMAP_SplitByGroup.pdf"), p_split, width = 16, height = 5.5)
ggsave(file.path(PLOT_DIR2, "MHP_ATAC_UMAP_SplitByGroup.png"), p_split, width = 16, height = 5.5, dpi = 300)


# ============= 段5：kNN 定量（Harmony 空间，聚类真正依据的空间）=============
# 137 个 Macrophages_Inf 各自 30 近邻中巨噬占比 vs 137 个随机普通巨噬的同一指标
rd <- getReducedDims(proj, reducedDims = "Harmony")
idx_inf <- which(proj$celltype_fine == "Macrophages_Inf")
knn <- ArchR:::.computeKNN(data = rd[-idx_inf, ], query = rd[idx_inf, ], k = 30)
nb  <- matrix(proj$celltype_fine[-idx_inf][knn], nrow = length(idx_inf))
frac_inf <- rowMeans(nb == "Macrophages")
cat("\n[Macrophages_Inf] 30 近邻中巨噬占比：\n"); print(summary(frac_inf))
cat("占比>=0.5 的细胞:", sum(frac_inf >= 0.5), "/", length(frac_inf), "\n")
cat("近邻类型总体分布：\n"); print(round(sort(table(as.vector(nb)) / length(nb), decreasing = TRUE), 3))

set.seed(42); idx_ctl <- sample(which(proj$celltype_fine == "Macrophages"), length(idx_inf))
knn_c <- ArchR:::.computeKNN(data = rd[-idx_ctl, ], query = rd[idx_ctl, ], k = 30)
nb_c  <- matrix(proj$celltype_fine[-idx_ctl][knn_c], nrow = length(idx_ctl))
frac_ctl <- rowMeans(matrix(nb_c %in% c("Macrophages", "Macrophages_Inf"), nrow = nrow(nb_c)))  # %in% 会拍平矩阵，必须重建维度
cat("\n[对照: 随机 137 个普通巨噬] 30 近邻中巨噬占比：\n"); print(summary(frac_ctl))
write.csv(data.frame(cell = rownames(rd)[idx_inf], frac_mac_neighbors = round(frac_inf, 3)),
          file.path(PLOT_DIR2, "MacInf_kNN_macrophage_fraction.csv"), row.names = FALSE)

cat("\n✅ MHP 04 完成 | 图在", PLOT_DIR2, "\n")


ct_col <- c("Macrophages"="#f69886","Neutrophils"="#d984b7","T_NK_cells"="#ef3f33","Plasma_cells"="#895da7",
            "B_cells"="#26a870","Epithelial_cells"="#d4c1ee","Endothelial_cells"="#a8d490","Fibroblasts"="#fdc182")
ct_order <- names(ct_col)

theme_no_axes <- theme_void() + theme(
  plot.title = element_text(face = "bold", hjust = 0.5, size = 20),
  legend.text = element_text(size = 14), legend.title = element_blank(),
  legend.key.size = unit(7, "mm"), aspect.ratio = 1, plot.margin = margin(10, 10, 10, 10))

# 坐标轴参数写到全局，切换 embedding 时重算一次
set_axis <- function(df) {
  x_min <<- min(df$UMAP_1); x_max <<- max(df$UMAP_1); y_min <<- min(df$UMAP_2); y_max <<- max(df$UMAP_2)
  x_range <<- x_max - x_min; y_range <<- y_max - y_min
  x_lower <<- x_min - x_range*0.08; x_upper <<- x_max + x_range*0.05
  y_lower <<- y_min - y_range*0.08; y_upper <<- y_max + y_range*0.05
  arrow_len_x <<- x_range*0.15; arrow_len_y <<- y_range*0.15
  origin_x <<- x_min - x_range*0.02; origin_y <<- y_min - y_range*0.02
  invisible(NULL)
}
add_L_axis <- function(p, ax = "UMAP") p +
  annotate("segment", x=origin_x, y=origin_y, xend=origin_x+arrow_len_x, yend=origin_y,
           arrow=arrow(length=unit(0.25,"cm"), type="closed"), linewidth=1.2) +
  annotate("segment", x=origin_x, y=origin_y, xend=origin_x, yend=origin_y+arrow_len_y,
           arrow=arrow(length=unit(0.25,"cm"), type="closed"), linewidth=1.2) +
  annotate("text", x=origin_x+arrow_len_x/2, y=origin_y-y_range*0.04, label=paste0(ax,"1"), size=5.5, fontface="bold") +
  annotate("text", x=origin_x-x_range*0.04, y=origin_y+arrow_len_y/2, label=paste0(ax,"2"), angle=90, size=5.5, fontface="bold")
# df 需含 UMAP_1/UMAP_2/celltype/group 四列（t-SNE 也沿用这两个列名，省得改函数）
mk_umap <- function(df, ttl = NULL, leg = TRUE, pt = 0.7, lab = FALSE) {
  p <- ggplot(df, aes(UMAP_1, UMAP_2, color = celltype)) + geom_point(size = pt, stroke = 0) +
    scale_color_manual(values = ct_col, drop = FALSE) +
    coord_fixed(1, xlim = c(x_lower, x_upper), ylim = c(y_lower, y_upper), clip = "off") + theme_no_axes
  if (lab) {
    cen <- df %>% group_by(celltype) %>% summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), .groups = "drop")
    p <- p + geom_label(data = cen, aes(label = celltype), color = "black", fill = alpha("white", 0.75),
                        size = 4.5, fontface = "bold", label.size = 0, label.padding = unit(0.15, "lines"), show.legend = FALSE)
  }
  if (!is.null(ttl)) p <- p + ggtitle(ttl)
  if (leg) p + theme(legend.position = "right") + guides(color = guide_legend(override.aes = list(size = 5))) else p + guides(color = "none")
}
plot_set <- function(df, tag, ax) {   # tag = "UMAP"/"TSNE"，一次出 单张 / 带标签 / 分组三合一
  set_axis(df)
  p1 <- add_L_axis(mk_umap(df), ax); p2 <- add_L_axis(mk_umap(df, lab = TRUE), ax)
  ggsave(file.path(PLOT_DIR, paste0("MHP_ATAC_", tag, "_CellType.png")), p1, width=9, height=7.5, dpi=300, bg="white")
  ggsave(file.path(PLOT_DIR, paste0("MHP_ATAC_", tag, "_CellType.pdf")), p1, width=9, height=7.5, bg="white")
  ggsave(file.path(PLOT_DIR, paste0("MHP_ATAC_", tag, "_CellType_labeled.png")), p2, width=9, height=7.5, dpi=300, bg="white")
  ggsave(file.path(PLOT_DIR, paste0("MHP_ATAC_", tag, "_CellType_labeled.pdf")), p2, width=9, height=7.5, bg="white")
  ps <- add_L_axis(mk_umap(df, "All", pt=0.6), ax) +
    mk_umap(df[df$group=="Mock",], sprintf("Mock (n = %d)", sum(df$group=="Mock")), FALSE, pt=0.6) +
    mk_umap(df[df$group=="Mhp",],  sprintf("Mhp (n = %d)",  sum(df$group=="Mhp")),  FALSE, pt=0.6) +
    plot_layout(ncol = 3, guides = "collect") & theme(legend.position = "right")
  ggsave(file.path(PLOT_DIR, paste0("MHP_ATAC_", tag, "_SplitByGroup.png")), ps, width=18, height=6, dpi=300, bg="white")
  ggsave(file.path(PLOT_DIR, paste0("MHP_ATAC_", tag, "_SplitByGroup.pdf")), ps, width=18, height=6, bg="white")
}

umap_df <- data.frame(UMAP_1 = proj@embeddings$UMAP$df[,1], UMAP_2 = proj@embeddings$UMAP$df[,2],
                      celltype = factor(proj$celltype, levels = ct_order), group = factor(proj$group, levels = c("Mock","Mhp")))
plot_set(umap_df, "UMAP", "UMAP")

if (!"TSNE" %in% names(proj@embeddings)) {
  proj <- addTSNE(ArchRProj = proj, reducedDims = "Harmony", name = "TSNE", perplexity = 30, force = TRUE)
  saveArchRProject(proj, OUT_DIR, load = FALSE)
}
tsne_df <- data.frame(UMAP_1 = proj@embeddings$TSNE$df[,1], UMAP_2 = proj@embeddings$TSNE$df[,2],
                      celltype = factor(proj$celltype, levels = ct_order), group = factor(proj$group, levels = c("Mock","Mhp")))
plot_set(tsne_df, "TSNE", "t-SNE")
set_axis(umap_df)   # 切回 UMAP 轴参数
cat("✅ 04 出图完成 →", PLOT_DIR, "\n")
