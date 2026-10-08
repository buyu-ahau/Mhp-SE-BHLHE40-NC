# ============================================================================
# 03_导出_全基因伪bulk_GeneScore_校正前UMAP_服务器.R — Fig. 1 与补充图 1–2 的输入：导出 scRNA-seq 伪 bulk 表达、snATAC-seq 基因活性分数与嵌入坐标
# 运行位置：服务器（conda 环境 Archr）；路径变量见 config/paths.R
# ============================================================================

source(file.path(Sys.getenv("MHP_CODE", "."), "config", "paths.R"))

suppressPackageStartupMessages({ library(qs); library(Seurat); library(Matrix) })
a <- qread(MHP_RNA_OBJ, nthreads = 8)
a <- NormalizeData(a, verbose = FALSE)

# ---- R1：全基因 × 11 类（counts 之和 / 表达比例 / DotPlot 口径均值）----
OUT <- file.path(MHP_ARCHR, "outputs/30_Fig1v2_RNA导出")
id2sym <- readRDS(file.path(MHP_ARCHR_ANNO, "gene_id2sym.rds"))
s2i <- tapply(names(id2sym), id2sym, paste, collapse = ";")                  # Symbol → Ensembl（一对多用分号）
g <- rownames(a); ens <- ifelse(grepl("^ENSSSCG", g), g, unname(s2i[g]))
cat("基因", length(g), "| 有 Ensembl:", sum(!is.na(ens)), "| 无对应:", sum(is.na(ens)), "| 一对多:", sum(grepl(";", ens)), "\n")
w <- function(m, fn) data.table::fwrite(data.frame(ensembl_id = ens, symbol = g, m, check.names = FALSE), file.path(OUT, fn))
cnt <- LayerData(a, assay = "RNA", layer = "counts"); dat <- LayerData(a, assay = "RNA", layer = "data")
f1 <- factor(a$celltype_new); M1 <- Matrix::t(Matrix::fac2sparse(f1)); n1 <- as.numeric(table(f1))    # 细胞 × 类 指示矩阵
pb1  <- as.matrix(cnt %*% M1)
pct1 <- sweep(as.matrix(((cnt > 0) * 1) %*% M1), 2, n1, "/")
avg1 <- sweep(as.matrix(expm1(dat) %*% M1), 2, n1, "/")
stopifnot(identical(colnames(pb1), levels(f1)), nrow(pb1) == 33029)
w(pb1, "RNA_pb_counts_byCelltype.csv.gz"); w(round(pct1, 5), "RNA_pct_byCelltype.csv.gz"); w(round(avg1, 5), "RNA_avgexp_byCelltype.csv.gz")
cpm1 <- sweep(pb1, 2, colSums(pb1), "/") * 1e6
print(round(cpm1[c("CD163", "CSF3R", "FLT3", "CD4", "CD8A", "NKG7", "CD79A", "JCHAIN", "SFTPC", "CDH5", "COL1A1"), ], 1))

# ---- R2：全基因 × 类×样本（66 组，备样本层面分析）----
f2 <- factor(paste(a$celltype_new, a$orig.ident, sep = "|")); M2 <- Matrix::t(Matrix::fac2sparse(f2)); n2 <- as.numeric(table(f2))
pb2  <- as.matrix(cnt %*% M2); pct2 <- sweep(as.matrix(((cnt > 0) * 1) %*% M2), 2, n2, "/")
cat("类×样本组数:", ncol(pb2), "| 最少细胞:", min(n2), levels(f2)[which.min(n2)], "\n")
w(pb2, "RNA_pb_counts_byCelltypeSample.csv.gz"); w(round(pct2, 5), "RNA_pct_byCelltypeSample.csv.gz")
data.table::fwrite(data.frame(group = c(levels(f1), levels(f2)), level = rep(c("celltype", "celltype|sample"), c(nlevels(f1), nlevels(f2))), n_cells = c(n1, n2)), file.path(OUT, "RNA_pb_ncells.csv"))

# ---- R3：每细胞 QC + Harmony 校正前 UMAP（Supp Fig 1C）+ PCA 前 50 维（本地可重跑 UMAP）----
cat("pca vs pca_RNA PC1 相关:", round(cor(Embeddings(a, "pca")[, 1], Embeddings(a, "pca_RNA")[, 1]), 4), "| pca 维数:", ncol(Embeddings(a, "pca")), "| harmony 维数:", ncol(Embeddings(a, "harmony_RNA")), "\n")
a <- RunUMAP(a, reduction = "pca", dims = 1:30, reduction.name = "umap_pre", reduction.key = "umapPre_", seed.use = 42, verbose = FALSE)
qc <- data.frame(cell = colnames(a), a@meta.data[, c("orig.ident", "Groups", "celltype_new", "nCount_RNA", "nFeature_RNA", "percent.mito", "percent.hgb", "Phase", "seurat_clusters")], Embeddings(a, "umap_pre")[, 1:2], row.names = NULL)
stopifnot(nrow(qc) == 35330, !anyNA(qc)); data.table::fwrite(qc, file.path(OUT, "RNA_cellQC_umapPre.csv.gz")); print(head(qc, 3))
npc <- min(50, ncol(Embeddings(a, "pca")))
data.table::fwrite(data.frame(cell = colnames(a), round(Embeddings(a, "pca")[, 1:npc], 5), row.names = NULL), file.path(OUT, "RNA_pca_embeddings.csv.gz"))
print(file.info(list.files(OUT, full.names = TRUE))[, "size", drop = FALSE])

# ======================== ATAC（conda Archr，新开一个 R）========================
# ---- A1：载入项目，看基因注释列 ----
suppressPackageStartupMessages({ library(ArchR); library(BSgenome.Sscrofa.Ensembl.v11.1) })
addArchRChrPrefix(chrPrefix = FALSE); addArchRThreads(threads = 32); set.seed(42)
proj <- loadArchRProject(file.path(MHP_ARCHR, "ArchR_filtered"), showLogo = FALSE)
OUT <- file.path(MHP_ARCHR, "outputs/31_Fig1v2_ATAC导出"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
cc <- getCellColData(proj); cat("细胞:", nrow(cc), "\n"); print(table(cc$celltype, cc$group))
gn <- as.data.frame(getGenes(proj)); cat("基因注释列:", colnames(gn), "| 基因数:", nrow(gn), "\n")
id2sym <- readRDS(file.path(MHP_ARCHR_ANNO, "gene_id2sym.rds"))
s2i <- tapply(names(id2sym), id2sym, paste, collapse = ";")
s2i <- tapply(gn$gene_id, gn$symbol, paste, collapse = ";"); cat("注释内 Symbol 数:", length(s2i), "| 一对多:", sum(grepl(";", s2i)), "\n")

# ---- A2：全基因 GeneScore × 8 类均值（Fig 1E 右半）----
se1 <- getGroupSE(proj, useMatrix = "GeneScoreMatrix", groupBy = "celltype", divideN = TRUE)
rd <- as.data.frame(rowData(se1)); cat("rowData 列:", colnames(rd), "| 基因:", nrow(rd), "| 重名:", sum(duplicated(rd$name)), "\n")
sym <- rd$name; ens <- ifelse(grepl("^ENSSSCG", sym), sym, unname(s2i[sym]))
wa <- function(m, fn) data.table::fwrite(data.frame(ensembl_id = ens, symbol = sym, chr = if (is.null(rd$seqnames)) NA else as.character(rd$seqnames), start = if (is.null(rd$start)) NA else rd$start, round(as.matrix(m), 5), check.names = FALSE), file.path(OUT, fn))
gs1 <- assay(se1); rownames(gs1) <- sym                                    # ★行名取自同一 se 的 rowData，不按外部顺序索引
mk8 <- c("CD163", "CSF3R", "CD8A", "CD79A", "MZB1", "SFTPC", "CDH5", "COL1A1")
print(data.frame(gene = mk8, top = colnames(gs1)[apply(gs1[mk8, ], 1, which.max)]))   # 预期：巨噬/中性粒/T_NK/B/浆/上皮/内皮/成纤维
cat("各类列和(应≈1e4):", round(colSums(gs1)), "\n")
wa(gs1, "ATAC_GeneScore_mean_byCelltype.csv.gz")

# ---- A3：组×类（16 列）与 样本×类（47 列）----
se2 <- getGroupSE(proj, useMatrix = "GeneScoreMatrix", groupBy = "GroupCellType", divideN = TRUE)
proj <- addCellColData(proj, data = paste(cc$Sample, cc$celltype, sep = "|"), name = "SampleCT", cells = rownames(cc), force = TRUE)   # ★只在内存加列，不要 saveArchRProject
se3 <- getGroupSE(proj, useMatrix = "GeneScoreMatrix", groupBy = "SampleCT", divideN = TRUE)
stopifnot(identical(rowData(se2)$name, sym), identical(rowData(se3)$name, sym))
wa(assay(se2), "ATAC_GeneScore_mean_byGroupCelltype.csv.gz"); wa(assay(se3), "ATAC_GeneScore_mean_bySampleCelltype.csv.gz")
tb <- list(table(cc$celltype), table(cc$GroupCellType), table(proj$SampleCT))
data.table::fwrite(data.frame(group = unlist(lapply(tb, names)), n_cells = unlist(lapply(tb, as.integer))), file.path(OUT, "ATAC_ncells.csv"))
cat("GroupCellType 列:", ncol(se2), "| SampleCT 列:", ncol(se3), "\n")

# ---- A4：GeneScore>0 细胞比例 × 8 类（可选，图谱 08 未用）----
m <- getMatrixFromProject(proj, useMatrix = "GeneScoreMatrix"); X <- assay(m); stopifnot(identical(rowData(m)$name, sym))
ct <- cc[colnames(X), "celltype"]; stopifnot(!anyNA(ct))                   # ★getMatrixFromProject 列序 ≠ cellColData 行序，按细胞名取
Mc <- Matrix::t(Matrix::fac2sparse(factor(ct))); pctA <- sweep(as.matrix(((X > 0) * 1) %*% Mc), 2, as.numeric(table(factor(ct))), "/")
wa(pctA, "ATAC_GeneScore_pct_byCelltype.csv.gz"); rm(m, X); invisible(gc())

# ---- A5：全部降维坐标（带细胞名）+ 文件清单 ----
cat("embeddings:", names(proj@embeddings), "\n")
emb <- do.call(cbind, lapply(names(proj@embeddings), function(k) setNames(getEmbedding(proj, embedding = k, returnDF = TRUE)[rownames(cc), 1:2], paste0(k, "_", 1:2))))
data.table::fwrite(data.frame(cell = rownames(cc), Sample = cc$Sample, group = cc$group, celltype = cc$celltype, celltype_fine = cc$celltype_fine, emb, row.names = NULL), file.path(OUT, "ATAC_embeddings_all.csv.gz"))
print(file.info(list.files(OUT, full.names = TRUE))[, "size", drop = FALSE])
