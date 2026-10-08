# ============================================================================
# 01_createArrows.R — snATAC-seq：由 fragments 文件建 Arrow（质控过滤）并创建 ArchR 项目
# 运行位置：服务器（conda 环境 Archr）；路径变量见 config/paths.R
# ============================================================================

source(file.path(Sys.getenv("MHP_CODE", "."), "config", "paths.R"))

# ============= 段A：环境与共享注释 =============
suppressPackageStartupMessages({
  library(ArchR)
  library(BSgenome.Sscrofa.Ensembl.v11.1)
})
addArchRChrPrefix(chrPrefix = FALSE)   # Ensembl 命名，务必 FALSE
addArchRThreads(threads = 12)          # 建 Arrow 阶段不宜过高（内存）
set.seed(42)

FRAG_DIR <- MHP_FRAG
PROJ_DIR <- MHP_ARCHR
SHARED   <- MHP_ARCHR_ANNO
QC_DIR   <- file.path(PROJ_DIR, "outputs", "00_QC")

dir.create(file.path(PROJ_DIR, "arrows"), recursive = TRUE, showWarnings = FALSE)
dir.create(QC_DIR, recursive = TRUE, showWarnings = FALSE)
setwd(file.path(PROJ_DIR, "arrows"))   # Arrow 就地生成，避免默认写到当前目录

genomeAnno <- readRDS(file.path(SHARED, "genomeAnnotation.rds"))
geneAnno   <- readRDS(file.path(SHARED, "geneAnnotation_symbol.rds"))

cat("染色体:", paste(seqlevels(genomeAnno$chromSizes), collapse = ","), "\n")
cat("基因数:", length(geneAnno$genes), "| TSS:", length(geneAnno$TSS), "\n")


# ============= 段B：样本表 =============
samples <- c("Mock_1", "Mock_2", "Mock_3", "Mhp_1", "Mhp_2", "Mhp_3")
frags <- setNames(file.path(FRAG_DIR, paste0(samples, ".fragments.tsv.gz")), samples)

cat("全部存在:", all(file.exists(frags)), "\n")
print(data.frame(sample = samples,
                 size_GB = round(file.size(frags) / 1e9, 2),
                 tbi = file.exists(paste0(frags, ".tbi"))))


# ============= 段C：建 Arrow（最慢，数小时）=============
ArrowFiles <- createArrowFiles(
  inputFiles      = frags,
  sampleNames     = names(frags),
  geneAnnotation  = geneAnno,
  genomeAnnotation = genomeAnno,
  minTSS          = 2,        # 宽松，后续再收紧
  minFrags        = 500,      # 宽松
  maxFrags        = 1e+05,
  addTileMat      = TRUE,
  addGeneScoreMat = TRUE,
  excludeChr      = c("MT"),
  threads         = 12,
  force           = TRUE
)
saveRDS(ArrowFiles, file.path(PROJ_DIR, "ArrowFiles.rds"))
cat("Arrow 完成:", length(ArrowFiles), "个\n")


# ============= 段D：双胞打分 + 建项目 =============
doubScores <- addDoubletScores(
  input = ArrowFiles, k = 10, knnMethod = "UMAP", LSIMethod = 1, threads = 12
)

proj <- ArchRProject(
  ArrowFiles       = ArrowFiles,
  outputDirectory  = file.path(PROJ_DIR, "ArchR"),
  geneAnnotation   = geneAnno,
  genomeAnnotation = genomeAnno,
  copyArrows       = TRUE
)

proj$group <- factor(ifelse(grepl("^Mock", proj$Sample), "Mock", "Mhp"),
                     levels = c("Mock", "Mhp"))
print(table(proj$Sample, proj$group))
cat("总细胞数:", nCells(proj), "\n")
saveArchRProject(proj, load = FALSE)


# ============= 段E：QC 分布，据此定真实阈值 =============
df <- getCellColData(proj, select = c("Sample", "group", "TSSEnrichment", "nFrags", "DoubletScore"))
df <- as.data.frame(df)

cat("\n--- TSSEnrichment 按样本 ---\n")
print(round(do.call(rbind, tapply(df$TSSEnrichment, df$Sample, summary)), 2))
cat("\n--- nFrags 按样本 ---\n")
print(round(do.call(rbind, tapply(df$nFrags, df$Sample, summary)), 0))
cat("\n--- 若按 TSS>=4 & nFrags>=1000 过滤，各样本保留数 ---\n")
print(table(df$Sample[df$TSSEnrichment >= 4 & df$nFrags >= 1000]))

write.csv(df, file.path(QC_DIR, "cellColData_stage0.csv"), row.names = TRUE)

p1 <- plotGroups(proj, groupBy = "Sample", colorBy = "cellColData",
                 name = "TSSEnrichment", plotAs = "violin", alpha = 0.4)
p2 <- plotGroups(proj, groupBy = "Sample", colorBy = "cellColData",
                 name = "log10(nFrags)", plotAs = "violin", alpha = 0.4)
p3 <- plotFragmentSizes(ArchRProj = proj)
p4 <- plotTSSEnrichment(ArchRProj = proj)
plotPDF(p1, p2, p3, p4, name = "stage0_QC", ArchRProj = proj, addDOC = FALSE,
        width = 8, height = 5)

cat("\n完成。QC 图在 ArchR/Plots/stage0_QC.pdf，明细在", QC_DIR, "\n")
