# ============================================================================
# config/paths.R — 全部 R 脚本共用的路径变量（各脚本开头 source 本文件）
# 用法：按本机 / 服务器的实际位置设置同名环境变量（推荐），或直接修改下面的默认值
# ============================================================================
.mhp_path <- function(var, default) { v <- Sys.getenv(var, unset = ""); if (nzchar(v)) v else default }

MHP_CODE       <- normalizePath(.mhp_path("MHP_CODE", "."), mustWork = TRUE)              # 本仓库根目录
MHP_RES        <- .mhp_path("MHP_RES",        file.path(MHP_CODE, "results"))             # 本地工作目录：输入表、面板与终稿的读写根目录
MHP_ARCHR      <- .mhp_path("MHP_ARCHR",      "/path/to/ArchR_MHP")                       # 服务器：ArchR 项目目录（arrows/、ArchR_filtered/、outputs/）
MHP_FRAG       <- .mhp_path("MHP_FRAG",       "/path/to/fragments")                       # 服务器：{Mock,Mhp}_{1,2,3}.fragments.tsv.gz 所在目录
MHP_ARCHR_ANNO <- .mhp_path("MHP_ARCHR_ANNO", "/path/to/ArchR_annotation")                # 服务器：Sscrofa11.1 的 ArchR 注释（genomeAnnotation.rds、geneAnnotation_symbol.rds、gene_id2sym.rds）
MHP_RNA_OBJ    <- .mhp_path("MHP_RNA_OBJ",    "/path/to/rna_final_annotated_symbol.qs")   # 服务器：已注释的 scRNA-seq Seurat 对象（qs 格式）
MHP_RAW        <- .mhp_path("MHP_RAW",        "/path/to/bulk_RNA_tracks")                 # bulk RNA-seq 轨迹文件目录（04.RNA/{Con,Mhp}{1,2,3}_RNA.bw）
MHP_ROSE       <- .mhp_path("MHP_ROSE",       file.path(MHP_RES, "6、SE/ROSE"))           # ROSE 超级增强子结果（{Con,Mhp}_SuperEnhancers.table.txt / .bed / _ENHANCER_TO_GENE.txt 等）
MHP_TE         <- .mhp_path("MHP_TE",         file.path(MHP_ROSE, "TE"))                  # ROSE 全部增强子的基因注释（{Con,Mhp}_AllEnhancers_ENHANCER_TO_GENE.txt）
MHP_HIC        <- .mhp_path("MHP_HIC",        file.path(MHP_RES, "7、HI-C结果/HiC_result")) # Hi-C 结果根目录（single_result/02.align/*.cool、single_result/05.loop/*.bedpe）
