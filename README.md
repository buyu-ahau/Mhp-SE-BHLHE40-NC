# Mhp-SE-BHLHE40

Analysis code for **"An infection-augmented super-enhancer element sustains BHLHE40-dependent responses of porcine macrophages to *Mycoplasma hyopneumoniae*"** (Bu et al.).

This repository holds the custom analysis and figure code for the six main figures, one folder per figure. Each folder contains the core scripts for that figure; preview and layout-only scripts are not included. Read alignment, peak calling, super-enhancer calling (ROSE), Hi-C contact-matrix generation and loop calling used standard tools with the parameters given in the Methods.

## Data

Raw sequencing data (scRNA-seq, snATAC-seq, bulk RNA-seq, histone and BHLHE40 CUT&Tag, and Hi-C) are available from the NCBI Sequence Read Archive under BioProject [PRJNA1475599](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1475599).

## Contents

Scripts marked † run on a Linux server (ArchR or pySCENIC environment); all other scripts run locally in R.

| Folder | Figure | Scripts |
|---|---|---|
| `01_Fig1_atlas` | Fig. 1, scRNA-seq and snATAC-seq atlas | `01_createArrows.R`† (Arrow files, quality control, ArchR project); `04_assign_celltype_and_plot.R`† (snATAC-seq annotation); `03_导出_…_服务器.R`† (export of embeddings, gene scores and pseudobulk expression); `04_Fig1CD_…` (Fig. 1b,d); `05_Fig1E_…` (Fig. 1e,f); `06_Fig1F_…` (Fig. 1g) |
| `02_Fig2_TF_programs` | Fig. 2, transcription-factor programs | `02a_…`† and `02c_…`† (pySCENIC GRN inference, cisTarget pruning and AUCell); `05_SCENIC_B_…` (Fig. 2b); `06_SCENIC_D_…` (Fig. 2e); `07_infection_DA.R`† (differential gene scores, peaks and chromVAR deviations); `07e_motif_jaspar_bgd.R`† (motif enrichment in differentially accessible regions, Fig. 2c) |
| `03_Fig3_BHLHE40_knockdown` | Fig. 3, BHLHE40 knockdown | `01_…` (quantitative panels from raw replicate values); `04_…` (flow cytometry plots from .fcs files, Fig. 3h,i); `11_…` (Supplementary Fig. 5) |
| `04_Fig4_super_enhancers` | Fig. 4, super-enhancers | `00b_…` (pairing of Con and Mhp super-enhancers); `02_d_…` (Fig. 4d); `03_e_…` (Fig. 4e); `04_b_…` (bulk RNA-seq DESeq2, Fig. 4b); `08_g_…` (cell-type attribution by snATAC-seq, Fig. 4g) |
| `05_Fig5_E2_locus` | Fig. 5, *BHLHE40* locus and enhancer E2 | `01_a_…` (Hi-C loop classes, Fig. 5a); `07_4C_…` (virtual 4C, Fig. 5b); `02_b_…` (locus tracks, Fig. 5c); `05_e_…` (peak-to-gene peaks in E2, Fig. 5f); `09_g_…` (BHLHE40 CUT&Tag at E2, Fig. 5h) |
| `06_Fig6_E2_CRISPRi` | Fig. 6, E2 CRISPRi | `04_…` (CCK-8, Fig. 6d); `05_…` (apoptosis rates from .fcs files, Fig. 6e); `09_…` and `10_…` (Transwell image segmentation and statistics, Fig. 6f); helpers `_fig6_common_本地.R` and `_transwell_seg_本地.R` |
| `shared` | All | Color palettes, ggplot2 theme and figure-export functions |
| `config` | All | Path configuration (`paths.R`) |

Panel letters in some file names refer to an earlier figure layout; the header of each script gives the current panel.

## Usage

1. Set the data locations in `config/paths.R`, or export environment variables with the same names:
   - `MHP_RES`: local working directory for input tables and figure output (default: `results/` in the repository).
   - Sequencing-provider outputs read by the Fig. 4 and Fig. 5 scripts: `MHP_ROSE` (ROSE super-enhancer calls), `MHP_TE` (ROSE enhancer-to-gene tables for all enhancers) and `MHP_HIC` (Hi-C contact matrices and loop calls). Defaults are subfolders of `MHP_RES`.
   - Server scripts: `MHP_ARCHR` (ArchR project), `MHP_FRAG` (fragment files), `MHP_ARCHR_ANNO` (ArchR genome and gene annotation for Sscrofa11.1), `MHP_RNA_OBJ` (annotated scRNA-seq Seurat object), `MHP_RAW` (bulk RNA-seq bigWig tracks), and `MHP_SCENIC` and `MHP_SCENIC_REF` (pySCENIC working directory and cisTarget databases).
2. Run scripts from the repository root, or set `MHP_CODE` to the repository path. Each R script begins with `source(file.path(Sys.getenv("MHP_CODE", "."), "config", "paths.R"))`.
3. Plotting scripts write 300-dpi PNG previews by default; set `MHP_FINAL=1` to write the final SVG and 900-dpi TIFF files.

Paths under `MHP_RES` keep the directory layout of the original project. File names and code comments are in Chinese.

## Software

R 4.6.1; Seurat v5 (scRNA-seq) and v4.4.0 (ArchR integration); ArchR 1.0.3; chromVAR 1.30.1; pySCENIC 0.12.1; DESeq2 1.52.0; edgeR 4.10.3; ChIPseeker 1.48.0; clusterProfiler 4.20.0; flowCore 2.24.0; cooler 0.10.2; MACS2 2.2.9.1; Cell Ranger 7.1.0; Cell Ranger ATAC 2.1.0. See the Methods for details.

## License

MIT License; see `LICENSE`.
