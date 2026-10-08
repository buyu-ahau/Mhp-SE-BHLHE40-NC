#!/bin/bash
# ============================================================================
# 02c_pySCENIC_CTX_AUCell_新服务器.sh — pySCENIC 第 2、3 步：cisTarget 基序剪枝（hg38 v10 排序库）与 AUCell 调控子活性打分
# 运行位置：服务器（conda 环境 pyscenic）；需先设置环境变量 MHP_SCENIC 与 MHP_SCENIC_REF
# 输入：02a 的输出 grn_output.tsv 与 MHP_scenic_input.loom；输出 output/d1_reg.csv（调控子）与 output/d1_auc.loom（AUCell）
# ============================================================================
D="${MHP_SCENIC:?请先设置 MHP_SCENIC（pySCENIC 工作目录，即 02a 的输出目录）}"
REF="${MHP_SCENIC_REF:?请先设置 MHP_SCENIC_REF（cisTarget 排序库与注释表所在目录）}"
O="$D/output"
mkdir -p "$O" && cd "$D" || exit 1

pyscenic ctx grn_output.tsv \
    "$REF/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather" \
    "$REF/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather" \
    --annotations_fname "$REF/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl" \
    --expression_mtx_fname MHP_scenic_input.loom \
    --num_workers 96 --output "$O/d1_reg.csv" \
&& pyscenic aucell MHP_scenic_input.loom "$O/d1_reg.csv" \
    --output "$O/d1_auc.loom" --num_workers 96 --seed 42

ls -lh "$O"
