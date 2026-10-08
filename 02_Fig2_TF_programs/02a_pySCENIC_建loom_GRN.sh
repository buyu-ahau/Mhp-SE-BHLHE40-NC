#!/bin/bash
# ============================================================================
# 02a_pySCENIC_建loom_GRN.sh — Fig. 2a、2b、2e 的上游：pySCENIC 第 1 步，由表达矩阵建 loom 并用 GRNBoost2 推断 TF–靶基因共表达网络
# 运行位置：服务器（conda 环境 pyscenic）；需先设置环境变量 MHP_SCENIC 与 MHP_SCENIC_REF
# 输入：MHP_counts_human.csv（细胞 × 人同源基因计数矩阵）；输出 MHP_scenic_input.loom 与 grn_output.tsv
# ============================================================================
REF="${MHP_SCENIC_REF:?请先设置 MHP_SCENIC_REF（allTFs_hg38.txt 所在目录）}"
DIR="${MHP_SCENIC:?请先设置 MHP_SCENIC（pySCENIC 工作目录，含 MHP_counts_human.csv）}"
ORG="MHP"
cd "$DIR" || { echo "[!] $DIR 不存在"; exit 1; }

echo "== [$ORG] 建 loom =="
python3 -c "
import pandas as pd, loompy, numpy as np
df = pd.read_csv('${ORG}_counts_human.csv', index_col=0)
print('矩阵(细胞×基因):', df.shape)
loompy.create('${ORG}_scenic_input.loom', df.values.T,
              row_attrs={'Gene': np.array(df.columns.tolist())},
              col_attrs={'CellID': np.array(df.index.tolist())})
print('loom OK')
"

echo "== [$ORG] Step1/3 GRN (GRNBoost2) =="
pyscenic grn "$DIR/${ORG}_scenic_input.loom" "$REF/allTFs_hg38.txt" \
    -o "$DIR/grn_output.tsv" --method grnboost2 --num_workers 48 --seed 42

echo "[DONE] GRN"; ls -lh "$DIR"
