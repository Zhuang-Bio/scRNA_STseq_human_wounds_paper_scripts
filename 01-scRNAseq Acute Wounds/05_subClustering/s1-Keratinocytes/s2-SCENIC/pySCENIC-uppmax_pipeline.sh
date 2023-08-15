#!/bin/bash -l
#SBATCH -A snic2021-22-701
#SBATCH -p node -n 20
#SBATCH -t 5-00:00:00
#SBATCH -J allmmpyscenic
#SBATCH -e allmmpyscenic_SLURM_Job_id%j.sderr
#SBATCH -o allmmpyscenic_SLURM_Job_id%j.sdout
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zhuang.liu@ki.se

cd /crex/proj/snic2021-23-156/proj_10XscRNAseq/proj_human/pySCENIC

conda activate pyscenic

arboreto_with_multiprocessing.py \
hswound_pySCENIC.loom \
hs_hgnc_tfs.txt \
--method grnboost2 \
--output hswound_adj_all.tsv \
--num_workers 20 \
--seed 777

pyscenic ctx \
hswound_adj_all.tsv \
hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather \
--annotations_fname motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
--expression_mtx_fname hswound_pySCENIC.loom \
--output hswound_reg.csv \
--mode "dask_multiprocessing" \
--mask_dropouts \
--num_workers 20

pyscenic aucell \
hswound_pySCENIC.loom \
hswound_reg.csv \
--output hswound_SCENIC_AUC.loom \
--num_workers 20

