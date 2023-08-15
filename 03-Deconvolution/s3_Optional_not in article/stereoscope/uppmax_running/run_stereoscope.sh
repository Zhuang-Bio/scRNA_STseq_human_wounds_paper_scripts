#!/bin/bash -l
#SBATCH -A snic2021-22-701
#SBATCH -p node -n 19
#SBATCH -t 10-00:00:00
#SBATCH -J stereo
#SBATCH -e stereo_SLURM_Job_id%j.sderr.txt
#SBATCH -o stereo_SLURM_Job_id%j.sdout.txt
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zhuang.liu@ki.se

cd /proj/snic2021-23-156/stereoscope_deconv 

conda activate scvi-env

python stereoscope_run.py

