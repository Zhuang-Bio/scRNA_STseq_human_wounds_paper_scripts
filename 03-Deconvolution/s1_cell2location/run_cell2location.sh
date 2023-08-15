#!/bin/bash -l
#SBATCH -A snic2021-22-701
#SBATCH -p node -n 18
#SBATCH -t 10-00:00:00
#SBATCH -J n5
#SBATCH -e n5_SLURM_Job_id%j.sderr.txt
#SBATCH -o n5_SLURM_Job_id%j.sdout.txt
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zhuang.liu@ki.se

cd /proj/snic2021-23-156/cell2location_0704

conda activate cell2loc_scvi0418

python cell2location_run.py

