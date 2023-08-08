#! /bin/bash -l
#SBATCH -A sens2020010
#SBATCH -p core -n 16
#SBATCH -t 5-00:00:00
#SBATCH -J hs38CR2
#SBATCH -e hs38CR2.SLURM_Job_id=%j.sderr
#SBATCH -o hs38CR2.SLURM_Job_id=%j.sdout
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zhuang.liu@ki.se

module load bioinfo-tools cellranger/5.0.1

cd /proj/sens2020010/proj_10X_woundhealing/01_cellranger

JOBID="PWH27D0"
SAMPLE_IDS="P20852_1002"
TRANSCRIPTOME="/proj/sens2020010/ref10X_GRCh38_2021Jul/GRCh38"
FASTQS="/proj/sens2020010/proj_10X_woundhealing/P20852/P20852_1002/02-FASTQ/210531_A00621_0413_BHVFK2DSXY"

cellranger count \
--id=$JOBID \
--transcriptome=$TRANSCRIPTOME \
--fastqs=$FASTQS \
--sample=$SAMPLE_IDS \
--localcores=16 \
--localmem=90 \
--expect-cells=10000

/proj/sens2020010/proj_10X_woundhealing/P20852/P20852_1002/02-FASTQ/210531_A00621_0413_BHVFK2DSXY/P20852_1002_S10_L002_R1_001.fastq.gz \
/proj/sens2020010/proj_10X_woundhealing/P20852/P20852_1002/02-FASTQ/210531_A00621_0413_BHVFK2DSXY/P20852_1002_S10_L002_R2_001.fastq.gz