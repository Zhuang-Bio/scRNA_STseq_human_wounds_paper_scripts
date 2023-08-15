#! /bin/bash -l
#SBATCH -A sens2020010
#SBATCH -p core -n 16
#SBATCH -t 5-00:00:00
#SBATCH -J P20063_102
#SBATCH -e P20063_102.SLURM_Job_id=%j.sderr
#SBATCH -o P20063_102.SLURM_Job_id=%j.sdout
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zhuang.liu@ki.se

module load bioinfo-tools spaceranger/1.2.0

cd /proj/sens2020010/proj_STseqData/Donor_2/spaceranger-latestRef-manualAlignment

SAMPLE_IDS="P20063_102"
TRANSCRIPTOME="/proj/sens2020010/ref10X_GRCh38_2021Jul/GRCh38"
FASTQS="/proj/sens2020010/proj_STseqData/Donor_2/P20063_102/02-FASTQ/210304_A00187_0447_BH3L7JDRXY/"
IMAGES="/proj/sens2020010/proj_STseqData/Donor_2/01-Images/210204_P20063_V10B01-24_LP1_HE_AM_B1-Spots/210204_P20063_V10B01-24_LP1_HE_AM_B1-Spot000001.jpg"

spaceranger count \
--id=$SAMPLE_IDS \
--transcriptome=$TRANSCRIPTOME \
--fastqs=$FASTQS \
--sample=$SAMPLE_IDS \
--image=$IMAGES \
--slide=V10B01-024 \
--area=B1 \
--slidefile=/proj/sens2020010/proj_STseqData/Donor_2/00-slidefile/V10B01-024.gpr \
--loupe-alignment=/proj/sens2020010/proj_STseqData/manual_alignment_files/V10B01-024-B1.json \
--localcores=12 \
--localmem=100
