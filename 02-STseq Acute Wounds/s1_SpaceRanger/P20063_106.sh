#! /bin/bash -l
#SBATCH -A sens2020010
#SBATCH -p core -n 16
#SBATCH -t 5-00:00:00
#SBATCH -J P20063_106
#SBATCH -e P20063_106.SLURM_Job_id=%j.sderr
#SBATCH -o P20063_106.SLURM_Job_id=%j.sdout
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zhuang.liu@ki.se

module load bioinfo-tools spaceranger/1.2.0

cd /proj/sens2020010/proj_STseqData/Donor_3_4/P20063/spaceranger-latestRef-manualAlignment

SAMPLE_IDS="P20063_106"
TRANSCRIPTOME="/proj/sens2020010/ref10X_GRCh38_2021Jul/GRCh38"
FASTQS="/proj/sens2020010/proj_STseqData/Donor_3_4/P20063/P20063_106/02-FASTQ/210531_A00621_0413_BHVFK2DSXY/"
IMAGES="/proj/sens2020010/proj_STseqData/Donor_3_4/P20063/01-Images-V10B01-025-026/210311_P20063_106_V10B01-025_PWH23_D1_B1-Spots/210311_P20063_106_V10B01-025_PWH23_D1_B1-Spot000001.jpg"

spaceranger count \
--id=$SAMPLE_IDS \
--transcriptome=$TRANSCRIPTOME \
--fastqs=$FASTQS \
--sample=$SAMPLE_IDS \
--image=$IMAGES \
--slide=V10B01-025 \
--area=B1 \
--slidefile=/proj/sens2020010/proj_STseqData/Donor_3_4/P20063/00-slidefile/V10B01-025.gpr \
--loupe-alignment=/proj/sens2020010/proj_STseqData/manual_alignment_files/V10B01-025-B1.json \
--localcores=12 \
--localmem=100
