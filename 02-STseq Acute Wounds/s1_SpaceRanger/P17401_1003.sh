#! /bin/bash -l
#SBATCH -A sens2020010
#SBATCH -p core -n 16
#SBATCH -t 5-00:00:00
#SBATCH -J P17401_1003
#SBATCH -e P17401_1003.SLURM_Job_id=%j.sderr
#SBATCH -o P17401_1003.SLURM_Job_id=%j.sdout
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zhuang.liu@ki.se

module load bioinfo-tools spaceranger/1.2.0

cd /proj/sens2020010/proj_STseqData/Donor_1/01-Results/spaceranger-latestRef-manualAlignment

SAMPLE_IDS="P17401_1003"
TRANSCRIPTOME="/proj/sens2020010/ref10X_GRCh38_2021Jul/GRCh38"
FASTQS="/proj/sens2020010/proj_STseqData/Donor_1/P17401_1003/02-FASTQ/201202_A00187_0394_BHWF2YDRXX/"
IMAGES="/proj/sens2020010/proj_STseqData/Donor_1/01-Results/00-images/201006_P17401_V10J29-111_AM.C1-Spots/201006_P17401_V10J29-111_AM.C1-Spot000001.jpg"

spaceranger count \
--id=$SAMPLE_IDS \
--transcriptome=$TRANSCRIPTOME \
--fastqs=$FASTQS \
--sample=$SAMPLE_IDS \
--image=$IMAGES \
--slide=V10J29-111 \
--area=C1 \
--slidefile=/proj/sens2020010/proj_STseqData/Donor_1/01-Results/00-slidefile/V10J29-111.gpr \
--loupe-alignment=/proj/sens2020010/proj_STseqData/manual_alignment_files/V10J29-111-C1.json \
--localcores=12 \
--localmem=100
