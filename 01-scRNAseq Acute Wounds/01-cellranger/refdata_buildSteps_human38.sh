#! /bin/bash -l
#SBATCH -A sens2020010
#SBATCH -p node -n 16
#SBATCH -t 5-00:00:00
#SBATCH -J mus39ref
#SBATCH -e musGRCm39.SLURM_Job_id=%j.sderr
#SBATCH -o musGRCm39.SLURM_Job_id=%j.sdout
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zhuang.liu@ki.se

module load bioinfo-tools cellranger/5.0.1

# Modify sequence headers in the Ensembl FASTA to match the file
# "GRCh38.primary_assembly.genome.fa" from GENCODE. Unplaced and unlocalized
# sequences such as "KI270728.1" have the same names in both versions.
#
# Input FASTA:
#   >1 dna:chromosome chromosome:GRCh38:1:1:248956422:1 REF
#
# Output FASTA:
#   >chr1 1
# sed commands:
# 1. Replace metadata after space with original contig name, as in GENCODE
# 2. Add "chr" to names of autosomes and sex chromosomes
# 3. Handle the mitochrondrial chromosome

cat Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    | sed -E 's/^>(\S+).*/>\1 \1/' \
    | sed -E 's/^>([0-9]+|[XY]) />chr\1 /' \
    | sed -E 's/^>MT />chrM /' \
    > Homo_sapiens.GRCh38.dna.primary_assembly_modified.fa

#Check the head file
grep -E ">" Homo_sapiens.GRCh38.dna.primary_assembly_modified.fa | less


# Remove version suffix from transcript, gene, and exon IDs in order to match
# previous Cell Ranger reference packages
#
# Input GTF:
#     ... gene_id "ENSG00000223972.5"; ...
# Output GTF:
#     ... gene_id "ENSG00000223972"; gene_version "5"; ...
# Pattern matches Ensembl gene, transcript, and exon IDs for human or mouse:
ID="(ENS(MUS)?[GTE][0-9]+)\.([0-9]+)"
cat gencode.v38.primary_assembly.annotation.gtf \
    | sed -E 's/gene_id "'"$ID"'";/gene_id "\1"; gene_version "\3";/' \
    | sed -E 's/transcript_id "'"$ID"'";/transcript_id "\1"; transcript_version "\3";/' \
    | sed -E 's/exon_id "'"$ID"'";/exon_id "\1"; exon_version "\3";/' \
    > gencode.v38.primary_assembly.annotation_modified.gtf


# Define string patterns for GTF tags
# NOTES:
# - Since GENCODE release 31/M22 (Ensembl 97), the "lincRNA" and "antisense"
#   biotypes are part of a more generic "lncRNA" biotype.
# - These filters are relevant only to GTF files from GENCODE. The GTFs from
#   Ensembl release 98 have the following differences:
#   - The names "gene_biotype" and "transcript_biotype" are used instead of
#     "gene_type" and "transcript_type".
#   - Readthrough transcripts are present but are not marked with the
#     "readthrough_transcript" tag.
#   - Only the X chromosome versions of genes in the pseudoautosomal regions
#     are present, so there is no "PAR" tag.
BIOTYPE_PATTERN=\
"(protein_coding|lncRNA|\
IG_C_gene|IG_D_gene|IG_J_gene|IG_LV_gene|IG_V_gene|\
IG_V_pseudogene|IG_J_pseudogene|IG_C_pseudogene|\
TR_C_gene|TR_D_gene|TR_J_gene|TR_V_gene|\
TR_V_pseudogene|TR_J_pseudogene)"
GENE_PATTERN="gene_type \"${BIOTYPE_PATTERN}\""
TX_PATTERN="transcript_type \"${BIOTYPE_PATTERN}\""
READTHROUGH_PATTERN="tag \"readthrough_transcript\""
PAR_PATTERN="tag \"PAR\""


# Construct the gene ID allowlist. We filter the list of all transcripts
# based on these criteria:
#   - allowable gene_type (biotype)
#   - allowable transcript_type (biotype)
#   - no "PAR" tag (only present for Y chromosome PAR)
#   - no "readthrough_transcript" tag
# We then collect the list of gene IDs that have at least one associated
# transcript passing the filters.
cat gencode.v38.primary_assembly.annotation_modified.gtf \
    | awk '$3 == "transcript"' \
    | grep -E "$GENE_PATTERN" \
    | grep -E "$TX_PATTERN" \
    | grep -Ev "$READTHROUGH_PATTERN" \
    | grep -Ev "$PAR_PATTERN" \
    | sed -E 's/.*(gene_id "[^"]+").*/\1/' \
    | sort \
    | uniq \
    > gencode.v38.primary_assembly.annotation_modified_geneallows.gtf


# Filter the GTF file based on the gene allowlist
# Copy header lines beginning with "#"
grep -E "^#" gencode.v38.primary_assembly.annotation_modified.gtf > gencode.v38.primary_assembly.annotation_modified_final.gtf
# Filter to the gene allowlist
grep -Ff gencode.v38.primary_assembly.annotation_modified_geneallows.gtf gencode.v38.primary_assembly.annotation_modified.gtf \
    >> gencode.v38.primary_assembly.annotation_modified_final.gtf

#move prepared filed into a folder
mkdir prepared
mv gencode.v38.primary_assembly.annotation.gtf \
    gencode.v38.primary_assembly.annotation_modified_geneallows.gtf \
    gencode.v38.primary_assembly.annotation_modified.gtf \
    Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    prepared

# Create reference package
cd /proj/sens2020010/ref10X_GRCh38_2021Jul
cellranger mkref \
    --ref-version=2021Jul \
    --genome=GRCh38 \
    --fasta=Homo_sapiens.GRCh38.dna.primary_assembly_modified.fa \
    --genes=gencode.v38.primary_assembly.annotation_modified_final.gtf
