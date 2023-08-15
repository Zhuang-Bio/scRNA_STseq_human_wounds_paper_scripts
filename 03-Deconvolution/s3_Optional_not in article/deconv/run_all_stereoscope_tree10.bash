# Run Stereoscope with all samples.


SAMPLES="Donor1_Wound1 Donor4_Skin Donor1_Wound7 Donor2_Wound30 Donor1_Skin Donor2_Skin Donor3_Skin Donor2_Wound1 Donor3_Wound1 Donor4_Wound1 Donor1_Wound30 Donor3_Wound30 Donor4_Wound30 Donor2_Wound7 Donor3_Wound7  Donor4_Wound7"
#SAMPLES="Donor1_Wound30 Donor3_Wound30 Donor4_Wound30 Donor2_Wound7 Donor3_Wound7  Donor4_Wound7"

# do before running bash script instead?
conda activate scanpy_scvi

SC_LIST="../../deconv/inputs/sc_data/degs_tree_top10.txt"
SC_DATA="../../deconv/inputs/sc_data/s1_subsampled_no_basal4.h5ad"
OUTDIR="../../deconv/results/stereoscope10"
mkdir -p $OUTDIR

python stereoscope_create_sc_ref.py -i $SC_DATA -o $OUTDIR/sc_ref -a cl.annot -g $SC_LIST

for sample in $SAMPLES
do
    echo $sample
    python stereoscope_run_prediction.py -i ../../deconv/inputs/st_data/ -o $OUTDIR/ -r $OUTDIR/sc_ref -s $sample
done

