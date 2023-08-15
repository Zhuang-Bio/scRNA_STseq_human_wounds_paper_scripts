



#SAMPLES="Donor1_Wound1 Donor4_Skin Donor1_Wound7 Donor2_Wound30" 

SAMPLES="Donor1_Skin Donor2_Skin Donor3_Skin Donor2_Wound1 Donor3_Wound1 Donor4_Wound1 Donor1_Wound30 Donor3_Wound30 Donor4_Wound30 Donor2_Wound7 Donor3_Wound7  Donor4_Wound7"


conda activate scanpy_scvi

#python stereoscope_create_sc_ref.py -i ../../deconv/inputs/sc_data/s1_subsampled.h5ad -o ../../deconv/results/stereoscope/sc_ref -a cl.annot -g ../../deconv/inputs/sc_data/degs_fc0.5_pval0.01.txt

for sample in $SAMPLES
do
    echo $sample
    python stereoscope_run_prediction.py -i ../../deconv/inputs/st_data/ -o ../../deconv/results/stereoscope/ -r ../../deconv/results/stereoscope/sc_ref -s $sample
done
