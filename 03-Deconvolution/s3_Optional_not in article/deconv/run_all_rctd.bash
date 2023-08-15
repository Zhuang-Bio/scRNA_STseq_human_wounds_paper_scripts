



#SAMPLES="Donor1_Wound1 Donor4_Skin Donor1_Wound7 Donor2_Wound30" 

SAMPLES="Donor1_Skin Donor2_Skin Donor3_Skin Donor2_Wound1 Donor3_Wound1 Donor4_Wound1 Donor1_Wound30 Donor3_Wound30 Donor4_Wound30 Donor2_Wound7 Donor3_Wound7  Donor4_Wound7"


conda activate R_rctd

#Rscript rctd_create_sc_ref.R  -i ../../deconv/inputs/sc_data/s1_subsampled -o ../../deconv/results/rctd/sc_ref 


for sample in $SAMPLES
do
    echo $sample
    Rscript rctd_run_prediction.R -i ../../deconv/inputs/st_data/ -o ../../deconv/results/rctd/ -r ../../deconv/results/rctd/sc_ref -s $sample
done
