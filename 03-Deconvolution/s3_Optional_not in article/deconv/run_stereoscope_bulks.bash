# Run Stereoscope with bulk samples

# do before running bash script instead?
conda activate scanpy_scvi

SC_LIST="../../deconv/inputs/bulks/bulk_intersection_genes.txt"
SC_DATA="../../deconv/inputs/sc_data/s1_subsampled_no_basal4.h5ad"
OUTDIR="../../deconv/results/bulks"
mkdir -p $OUTDIR


#python stereoscope_create_sc_ref.py -i $SC_DATA -o $OUTDIR/sc_ref -a cl.annot -g $SC_LIST

python stereoscope_run_bulks.py -i ../../deconv/inputs/bulks/bulk_data.h5ad -o $OUTDIR/ -r $OUTDIR/sc_ref 

