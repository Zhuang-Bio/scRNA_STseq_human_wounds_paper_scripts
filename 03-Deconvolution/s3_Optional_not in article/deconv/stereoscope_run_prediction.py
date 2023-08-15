import scvi
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import os, sys

from scvi.data import register_tensor_from_anndata
from scvi.external import RNAStereoscope, SpatialStereoscope

import warnings
warnings.filterwarnings('ignore')

import argparse
parser = argparse.ArgumentParser()


parser.add_argument("-i","--indir", help="input directory", type=str)
parser.add_argument("-o","--outdir", help="output directory", type=str)
parser.add_argument("-s","--sample", help="ST sample name", type=str)
parser.add_argument("-r","--sc_ref", help="saved single cell model directory", type=str)


args = parser.parse_args()
print(args)

#python stereoscope_run_prediction.py -i ../../deconv/inputs/st_data/ -o ../../deconv/results/stereoscope/ -r ../../deconv/results/stereoscope/sc_ref -s Donor3_Wound7

sample = args.sample
print("Running stereoscope for sample "+sample)


# Setup folders
if not os.path.exists(args.outdir):
    os.makedirs(args.outdir)

    
outdir = os.path.join(args.outdir,sample)
if not os.path.exists(outdir):
    os.makedirs(outdir)

indir = os.path.join(args.indir, sample)
if not os.path.exists(indir):
    sys.exit("ERROR: No input folder: " + indir)

raw_file = os.path.join(indir, "rawdata_path.csv")
if not os.path.exists(raw_file):
     sys.exit("ERROR: No rawdata path file: " + raw_file)
     
fh = open(os.path.join(indir, "rawdata_path.csv"),"r+")
raw_path = fh.read()
fh.close()
raw_path = raw_path.strip()

if not os.path.exists(raw_path):
     sys.exit("ERROR: No rawdata file: " + raw_path)


barcodes_file =  os.path.join(indir, "barcodes.csv")
if not os.path.exists(barcodes_file):
     sys.exit("ERROR: No barcodes file: " + barcodes_file)


    

# Load ST data    
print("Reading ST data from :" + raw_path)   
adata_st = sc.read_visium(raw_path, library_id = sample)
adata_st.var_names_make_unique()

# read barcodes and filter spots
barcodes = pd.read_csv(barcodes_file, header = None)
adata_st = adata_st[adata_st.obs.index.isin(barcodes[0]),:]


# load sc model
adata_path =  os.path.join(args.sc_ref,"adata.h5ad")
sc_adata = sc.read_h5ad(adata_path)
print("Loaded SC data from file :" + adata_path)
model_dir = os.path.join(args.sc_ref,"model")
sc_model = RNAStereoscope.load(model_dir, sc_adata)
print("Loaded RNA model from file :" + model_dir)



# subst for same genes as sc data
genes = sc_model.adata.var_names
adata_st.layers["counts"] = adata_st.X.copy()
adata_st = adata_st[:, genes].copy()

# create object and train 
SpatialStereoscope.setup_anndata(adata_st, layer="counts")
spatial_model = SpatialStereoscope.from_rna_model(adata_st, sc_model)
spatial_model.train(max_epochs = 10000)
#spatial_model.history["elbo_train"][10:].plot()
spatial_model.save(os.path.join(outdir,"stmodel"), overwrite = True)

# write results file
adata_st.obsm["deconvolution"] = spatial_model.get_proportions()
adata_st.obsm["deconvolution"].to_csv(outdir + "/proportions.csv")
