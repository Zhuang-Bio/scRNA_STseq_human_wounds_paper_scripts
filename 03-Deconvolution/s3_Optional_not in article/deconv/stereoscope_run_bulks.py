import scvi
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import scanpy as sc
import os, sys


from scvi.data import register_tensor_from_anndata
from scvi.external import RNAStereoscope, SpatialStereoscope

import warnings
warnings.filterwarnings('ignore')

import argparse
parser = argparse.ArgumentParser()


parser.add_argument("-i","--infile", help="input directory", type=str)
parser.add_argument("-o","--outdir", help="output directory", type=str)
parser.add_argument("-r","--sc_ref", help="saved single cell model directory", type=str)


args = parser.parse_args()
print(args)

#python stereoscope_run_bulks.py -i ../../deconv/inputs/bulks/bulk_data.h5ad -o ../../deconv/results/bulks/ -r ../../deconv/results/bulks/sc_ref 

print("Running stereoscope for file "+args.infile)

#outdir = "../../deconv/results/bulks/"
#infile = "../../deconv/inputs/bulks/bulk_data.h5ad"
#ref = "../../deconv/results/bulks/sc_ref"
#adata_path =  os.path.join(ref,"adata.h5ad")
#model_dir = os.path.join(ref,"model")
#adata_st = sc.read_h5ad(infile)


# Setup folders
outdir = args.outdir
if not os.path.exists(outdir):
    os.makedirs(outdir)

# read count matrix
adata_st = sc.read_h5ad(args.infile)
adata_st = adata_st.raw.to_adata()


# load sc model
adata_path =  os.path.join(args.sc_ref,"adata.h5ad")
sc_adata = sc.read_h5ad(adata_path)
print("Loaded SC data from file :" + adata_path)
model_dir = os.path.join(args.sc_ref,"model")
sc_model = RNAStereoscope.load(model_dir, sc_adata)
print("Loaded RNA model from file :" + model_dir)



# subst for same genes as sc data
genes = sc_model.adata.var_names
genes = genes.intersection(adata_st.var.index)
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


#ValueError: Value is not broadcastable with batch_shape+event_shape: torch.Size([68, 1233]) vs torch.Size([68, 1406]).
#ValueError: Expected value argument (Tensor of shape (68, 4555)) to be within the support (IntegerGreaterThan(lower_bound=0)) of the distribution NegativeBinomial(total_count: torch.Size([68, 4555]), logits: torch.Size([68, 4555])), but found invalid values:
