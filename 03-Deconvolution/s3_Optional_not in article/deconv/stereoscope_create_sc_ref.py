import scvi
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import os

from scvi.data import register_tensor_from_anndata
from scvi.external import RNAStereoscope, SpatialStereoscope

import warnings
warnings.filterwarnings('ignore')

import argparse
parser = argparse.ArgumentParser()


parser.add_argument("-i","--input", help="single cell anndata object", type=str)
parser.add_argument("-o","--outdir", help="output directory", type=str)
parser.add_argument("-a","--annot_column", help="column in metadata for celltype annotation", type=str)
parser.add_argument("-g","--gene_list", help="csv file with genes to use for training", type=str)

args = parser.parse_args()
print(args)

#python stereoscope_create_sc_ref.py -i ../../deconv/inputs/sc_data/s1_subsampled.h5ad -o ../../deconv/results/stereoscope/sc_ref -a cl.annot -g ../../deconv/inputs/sc_data/degs_fc0.5_pval0.01.txt

if not os.path.exists(args.outdir):
    os.makedirs(args.outdir)

# read in the sc data
adata_sc = sc.read_h5ad(args.input)

# OBS! The raw matrix has numericals instead of gene names as var.index, also make sure column name is correct!
adata_sc.raw.var.rename(columns = {'_index':'features'}, inplace = True)
adata_sc.raw.var.index = adata_sc.var.index
adata_sc = adata_sc.raw.to_adata()

# read genes
genes = pd.read_csv(args.gene_list, header=None)
genes = list(genes[0])


# add counts layer
adata_sc.layers["counts"] = adata_sc.X.copy()

# subset for the genes
sc_adata = adata_sc[:, genes].copy()

# create stereoscope object
RNAStereoscope.setup_anndata(sc_adata, layer = "counts", labels_key = args.annot_column)


sc_model = RNAStereoscope(sc_adata)
sc_model.train(max_epochs = 100)
sc_model.history["elbo_train"][10:].plot()

model_dir = os.path.join(args.outdir,"model")
sc_model.save(model_dir, overwrite=True)
adata_path =  os.path.join(args.outdir,"adata.h5ad")
sc_model.adata.write_h5ad(adata_path)


