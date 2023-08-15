#!/usr/bin/python

import scvi
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
import os

#from scvi.data import register_tensor_from_anndata
from scvi.external import RNAStereoscope, SpatialStereoscope

results_folder = '/proj/snic2021-23-156/stereoscope_deconv/'
# reload sc model
ref_run_name = f"{results_folder}/regression_model"

adata_path = os.path.join(f"{ref_run_name}","sc.h5ad")
adata_sc = sc.read_h5ad(adata_path)
print("Loaded SC data from file :" + adata_path)
sc_model = RNAStereoscope.load(f"{ref_run_name}", adata_sc)
print("Loaded RNA model from file :" + ref_run_name)

# ST data (raw counts)
adata_vis = sc.read_h5ad(results_folder + "Seurat_STseq_integrated.h5ad")

# OBS! The raw matrix has numericals instead of gene names as var.index, also make sure column name is correct!
adata_vis.var = adata_vis.var.drop(columns = ['_index'])
adata_vis.var.index = adata_vis.var['features']

# subst for same genes as sc data
genes = sc_model.adata.var_names
adata_vis.layers["counts"] = adata_vis.X.copy()
adata_vis = adata_vis[:, genes].copy()

# create object and train 
SpatialStereoscope.setup_anndata(adata_vis, layer="counts")
spatial_model = SpatialStereoscope.from_rna_model(adata_vis, sc_model)
spatial_model.train(max_epochs = 10000)
spatial_model.history["elbo_train"][10:].plot()

spatial_model.save(f"{results_folder}/stereo_map", overwrite = True)

# write results file
adata_vis.obsm["deconvolution"] = spatial_model.get_proportions()
adata_vis.obsm["deconvolution"].to_csv(results_folder + "/stereo_celltype_proportions.csv")

