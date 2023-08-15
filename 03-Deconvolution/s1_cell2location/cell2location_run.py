#!/usr/bin/python

import sys
import scanpy as sc
import anndata
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

import cell2location
import scvi

from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42 # enables correct plotting of text
import seaborn as sns

# silence scanpy that prints a lot of warnings
import warnings

from cell2location.models import RegressionModel

results_folder = '/proj/snic2021-23-156/cell2location_0704'

adata_vis=sc.read("/proj/snic2021-23-156/cell2location_0704/ST-seq_scanpy0519.h5ad")
adata_sc=sc.read("/proj/snic2021-23-156/cell2location_0704/scRNAseq_allcleanCells_Top200DEgenes.h5ad")

# prepare anndata for the regression model
RegressionModel.setup_anndata(adata=adata_sc,
                              # 10X reaction / sample / batch
                              batch_key='orig.ident',
                              # cell type, covariate used for constructing signatures
                              labels_key='newCellTypes',
                              # multiplicative technical effects (platform, 3' vs 5', donor effect)
                              categorical_covariate_keys=None,
                              continuous_covariate_keys=None
                              )

mod = RegressionModel(adata_sc)
mod.view_anndata_setup(adata_sc)

# Use all data for training (validation not implemented yet, train_size=1)
mod.train(max_epochs=250, # If None, defaults to np.min([round((20000 / n_cells) * 400), 400])
          batch_size=2500, 
          train_size=1, # proportion of cells in the training set (for cross-validation)
          lr=0.002) #max_epochs=250

mod.history["elbo_train"].iloc[:].plot()
plt.show()

# In this section, we export the estimated cell abundance (summary of the posterior distribution).
adata_sc = mod.export_posterior(
    adata_sc, sample_kwargs={'num_samples': 1000, 'batch_size': 2500}
)

ref_run_name2 = f'/proj/snic2021-23-156/cell2location_0704/regressionMod'
# Save model
mod.save(f"{ref_run_name2}", overwrite=True)

# Save anndata object with results
adata_file = f"{ref_run_name2}/sc.h5ad"
adata_sc.write(adata_file)

mod.plot_QC()

# export estimated expression in each cluster
if 'means_per_cluster_mu_fg' in adata_sc.varm.keys():
    inf_aver = adata_sc.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_sc.uns['mod']['factor_names']]].copy()
else:
    inf_aver = adata_sc.var[[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_sc.uns['mod']['factor_names']]].copy()
inf_aver.columns = adata_sc.uns['mod']['factor_names']
print(len(inf_aver))
inf_aver.iloc[0:5, :]

pd.DataFrame(inf_aver).to_csv(f"{ref_run_name2}/infer_signatures_sc_filter_0704.csv")

# find shared genes and subset both anndata and reference signatures
intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
print(len(intersect))
adata_vis = adata_vis[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()

# prepare anndata for cell2location model
cell2location.models.Cell2location.setup_anndata(adata=adata_vis, batch_key="library_id")

# create and train the model
mod = cell2location.models.Cell2location(
    adata_vis, 
    cell_state_df=inf_aver,
    # the expected average cell abundance: tissue-dependent
    # hyper-prior which can be estimated from paired histology:
    N_cells_per_location=5,
    # hyperparameter controlling normalisation of
    # within-experiment variation in RNA detection (using default here):
    detection_alpha=20
)

mod.view_anndata_setup(adata_vis)

mod.train(max_epochs=30000,
          # train using full data (batch_size=None)
          batch_size=None,
          # use all data points in training because
          # we need to estimate cell abundance at all locations
          train_size=1,
          use_gpu=False)

# In this section, we export the estimated cell abundance (summary of the posterior distribution).
adata_vis = mod.export_posterior(
    adata_vis, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'use_gpu': False}
)

# Save model
mod.save(f"{results_folder}", overwrite=True)

# Save anndata object with results
adata_file = f"{results_folder}/st_deconv.h5ad"
adata_vis.write(adata_file)


