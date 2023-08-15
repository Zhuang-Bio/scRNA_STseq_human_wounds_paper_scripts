# Use SeuratDisk package to convert the seurat objects to a h5ad object to load into scanpy/squidpy.

library(Seurat)
library(SeuratDisk)
library(hdf5r)



# first the scRNAseq data.
outdir = "/Users/zhuliu/Desktop/scRNA_STseq/proj_10X_woundhealing/03_results/01-Seurat-PreAnalysis/02_Seurat_BatchCorrection/s1_Seurat_allSample_harmony"
savefile = file.path(outdir,"s1_cleanSeurat_harmony_allSamples_clusters.rds")
sc.data = readRDS(file = savefile)


# add new celltype annotation
anno <- readxl::read_xlsx(file.path(outdir,"s5_cell_types_for_each_cluster_manually_updated.xlsx"), sheet = 1)
ct = plyr::mapvalues(x=Idents(sc.data), from=c(anno$Cluster), to=c(anno$CellType))
sc.data$Celltype = as.character(ct)

subtype = plyr::mapvalues(x=Idents(sc.data), from=c(anno$Cluster), to=c(anno$Detailed_info))
sc.data$Subtype = as.character(subtype)

type = plyr::mapvalues(x=Idents(sc.data), from=c(anno$Cluster), to=c(anno$MainCellTypes))
sc.data$MainCelltype = as.character(type)


# OBS! All factor columns are converted to numbers in scanpy, fix!
cl = sapply(sc.data@meta.data, is.factor)

for (col in names(cl)[cl]){
  print(col)
  sc.data[[col]] = as.character(sc.data[[col]][,1])
}


# to keep the reductions even if removing SCT assay, change the assay.used
sc.data@reductions$harmony@assay.used = "RNA"
sc.data@reductions$umap@assay.used = "RNA"
sc.data@reductions$pca@assay.used = "RNA"

sc.data@active.assay = "RNA"
sc.data = DietSeurat(sc.data, assays = "RNA", dimreducs = c("pca", "umap", "harmony"))

sc.data@assays$RNA@counts[1:50,1:10]
dim(sc.data@assays$RNA@counts)
#25778 12259

#cover the data slot with raw counts, otherwise the normalized data will be the scanpy.anndata.X
#or you can change it in scanpy pipeline: anndata=anndata.raw.to_adata(), get the raw counts as anndata.X
sc.data@assays$RNA@data = sc.data@assays$RNA@counts

SaveH5Seurat(sc.data, filename = file.path(outdir,"s1_cleanSeurat_harmony_allSamples_clusters.h5Seurat"), overwrite = T)
Convert(file.path(outdir,"s1_cleanSeurat_harmony_allSamples_clusters.h5Seurat"), dest = "h5ad", overwrite = T)

rm(sc.data)
file.remove(file.path(outdir,"s1_cleanSeurat_harmony_allSamples_clusters.h5Seurat"))
gc()


do.st = T

if(do.st){
  # then the ST data
  outdir = "/Users/zhuliu/Desktop/scRNA_STseq/proj_spatialTranseq/03_results/01-SpatialAnalysis"
  savefile = file.path(outdir,"Seurat_STseq_integrated.rds")
  st.data = readRDS(file = savefile)
  # OBS! All factor columns are converted to numbers in scanpy, fix!
  cl = sapply(st.data@meta.data, is.factor)
  
  for (col in names(cl)[cl]){
    print(col)
    st.data[[col]] = as.character(st.data[[col]][,1])
  }
  # to keep the reductions even if removing SCT assay, change the assay.used
  st.data@reductions$harmony@assay.used = "Spatial"
  st.data@reductions$umap@assay.used = "Spatial"
  st.data@reductions$pca@assay.used = "Spatial"
  
  st.data@active.assay = "Spatial"
  st.data@assays$Spatial@data = st.data@assays$Spatial@counts
  st.data = DietSeurat(st.data, assays = "Spatial", dimreducs = c("pca", "umap", "harmony"))
  
  SaveH5Seurat(st.data, filename = file.path(outdir,"Seurat_STseq_integrated.h5Seurat"))
  Convert(file.path(outdir,"Seurat_STseq_integrated.h5Seurat"), dest = "h5ad")
  
  file.remove(file.path(outdir,"Seurat_STseq_integrated.h5Seurat"))
}


wound <- readRDS("/Users/zhuliu/Desktop/scRNA_STseq/proj_10X_woundhealing/03_results/01-Seurat-PreAnalysis/02_Seurat_BatchCorrection/s3_Seurat_subFibroblasts/s1_clean_fibroblast.rds")
wound@assays$RNA@counts[1:50,1:10]
dim(wound@assays$RNA@counts)
#25778 12259

