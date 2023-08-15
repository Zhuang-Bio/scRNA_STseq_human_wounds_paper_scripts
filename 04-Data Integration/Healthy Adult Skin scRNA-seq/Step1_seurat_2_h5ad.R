#########################
# Use SeuratDisk package 
# to convert the seurat objects to a h5ad object to load into scanpy/squidpy.
library(Seurat)
library(SeuratDisk)
library(hdf5r)
library(tidyverse)

# output dir
outdir = "/Users/zhuliu/Desktop/sc st/compareStudies/NingGroup_woundsData"


###########################
##---- Main clusters ----##
# load the scRNAseq data (Main clusters)
sc.data = readRDS(file = "/Users/zhuliu/Desktop/sc st/plotting scripts/Figure 1/allcombined_wounds_newAnnotation.rds")
mt <- sc.data@meta.data %>% rownames_to_column(var = "barcodeID")
mt <- mt %>% mutate(barcode = barcodeID) %>% select(-c(3,4,9,18:22))

mt_f <- mt %>% select(1, 20, 2:10, 19, 11:14, 18, everything()) %>% column_to_rownames(var = "barcodeID") %>% 
  rename("SampleID" = "orig.ident",
         "ClusterID" = "newAnnoNum",
         "CellType_fullname" = "CellTypes",
         "CellType" = "newCellTypes",
         "BroadCellType" = "newMainCellTypes")
sc.data@meta.data <- mt_f

DimPlot(sc.data, group.by = "CellType")

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

# covert to h5ad file
SaveH5Seurat(sc.data, filename = file.path(outdir,"HumanWounds_scRNA_MainClusters.h5Seurat"), overwrite = T)
Convert(file.path(outdir,"HumanWounds_scRNA_MainClusters.h5Seurat"), dest = "h5ad", overwrite = T)

# remove unneeded file
rm(sc.data)
file.remove(file.path(outdir,"HumanWounds_scRNA_MainClusters.h5Seurat"))
gc()



####################################
##---- myeloid/lymphoid cells ----##
# load the scRNAseq data
sc.data = lym_sub
mt <- sc.data@meta.data %>% rownames_to_column(var = "barcodeID")
mt <- mt %>% mutate(barcode = barcodeID) %>% select(-c(3,4,9,19:22, 24:25))

colnames(mt)
mt_f <- mt %>% select(1, 18, everything()) %>% column_to_rownames(var = "barcodeID") %>% 
  rename("SampleID" = "orig.ident",
         "ClusterID" = "SCT_snn_res.0.5",
         "CellType_sameMain" = "MainCellTypes",
         "CellType" = "upCellTypes") %>% select(1:14, 16, everything())
sc.data@meta.data <- mt_f

DimPlot(sc.data, group.by = "CellType")

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

# covert to h5ad file
SaveH5Seurat(sc.data, filename = file.path(outdir,"HumanWounds_scRNA_LymphoidCells.h5Seurat"), overwrite = T)
Convert(file.path(outdir,"HumanWounds_scRNA_LymphoidCells.h5Seurat"), dest = "h5ad", overwrite = T)

# remove unneeded file
rm(sc.data)
file.remove(file.path(outdir,"HumanWounds_scRNA_LymphoidCells.h5Seurat"))
gc()
