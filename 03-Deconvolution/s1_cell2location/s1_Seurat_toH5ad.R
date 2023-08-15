suppressPackageStartupMessages(require(Matrix))
suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(Seurat))
suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(SeuratDisk))
suppressPackageStartupMessages(require(hdf5r))
library(tidyverse)


## Prepare and read the scRNA-seq data
# load the clean rds object from the step 1 (Figure 1)
sc.data <- readRDS("/Users/zhuliu/Desktop/sc st/plotting scripts/Figure 1/allcombined_wounds_newAnnotation.rds")
DimPlot(sc.data, group.by = "newCellTypes", label = T) + NoAxes()
DimPlot(sc.data, group.by = "newMainCellTypes", label = T) + NoAxes()


# OBS! All factor columns are converted to numbers in scanpy, fix!
cl = sapply(sc.data@meta.data, is.factor)

for (col in names(cl)[cl]){
  print(col)
  sc.data[[col]] = as.character(sc.data[[col]][,1])
}

sc.data$newMainCellTypes <- gsub(" ", "_",sc.data$newMainCellTypes)
sc.data$newMainCellTypes <- gsub("\\&", "_",sc.data$newMainCellTypes)
sc.data$newMainCellTypes <- gsub("__", "_",sc.data$newMainCellTypes)
sc.data$newMainCellTypes <- gsub("__", "_",sc.data$newMainCellTypes)
table(sc.data$newMainCellTypes)


# load the sub cell type annotation
mt1 <- data.table::fread("/Users/zhuliu/Desktop/sc st/plotting scripts/newsubCluster_annotation/MEL_Schwann.txt") %>% setNames(c("barcode", "upCellTypes"))
mt2 <- data.table::fread("/Users/zhuliu/Desktop/sc st/plotting scripts/newsubCluster_annotation/newFibroblast_metadata.txt")
mt3 <- data.table::fread("/Users/zhuliu/Desktop/sc st/plotting scripts/newsubCluster_annotation/newKeratinocyte_metadata.txt")
mt4 <- data.table::fread("/Users/zhuliu/Desktop/sc st/plotting scripts/newsubCluster_annotation/newMyeloid_metadata.txt")
mt5 <- data.table::fread("/Users/zhuliu/Desktop/sc st/plotting scripts/newsubCluster_annotation/newOther_metadata.txt")
mt6 <- data.table::fread("/Users/zhuliu/Desktop/sc st/plotting scripts/newsubCluster_annotation/newLymphoid_metadata.txt")
mt_all <- rbind(mt1, mt2, mt3, mt4, mt5, mt6)
table(mt_all$upCellTypes)

mt_sc <- sc.data@meta.data %>% tibble::rownames_to_column(var = "barcode") %>% left_join(., mt_all) %>% select(-18:-23, -26)
mt_sc$upCellTypes <- ifelse(is.na(mt_sc$upCellTypes), mt_sc$newCellTypes, mt_sc$upCellTypes)


mt_sc$upCellTypes <- gsub("-", "_", mt_sc$upCellTypes)
mt_sc$upCellTypes <- gsub("\\(|\\)", "_", mt_sc$upCellTypes)
mt_sc$upCellTypes <- gsub("\\+", "_", mt_sc$upCellTypes)
mt_sc$upCellTypes <- gsub("_$", "", mt_sc$upCellTypes)
mt_sc$upCellTypes <- gsub("_$", "", mt_sc$upCellTypes)
mt_sc$upCellTypes <- gsub("/", "_", mt_sc$upCellTypes)
mt_sc$upCellTypes <- gsub("Spi_II", "Spi_II_a", mt_sc$upCellTypes)
mt_sc$upCellTypes <- gsub("Spi_II_a_a", "Spi_II_a", mt_sc$upCellTypes)
mt_sc$upCellTypes <- gsub("Spi_II_a_b|Spi_II_aI", "Spi_II_b", mt_sc$upCellTypes)
mt_sc$upCellTypes <- gsub("M1", "Mac1", mt_sc$upCellTypes)
mt_sc$upCellTypes <- gsub("M2", "Mac2", mt_sc$upCellTypes)
mt_sc$upCellTypes <- gsub("FB_II_APOD_ITMac2A", "FB_II_APOD_ITM2A", mt_sc$upCellTypes)
mt_sc$upCellTypes <- gsub("Mac_mig", "Mac3", mt_sc$upCellTypes)
table(mt_sc$upCellTypes)
identical(colnames(sc.data), mt_sc$barcode)

mt_sc_f <- mt_sc %>% column_to_rownames(var = "barcode") %>% select(1:8, 19, everything())

sc.data@meta.data <- mt_sc_f
# to keep the reductions even if removing SCT assay, change the assay.used
sc.data@reductions$harmony@assay.used = "RNA"
sc.data@reductions$umap@assay.used = "RNA"
sc.data@reductions$pca@assay.used = "RNA"

sc.data <- DietSeurat(sc.data, assays = "RNA", dimreducs = c("pca", "umap", "harmony"))

DimPlot(sc.data, group.by = "upCellTypes", label = T) + NoAxes()


outdir = "/Users/zhuliu/Desktop/scRNA_STseq/proj_spatialTranseq/03_results/deconvolution/cell2location_1118_allsubCT"
# Save as AnnData
SaveH5Seurat(sc.data, filename = file.path(outdir, "scRNAseq_allcleanCells_tissueMAP.h5seurat"),
             overwrite = T)
Convert(file.path(outdir, "scRNAseq_allcleanCells_tissueMAP.h5seurat"), dest = "h5ad", overwrite = T)

file.remove(file.path(outdir, "scRNAseq_allcleanCells_tissueMAP.h5seurat"))

