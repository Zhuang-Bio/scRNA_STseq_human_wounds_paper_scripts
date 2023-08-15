setwd("/Users/zhuliu/Desktop/scRNA_STseq/proj_10X_woundhealing/03_results/01-Seurat-PreAnalysis/02_Seurat_BatchCorrection/s1_Seurat_allSample_harmony/s8_Cell_Cell_L_R_analysis_upCellTypes/CellPhoneDB")

library(Seurat)
library(SeuratObject)
library(tidyverse)
library(Matrix)

# load the seurat object
hswounds <- readRDS("/Users/zhuliu/Desktop/sc st/plotting scripts/Figure 1/allcombined_wounds_newAnnotation.rds")
hswounds$newCellTypes <- as.character(hswounds$newCellTypes)
str(hswounds$newCellTypes)

# load the clean metadata
mt <- data.table::fread("/Users/zhuliu/Desktop/sc st/plotting scripts/newsubCluster_annotation/Subclusters.txt")
hswounds <- hswounds[, mt$barcode]

identical(colnames(hswounds), mt$barcode)
hswounds$CellType <- mt$CellType

hswounds <- DietSeurat(hswounds, counts = TRUE, assays = "RNA")
hswounds <- NormalizeData(hswounds)

hswounds_sub <- subset(hswounds, subset = Condition == "Skin")

# Write gene expression in mtx format
# Save normalised counts - NOT scaled!
writeMM(hswounds_sub@assays$RNA@data, file = 'countsData_Skin/matrix.mtx')
# save gene and cell names
write(x = rownames(hswounds_sub@assays$RNA@data), file = "countsData_Skin/features.tsv")
write(x = colnames(hswounds_sub@assays$RNA@data), file = "countsData_Skin/barcodes.tsv")

# Generate your meta data
hswounds_sub@meta.data$Cell = rownames(hswounds_sub@meta.data)
df = hswounds_sub@meta.data[, c('Cell', 'CellType')]
write.table(df, file ='hswounds_metadata_Skin.tsv', sep = '\t', quote = F, row.names = F)


hswounds_sub <- subset(hswounds, subset = Condition == "Wound1")

# Write gene expression in mtx format
# Save normalised counts - NOT scaled!
writeMM(hswounds_sub@assays$RNA@data, file = 'countsData_Wound1/matrix.mtx')
# save gene and cell names
write(x = rownames(hswounds_sub@assays$RNA@data), file = "countsData_Wound1/features.tsv")
write(x = colnames(hswounds_sub@assays$RNA@data), file = "countsData_Wound1/barcodes.tsv")

# Generate your meta data
hswounds_sub@meta.data$Cell = rownames(hswounds_sub@meta.data)
df = hswounds_sub@meta.data[, c('Cell', 'CellType')]
write.table(df, file ='hswounds_metadata_Wound1.tsv', sep = '\t', quote = F, row.names = F)



hswounds_sub <- subset(hswounds, subset = Condition == "Wound7")

# Write gene expression in mtx format
# Save normalised counts - NOT scaled!
writeMM(hswounds_sub@assays$RNA@data, file = 'countsData_Wound7/matrix.mtx')
# save gene and cell names
write(x = rownames(hswounds_sub@assays$RNA@data), file = "countsData_Wound7/features.tsv")
write(x = colnames(hswounds_sub@assays$RNA@data), file = "countsData_Wound7/barcodes.tsv")


# Generate your meta data
hswounds_sub@meta.data$Cell = rownames(hswounds_sub@meta.data)
df = hswounds_sub@meta.data[, c('Cell', 'CellType')]
write.table(df, file ='hswounds_metadata_Wound7.tsv', sep = '\t', quote = F, row.names = F)



hswounds_sub <- subset(hswounds, subset = Condition == "Wound30")

# Write gene expression in mtx format
# Save normalised counts - NOT scaled!
writeMM(hswounds_sub@assays$RNA@data, file = 'countsData_Wound30/matrix.mtx')
# save gene and cell names
write(x = rownames(hswounds_sub@assays$RNA@data), file = "countsData_Wound30/features.tsv")
write(x = colnames(hswounds_sub@assays$RNA@data), file = "countsData_Wound30/barcodes.tsv")


# Generate your meta data
hswounds_sub@meta.data$Cell = rownames(hswounds_sub@meta.data)
df = hswounds_sub@meta.data[, c('Cell', 'CellType')]
write.table(df, file ='hswounds_metadata_Wound30.tsv', sep = '\t', quote = F, row.names = F)


###############
# Run the cellphoneDB
cellphonedb method statistical_analysis \
/Users/zhuliu/Desktop/scRNA_STseq/proj_10X_woundhealing/03_results/01-Seurat-PreAnalysis/02_Seurat_BatchCorrection/s1_Seurat_allSample_harmony/s8_Cell_Cell_L_R_analysis_upCellTypes/CellPhoneDB/hswounds_metadata.tsv \
/Users/zhuliu/Desktop/scRNA_STseq/proj_10X_woundhealing/03_results/01-Seurat-PreAnalysis/02_Seurat_BatchCorrection/s1_Seurat_allSample_harmony/s8_Cell_Cell_L_R_analysis_upCellTypes/CellPhoneDB/countsData \
--counts-data=hgnc_symbol \
--threshold 0.1 \
--output-path=custom_folder



###############
# DFU project
seu_dfu <- readRDS("/Users/zhuliu/Desktop/sc st/compareStudies/Diabetic foot ulcers/step0_rawDFUdata.rds")
mt <- data.table::fread("metadata/DFU_metadata.txt")
seu_dfu <- seu_dfu[, mt$barcode]

identical(colnames(seu_dfu), mt$barcode)
seu_dfu$CellType <- mt$PredictedCellType

seu_dfu <- DietSeurat(seu_dfu, counts = TRUE, assays = "RNA")
seu_dfu <- NormalizeData(seu_dfu)

table(seu_dfu$Condition)

seu_dfu_sub <- subset(seu_dfu, subset = Condition == "H")
# Write gene expression in mtx format
# Save normalised counts - NOT scaled!
writeMM(seu_dfu_sub@assays$RNA@data, file = 'countData_healthy/matrix.mtx')
# save gene and cell names
write(x = rownames(seu_dfu_sub@assays$RNA@data), file = "countData_healthy/features.tsv")
write(x = colnames(seu_dfu_sub@assays$RNA@data), file = "countData_healthy/barcodes.tsv")

# Generate your meta data
seu_dfu_sub@meta.data$Cell = rownames(seu_dfu_sub@meta.data)
df = seu_dfu_sub@meta.data[, c('Cell', 'CellType')]
write.table(df, file ='DFU_metadata_healthy.tsv', sep = '\t', quote = F, row.names = F)


seu_dfu_sub <- subset(seu_dfu, subset = Condition == "DFU_H")
# Write gene expression in mtx format
# Save normalised counts - NOT scaled!
writeMM(seu_dfu_sub@assays$RNA@data, file = 'countData_DFU_H/matrix.mtx')
# save gene and cell names
write(x = rownames(seu_dfu_sub@assays$RNA@data), file = "countData_DFU_H/features.tsv")
write(x = colnames(seu_dfu_sub@assays$RNA@data), file = "countData_DFU_H/barcodes.tsv")

# Generate your meta data
seu_dfu_sub@meta.data$Cell = rownames(seu_dfu_sub@meta.data)
df = seu_dfu_sub@meta.data[, c('Cell', 'CellType')]
write.table(df, file ='DFU_metadata_DFU_H.tsv', sep = '\t', quote = F, row.names = F)


seu_dfu_sub <- subset(seu_dfu, subset = Condition == "DFU_NH")
# Write gene expression in mtx format
# Save normalised counts - NOT scaled!
writeMM(seu_dfu_sub@assays$RNA@data, file = 'countData_DFU_NH/matrix.mtx')
# save gene and cell names
write(x = rownames(seu_dfu_sub@assays$RNA@data), file = "countData_DFU_NH/features.tsv")
write(x = colnames(seu_dfu_sub@assays$RNA@data), file = "countData_DFU_NH/barcodes.tsv")

# Generate your meta data
seu_dfu_sub@meta.data$Cell = rownames(seu_dfu_sub@meta.data)
df = seu_dfu_sub@meta.data[, c('Cell', 'CellType')]
write.table(df, file ='DFU_metadata_DFU_NH.tsv', sep = '\t', quote = F, row.names = F)

