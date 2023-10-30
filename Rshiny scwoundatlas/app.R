# check the installed packages
check.packages <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])] #Check if packages has been installed; installed.packages() return a matrix package names, library paths and version numbers.
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
}

packages <- c("data.table", "Matrix", "hdf5r", "reticulate", "ggplot2", 
              "gridExtra", "glue", "readr", "RColorBrewer", "R.utils", "Seurat",
              "shiny", "shinyhelper", "data.table", "Matrix", "DT", "hdf5r", 
              "reticulate", "ggplot2", "gridExtra", "magrittr", "ggdendro")


check.packages(packages)

# run the app
shiny::runApp()


### Change the file formats and defaults
#sc1conf = readRDS("sc1conf.rds")
#sc1conf$ID[1] = "SampleID"
#sc1conf=sc1conf[-c(2,4),]
#ct.cols <- paste("#d94701", "#fdae61", "#fd8d3c", "#fdbe85", #Basal clusters
#                 "#33A02C", "#72BF5A", "#B2DF8A", #Spinous clusters
#                 "#f768a1", "#d4b9da", #Granular, Hair follicle
#                 "#737373", #MEL
#                 "#0570b0", "#3690c0", "#92c5de", "#d1e5f0", #Fibroblast clusters
#                 "#c0a700", "#1a9850", "#fb9a99", "#8d4720",  # Schwann,PCvSMC,LE,VE
#                 "#35978f", "#41b6c4", "#80cdc1","#df65b0", #"NK-cell", "Th", "Plasma_Bcell", "Mast-cell", 
#                 "#dd3497", "#807dba","#6a3d9a","#9e9ac8", "#b15928" #"Mono-Mac", "cDC1", "cDC2", "DC3", "LC"
#                 , sep = "|", collapse = "|")
#sc1conf$fCL[5] <- ct.cols
#
#main.col <- paste("#ce8b1a", "#737373", "#d4b9da", "#3690c0", "#1a9850", 
#              "#fb9a99", "#c0a700", "#df65b0", "#807dba", "#35978f"
#              , sep = "|", collapse = "|")
#sc1conf$fCL[6] <- main.col
#saveRDS(sc1conf, "sc1conf.rds")
#
#
#sc1def  = readRDS("sc1def.rds")
#sc1def$meta1 <- "CellTypes"
#sc1def$meta2 <- "Condition"
#gene1 <- "FOSL1";names(gene1) <- "FOSL1"
#gene2 <- "CXCL1";names(gene2) <- "CXCL1"
#sc1def$gene1 <- gene1
#sc1def$gene2 <- gene2
#genes <- c("MMP1", "MMP3", "AREG", "TNFRSF12A", "FGFBP1", "S100A2", "S100A10",
#           "KRT6A", "KRT6B", "KRT17", "KRT16", "S100A7", "S100A8", "S100A9")
#names(genes) <- c("MMP1", "MMP3", "AREG", "TNFRSF12A", "FGFBP1", "S100A2", "S100A10",
#                  "KRT6A", "KRT6B", "KRT17", "KRT16", "S100A7", "S100A8", "S100A9")
#sc1def$genes <- genes
#sc1def$grp2 <- "Condition"
#sc1def$grp1 <- "CellTypes"
#saveRDS(sc1def, "sc1def.rds")
#
#sc1gene = readRDS("sc1gene.rds")
#
#
#sc1meta = readRDS("sc1meta.rds")
#sc1meta=sc1meta[,-c(3,5)]
#colnames(sc1meta) = c("CellID", "SampleID", "Gender", "Condition", "seurat_clusters", "CellTypes",
#                      "MainCellTypes", "UMAP_1", "UMAP_2")
#saveRDS(sc1meta, "sc1meta.rds")
