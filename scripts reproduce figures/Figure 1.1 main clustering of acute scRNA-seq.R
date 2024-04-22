library(Seurat)
library(tidyverse)
library(magrittr)
library(ggplot2)
library(patchwork)
library(circlize)

setwd("")
rm(list = ls())

############################################################
##Load the final clean rds object (human skin wound healing)
#This object was filtered after subclustering analysis since there were still some doublets existing in the subclusters
hswound.combined.sct <- readRDS("allcombined_wounds_newAnnotation.rds")

# The Bas-II and HF were excluded from subclustering analysis of keratinocytes
fact_lev <- c("Bas-I", "Bas-II", "Bas-prolif", "Bas-mig", 
              "Spi-I", "Spi-II", "Spi-mig", 
              "Gra-I", "HF", 
              "MEL", 
              "FB-I", "FB-II", "FB-III", "FB-prolif", 
              "Schwann", "PC-vSMC", "LE", "VE", 
              "NK-cell", "Th", "Plasma_Bcell", "Mast-cell", 
              "Mono-Mac", "cDC1", "cDC2", "DC3", "LC")

ct.cols <- c("#d94701", "#fdae61", "#fd8d3c", "#fdbe85", #Basal clusters
             "#33A02C", "#72BF5A", "#B2DF8A", #Spinous clusters
             "#f768a1", "#d4b9da", #Granular, Hair follicle
             "#737373", #MEL
             "#0570b0", "#3690c0", "#92c5de", "#d1e5f0", #Fibroblast clusters
             "#c0a700", "#1a9850", "#fb9a99", "#8d4720",  # Schwann,PCvSMC,LE,VE
             "#35978f", "#41b6c4", "#80cdc1","#df65b0", #"NK-cell", "Th", "Plasma_Bcell", "Mast-cell", 
             "#dd3497", "#807dba","#6a3d9a","#9e9ac8", "#b15928" #"Mono-Mac", "cDC1", "cDC2", "DC3", "LC"
)
# add the color names using new cell types
names(ct.cols) <- fact_lev

##########################################
##---- Figure 1 UMAP all cell types ----##
##########################################
p_umap.clean <- (DimPlot(object = hswound.combined.sct, 
         reduction = "umap", 
         label = F, 
         group.by = "newCellTypes", cols = ct.cols) + NoAxes() + ggtitle("")) + 
  (DimPlot(object = hswound.combined.sct, 
           reduction = "umap", 
           label = T, 
           group.by = "newAnnoNum", cols = ct.cols, label.size = 2) + NoAxes() + ggtitle(""))
p_umap.clean

pdf("Fig1_mainUMAP_clean.pdf", useDingbats = F, width = 6, height = 4)
p_umap.clean
dev.off()

fact_lev_main <- c("Keratinocyte", "Melanocyte", "Hair follicle", "Fibroblast", "Pericyte & Smooth muscle", 
              "Endothelial", "Schwann", "Mast", 
              "Myeloid", "Lymphoid")

hswound.combined.sct$newMainCellTypes <- factor(hswound.combined.sct$newMainCellTypes, levels = fact_lev_main)
main.col <- c("#ce8b1a", "#737373", "#d4b9da", "#3690c0", "#1a9850", 
              "#fb9a99", "#c0a700", "#df65b0", "#807dba", "#35978f")
names(main.col) <- fact_lev_main
p_umap.clean.main <- (DimPlot(object = hswound.combined.sct, 
                         reduction = "umap", 
                         label = F, 
                         group.by = "newMainCellTypes", label.size = 3, cols = main.col) + NoAxes() + ggtitle(""))

#pdf("Fig1_mainUMAP_clean_combined.pdf", useDingbats = F, width = 7, height = 4)
p_umap.clean.main
#dev.off()


#############################################
##---- Figure 1 markers all cell types ----##
#############################################
# This step can set the active.ident cell names into defined number labels
hswound.combined.sct$newCellTypes <- factor(hswound.combined.sct$newCellTypes, levels = rev(fact_lev))
hswound.combined.sct$newCellTypes <- factor(hswound.combined.sct$newCellTypes, levels = fact_lev)
hswound.combined.sct <- SetIdent(hswound.combined.sct, value = hswound.combined.sct@meta.data$newCellTypes)

canonical_markers <- c("KRT5", "KRT15", "MT1G", #"COL17A1",#Bas-I, Bas-II
                       "CENPF", "MKI67", "MMP3", "MMP1", #Bas-prolif, Bas-mig
                       "KRT10", "MT1X", "DMKN", #"LGALS7B", #Spi-I and Spi-II
                       "KRT16", "KRT6A", "FLG", #"LORICRIN", #Spi-mig, Gra-I
                       "SOX9", "KRT17", "ADGRL3", #HF
                       "TYRP1", "PMEL", #MEL
                       "COL1A1", "POSTN", "ADAM12", "CCL19", "CXCL12", "INHBA", #FB-I to IV #"PTGDS", 
                       "SOX10", "SOX2",
                       "ACTA2", "MYH11", "TAGLN", #pc-vsmc
                       "PECAM1", "CCL21", "TFF3", "SELE", "ACKR1",#"VWF",  #LE, VE
                       "NKG7", "GNLY","CD3D", #"TRAC", #NK-T, Th
                       "JCHAIN", "MS4A1", "TPSAB1", "TPSB2", #Bcell, mast-cell
                       "CD163", "CD14", "CLEC9A", "XCR1", "CD1C", "IL1R2", "LAMP3", "FSCN1","CD207", "CD1A")

plot_marker <- DotPlot(hswound.combined.sct, features = canonical_markers, 
        group.by = "newCellTypes", cols = c("white", "#cb181d"), 
        dot.scale = 4, col.min = 0, dot.min = 0.1
) + theme(axis.text.x = element_text(angle = 90),
          panel.border = element_rect(colour = "black", fill=NA)) + labs(x="", y="") #coord_flip() + 


pdf("Fig1_markers_new_vertical.pdf", useDingbats = F, width = 12, height = 6)
plot_marker
dev.off()


################################################
##---- Figure 1 dendrogram all cell types ----##
################################################
DefaultAssay(hswound.combined.sct) <- "SCT"
hswound.combined.sct$newCellTypes <- factor(hswound.combined.sct$newCellTypes, levels = fact_lev)
hswound.combined.sct <- SetIdent(hswound.combined.sct, value = hswound.combined.sct@meta.data$newCellTypes)

hswound.combined.sct <- BuildClusterTree(object = hswound.combined.sct, assay = "SCT", slot = "scale.data")
# pull the tree
hswound.clustertree <- Tool(object = hswound.combined.sct, slot = "BuildClusterTree")
# plot the tree
suppressMessages(require(ape))
{plot.phylo(x = hswound.clustertree, direction = "rightwards", use.edge.length = TRUE, font = 4, tip.color = ct.cols, edge.width = 2, label.offset = 5, no.margin = TRUE)
  tiplabels(pch = 19, col = ct.cols, adj = 3, cex = 2)}

pdf("Fig1_dendrogram.pdf", useDingbats = F, width = 4, height = 4)
dev.off()


#################################################
##---- Figure 1 Supplementary figures (QC) ----##
#################################################
# This step can set the active.ident cell names into defined number labels
hswound.combined.sct$newCellTypes <- factor(hswound.combined.sct$newCellTypes, levels = fact_lev)
hswound.combined.sct <- SetIdent(hswound.combined.sct, value = hswound.combined.sct@meta.data$newCellTypes)

# Determine metrics to plot present in hswound.combined.sct@meta.data
metrics <-  c("nCount_RNA", "nFeature_RNA", "percent.mt", "doublet_scores")

VlnPlot(hswound.combined.sct, features = metrics, group.by='newCellTypes', pt.size = 0, ncol = 2, cols = ct.cols)

pdf("Fig1_supple_metrics_vlnplot.pdf", useDingbats = F, width = 12, height = 6)
VlnPlot(hswound.combined.sct, features = metrics, group.by='newCellTypes', pt.size = 0, ncol = 2, cols = ct.cols)
dev.off()

p_umap.harm2 <- (DimPlot(object = hswound.combined.sct, 
                         reduction = "umap", 
                         label = F, 
                         group.by = "Condition") + ggtitle(""))

pdf("Fig1_supple_UMAP_harmony_donor_condition.pdf", useDingbats = F, width = 14, height = 6)
p_umap.harm1 + p_umap.harm2
dev.off()


 ###########################################################
##---- Figure 1 Supple umap without batch correction ----##
###########################################################
library(Seurat)
library(SeuratObject)
library(ggplot2)

# load the rds object from Uppmax without-batch-correction fold
hswound.combined.sct <- readRDS("/Users/zhuliu/Downloads/s1_cleanSeurat_withoutharmony_allSamples_clusters.rds")
hswound.combined.sct$orig.ident <- factor(x = hswound.combined.sct$orig.ident, levels =
                                        c("PWH26D0", "PWH26D1", "PWH26D7","PWH26D30",
                                          "PWH27D0", "PWH27D1", "PWH27D7","PWH27D30",
                                          "PWH28D0", "PWH28D1", "PWH28D7","PWH28D30"))
hswound.combined.sct$Condition <- factor(x = hswound.combined.sct$Condition, levels =
                                            c("Skin", "Wound1", "Wound7", "Wound30"))

p_umap.withoutharm2 <- (DimPlot(object = hswound.combined.sct, 
                        reduction = "umap", 
                        label = F, 
                        group.by = "Condition") + ggtitle(""))

pdf("Fig1_supple_UMAP_withountHarmony_donor_condition.pdf", useDingbats = F, width = 14, height = 6)
p_umap.withoutharm1 + p_umap.withoutharm2
dev.off()

