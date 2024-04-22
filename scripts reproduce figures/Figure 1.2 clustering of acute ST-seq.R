library(Seurat)
library(tidyverse)
library(magrittr)
library(ggplot2)
library(patchwork)
library(circlize)

setwd("")
rm(list = ls())

##############################################
# Fig 1 supple Dim plots of batch correction
# the plots without batch corrections drawn by the Rmd script
hswound.STcom <- readRDS("Seurat_STseq_integrate_update.rds")

DimPlot(hswound.STcom, label = T, label.size = 10) + NoAxes()

p1 = DimPlot(hswound.STcom, group.by = "Donor") + NoAxes()
p2 = DimPlot(hswound.STcom, group.by = "Condition") + NoAxes()
p3 = DimPlot(hswound.STcom, group.by = "Sample_name") + NoAxes()
#pdf("st_Fig_1_Harmony.pdf", useDingbats = F, width = 12, height = 8)
wrap_plots(list(p1, p2, p3), ncol = 2)
#dev.off()

mg <- data.table::fread("Seurat_STseq_integrated_markerGene.txt")

# re-annotate the cell clusters
fac_levs <- c("Basal_epi", "Suprabasal_epi", "Wound_edges",
              "Hair_follicle", "Papillary_dermis", 
              "Sebaceous_gland", "Sweat_gland", "Sweatgland_FB",
              "FB1", "FB2", 
              "Dermis1", "Dermis2","Smooth_muscle", 
              "Immune",  "Immune_endo", "LE", "Mast_cell")

current_cluster_ids <- c(10, 8, 5,
                         16, 9, 13, 14, 4,
                         7, 3,
                         1, 15, 11, 
                         6, 2, 12, 17) # List of current cluster IDs
current_cluster_ids <- seq(1:17)
tmp <- plyr::mapvalues(x=mg$cluster, from=current_cluster_ids, to=fac_levs)
mg$clusterAnnotation <- tmp
mg_f <- mg %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC) %>% arrange(cluster)

tmp_mg <- plyr::mapvalues(x=mg_f$clusterAnnotation, from=fac_levs, to=current_cluster_ids)
mg_f$newcluster <- tmp_mg

data.table::fwrite(mg_f, "Markers_top50maers.txt", sep = "\t")

tmp <- plyr::mapvalues(x=hswound.STcom$SCT_snn_res.0.5, from=current_cluster_ids, to=fac_levs)
hswound.STcom$AnnoType <- tmp
hswound.STcom$AnnoType <- factor(hswound.STcom$AnnoType, levels = fac_levs)
hswound.STcom$AnnoTypeNum <- paste0(hswound.STcom$SCT_snn_res.0.5, " ", tmp)
DimPlot(hswound.STcom, label = T, group.by = "AnnoType") + NoAxes()
DimPlot(hswound.STcom, label = T, group.by = "AnnoTypeNum") + NoAxes()

########################
# Fig 1 UMAP of clusters
# define the colors for clusters
ct.cols <- c("#33A02C", "#fdbe85", "#cb181d",
             "#4c0097", "#d4b9da",
             "#54278f", "#6a51a3", "#807dba",
             "#056aa6", "#0570b0",
             "#3690c0", "#74a9cf", "#99d8c9", 
             "#bf057a", "#dd3497", "#f768a1", "#fb6a4a"
             )

names(ct.cols) <- fac_levs

#pdf("ST_Fig1_UMAP_clusters.pdf", useDingbats = F, width = 6, height = 5)
DimPlot(hswound.STcom, label = F, label.size = 5, cols = ct.cols, group.by = "AnnoType") + NoAxes()
#dev.off()


################################################
##---- Figure 1 dendrogram all cell types ----##
################################################
DefaultAssay(hswound.STcom) <- "SCT"
hswound.STcom$AnnoType <- factor(hswound.STcom$AnnoType, levels = fac_levs)
hswound.STcom <- SetIdent(hswound.STcom, value = hswound.STcom@meta.data$AnnoType)

hswound.STcom <- BuildClusterTree(object = hswound.STcom, assay = "SCT", slot = "scale.data")
# pull the tree
hswound.clustertree <- Tool(object = hswound.STcom, slot = "BuildClusterTree")
# plot the tree
suppressMessages(require(ape))
{plot.phylo(x = hswound.clustertree, direction = "rightwards", use.edge.length = TRUE, font = 4, tip.color = ct.cols, edge.width = 2, label.offset = 5, no.margin = TRUE)
  tiplabels(pch = 19, col = ct.cols, adj = 3, cex = 2)}

#pdf("Fig1_ST_dendrogram.pdf", useDingbats = F, width = 4, height = 4)
#dev.off()


##########################################
# Fig 1 Projecting clusters onto HE images
img.names <- unique(hswound.STcom@meta.data$Sample_name) %>% as.character()
img.names.plot <- img.names[grep("Donor2", img.names)]
plot_list <- list()
ratio_list <- list()
for (i in seq_along(img.names.plot)) {
  coord <- GetTissueCoordinates(object = hswound.STcom, image = img.names.plot[i])
  # calculate the aspect ratio of rows to columns
  myratio <- (max(coord$imagerow) - min(coord$imagerow)) / (max(coord$imagecol) - min(coord$imagecol))
  ratio_list[[i]] <- myratio
  plot_list[[i]] <- SpatialDimPlot(hswound.STcom, #group.by = "seurat_clusters",
                                 images = img.names.plot[i], stroke = 0, 
                                 pt.size.factor = 1.8*myratio, 
                                 combine = TRUE,
                                 alpha = c(1, 1), 
                                 image.alpha = 0, crop = TRUE, 
                                 group.by = "AnnoType",
                                 cols = ct.cols) +
    theme(aspect.ratio = myratio) 
}
# replace the wound1 plot
coord <- GetTissueCoordinates(object = hswound.STcom, image = img.names.plot[2])
# calculate the aspect ratio of rows to columns
myratio <- (max(coord$imagerow) - min(coord$imagerow)) / (max(coord$imagecol) - min(coord$imagecol))

plot_list[[2]] = SpatialDimPlot(hswound.STcom, #group.by = "seurat_clusters",
                                images = img.names.plot[2], stroke = 0, 
                                pt.size.factor = 1.5*myratio, 
                                combine = TRUE,
                                alpha = c(1, 1), 
                                image.alpha = 0, crop = TRUE, 
                                group.by = "AnnoType",
                                cols = ct.cols) +
  theme(aspect.ratio = myratio) 
# apply theme to each feature plot
#for(i in 1:length(plot_list)) {
#  plot_list[[i]] <- plot_list[[i]] + ggtitle(img.names.plot[i])
#}   
#cowplot::plot_grid(plotlist = plot_list, ncol = 4)
pdf("ST_Fig_1_donor2_cluster_projection_up.pdf", useDingbats = F, width = 8, height = 4)
patchwork::wrap_plots(plot_list, ncol = 2) & theme(legend.position = 'none')
dev.off()


############################
# Fig 1 Plot of marker genes
markers <- data.table::fread("Seurat_STseq_integrated_markerGene.txt") %>% 
  mutate(diff.pct = pct.1 - pct.2) %>% filter(gene_biotype == "protein_coding")
markers$cluster <- as.factor(markers$cluster)
tmp_1 <- plyr::mapvalues(x=markers$cluster, from=current_cluster_ids, to=fac_levs)
markers$AnnoType <- tmp_1
markers <- markers %>% select(13,12, everything()) 

# check the top marker genes of different clusters
top.markers <- c("KRT15", "KRT5", "COL17A1", "FGFR3", "SYT8",
                 "KRT1", "KRT2", "LORICRIN", "FLG", "SBSN",
                 "KRT6A", "KRT6B", "KRT6C", "KRT16", "KRT17", "S100A8", "S100A9",
                 "ITGB4", "HR", "FZD7", "KRT75",
                 "DSG1", "PERP",
                 "MGST1", "FADS2", "FAR2", "KRT79",
                 "KRT77", "MMP7", "DEFB1", "SCGB1D2",
                 "DCD", "SCGB2A2", "MUCL1", "PIP",
                 "ADAM12", "ASPN", "POSTN", "COL12A1",
                 "MMP2", "FBLN1", "CCN5", "TNXB", "PI16",
                 "SAMHD1", "CCDC80", 
                 "APOD", "IGFBP6",
                 "MYH11", "ACTA2", "TAGLN", "MYL9",
                 "FABP4", "LYZ", "CD36", "CD68", "TYROBP", "CD163",
                 "HLA-DPA1", "HLA-DRA", "PECAM1", "VWF", "CD74",
                 "CCL21", "TFF3", "MMRN1", "LYVE1",
                 "TPSB2", "TPSAB1", "MS4A2"
                 )

p = DotPlot(hswound.STcom, features = rev(top.markers), assay = "Spatial", 
            cols = c("white", "red"),  dot.min = 0.1, col.min = 0, 
            group.by = "AnnoType", dot.scale = 4) + coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 90),
        axis.text.y = element_text(size = 12),
        panel.border = element_rect(colour = "black", fill=NA))

#pdf("ST_clusters_markers_dotplot_v.pdf", useDingbats = F, width = 8, height = 14)
print(p)
#dev.off()

genes <- c("KRT15", "KRT1", "IVL", "FLG", "PMEL", "COL1A2", "ADAM12", "SFRP2", "SFRP4", 
           "CD3D", "CD3E", "LYZ", "IL1B", "CD68", "CD163", "C1QB", "MKI67", "CD14", "CD4",
           "CD8A", "CD8B", "CD207", "CD79A", "JCHAIN", "HEY1", "SEMA3G", "MYH11", "COL6A3")
FeaturePlot(hswound.STcom, features = genes[1:14], cols = c("grey90", "red"), ncol = 5)
FeaturePlot(hswound.STcom, features = genes[15:28], cols = c("grey90", "red"), ncol = 5)

#Check the expression of neutrophil marker genes
FeaturePlot(hswound.STcom, features = c('CSF3R', 'FPR1', 'FCGR3B', 'NAMPT', 'MNDA', 
                                        'G0S2', 'CMTM2', 'CXCR2', 'PROK2'), 
            cols = c("grey90", "red"),
            order = T)

FeaturePlot(hswound.STcom, features = c('EGF'), 
            cols = c("grey90", "red"))
#########################################################
# Fig 1 supple Violin plots of metrics of quality control
p_vln <- VlnPlot(
  hswound.STcom,
  features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"),
  #group.by = "Sample_name", cols = sm.color
  ncol = 2, pt.size=0, cols = ct.cols)

pdf("st_Fig_1_metrics_clusters_Vlnplots.pdf", useDingbats = F, width = 12, height = 7)
p_vln
dev.off()


###########################################################################
# Fig 1 Calculate the scores of Bas & Spi migratory clusters onto HE images
# read the marker genes from scRNA-seq data
markers <- data.table::fread("/Top200_marker_genes_cluster.txt") %>% 
  group_by(cluster) %>% arrange(cluster, desc(avg_log2FC)) %>% ungroup()

basmig_gene_hs <- list(Bas_mig= markers %>% dplyr::filter(cluster == "Bas_III") %>% 
                         top_n(50, pct.dif) %>% top_n(30, avg_log2FC) %>% pull(gene))

spimig_gene_hs <- list(Spi_mig= markers %>% dplyr::filter(cluster == "Spi_V") %>% 
                         top_n(50, pct.dif) %>% top_n(30, avg_log2FC) %>% pull(gene))

basprolif_gene_hs <- list(bas_prolif= markers %>% dplyr::filter(cluster == "Bas_II") %>% 
                          top_n(100, pct.dif) %>% top_n(50, avg_log2FC) %>% pull(gene))

hswound.STcom <- AddModuleScore(object = hswound.STcom, features = basmig_gene_hs, name = "Bas_mig_Score", assay = "Spatial") 
hswound.STcom <- AddModuleScore(object = hswound.STcom, features = spimig_gene_hs, name = "Spi_mig_Score", assay = "Spatial") 
hswound.STcom <- AddModuleScore(object = hswound.STcom, features = basprolif_gene_hs, name = "Bas_prolif_Score", assay = "Spatial") 


p4 <- SpatialFeaturePlot(hswound.STcom, features = c("Bas_mig_Score1", "Spi_mig_Score1"),
                         images = "Donor2_Wound30", #"Donor2_Skin", "Donor2_Wound1", "Donor2_Wound7", "Donor2_Wound30"
                         stroke = 0, combine = T, pt.size.factor = 1.8, 
                         alpha = c(1, 1), image.alpha = 0.6, crop = TRUE,
                         min.cutoff = 0) +
  theme(aspect.ratio = myratio) 
  scale_fill_gradientn(limits = c(0,5), breaks = c(0,1,2,3,4,5), 
                       colours = rev(RColorBrewer::brewer.pal(11, "Spectral"))) 
p4
#pdf("ST_Fig_1_donor2_cluster_bas_spi_Scores.pdf", useDingbats = F, width = 8, height = 8)
patchwork::wrap_plots(list(p1, p2, p3, p4), ncol = 2) & theme(legend.position = 'top',
                                                              legend.text = element_text(size = 6),
                                                              legend.title = element_text(size = 8))
#dev.off()


###########################################################
# Fig 1 Projecting deconvolution results onto UMAP clusters
deconv <- read_csv("/deconvolution/cell2location_0706/st_deconv_celltypes.csv")
colnames(deconv) <- gsub("q05cell_abundance_w_sf_", "", colnames(deconv))

mt <- hswound.STcom@meta.data %>% rownames_to_column(var = "barcode") %>% left_join(., deconv, by=c("barcode" = "barcode")) %>% 
  column_to_rownames(var = "barcode")
hswound.STcom@meta.data <- mt

#cc = scale_color_gradientn(colors = c("grey", "yellow", "red", "black"))
small.leg <- theme(legend.text = element_text(size = 10), legend.key.width = unit(0.2, "cm"))

p1 <- FeaturePlot(hswound.STcom, features = "Bas-mig", slot = "data", pt.size = 0.01,
            order = T, reduction = "umap") + NoAxes() + small.leg + 
  scale_color_gradientn(colors = c("grey", "yellow", "red", "black"), limits = c(0, 10)) +
  theme(plot.title = element_text(size = 12))
p2 <- FeaturePlot(hswound.STcom, features = "Spi-mig", slot = "data", pt.size = 0.01,
                  order = T, reduction = "umap") + NoAxes() + small.leg + 
  scale_color_gradientn(colors = c("grey", "yellow", "red", "black"), limits = c(0, 10)) +
  theme(plot.title = element_text(size = 12))

#pdf("ST_Fig_1_donor2_deconv_projection.pdf", useDingbats = F, width = 4, height = 6)
p1 / p2 + plot_layout(guides = "collect")
#dev.off()


#####################
# fibroblast scores
fib.marker <- c('COL1A1', 'COL1A2', 'LUM')

markers <- data.table::fread("s1_clean_fibroblast_mg_MianCellTypes.txt") %>% 
  mutate(pct.dif = pct.1 - pct.2)

fb1 <- list(gene1 = unique(c(markers %>% dplyr::filter(cluster == "FB-I") %>% 
            top_n(100, pct.dif) %>% top_n(50, avg_log2FC) %>% pull(gene), fib.marker)))

fb2 <- list(gene2 = unique(c(markers %>% dplyr::filter(cluster == "FB-II") %>% 
            top_n(100, pct.dif) %>% top_n(50, avg_log2FC) %>% pull(gene), fib.marker)))

fb3 <- list(gene3 = unique(c(markers %>% dplyr::filter(cluster == "FB-III") %>% 
            top_n(100, pct.dif) %>% top_n(50, avg_log2FC) %>% pull(gene), fib.marker)))

fb4 <- list(gene4 = unique(c(markers %>% dplyr::filter(cluster == "FB-IV") %>% 
            top_n(100, pct.dif) %>% top_n(50, avg_log2FC) %>% pull(gene), fib.marker)))

hswound.STcom <- AddModuleScore(object = hswound.STcom, features = fb1, name = "FB_I_Score", assay = "Spatial") 
hswound.STcom <- AddModuleScore(object = hswound.STcom, features = fb2, name = "FB_II_Score", assay = "Spatial") 
hswound.STcom <- AddModuleScore(object = hswound.STcom, features = fb3, name = "FB_III_Score", assay = "Spatial") 
hswound.STcom <- AddModuleScore(object = hswound.STcom, features = fb4, name = "FB_prolif_Score", assay = "Spatial") 

images.names <- unique(hswound.STcom$Sample_name) %>% as.character()

images.io <- images.names[grep("Donor4", images.names)]
plot_st <- list()
for (i in seq_along(images.io)) {
  plot_st[[i]] <- SpatialFeaturePlot(hswound.STcom, features = c("FB_I_Score1", "FB_II_Score1", "FB_III_Score1", "FB_prolif_Score1"),
                                     images = images.io[i], #"Donor2_Skin", "Donor2_Wound1", "Donor2_Wound7", "Donor2_Wound30"
                                     stroke = 0, combine = T, pt.size.factor = 1.2, 
                                     alpha = c(1, 1), image.alpha = 0.6, crop = TRUE,
                                     min.cutoff = 0)
    
}

#pdf("ST_Fig_1_donor2_cluster_bas_spi_Scores.pdf", useDingbats = F, width = 8, height = 8)
patchwork::wrap_plots(plot_st, ncol = 4) & theme(legend.position = 'top',
                                                              legend.text = element_text(size = 6),
                                                              legend.title = element_text(size = 8))
#dev.off()


###########################################
# Figure 1f NMF analysis
###########################################
# Re-draw the spatial feature plots of niches
# load the niche data
library(Seurat)
library(tidyverse)
library(scCustomize)
library(RColorBrewer)
library(viridis)

nmf_n15 <- readxl::read_xlsx("microenvironment_NMF_n15.xlsx", sheet = 1)
colnames(nmf_n15)
newFactor = paste0("Niche",c(1,2,3,4,5,6,
                             7,8,9,10,11,
                             12,13,14,15))
rawFactor = paste0("mean_nUMI_factorsfact_",c(0,1,2,3,4,5,
                                              6,14,12,13,8,
                                              11,9,10,7))
# rename the factor, keep it same as in Figure 1f
tmp_match <- plyr::mapvalues(colnames(nmf_n15)[4:18], from = rawFactor, to = newFactor)
colnames(nmf_n15)[4:18] <- tmp_match

# load the spatial data
st_seu <- readRDS("/Users/zhuliu/Desktop/sc st/ST plotting scripts/Seurat_STseq_integrate_update.rds")
mt <- st_seu@meta.data %>% rownames_to_column(var = "barcode")
mt <- mt %>% left_join(., nmf_n15[, -c(2:3)], by=c("barcode" = "SpotIndex" ))
identical(mt$barcode, nmf_n15$SpotIndex)

# modified the metadata of spatial seurat object
mt_f <- mt %>% column_to_rownames(var = "barcode")
st_seu@meta.data <- mt_f


########################
# plot interested niches
img.names <- unique(st_seu@meta.data$Sample_name) %>% as.character();img.names
# choose which donor
plot_donors <- img.names[grep("Donor4_", img.names)];plot_donors

ratio_list <- list()
for (i in seq_along(plot_donors)) {
  coord <- GetTissueCoordinates(object = st_seu, image = plot_donors[i])
  # calculate the aspect ratio of rows to columns
  myratio <- (max(coord$imagerow) - min(coord$imagerow)) / (max(coord$imagecol) - min(coord$imagecol))
  ratio_list[[i]] <- myratio
}
ratio_list

colours = magma(10)
plot_list <- list()
for (i in seq_along(plot_donors)) {
  plot_list[[i]] <- SpatialFeaturePlot(st_seu, 
                                       features = "Niche1",
                                       images = plot_donors[i], 
                                       stroke = 0, combine = T, 
                                       pt.size.factor = 1.2, 
                                       #pt.size.factor = ratio_list[[i]], 
                                       alpha = c(1, 1), 
                                       max.cutoff = "q98",
                                       image.alpha = 0, crop = FALSE) +
    #theme(aspect.ratio = ratio_list[[i]]) + 
    #scale_fill_gradientn(limits=c(0, 1), colours=magma(10), na.value = "#FCFDBFFF") # niche13
    #scale_fill_gradientn(limits=c(0, 5), colours=magma(10), na.value = "#FCFDBFFF") # niche8
    #scale_fill_gradientn(limits=c(0, 3), colours=magma(10), na.value = "#FCFDBFFF") # niche7
    scale_fill_gradientn(limits=c(0, 10), colours=magma(10), na.value = "#FCFDBFFF") # niche10
}
patchwork::wrap_plots(plot_list, ncol = 2)

#pdf("ST_Donor4_genePlot_Niche10_defaultColor.pdf", useDingbats = F, width = 8, height = 8)
pdf("ST_Donor4_genePlot_Niche1_custColor.pdf", useDingbats = F, width = 18, height = 18)
patchwork::wrap_plots(plot_list, ncol = 4)
dev.off()


########################################################
# SpatialfeaturePlot of gene expression in spatial data
########################################################
hswound.STcom <- readRDS("Seurat_STseq_integrate_update.rds")
img.names <- unique(hswound.STcom@meta.data$Sample_name) %>% as.character()

plot_donors <- img.names[grep("Donor2_", img.names)]

ratio_list <- list()
for (i in seq_along(plot_donors)) {
  coord <- GetTissueCoordinates(object = hswound.STcom, image = plot_donors[i])
  # calculate the aspect ratio of rows to columns
  myratio <- (max(coord$imagerow) - min(coord$imagerow)) / (max(coord$imagecol) - min(coord$imagecol))
  ratio_list[[i]] <- myratio
}
ratio_list
RColorBrewer::display.brewer.all()
colours = rev(RColorBrewer::brewer.pal(11, "Spectral"))
plot_list <- list()
for (i in seq_along(plot_donors)) {
  plot_list[[i]] <- SpatialFeaturePlot(hswound.STcom, 
                                       features = "FOSL1",
                                       images = plot_donors[i], 
                                       stroke = 0, combine = T, 
                                       pt.size.factor = 1.8*ratio_list[[i]], 
                                       alpha = c(1, 1), 
                                       image.alpha = 0.6, crop = TRUE) + 
    theme(aspect.ratio = ratio_list[[i]]) #+
  #scale_fill_gradientn(limits = c(0,2), breaks = c(0,0.5,1,1.5,2))
}
patchwork::wrap_plots(plot_list, ncol = 2)

plot_list[3]

patchwork::wrap_plots(plot_list, ncol = 2)
pdf("ST_Donor2_genePlot_FOSL1_HE.pdf", useDingbats = F, width = 8, height = 8)
patchwork::wrap_plots(plot_list, ncol = 4)
dev.off()


# replace the skin/wound1 plot
coord <- GetTissueCoordinates(object = hswound.STcom, image = plot_donors[2])
# calculate the aspect ratio of rows to columns
myratio <- (max(coord$imagerow) - min(coord$imagerow)) / (max(coord$imagecol) - min(coord$imagecol))

plot_list[[2]] = SpatialFeaturePlot(hswound.STcom, 
                                    features = "FOSL1",
                                    images = plot_donors[2], # change the Donor
                                    stroke = 0, combine = T, 
                                    pt.size.factor = 1.6*myratio, 
                                    alpha = c(1, 1), 
                                    image.alpha = 0, crop = TRUE) + 
  theme(aspect.ratio = myratio)


pdf("ST_Donor2_genePlot_FOSL1_STexp2.pdf", useDingbats = F, width = 8, height = 8)
patchwork::wrap_plots(plot_list, ncol = 2)
dev.off()
