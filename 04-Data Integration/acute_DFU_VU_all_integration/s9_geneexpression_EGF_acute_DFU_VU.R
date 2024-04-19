library(tidyverse)
library(ggplot2)
library(patchwork)
library(ggpubr)

library(Seurat)
library(scCustomize)
library(scRNAtoolVis)

# load the data

inteData <- readRDS("00-InteAllMainCellType/allacuteWound_DFU_VU_integrated_reduced.rds")

egf <- c("TGFA", "AREG", "HBEGF", "EREG", "EGFR", "ERBB2")

#my_comparisons <- list( c("Skin", "Wound1"), c("Skin", "Wound7"), c("Skin", "Wound30"))
#for (i in seq_along(egf)) {
#  print((VlnPlot_scCustom(seurat_object = inteData,
#                          features = egf[i], 
#                          group.by = "Condition",
#                          plot_median = TRUE,y.max = 8) & NoLegend()) +
#          stat_compare_means(comparisons = my_comparisons))
#}

# AverageExpression OR AggregateExpression
gene_exp <- AggregateExpression(inteData, features = egf, 
                                group.by = c("Project", "Condition"),
                                assays = "RNA", normalization.method = "LogNormalize", scale.factor = 10000)$RNA
gene_exp_f <- gene_exp %>% as.data.frame() %>% 
  mutate(Wound1_FC = AcuteWound_Wound1/AcuteWound_Skin, 
         Wound7_FC= AcuteWound_Wound7/AcuteWound_Skin, 
         Wound30_FC=AcuteWound_Wound30/AcuteWound_Skin,
         DFU_H_FC = DFU_DFU_H/DFU_H, DFU_NH_FC= DFU_DFU_NH/DFU_H,
         VU_FC = VU_VU/VU_NS) %>% select(10:15) %>% 
  rownames_to_column(var = "ID") %>% 
  pivot_longer(!ID, names_to = "Group", values_to = "FoldChange") %>% 
  mutate(Condition = rep(c("Wound1", "Wound7", "Wound30", "DFU_H", "DFU_NH", "VU") , length(egf)))


gene_exp_f$Condition <- factor(gene_exp_f$Condition, levels = c("Wound1", "Wound7", "Wound30", "DFU_H", "DFU_NH", "VU"))
gene_exp_f$ID <- factor(gene_exp_f$ID, levels = rev(egf))

# Do the scale or not
## Setting the color to blue while values are lower than 1 (Down-regulated)
## red colors while values are higher than 1 (Up-regulated) compared to normal skin
ggplot(data = gene_exp_f, aes(x=ID, y = Condition, 
                              color = FoldChange)) + 
  geom_point(size=3) +
  scale_color_gradient2(limits = c(0,2), low = "blue", mid = "white", high = "red", midpoint = 1, na.value = "red") +
  #scale_size(range = c(0,8), limits = c(-2,2))+
  theme_bw() + 
  ylab("") + 
  xlab("") +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1)) +
  coord_flip()

#pdf("Allinte_EGF_signaling_DFU_VU_otherLigand.pdf", useDingbats = F, width = 3, height = 3)
#dev.off()



#####################################################################################
# Try to re-annotate the integrated main clusters (Integrated acute wound, DFU, VU) #
# Since we filter out the HFs, Mast cells, Melanocytes in downstream analysis
# thus, we first annotated these clusters based on the res1.2, 
# then, we annotated the clusters based on subcluster annotations
#####################################################################################
# load the annotated metadata (newCellType, mainCellType)
inteData <- readRDS("00-InteAllMainCellType/allacuteWound_DFU_VU_integrated_reduced.rds")

# remove the main cell type as Other
inteData <- subset(inteData, subset = mainCellType != "Other")
inteData@meta.data <- droplevels(inteData@meta.data)

Idents(inteData) <- inteData$mainCellType
Cluster_Highlight_Plot(seurat_object = inteData, 
                       cluster_name = c("HF", "Bas-mig", "Spi-mig"), 
                       highlight_color = c("skyblue", "pink","forestgreen"),
                       background_color = "lightgray")

# set the factor for cell type and plot
fact_lev <- c("Bas-I", "Bas-prolif", "Bas-mig", 
              "Spi-I", "Spi-II", "Spi-mig", 
              "Gra", "HF", 
              "MEL", 
              "FB-I", "FB-II", "FB-III", "FB-prolif", 
              "PC_vSMC", "LE", "VE", 
              "NK", "Tcell", "Plasma_Bcell", "Mast-cell", 
              "Mac_inf", "Mac1", "Mac2", "Mac3",
              "cDC1", "cDC2", "DC3", "LC"#, "Other"
              )

ct.cols <- c("#d94701", "#fd8d3c", "#fdbe85", #Basal clusters
             "#33A02C", "#72BF5A", "#B2DF8A", #Spinous clusters
             "#f768a1", "#d4b9da", #Granular, Hair follicle
             "#737373", #MEL
             "#0570b0", "#3690c0", "#92c5de", "#d1e5f0", #Fibroblast clusters
             "#1a9850", "#fb9a99", "#8d4720",  #PCvSMC,LE,VE
             "#35978f", "#41b6c4", "#80cdc1","#df65b0", #"NK-cell", "Th", "Plasma_Bcell", "Mast-cell", 
             "#dd3497", "#CF384D", "#FA9C58", "#93d741", # Mono-Mac"
             "#807dba","#6a3d9a","#9e9ac8", "#b15928"#, #"cDC1", "cDC2", "DC3", "LC"
             #"grey90" # Other
)
# add the color names using new cell types
names(ct.cols) <- fact_lev

inteData$mainCellType <- factor(inteData$mainCellType, levels = fact_lev)
inteData$Condition %>% table()

DimPlot(inteData, group.by = "mainCellType", label = T, cols = ct.cols) & NoLegend()
DimPlot(inteData, group.by = "mainCellType", label = T, cols = ct.cols, split.by = "Project") & NoLegend()

# EGF signaling genes
egf <- c("TGFA", "AREG", "HBEGF", "EREG", "EGFR", "ERBB2")
#hgf <- c("HGF", "MET")
#cxcl1 <- c("CXCL1", "ACKR1", "CXCR1", "CXCR2")

Idents(inteData) <- inteData$Project

p1 <- Stacked_VlnPlot(seurat_object = subset(inteData, idents = "AcuteWound"), 
                features = egf, 
                x_lab_rotate = TRUE,pt.size = 0,
                colors_use = c("#fdc086","#FC4E07", "#4daf4a", "#377eb8", "#fdc086",'#a1d76a','#e9a3c9', "#fdc086","#ff7f00"), 
                group.by = "mainCellType", 
                split.by = "Condition", plot_legend = TRUE)
p2 <- Stacked_VlnPlot(seurat_object = subset(inteData, idents = "DFU"), 
                features = egf, 
                x_lab_rotate = TRUE,pt.size = 0,
                colors_use = c("#fdc086","#FC4E07", "#4daf4a", "#377eb8", "#fdc086",'#a1d76a','#e9a3c9', "#fdc086","#ff7f00"), 
                group.by = "mainCellType", 
                split.by = "Condition", plot_legend = TRUE)
p3 <- Stacked_VlnPlot(seurat_object = subset(inteData, idents = "VU"), 
                features = egf, 
                x_lab_rotate = TRUE,pt.size = 0,
                colors_use = c("#fdc086","#FC4E07", "#4daf4a", "#377eb8", "#fdc086",'#a1d76a','#e9a3c9', "#fdc086","#ff7f00"), 
                group.by = "mainCellType", 
                split.by = "Condition", plot_legend = TRUE)

wrap_plots(list(p1, p2, p3), ncol = 1)
#pdf("Allinte_EGF_signaling_DFU_VU_otherLigand_exp.pdf", useDingbats = F, width = 16, height = 16)
#dev.off()

p1 <- jjDotPlot(object = inteData, 
                gene = egf,
                xtree = FALSE, ytree = FALSE,
                id="Condition",
                cluster.order = rev(levels(inteData$Condition)),
                rescale = T,
                rescale.min = 0,
                rescale.max = 1) + 
  scale_size(breaks = c(0, 25, 50, 75, 100), range = c(1,4)) 

p2 <- jjDotPlot(object = inteData, 
                gene = egf,
                xtree = FALSE, ytree = FALSE,
                id="mainCellType",
                cluster.order = rev(levels(inteData$mainCellType)),
                rescale = T,
                rescale.min = 0,
                rescale.max = 1) + 
  scale_size(breaks = c(0, 25, 50, 75, 100), range = c(1,4))

#pdf("Allinte_CellType_gene_expression_egf.pdf", useDingbats = F, width = 14, height = 10)
p1/p2
#dev.off()

#pdf("Allinte_CellType_gene_expression_egf_splitCondition.pdf", useDingbats = F, width = 14, height = 40)
jjDotPlot(object = inteData, 
          gene = egf,
          xtree = FALSE, ytree = FALSE,
          id="mainCellType",
          split.by = "Condition",
          cluster.order = unlist(lapply(rev(levels(inteData$mainCellType)), FUN = function(x) paste0(x, " (",rev(levels(inteData$Condition)), ")"))),
          split.by.aesGroup = T,
          rescale = T,
          rescale.min = 0,
          rescale.max = 1) + 
  scale_size(breaks = c(0, 25, 50, 75, 100), range = c(1,4))
#dev.off()



#####################################################################
# collagen, MMP, TIMP expression in acute wounds and chronic wounds #
#####################################################################
inteData
collagens <- rownames(inteData)[grep("^COL", rownames(inteData))]
mmps <- rownames(inteData)[grep("^MMP", rownames(inteData))]
timp <- rownames(inteData)[grep("^TIMP", rownames(inteData))]
bmps <- rownames(inteData)[grep("^BMP", rownames(inteData))]
serps <- rownames(inteData)[grep("^SERP", rownames(inteData))]
cxcls <- rownames(inteData)[grep("^CXC", rownames(inteData))]
ccls <- rownames(inteData)[grep("^CCL", rownames(inteData))]

Idents(inteData) <- inteData$mainCellType
jjDotPlot(object = inteData, 
          gene = c(mmps, timp, bmps),
          xtree = FALSE, 
          ytree = FALSE,
          id="Condition",
          cluster.order = rev(levels(inteData$Condition)),
          rescale = T,
          rescale.min = 0,
          rescale.max = 1)

jjDotPlot(object = subset(inteData, idents = c("FB-I", "FB-II", "FB-III", "FB-prolif")), 
          gene = c(mmps, timp, bmps),
          xtree = FALSE, 
          ytree = FALSE,
          id="Condition",
          cluster.order = rev(levels(inteData$Condition)),
          rescale = T,
          rescale.min = 0,
          rescale.max = 1) #+ 
  #scale_size(breaks = c(0, 25, 50, 75, 100), range = c(1,4)) 

#pdf("Allinte_CellType_gene_expression_egf_splitCondition.pdf", useDingbats = F, width = 14, height = 40)
jjDotPlot(object = subset(inteData, idents = c("FB-I")), #"FB-II", "FB-III", "FB-prolif")), 
          gene = c(mmps, timp, bmps),
          xtree = FALSE, ytree = FALSE,
          id="mainCellType",
          split.by = "Condition",
          cluster.order = unlist(lapply(rev(levels(inteData$mainCellType)), FUN = function(x) paste0(x, " (",rev(levels(inteData$Condition)), ")"))),
          split.by.aesGroup = T,
          rescale = T,
          rescale.min = 0,
          rescale.max = 1) + 
  scale_size(breaks = c(0, 25, 50, 75, 100), range = c(1,4))
#dev.off()

acutewound <- subset(inteData, subset = Project == "AcuteWound")
#pdf("Allinte_CellType_gene_expression_egf_splitCondition.pdf", useDingbats = F, width = 14, height = 40)
jjDotPlot(object = subset(acutewound, idents = c("FB-I")), #"FB-II", "FB-III", "FB-prolif")), 
          gene = c(mmps, timp, bmps),
          xtree = FALSE, ytree = FALSE,
          id="mainCellType",
          split.by = "Condition",
          cluster.order = unlist(lapply(rev(levels(acutewound$mainCellType)), FUN = function(x) paste0(x, " (",rev(levels(inteData$Condition)), ")"))),
          split.by.aesGroup = T,
          rescale = T,
          rescale.min = 0,
          rescale.max = 1) + 
  scale_size(breaks = c(0, 25, 50, 75, 100), range = c(1,4))
#dev.off()