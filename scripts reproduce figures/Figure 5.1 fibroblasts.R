library(Seurat)
library(scCustomize)
library(tidyverse)
library(patchwork)
library(RColorBrewer)
library(viridis)
library(ComplexHeatmap)
library(circlize)

setwd("./")

rm(list = ls())

fib_sub <- readRDS("s1_clean_fibroblast.rds")

# re-annotate the cell clusters
fac_levs <- c("FB-I(POSTN+COL11A1+)", "FB-I(POSTN+MMP11+)", 
              "FB-I(POSTN+COL4A1+)", "FB-I(SFRP4+COMP+)", 
              "FB-II(APOD+ITM2A+)", "FB-II(APOE+CCL19+)",
              "FB-III(ELN+LEPR+)", 
              "FB-prolif",
              "FB(SFRP1+CRABP1+)", "FB(ELN+SFRP4+)")
current_cluster_ids <- c(0,2,8,4,
                         3,5,1,6,9,7) # List of current cluster IDs
tmp <- plyr::mapvalues(x=fib_sub$SCT_snn_res.0.5, from=current_cluster_ids, to=fac_levs)
fib_sub$upCellTypes <- tmp
fib_sub$upCellTypes <- factor(fib_sub$upCellTypes, levels = fac_levs)

# set the colors for the fibroblasts
ct.cols <- c("#57a9e2", "#31c783", "#c6cc85", "#c072fd", "#e36767", 
               "#aa6d60", "#ff9a41", "#eca1d5", "#31d7e8", "#b0b0b0")
names(ct.cols) <- fac_levs

#pdf("Fig5_UMAP_fb_PI16.pdf", useDingbats = F, width = 5, height = 4)
FeaturePlot(fib_sub, features = "PI16")
dev.off()


################################################
####---- Fig 5 UMAP plot of fibroblasts ----####
################################################
p.umap.fb <- DimPlot(fib_sub, group.by = "upCellTypes", label = F, 
                     cols = ct.cols, split.by = "Condition", ncol = 1) +
  NoAxes() + ggtitle("")
#pdf("Fig5_UMAP_fb_split.pdf", useDingbats = F, width = 14, height = 4)
p.umap.fb
#dev.off()

DimPlot(fib_sub, label = F, group.by = "upCellTypes", cols = ct.cols)

###############################################
####---- Fig 5 dotplot of marker genes ----####
###############################################
# marker genes for subclusters of fibroblasts
genes <- c("POSTN", "ADAM12", "ASPN", "COL16A1", "COL11A1", 
           "MMP11", "ISG15", "IER5", "IFI6",
           "COL4A1", "COL4A2", 
           "IGFBP6", "COMP", "EGFL6", 
           "APOD", "FGF7", "PLA2G2A", "ITM2A", 
           "C3", "APOE", "CCL19", "CD74", "COL6A5", 
           "LEPR", "PLPP3", "CD9", "SLIT3", "DPP4", 
           "PCLAF", "MKI67",
           "SFRP1", "CRABP1", "TNN", "COCH", "G0S2",
           "ELN", "AKR1C1", "GPRC5A", "SFRP4")
fib_sub$upCellTypes <- factor(fib_sub$upCellTypes, levels = rev(fac_levs))
plot_marker <- DotPlot(fib_sub, features = genes, 
                       group.by = "upCellTypes", cols = c("white", "#cb181d"), 
                       dot.scale = 4, col.min = 0, dot.min = 0.1
) + theme(axis.text.x = element_text(angle = 90, hjust = 1),
          panel.border = element_rect(colour = "black", fill=NA)) + labs(x="", y="")
plot_marker

pdf("Fig5_Dotplot_markers.pdf", useDingbats = F, width = 10, height = 4)
plot_marker
dev.off()


################################################
####---- Fig 5 cell proportion analysis ----####
################################################
# step 1. 12 individuals divided into cell types and their numbers
step1 <- table(fib_sub$orig.ident, fib_sub$upCellTypes) %>% as.data.frame()
step1$sepa <- "Dermis"

# step 2. read the cell numbers of dermis and epidermis
step2 <- data.table::fread("all_epidermis_dermis_perIndi_220613.txt")

# step 3. calculate the proportion of each cell type for each individual
step3 <- step1 %>% left_join(., step2, by=c("sepa"="sepa", "Var1"="Var1")) %>% distinct() %>% 
  mutate(Prop=Freq/a_sum)

# step 4. calculate the total normalized proportions of each cell type per condition
df.group <- step3 %>% mutate(Sample = gsub("PWH..", "", Var1)) %>% group_by(Sample, Var2) %>% summarise(Freq=sum(Prop)) %>% 
  ungroup() %>% group_by(Sample) %>% mutate(Freq_new = Freq/sum(Freq), lbl = scales::percent(Freq_new)) %>% ungroup() %>% 
  rename("Cluster" = "Var2")

df.group$Sample <- gsub("D0", "Skin", df.group$Sample)
df.group$Sample <- gsub("D1", "Wound1", df.group$Sample)
df.group$Sample <- gsub("D7", "Wound7", df.group$Sample)
df.group$Sample <- gsub("D30", "Wound30", df.group$Sample)
df.group$Sample <- factor(df.group$Sample, levels = c("Skin", "Wound1", "Wound7", "Wound30"))
df.group$Cluster <- factor(df.group$Cluster, levels = fac_levs)

p_cp_5 <- ggplot(df.group, aes(x = Sample, y = Freq_new, fill = Cluster)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = ct.cols) + 
  xlab('') +
  scale_y_continuous(breaks = seq(0, 1, .2), 
                     expand = c(0, 0.01),
                     labels = c("0%", "20%", "40%", "60%", "80%", "100%"),
                     name = 'Percentage') +
  #geom_text(aes(label = lbl), 
  #          size = 4, 
  #          position = position_stack(vjust = 0.5)) +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    axis.text.x = element_text(size = 14, color = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.y = element_text(size = 16, color = "black"),
    legend.position = "right",
    legend.text = element_text(size = 12, color = "black")
  ) 
p_cp_5

#pdf("Fig5_cellproportion_fb.pdf", useDingbats = F, width = 5, height = 4)
#p_cp_5
#dev.off()


#########################################################
# run the statistics of cell proportions compared to Skin
df.group.out <- step3 %>% mutate(Group = gsub("PWH..", "", Var1)) %>% select(-4)
data.table::fwrite(df.group.out, file = "AcuteWound_normCellProportion_FB.txt", sep = "\t")

# summary the MKI67 positive cell numbers
genes_all <- FetchData(fib_sub, vars = c("Condition", "orig.ident", "upCellTypes",
                                         "MKI67"), slot = "data")
genes_all_f <- genes_all %>% filter(MKI67 > 0)
table(genes_all_f$upCellTypes, genes_all_f$Condition)

sm.tol <- step3 %>% select(1,5) %>% distinct()
mki67 <- table(genes_all_f$orig.ident) %>% as.data.frame() %>% left_join(., sm.tol, by=c("Var1" = "Var1")) %>% 
  mutate(Prop=Freq/a_sum, CellType="MKI67positive", Sample = Var1) %>% 
  mutate(Sample = gsub("^PWH[0-9]{2}", "", Sample)) %>%
  dplyr::select(1,5,everything())
colnames(mki67) <- colnames(df.group.out)
com <- rbind(df.group.out, mki67)
data.table::fwrite(com, file = "AcuteWound_normCellProportion_FB_mki67positive.txt", sep = "\t")
df.group.out <- com

# make a function to run all the test using the
# quasibinomial and normal Wilcoxon test
ct <- unique(df.group.out$Var2) %>% as.character();ct

stat_sig <- list()
for (i in seq_along(ct)){
  ct_sel <- ct[i]
  ct_por <- df.group.out %>% dplyr::filter(Var2 == ct_sel) %>% arrange(Group)
  
  stats_all <- function(GroupInfo = NULL){
    ct_por_gr <- ct_por %>% dplyr::filter(Group %in% GroupInfo)
    # check the normality of the variables
    print(paste0("Shapiro.test on ", GroupInfo[1]))
    ct_por %>% dplyr::filter(Group %in% GroupInfo[1]) %>% pull(Prop) %>% shapiro.test() %>% print()
    print(paste0("Shapiro.test on ", GroupInfo[2]))
    ct_por %>% dplyr::filter(Group %in% GroupInfo[2]) %>% pull(Prop) %>% shapiro.test() %>% print()
    
    # paired samples Wilcoxon test
    if(TRUE){
      res <- t.test(Prop ~ Group, data = ct_por_gr)
    }else{
      res <- t.test(Prop ~ Group, data = ct_por_gr)
    }
    print(res)
    wilcoxon_paired = res$p.value
    # quasibinomial test
    test.quasi = glm(formula = Prop ~ Group, data = ct_por_gr, family=quasibinomial)
    print(summary(test.quasi))
    quasibinomial=anova(test.quasi, test = "LRT")$`Pr(>Chi)`[2]
    return(c(wilcoxon_paired, quasibinomial))
  }
  # Test of D1 versus D0
  res_all_1=stats_all(GroupInfo = c("D0", "D1"))
  res_all_2=stats_all(GroupInfo = c("D0", "D7"))
  res_all_3=stats_all(GroupInfo = c("D0", "D30"))
  stat_sig[[ct_sel]] <- append(stat_sig[[ct_sel]], res_all_1)
  stat_sig[[ct_sel]] <- append(stat_sig[[ct_sel]], res_all_2)
  stat_sig[[ct_sel]] <- append(stat_sig[[ct_sel]], res_all_3)
}
# combine all the data
stat_sig_all <- do.call("rbind", stat_sig) %>% as.data.frame() %>% 
  rownames_to_column(var = "CellType") %>% 
  setNames(c("CellType", "WilcoPair_D1D0", "QuasiBino_D1D0", "WilcoPair_D7D0", "QuasiBino_D7D0", "WilcoPair_D30D0", "QuasiBino_D30D0"))
data.table::fwrite(stat_sig_all, "Fig_AcuteWound_normCellProportion_FB_addSig.txt", sep = "\t")



################################################
##---- Figure 5 dendrogram all cell types ----##
################################################
DefaultAssay(fib_sub) <- "SCT"
fib_sub <- SetIdent(fib_sub, value = fib_sub@meta.data$upCellTypes)

fib_sub <- BuildClusterTree(object = fib_sub, assay = "SCT", slot = "scale.data")
# pull the tree
hswound.clustertree <- Tool(object = fib_sub, slot = "BuildClusterTree")
# plot the tree
suppressMessages(require(ape))
{plot.phylo(x = hswound.clustertree, direction = "rightwards", use.edge.length = TRUE, font = 4, tip.color = ct.cols, edge.width = 2, label.offset = 5, no.margin = TRUE)
  tiplabels(pch = 19, col = ct.cols, adj = 3, cex = 2)}

#pdf("Fig5_dendrogram.pdf", useDingbats = F, width = 8, height = 4)
#dev.off()


########################################
##---- Figure 5 RNA velocity plot ----##
########################################
velocity <- data.table::fread("FB_RNA_velocity_initial.txt", header = T)
identical(colnames(fib_sub), velocity$barcode)
colnames(velocity)[2] <- 'Initial_states'
velocity2 <- data.table::fread("FB_RNA_velocity_terminal.txt", header = T)
identical(colnames(fib_sub), velocity2$barcode)
colnames(velocity2)[2] <- 'Terminal_states'

fib_sub$Initial_states <- velocity$Initial_states
fib_sub$Terminal_states <- velocity2$Terminal_states
p1 <- FeaturePlot(fib_sub, features = c("Terminal_states"),
            min.cutoff = 0, cols = c("#400080","#00d100","yellow"), max.cutoff = "q97", pt.size = 0.2) + NoAxes()
p2 <- FeaturePlot(fib_sub, features = c("Initial_states"),
                  min.cutoff = 0, cols = c("#400080","#00d100","yellow"), max.cutoff = "q97", pt.size = 0.2) + NoAxes()

#pdf("Fig5_states2.pdf", useDingbats = F, width = 5, height = 10)
p1 + p2
#dev.off()

FeaturePlot(fib_sub, features = c("Terminal_states"))
FeaturePlot(fib_sub, features = c("Initial_states"))


#################################################################
##---- Figure 5 marker genes overlap with previous studies ----##
#################################################################
# Previous study divide fibroblasts into four clusters:
# reticular (WISP2, SLPI, MFAP5, TSPAN8)
# papillary (APCDD1, ID1, WIF1, COL18A1, PTDGS)
# mesenchymal (ASPN, POSTN, GPC3, TNN, SFRP1)
# proinflammatory (CCL19, APOE, CXCL2, CXCL3, EFEMP1)
# proliferating (MKI67)
canonical_markers <- c("ASPN", "APOE", "LEPR", "PCLAF",
                       "POSTN", "CCL2", "COL13A1", "H1-5",
                       "COL11A1", "C3", "F13A1", "MKI67")

plot_list <- FeaturePlot(
  fib_sub,
  features=unlist(canonical_markers),
  combine=FALSE, cols = c("gray","red")
)

# apply theme to each feature plot
for(i in 1:length(plot_list)){
  plot_list[[i]] <- plot_list[[i]] + NoLegend() + NoAxes()
}

wrap_plots(plot_list, ncol = 4)


FeaturePlot(
  fib_sub,
  features="ACTB",
  combine=TRUE, cols = c("gray","red"), split.by = "Condition"
)


# re-annotate the cell clusters with only four clusters
fac_levs <- c("FB-I(mesenc)", "FB-I(mesenc)", 
              "FB-I(mesenc)", "FB-I(mesenc)", 
              "FB-II(proinf)", "FB-II(proinf)",
              "FB-III(papill)", 
              "FB-prolif",
              "FB-I(mesenc)", "FB-I(mesenc)")
current_cluster_ids <- c(0,2,8,4,
                         3,5,1,6,9,7) # List of current cluster IDs
tmp <- plyr::mapvalues(x=fib_sub$SCT_snn_res.0.5, from=current_cluster_ids, to=fac_levs)
fib_sub$comCellTypes <- tmp
fib_sub$comCellTypes <- factor(fib_sub$comCellTypes, levels = c("FB-I(mesenc)", 
                                                                "FB-II(proinf)", 
                                                                "FB-III(papill)", 
                                                                "FB-prolif"))

# set the colors for the fibroblasts
ct.cols <- c("#57a9e2", "#e36767", 
               "#ff9a41", "#eca1d5")
names(ct.cols) <- c("FB-I(mesenc)", 
                    "FB-II(proinf)", 
                    "FB-III(papill)", 
                    "FB-prolif")


#########################################################
####---- Fig 5 UMAP plot of combined fibroblasts ----####
#########################################################
p.umap.fb<- DimPlot(fib_sub, group.by = "comCellTypes", label = T, cols = ct.cols) +
  NoAxes() + ggtitle("")
p.umap.fb

#pdf("Fig5_UMAP_fb_comb.pdf", useDingbats = F, width = 5, height = 4)
p.umap.fb
#dev.off()


# markers
markers <- data.table::fread("s1_clean_fibroblast_mg_res05.txt")
top200markers <- markers %>% group_by(cluster) %>% top_n(n = 200, wt = avg_log2FC) %>% arrange(cluster)
data.table::fwrite(top200markers, "Markers_top200_subFB.txt", sep = "\t")

