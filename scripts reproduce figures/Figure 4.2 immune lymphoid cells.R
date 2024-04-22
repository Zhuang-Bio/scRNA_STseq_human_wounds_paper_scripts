library(Seurat)
library(scCustomize)
library(SeuratObject)
library(tidyverse)
library(magrittr)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
#library(viridis)
#library(ComplexHeatmap)
#library(circlize)

setwd("./")

rm(list = ls())

lym_sub <- readRDS( "s1_clean_immu_lymphoid_newUMAP.rds")

# re-annotate the cell clusters
fac_levs <- c("Treg", "Th", "ILCs", "Tc", 
              "ILC1/NK", "NK", "Ttol", 
              "Plasma", "Bcell")

current_cluster_ids <- c(5, 0, 1, 7, 
                         3, 2, 4,
                         6, 8) # List of current cluster IDs
tmp <- plyr::mapvalues(x=lym_sub$SCT_snn_res.0.5, from=current_cluster_ids, to=fac_levs)

lym_sub$upCellTypes <- tmp
lym_sub$upCellTypes <- factor(lym_sub$upCellTypes, levels = fac_levs)
Idents(lym_sub) <- lym_sub$upCellTypes

# use new colors
ct.cols <- c('#807dba','#9970ab','#810f7c', '#41ae76', 
             '#74add1', '#8dd3c7','#76d72f',
             '#f768a1', '#fdb462'
)

names(ct.cols) <- fac_levs

#############################################
####---- Fig 7 UMAP plot of myeloids ----####
#############################################
p.umap <- DimPlot(object = lym_sub, 
                  reduction = "umap", 
                  label = F, group.by = "upCellTypes", cols = ct.cols) +
  NoAxes() + ggtitle("")
p.umap

#pdf("Fig7_UMAP_lymphoid.pdf", useDingbats = F, width = 5, height = 4)
p.umap
#dev.off()

#pdf("Fig7_UMAP_split_lymphoid.pdf", useDingbats = F, width = 12, height = 4)
DimPlot(object = lym_sub, 
        reduction = "umap", 
        label = F, group.by = "upCellTypes", cols = ct.cols, split.by = "Condition", ncol = 1) +
  NoAxes() + ggtitle("")
#dev.off()


###############################################
####---- Fig 7 dotplot of marker genes ----####
###############################################
# choose the top 5-8 representative markers
top_repre_markers <- c("CD3D", "CD3G", "KLRB1", "KLRD1", "CD4", "CD40LG", "CD8A", "CD8B")
DotPlot(lym_sub, features = top_repre_markers,
        group.by = "upCellTypes", cols = c("white", "#cb181d"), 
        dot.scale = 4, scale = T, dot.min = 0.05 #col.min = 0.1, col.max=2.5 #col.min = -2.5
) + theme(axis.text.x = element_text(angle = 90, hjust = 1),
          panel.border = element_rect(colour = "black", fill=NA)) + labs(x="", y="")
#pdf("Fig_Dotplot_markers_lymphoid2.pdf", useDingbats = F, width = 4.5, height = 3.5)
#dev.off()


top_repre_markers <- c("CD3D", "CD3G", "KLRB1", "KLRD1", "CD4", "CD40LG", "CD8A", "CD8B",
                       "TIGIT", "BATF", "FOXP3", "CTLA4", "CORO1B", #cluster 5: Treg cells
                       "LDHB", "IL7R", "AP3M2", "GIMAP7", "KLF2", #cluster 0: Th cells
                       "AHR", "CCR6", "PTGER4", "ANKRD28", "LPAR6", #cluster 1: ILC cells
                       "FXYD2", "TRGC2", "KLRC3", "KLRC2","PDE4A", #cluster 7: Tc cells
                       "XCL1", "GNLY","XCL2", "FCER1G", #cluster 3: ILC1/NK cells
                       "NKG7", "GZMA","CRTAM", "GZMK", "TNFRSF9", #cluster 2: NK cells
                       "HSPA1B", "DNAJB1", "JUN", "FOS", "NR4A1", #cluster 4: Ttor cells
                       "PTGDS", "JCHAIN", "IL3RA", "CCR7", "CXCL8", #cluster 6: plasma cells
                       "IGKC", "MS4A1", "CD79A", "BANK1", "IGHM" #cluster 8: B cells
)

lym_sub$upCellTypes <- factor(lym_sub$upCellTypes, levels = rev(fac_levs))

plot_marker <- DotPlot(lym_sub, features = top_repre_markers, 
                       group.by = "upCellTypes", cols = c("white", "#cb181d"), 
                       dot.scale = 4,col.min = 0, dot.min = 0.05
) + theme(axis.text.x = element_text(angle = 90, hjust = 1),
          panel.border = element_rect(colour = "black", fill=NA)) + labs(x="", y="")
pdf("Fig_Dotplot_markers_lymphoid.pdf", useDingbats = F, width = 11, height = 3.5)
plot_marker
dev.off()

top_repre_markers <- c("TIGIT", "BATF", "FOXP3", #cluster 5: Treg cells
                       "LDHB", "GIMAP7", "KLF2", #cluster 0: Th cells
                       "AHR", "CCR6", "PTGER4", #cluster 1: ILC cells
                       "TRGC2", "KLRC3", "KLRC2", #cluster 7: Tc cells
                       "XCL1","XCL2", "FCER1G", #cluster 3: ILC1/NK cells
                       "NKG7", "GZMA", "GZMK", #cluster 2: NK cells
                       "HSPA1B", "DNAJB1", "NR4A1", #cluster 4: Ttor cells
                       "PTGDS", "JCHAIN", "IL3RA", #cluster 6: plasma cells
                       "MS4A1", "CD79A", "IGHM" #cluster 8: B cells
)
p2 <- jjDotPlot(object = lym_sub, 
                gene = top_repre_markers,
                xtree = FALSE, ytree = FALSE,
                id="upCellType",
                cluster.order = rev(fac_levs),
                rescale = T,
                rescale.min = 0,
                rescale.max = 1) + scale_size(breaks = c(0, 25, 50, 75, 100), range = c(1,4))

pdf("Lymphoid_marker_gene_expression_reduce.pdf", useDingbats = F, width = 8, height = 12)
p2
dev.off()

# markers
markers <- data.table::fread("s1_clean_immu_lymphoid_mg_res0.5.txt")
top200markers <- markers %>% group_by(cluster) %>% top_n(n = 200, wt = avg_log2FC) %>% arrange(cluster)
data.table::fwrite(top200markers, "Markers_top200_subLymphoid.txt", sep = "\t")


################################################
####---- Fig 7 cell proportion analysis ----####
################################################
# step 1. 12 individuals divided into cell types and their numbers
lym_sub$upCellTypes <- factor(lym_sub$upCellTypes, levels = fac_levs)
step1 <- table(lym_sub$orig.ident, lym_sub$upCellTypes) %>% as.data.frame()
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

#pdf("Fig7_cellproportion_lymphoid.pdf", useDingbats = F, width = 5, height = 5)
p_cp_5
#dev.off()


#########################################################
# run the statistics of cell proportions compared to Skin
df.group.out <- step3 %>% mutate(Group = gsub("PWH..", "", Var1)) 
data.table::fwrite(df.group.out, file = "AcuteWound_normCellProportion_Lymphoid.txt", sep="\t")
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
data.table::fwrite(stat_sig_all, "Figure 5-2_AcuteWound_normCellProportion_Lymphoid_addSig.txt", sep = "\t")

