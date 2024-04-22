library(Seurat)
library(SeuratObject)
library(tidyverse)
library(magrittr)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(viridis)
library(ComplexHeatmap)
library(circlize)

setwd(".,")

rm(list = ls())

mye_sub <- readRDS("s1_clean_immu_myeloid.rds")

# re-annotate the cell clusters
fac_levs <- c("Mac_inf", "Mac1","Mac2", "Mac3", 
              "pDC", "cDC1", "cDC2", "DC3", 
              "LC", "Apoptotic", "Cycling")

current_cluster_ids <- c(2, 3, 1, 5,
                         9, 4, 0, 6, 
                         8, 7, 10) # List of current cluster IDs
tmp <- plyr::mapvalues(x=mye_sub$SCT_snn_res.0.5, from=current_cluster_ids, to=fac_levs)

mye_sub$upCellTypes <- tmp
mye_sub$upCellTypes <- factor(mye_sub$upCellTypes, levels = fac_levs)
Idents(mye_sub) <- mye_sub$upCellTypes
#saveRDS(mye_sub, file = "myeloid_newAnnotation.rds")

# use new colors
ct.cols <- c('#f1b6da','#df65b0', '#41ab5d', '#addd8e',
             "#810f7c", "#807dba", "#9970ab", "#b2abd2", 
             "#bf812d", "#fee0b6", "#f1a340"
)
names(ct.cols) <- fac_levs


library(scRNAtoolVis)
library(scCustomize)
FeaturePlot_scCustom(seurat_object = mye_sub, 
                     features = c("APOE", "IL1B","DAB2", "MMP19"),
                     na_cutoff = 3,
                     #colors_use = c("grey90", "#F5AD32", "#1D078D"),
                     order = T) & NoAxes()
#pdf("Figure_mac_markers.pdf", useDingbats = F, width = 8, height = 8)
#dev.off()

#############################################
####---- Fig 6 UMAP plot of myeloids ----####
#############################################
p.umap <- DimPlot(object = mye_sub, 
                       reduction = "umap", 
                       label = F, group.by = "upCellTypes", cols = ct.cols) +
  NoAxes() + ggtitle("")
p.umap

#pdf("Fig_UMAP_myeloid.pdf", useDingbats = F, width = 5, height = 4)
p.umap
#dev.off()

#pdf("Fig_UMAP_split_myeloid.pdf", useDingbats = F, width = 12, height = 4)
DimPlot(object = mye_sub, 
        reduction = "umap", 
        label = F, group.by = "upCellTypes", cols = ct.cols, split.by = "Condition", ncol = 1) +
  NoAxes() + ggtitle("")
#dev.off()


# GOresults 
library(clusterProfiler)
library(org.Hs.eg.db)
markers <- data.table::fread("DEmarkers_immu_myeloid_mg_res0.5.txt")
markers$upCellTypes <- plyr::mapvalues(x=markers$cluster, from=current_cluster_ids, to=fac_levs)

markers_mac <- markers %>% filter(upCellTypes %in% c("Mac_inf", "Mac1", "Mac2", "Mac3")) %>% 
  group_by(upCellTypes) %>% top_n(200, avg_log2FC)
#data.table::fwrite(markers_mac, "goresults_inputmarkers.txt", sep = "\t")

markers_mac.list <- list(Mac_inf = markers_mac %>% filter(upCellTypes %in% "Mac_inf") %>% pull(entrezgene_id) %>% as.character(),
                         Mac1 = markers_mac %>% filter(upCellTypes %in% "Mac1") %>% pull(entrezgene_id) %>% as.character(),
                         Mac2 = markers_mac %>% filter(upCellTypes %in% "Mac2") %>% pull(entrezgene_id) %>% as.character(),
                         Mac3 = markers_mac %>% filter(upCellTypes %in% "Mac3") %>% pull(entrezgene_id) %>% as.character())

comSm.up.bp <- compareCluster(geneCluster = markers_mac.list, fun = enrichGO, 
                              OrgDb = org.Hs.eg.db, pvalueCutoff = 1, 
                              pAdjustMethod = "BH", qvalueCutoff  = 1,
                              minGSSize = 10,maxGSSize = 500,keyType = "ENTREZID")

comSm.up.bp <- setReadable(comSm.up.bp, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
outputdata <- comSm.up.bp %>% as.data.frame() %>% dplyr::filter(p.adjust < 0.05) %>% dplyr::arrange(Description, Cluster)

dotplot(comSm.up.bp, showCategory = 10, label_format = 100, size = "Count") + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12, color = "black", angle = 30, hjust = 1),
        axis.text.y = element_text(size = 12, color = "black"))


ego_BP_up <- enrichGO(gene = markers_mac.list$Mac_inf,
                        OrgDb         = org.Hs.eg.db,
                        ont           = "BP",
                        keyType = "ENTREZID",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.2,
                        minGSSize = 10,
                        maxGSSize = 500,
                        readable = T)
ego_BP_up <- as.data.frame(ego_BP_up)
#data.table::fwrite(ego_BP_up, file = paste0("GO_BPs_", "Mac_inf", ".txt"), sep = "\t")
dotplot(ego_BP_up, showCategory=10, label_format = 100) + 
    theme(
      panel.grid.minor = element_blank(),
      panel.border = element_rect(size = 0.5, color = "black"),
      axis.text.x = element_text(size = 14, color = "black"),
      axis.text.y = element_text(size = 14, color = "black"))
  

## GO results for final figure ##
go_inte <- readxl::read_xlsx("GOresults_ClusterProfiler_Web,GOresults_macrophages.xlsx", sheet = 10)
go_inte_f <- go_inte %>% select(2,5,9,10)

# reset the factor levels
go_inte_f$Group <- factor(go_inte_f$Group, levels = c("Mac_inf", "Mac1", "Mac2", "Mac3"))
GOs <- sort(unique(go_inte_f$description))[c(14, 7, 11, 23, 3, 18, 22, 5, 12, 16,
                                             2, 20, 15, 21, 
                                             19, 4, 13, 10, 17, 9,
                                             1, 8, 6)]
                                             
go_inte_f$description <- factor(go_inte_f$description, levels = rev(GOs))

p3 <- ggplot(data = go_inte_f, aes(x=Group, y = description, 
                            color = FDR, size = overlap)) + 
  geom_point() +
  scale_colour_gradient(limits = c(0,0.05),low = "red", high = "blue") +
  scale_size(range = c(1,6), limits = c(10,40))+
  theme_bw() + 
  ylab("") + 
  xlab("") +
  labs(color = "FDR",
       size = "Gene Count") +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1))
p3

pdf("Fig6_GOs_myeloid_mac.pdf", useDingbats = F, width = 5.5, height = 6)
p3
dev.off()

################################################
####---- Fig cell proportion analysis ----####
################################################
# step 1. 12 individuals divided into cell types and their numbers
step1 <- table(mye_sub$orig.ident, mye_sub$upCellTypes) %>% as.data.frame()
step1$sepa <- c(rep("Dermis", 12*8), rep("Epidemis", 12), rep("Dermis", 12*2))

# combine the mast cell into the myeloid cell type
step1 <- rbind(step1) #, mt_mast

# step 2. read the cell numbers of dermis and epidermis
step2 <- data.table::fread("all_epidermis_dermis_perIndi_220613.txt")

# step 3. calculate the proportion of each cell type for each individual
step3 <- step1 %>% left_join(., step2, by=c("sepa"="sepa", "Var1"="Var1")) %>% distinct() %>% 
  mutate(Prop=Freq/a_sum)

# step 4. calculate the total normalized proportions of each cell type per condition
df.group <- step3 %>% mutate(Sample = gsub("PWH..", "", Var1)) %>% group_by(Sample, Var2) %>% summarise(Freq=sum(Prop)) %>% 
  ungroup() %>% group_by(Sample) %>% mutate(Freq_new = Freq,sum(Freq), lbl = scales::percent(Freq_new)) %>% ungroup() %>% 
  dplyr::rename("Cluster" = "Var2")

df.group$Sample <- gsub("D0", "Skin", df.group$Sample)
df.group$Sample <- gsub("D1", "Wound1", df.group$Sample)
df.group$Sample <- gsub("D7", "Wound7", df.group$Sample)
df.group$Sample <- gsub("D30", "Wound30", df.group$Sample)
df.group$Sample <- factor(df.group$Sample, levels = c("Skin", "Wound1", "Wound7", "Wound30"))
fac_levs <- c(fac_levs) #
df.group$Cluster <- factor(df.group$Cluster, levels = c(fac_levs))
ct.cols=c(ct.cols, "grey");names(ct.cols)[12] <- "Mast-cell"


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

#pdf("Fig6_cellproportion_myeloid.pdf", useDingbats = F, width = 5, height = 8)
p_cp_5
#dev.off()

#########################################################
# run the statistics of cell proportions compared to Skin
df.group.out <- step3 %>% mutate(Group = gsub("PWH..", "", Var1)) 
data.table::fwrite(df.group.out, file = "AcuteWound_normCellProportion_Myeloid.txt", sep = "\t")
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
data.table::fwrite(stat_sig_all, "Figure_AcuteWound_normCellProportion_Myeloid_addSig.txt", sep = "\t")


###############################################
####---- Fig 2 dotplot of marker genes ----####
###############################################
top_repre_markers <- c('LYZ', 'HLA-DRA', "CD68","CD163", 
                       "APOE", "PHLDA1", "CXCL1", "CCL3", "CCL4",
                       "IL1B", "THBS1", "PTGS2", "EREG", "AQP9",
                       "DAB2", "C1QA", "C1QB", "MAF", "CCL13",
                       "MMP9", "MMP19", "SDS", "FGL2", "VEGFA",
                       "ACOT7", "LTB", "IFITM1", "TCF4","IGKC",#"CLEC4C", "NRP1",
                       "CLEC9A", "WDFY4", "CPNE3", "DNASE1L3", "CADM1",
                       "CD1C", "IL1R2", "CLEC10A", "CCR7",
                       "LAMP3", "CCL19", "FSCN1", "IL7R",
                       "FCGBP", "CD207", "CD1A", "CLDN1", 
                       "DNAJB1", "HSPA6", "HSPA1B", "HSPB1", "NR4A1",
                       "PCLAF", "H4C3", "H1-5", "HMGB2", "MKI67")

mye_sub$upCellTypes <- factor(mye_sub$upCellTypes, levels = rev(fac_levs))

plot_marker <- DotPlot(mye_sub, features = top_repre_markers, 
                       group.by = "upCellTypes", cols = c("white", "#cb181d"), 
                       dot.scale = 4, col.min = 0, dot.min = 0.05
) + theme(axis.text.x = element_text(angle = 90, hjust = 1),
          panel.border = element_rect(colour = "black", fill=NA)) + labs(x="", y="")
#pdf("Fig_Dotplot_markers_addPanmarker.pdf", useDingbats = F, width = 14, height = 4.5)
plot_marker
#dev.off()

library(scRNAtoolVis)
p1 <- jjDotPlot(object = mye_sub, 
                gene = top_repre_markers,
                xtree = FALSE, ytree = FALSE,
                id="Condition",
                cluster.order = rev(c("Skin","Wound1","Wound7","Wound30")),
                rescale = T,
                rescale.min = 0,
                rescale.max = 1) + scale_size(breaks = c(0, 25, 50, 75, 100), range = c(1,4))
p2 <- jjDotPlot(object = mye_sub, 
                gene = top_repre_markers,
                xtree = FALSE, ytree = FALSE,
                id="upCellType",
                cluster.order = rev(fac_levs),
                rescale = T,
                rescale.min = 0,
                rescale.max = 1) + scale_size(breaks = c(0, 25, 50, 75, 100), range = c(1,4))

pdf("Myeloid_marker_gene_expression.pdf", useDingbats = F, width = 14, height = 10)
p1/p2
dev.off()

# markers
markers <- data.table::fread("s1_clean_immu_myeloid_mg_res0.5.txt")
# re-annotate the cell clusters
fac_levs <- c("Mac_inf", "Mac1","Mac2", "Mac3", 
              "pDC", "cDC1", "cDC2", "DC3", 
              "LC", "Apoptotic", "Cycling")

current_cluster_ids <- c(2, 3, 1, 5,
                         9, 4, 0, 6, 
                         8, 7, 10) # List of current cluster IDs
tmp <- plyr::mapvalues(x=markers$cluster, from=current_cluster_ids, to=fac_levs)
markers$cluster <- tmp

top200markers <- markers %>% group_by(cluster) %>% top_n(n = 200, wt = avg_log2FC) %>% arrange(cluster)
data.table::fwrite(top200markers, "Markers_top200_subMyeloid.txt", sep = "\t")

