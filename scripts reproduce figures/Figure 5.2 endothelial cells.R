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

setwd("")

rm(list = ls())

other_sub <- readRDS("s1_clean_othercell.rds")

DimPlot(other_sub, group.by = "SCT_snn_res.0.5", label = T, label.size = 6)

# re-annotate the cell clusters
fac_levs <- c("SMC", "Pericytes", "LE",
              "VE_arteriole", "VE_capillary",
              "VE_venule1", "VE_venule2")

current_cluster_ids <- c(1, 3, 2, 
                         5, 6, 
                         0, 4) # List of current cluster IDs
tmp <- plyr::mapvalues(x=other_sub$SCT_snn_res.0.5, from=current_cluster_ids, to=fac_levs)

other_sub$upCellTypes <- tmp
other_sub$upCellTypes <- factor(other_sub$upCellTypes, levels = fac_levs)
Idents(other_sub) <- other_sub$upCellTypes

# use new colors
ct.cols <- c('#12a22d', '#ce8b1a', '#fe8Bac', 
             '#d93860', '#3ba4db', 
             '#807dba', '#d4b9da')

names(ct.cols) <- fac_levs

mt <- other_sub@meta.data %>% rownames_to_column(var = "barcode") %>% select(1, 30)
#data.table::fwrite(mt, "newOther_metadata.txt", sep = "\t")
#############################################
####---- Fig UMAP plot of  ----####
#############################################
p.umap <- DimPlot(object = other_sub, 
                  reduction = "umap", 
                  label = F, group.by = "upCellTypes", cols = ct.cols) +
  NoAxes() + ggtitle("")
p.umap

#pdf("Fig8_UMAP_other.pdf", useDingbats = F, width = 5, height = 4)
p.umap
#dev.off()

#pdf("Fig8_UMAP_split_other.pdf", useDingbats = F, width = 8, height = 6)
DimPlot(object = other_sub, 
        reduction = "umap", 
        label = F, group.by = "upCellTypes", cols = ct.cols, split.by = "Condition", ncol = 2) +
  NoAxes() + ggtitle("")
#dev.off()


###############################################
####---- Fig dotplot of marker genes ----####
###############################################
# choose the top 5-8 representative markers
Tpan <-  c("PECAM1", "VWF")
p3 <- FeaturePlot(
  other_sub, 
  order = T,
  features = Tpan,
  cols = c("gray","red"), combine = FALSE)
p2_list <- lapply(p3, function(x) {x + NoAxes() + theme(text = element_text(size = 10))})

#pdf("Fig8_other_panMarker.pdf", useDingbats = F, width = 12, height = 8)
cowplot::plot_grid(plotlist = p2_list, ncol = 1)
#dev.off()


top_repre_markers <- c("ACTA2", "TAGLN", "MYH11", "ADRA2A", "MUSTN1", # SMC
                       "COL6A3", "COL5A2", "THY1", "DCN", "GUCY1A2", # Pericytes
                       "CCL21", "TFF3", "LYVE1", "PDPN", # LE
                       "PECAM1","IGFBP3", "ARL15", "CXCL12", "HEY1", "SEMA3G", # VE_arteriole
                       "COL15A1", "VWA1", "RGCC", "H19", "APLN", # VE_capillary
                       "G0S2", "SELE", "ACKR1", "CSF3", "STC1", # VE_venule1
                       "CCL14", "AQP1", "ID1", "SOX18", "SNCG" # VE_venule2
)

other_sub$upCellTypes <- factor(other_sub$upCellTypes, levels = rev(fac_levs))

plot_marker <- DotPlot(other_sub, features = top_repre_markers, 
                       group.by = "upCellTypes", cols = c("white", "#cb181d"), 
                       dot.scale = 4, col.min = 0, dot.min = 0.1
) + theme(axis.text.x = element_text(angle = 90, hjust = 1),
          panel.border = element_rect(colour = "black", fill=NA)) + labs(x="", y="") 
#pdf("Fig_Dotplot_markers.pdf", useDingbats = F, width = 9, height = 3)
plot_marker
#dev.off()


################################################
####---- Fig 8 cell proportion analysis ----####
################################################
# step 1. 12 individuals divided into cell types and their numbers
other_sub$upCellTypes <- factor(other_sub$upCellTypes, levels = fac_levs)
step1 <- table(other_sub$orig.ident, other_sub$upCellTypes) %>% as.data.frame()
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

#pdf("Fig_cellproportion_other.pdf", useDingbats = F, width = 5, height = 5)
p_cp_5
#dev.off()

#########################################################
# run the statistics of cell proportions compared to Skin
df.group.out <- step3 %>% mutate(Group = gsub("PWH..", "", Var1)) 
data.table::fwrite(df.group.out, "AcuteWound_normCellProportion_Endo.txt", sep = "\t")
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
data.table::fwrite(stat_sig_all, "Figure_AcuteWound_normCellProportion_Endo_addSig.txt", sep = "\t")

# markers
markers <- data.table::fread("s1_clean_othercell_mg_res0.5.txt")
top200markers <- markers %>% group_by(cluster) %>% top_n(n = 200, wt = avg_log2FC) %>% arrange(cluster)
data.table::fwrite(top200markers, "Markers_top200_subEndo.txt", sep = "\t")


