library(tidyverse)
library(ggplot2)
library(patchwork)
library(ggbreak)
library(ggpubr)

setwd("./")


########################################
## Immune cell dynamics in wound healing
### read all the cell information (include epidermis and dermis)
acute_mt <- data.table::fread("/Users/zhuliu/Desktop/sc st/plotting scripts/Figure 1/allcombined_wounds_newAnnotation_metadata.txt")

# load the original cell number of epidermis and dermis
acute_mt_epider <- data.table::fread("/Users/zhuliu/Desktop/scRNA_STseq/proj_10X_woundhealing/03_results/01-Seurat-PreAnalysis/02_Seurat_BatchCorrection/s1_Seurat_allSample_harmony/all_epidermis_dermis_perIndi_220613.txt")
acute_mt_epider <- acute_mt_epider %>% filter(sepa == "Dermis") %>% select(-1) %>% setNames(c("Sample", "Total"))

allcell <- table(acute_mt$newMainCellTypes, acute_mt$orig.ident) %>% as.data.frame() %>% left_join(., acute_mt_epider, by=c("Var2" = "Sample")) %>% 
  mutate(pro=Freq / Total)

allmye <- allcell %>% filter(Var1 %in% c("Myeloid", "Lymphoid")) %>% 
  mutate(Group = gsub("PWH[0-9][0-9]", "", Var2)) %>% 
  group_by(Group, Var1) %>% summarise( avgProp = mean(pro)) %>% ungroup() %>% setNames(c("Group", "CellType", "Proportion"))
allmye$Group <- factor(allmye$Group, levels = c("D0", "D1", "D7", "D30"))

ggplot(data=allmye, aes(x=Group, y=Proportion, group=CellType, colour=CellType)) +
  geom_line()+
  geom_point(aes(shape=CellType)) +
  scale_color_manual(values=c("#E69F00", "purple"))

# Try to draw a subCellType proportion
acute_mt_sm <- table(acute_mt$orig.ident, acute_mt$newCellTypes) %>% as.data.frame() %>% setNames(c("Sample", "CellType", "Freq")) %>% 
  left_join(., acute_mt_epider, by=c("Sample" = "Sample")) %>% 
  mutate(pro=Freq / Total)

allimmune <- acute_mt_sm %>% filter(CellType %in% c("Mono-Mac", "cDC1", "cDC2", "DC3", "Mast-cell", "NK-cell", "Th", "Plasma_Bcell")) %>% 
  mutate(Group = gsub("PWH[0-9][0-9]", "", Sample)) %>% 
  group_by(Group, CellType) %>% summarise( avgProp = mean(pro)) %>% ungroup() %>% setNames(c("Group", "CellType", "Proportion"))
allimmune$Group <- factor(allimmune$Group, levels = c("D0", "D1", "D7", "D30"))

ggplot(data=allimmune, aes(x=Group, y=Proportion, group=CellType, colour=CellType)) +
  geom_line()+
  geom_point(aes(shape=CellType)) #+
  #scale_color_manual(values=c("#E69F00", "purple"))



###################################
####---- Manuscript figure ----####
# Try to draw a subCellType proportion from each subclustering analysis
### read all the cell information (include epidermis and dermis)
acute_mt <- data.table::fread("/Users/zhuliu/Desktop/sc st/plotting scripts/Figure 1/allcombined_wounds_newAnnotation_metadata.txt")
# load the original cell number of epidermis and dermis
acute_mt_epider <- data.table::fread("/Users/zhuliu/Desktop/scRNA_STseq/proj_10X_woundhealing/03_results/01-Seurat-PreAnalysis/02_Seurat_BatchCorrection/s1_Seurat_allSample_harmony/all_epidermis_dermis_perIndi_220613.txt")
acute_mt_epider <- acute_mt_epider %>% 
  #group_by(sepa,Var1) %>% summarise(Total = sum(a_sum)) %>% 
  #setNames(c("Sample", "Total"))
  #filter(sepa == "Dermis") %>% 
  #select(-1) %>% 
  setNames(c("Loca", "Sample", "Total")) # use this for separating dermis and epidermis
#acute_mt_cd45 <- data.table::fread("CD45_positiveCellNumber.txt") %>% setNames(c("Sample", "Total"))

neu <- data.table::fread("/Users/zhuliu/Desktop/sc st/plotting scripts/Figure 6/Neutrophils/neu_clustering_metadata.txt")
mye <- data.table::fread("/Users/zhuliu/Desktop/sc st/plotting scripts/Figure 6/myeloid_newAnnotation_metadata.txt")
lym <- data.table::fread("/Users/zhuliu/Desktop/sc st/plotting scripts/Figure 7/lymphoid_newAnnotation_metadata.txt")


neu_df <- table(neu$ID, neu$Condition) %>% as.data.frame() %>% mutate(Var2 = "Neu") %>% setNames(c("Sample", "CellType", "Freq")) %>% mutate(Type = "Myeloid", Loca = "Dermis")
mye_df <- table(mye$orig.ident, mye$upCellTypes) %>% as.data.frame() %>% setNames(c("Sample", "CellType", "Freq")) %>% mutate(Type = "Myeloid", Loca = "Dermis")
mye_df$Loca <- ifelse(mye_df$CellType == "LC", "Epidemis", mye_df$Loca)
lym_df <- table(lym$orig.ident, lym$upCellTypes) %>% as.data.frame() %>% setNames(c("Sample", "CellType", "Freq")) %>% mutate(Type = "Lymphoid", Loca = "Dermis")

allimmune_com <- rbind(neu_df, mye_df, lym_df) %>% 
  left_join(., acute_mt_epider, by=c("Sample", "Loca")) %>% # remove the Loca if we use total number of epidermis and dermis
  #left_join(., acute_mt_cd45, by=c("Sample")) %>% # normalize by CD45+ cell number
  mutate(pro=Freq / Total)
table(allimmune_com$CellType)

allimmune <- allimmune_com %>% filter(!CellType %in% c("Apoptotic", "Cycling")) %>% 
  mutate(Group = gsub("PWH[0-9][0-9]", "", Sample)) %>% 
  group_by(Group, CellType) %>% summarise( avgProp = mean(pro)) %>% ungroup() %>% setNames(c("Group", "CellType", "Proportion"))
allimmune$Group <- factor(allimmune$Group, levels = c("D0", "D1", "D7", "D30"))

require(ggrepel)
ggplot(data=allimmune, aes(x=Group, y=Proportion, group=CellType, colour=CellType)) +
  geom_line()+
  geom_point(aes(shape=CellType)) +
  geom_text_repel(
    aes(label = CellType), data = allimmune,
    fontface ="plain", color = "black", size = 3
  )
#pdf("all_immune_cells_dynamics.pdf", useDingbats = F, width = 10, height = 8)
#dev.off()
#scale_color_manual(values=c("#E69F00", "purple"))


tmp <- allimmune %>% group_by(CellType, Group) %>% summarise(avgProp = Proportion) %>% ungroup()

# separate the figure into four condition
skin <- c("Tc", "cDC2", "Mac2")
wound1 <- c("Mac_inf", "Mac1", "Neu", "DC3", "Treg", "Th", "ILCs")
wound7 <- c("Mac3", "pDC", "cDC1", "ILC1/NK", "NK", "Plasma", "Bcell")
wound30 <- c("LC")

allimmune <- allimmune %>% 
  mutate(condition = case_when(
    CellType %in% skin ~ "Skin",
    CellType %in% wound1 ~ "Wound1",
    CellType %in% wound7 ~ "Wound7",
    CellType %in% wound30 ~ "Wound30",
  )
)

allimmune$condition <- factor(allimmune$condition, levels = c("Skin", "Wound1", "Wound7", "Wound30"))
ggplot(data=allimmune, aes(x=Group, y=Proportion, group=CellType, colour=CellType)) +
  geom_line()+
  geom_point(aes(shape=CellType)) +
  geom_text_repel(
    aes(label = CellType), data = allimmune,
    fontface ="plain", color = "black", size = 3
  ) +
  theme_classic() + 
  facet_wrap(vars(condition))

pdf("all_immune_cells_dynamics_splitCondition4_cd45.pdf", useDingbats = F, width = 10, height = 10)
dev.off()


############################
# Remove the minor changes
## One-way anova
ct <- unique(allimmune_com$CellType)

sigtest <- list()
for (i in seq_along(ct)) {
  ct_tmp <- ct[i]
  df_tmp <- allimmune_com %>% dplyr::filter(CellType == ct_tmp) %>% distinct() %>% 
    mutate(Group = gsub("PWH[0-9][0-9]", "", Sample))
  # Compute the analysis of variance
  res.aov <- aov(pro~Group, data = df_tmp)
  # Summary of the analysis
  sigtest[[ct_tmp]] <- summary(res.aov)[[1]][["Pr(>F)"]][1]
}
df <- do.call("rbind", sigtest) %>% as.data.frame() %>% setNames(c("pvalue")) %>% rownames_to_column(var="CellType")
df.rm <- df %>% filter(pvalue>0.05)
#> df
#CellType       pvalue
#1        Neu 3.571894e-02
#2  Apoptotic 3.341339e-02
#3       cDC1 8.685537e-02
#4       cDC2 5.864249e-03
#5    Cycling 3.886038e-02
#6        DC3 1.758871e-07
#7         LC 4.532755e-01
#8    Mac_inf 9.253923e-03
#9       Mac1 2.446229e-04
#10      Mac2 1.946880e-02
#11      Mac3 8.635920e-02
#12       pDC 2.537743e-02
#13     Bcell 2.116552e-01
#14   ILC1/NK 6.116509e-02
#15      ILCs 5.956438e-03
#16        NK 1.790555e-01
#17    Plasma 7.822113e-03
#18        Tc 2.759176e-01
#19        Th 1.303164e-02
#20      Treg 3.795283e-01
#21      Ttol 9.001234e-01

allimmune <- allimmune_com %>% filter(!CellType %in% c("Apoptotic", "Cycling")) %>% 
  filter(!CellType %in% df.rm$CellType) %>% 
  mutate(Group = gsub("PWH[0-9][0-9]", "", Sample)) %>% 
  group_by(Group, CellType) %>% summarise( avgProp = mean(pro)) %>% ungroup() %>% setNames(c("Group", "CellType", "Proportion"))
allimmune$Group <- factor(allimmune$Group, levels = c("D0", "D1", "D7", "D30"))

require(ggrepel)
ggplot(data=allimmune, aes(x=Group, y=Proportion, group=CellType, colour=CellType)) +
  geom_line()+
  geom_point(aes(shape=CellType)) +
  geom_text_repel(
    aes(label = CellType), data = allimmune,
    fontface ="plain", color = "black", size = 3
  )
#pdf("all_immune_cells_dynamics_OneWayAnova.pdf", useDingbats = F, width = 10, height = 8)
#dev.off()
#scale_color_manual(values=c("#E69F00", "purple"))


# separate the figure into four condition
skingain <- c("cDC2", "Mac2", "ILCs")
allimmune$condition <- ifelse(allimmune$CellType %in% skingain, "Skin", "Wounds")

allimmune$condition <- factor(allimmune$condition, levels = c("Skin", "Wounds"))
ggplot(data=allimmune, aes(x=Group, y=Proportion, group=CellType, colour=CellType)) +
  geom_line()+
  geom_point(aes(shape=CellType)) +
  geom_text_repel(
    aes(label = CellType), data = allimmune,
    fontface ="plain", color = "black", size = 3
  ) +
  theme_classic() + 
  facet_wrap(vars(condition))

#pdf("all_immune_cells_dynamics_OneWayAnova2.pdf", useDingbats = F, width = 12, height = 4)
#dev.off()


##########################
# Final manuscript figure
# only keep the interested cell types
allimmune_f <- allimmune %>% filter(CellType %in% c("Neu", "Mac_inf", "Mac1", "DC3", "Th", "Plasma"))


migKC <- data.table::fread("/Users/zhuliu/Desktop/sc st/Versions_SC_ST_Project/V2_Figures_231017/01-Scripts/Figures plotting scripts/Fig2_3_cellProportion_exp/cellproportion_normalized_wounds.txt") %>% 
  filter(Var2 %in% c("Bas_mig", "Spi_mig")) %>% dplyr::group_by(Group, Var2) %>% summarise(Proportion = mean(Prop)) %>% 
  ungroup() %>% dplyr::mutate(Condition = "MigKC")

colnames(migKC) <- colnames(allimmune_f)

allimmune_com <- rbind(allimmune_f, migKC)
allimmune_com$Group <- factor(allimmune_com$Group, levels = c("D0", "D1", "D7", "D30"))

ggplot(data=allimmune_com, aes(x=Group, y=Proportion, group=CellType, colour=CellType)) +
  geom_line()+
  geom_point() +
  geom_text_repel(
    aes(label = CellType), data = allimmune_com,
    fontface ="plain", color = "black", size = 3
  ) +
  theme_classic() +
  facet_wrap(vars(condition))

allimmune_com1 <- allimmune_com %>% filter(condition == "Wounds")
ggplot(data=allimmune_com1, aes(x=Group, y=Proportion, group=CellType, colour=CellType)) +
  geom_line()+
  geom_point() +
  geom_text_repel(
    aes(label = CellType), data = allimmune_com1,
    fontface ="plain", color = "black", size = 3
  ) +
  scale_y_continuous(limits = c(0, 0.1), breaks = c(0, 0.02, 0.04, 0.06, 0.08, 0.1)) +
  theme_classic() 
pdf("Figure4j_all_immune_cells_dynamics1.pdf", useDingbats = F, width = 6, height = 4)
dev.off()

allimmune_com2 <- allimmune_com %>% filter(condition == "MigKC")
ggplot(data=allimmune_com2, aes(x=Group, y=Proportion, group=CellType, colour=CellType)) +
  geom_line()+
  geom_point() +
  geom_text_repel(
    aes(label = CellType), data = allimmune_com2,
    fontface ="plain", color = "black", size = 3
  ) +
  scale_y_continuous(limits = c(0, 0.25), breaks = c(0, 0.05, 0.1, 0.15, 0.2, 0.25)) +
  theme_classic() 
pdf("Figure4j_all_immune_cells_dynamics2.pdf", useDingbats = F, width = 6, height = 4)
dev.off()

