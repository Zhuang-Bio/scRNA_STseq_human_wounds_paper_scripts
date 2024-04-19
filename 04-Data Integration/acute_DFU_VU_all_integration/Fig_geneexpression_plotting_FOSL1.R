library(tidyverse)
library(ggplot2)
library(patchwork)
library(ggbreak)
library(ggpubr)

library(Seurat)
library(scCustomize)

setwd("./Fig2_3_cellProportion_exp")

######################################################################################
####---- Figure 2h Cell proportion of Bas/Spi KCs in acute and chronic wounds ----####
#### THE codes of extracting cell proportions of FOSL1 in different datasets were shown in the end. 
rm(list=ls())
wds <- data.table::fread("cellproportion_FOSL1_wounds.txt")
colnames(wds)[1:2] <- c("Sample", "CellType")
dfu <- data.table::fread("cellproportion_FOSL1_dfu.txt") %>% slice(1:188)
colnames(dfu) <- colnames(wds)
vu <- data.table::fread("cellproportion_FOSL1_VU_240118.txt")
vu$Group <- vu$Sample
#colnames(vu) <- colnames(wds)
#vu$CellType <- gsub("Basmig", "Bas_mig", vu$CellType)
#vu$CellType <- gsub("Spimig", "Spi_mig", vu$CellType)
#vu <- vu %>% filter(! Sample %in% c("NS23")) # compare if need to remove the NS23 sample

# combine all the data
alldf <- rbind(wds, dfu) %>% mutate(Prop = Exp/smTot) %>% 
  mutate(Group = Sample) %>% 
  mutate(Group = gsub("^PWH[0-9]{2}", "", Group)) %>%
  mutate(Group = gsub("H[0-9]$", "H", Group)) %>% 
  mutate(Group = gsub("NS[0-9]{2}$", "NS", Group)) %>% 
  mutate(Group = gsub("VU[0-9]$", "VU", Group))

# Bas_mig
## First calculate the mean value of skin
skin.avg <- alldf %>% 
  #filter(Group %in% c("D0", "D1", "D7", "D30")) %>% 
  filter(Group %in% c("H", "DFU_H", "DFU_NH")) %>% 
  #filter(Group %in% c("NS", "VU")) %>% 
  #filter(CellType %in% c("Bas_mig")) %>% 
  group_by(Group) %>% summarize(mean(Prop)) %>% ungroup();skin.avg

alldf.bas <- alldf %>% 
  #filter(Group %in% c("D0", "D1", "D7", "D30")) %>% 
  filter(Group %in% c("H", "DFU_H", "DFU_NH")) %>% 
  #filter(Group %in% c("NS", "VU")) %>% 
  #filter(CellType %in% c("Bas_mig")) %>% 
  mutate(NormProp = Prop /skin.avg$`mean(Prop)`[3])
#alldf.bas$Group <- factor(alldf.bas$Group, levels = c("D0", "D1", "D7", "D30"))
alldf.bas$Group <- factor(alldf.bas$Group, levels = c("H", "DFU_H", "DFU_NH"))
#alldf.bas$Group <- factor(alldf.bas$Group, levels = c("NS", "VU"))
alldf.bas.summary <- alldf.bas %>% group_by(Group) %>% 
  summarise(sd = sd(NormProp, na.rm = TRUE), NormProp = mean(NormProp));alldf.bas.summary

p1.vu <- ggplot(alldf.bas, aes(Group, log(NormProp+1))) +
  geom_col(data = alldf.bas.summary, fill = NA, color = "black") +
  #scale_y_break(c(30,240)) +
  geom_jitter( position = position_jitter(0.1), color = "black") + 
  #geom_errorbar(aes(ymin = NormProp-sd, ymax = NormProp+sd), 
  #              data = alldf.bas.summary, width = 0.2) +
  scale_y_continuous(limits = c(-0.1,6),breaks = c(0,1,2,3,4,5,6)) +
  ylab("log(Fold change of FOSL1+ cell proportions compared to skin)") +
  geom_hline(yintercept=log(1+1), linetype="dashed", color = "red") +
  theme_classic()
p1.vu
#pdf("CellProportion_Bas_mig_FOSL1.pdf", useDingbats = FALSE, width = 5, height = 3)
p1.wds + p1.dfu + p1.vu + plot_layout(widths = c(4,3,1.5))
#dev.off()

# add the significance of groups
#my_comparisons <- list(c("D1", "D0"), c("D7", "D0"), c("D30", "D0"))
#p1 + stat_compare_means(comparisons = my_comparisons, label = "p.signif")
skin.avg <- alldf %>% 
  #filter(Group %in% c("D0", "D1", "D7", "D30")) %>% 
  #filter(Group %in% c("H", "DFU_H", "DFU_NH")) %>% 
  filter(Group %in% c("NS", "VU")) %>% 
  filter(CellType %in% c("Bas_mig")) %>% group_by(Group) %>% summarize(mean(Prop))
skin.avg
wds.bas <- alldf %>% filter(CellType %in% c("Bas_mig")) %>% 
  mutate(NormProp = Prop /skin.avg$`mean(Prop)`[1])
#wds.bas$Group <- factor(wds.bas$Group, levels = c("D0", "D1", "D7", "D30"))
#wds.bas$Group <- factor(wds.bas$Group, levels = c("H", "DFU_H", "DFU_NH"))
wds.bas$Group <- factor(wds.bas$Group, levels = c("NS", "VU"))
# prepare the data for significance test
#df.prop <- wds.bas %>% filter(Group %in% c("DFU_NH", "DFU_H"))
#df.prop$Group <- factor(df.prop$Group, levels = c("DFU_NH", "DFU_H"))
test.sum <- aov(NormProp ~ Group, data = df.prop)
summary(test.sum)

t.test(NormProp ~ Group, data = df.prop)
wilcox.test(NormProp ~ Group, data = df.prop)

test.quasi = glm(formula = Prop ~ Group, data = df.prop, family=quasibinomial)
print(summary(test.quasi))
anova(test.quasi, test = "LRT")$`Pr(>Chi)`[2]


# all KC FOSL1+ expression

## only plot FOSL1 in all KCs in VU (20240118)
skin.avg <- vu %>% 
  filter(Group %in% c("NS", "VU")) %>% filter(! Var1 %in% c("NS23")) %>% 
  group_by(Group) %>% summarize(mean(Prop)) %>% ungroup();skin.avg

alldf.bas <- vu %>% 
  filter(Group %in% c("NS", "VU")) %>% filter(! Var1 %in% c("NS23")) %>% 
  mutate(NormProp = Prop /skin.avg$`mean(Prop)`[1]) %>% distinct(Group, NormProp)
alldf.bas$Group <- factor(alldf.bas$Group, levels = c("NS", "VU"))
alldf.bas.summary <- alldf.bas %>% group_by(Group) %>% 
  summarise(sd = sd(NormProp, na.rm = TRUE), NormProp = mean(NormProp));alldf.bas.summary

p1.vu <- ggplot(alldf.bas, aes(Group, log2(NormProp+1))) +
  geom_col(data = alldf.bas.summary, fill = NA, color = "black") +
  #scale_y_break(c(30,240)) +
  geom_jitter( position = position_jitter(0.1), color = "black") + 
  geom_errorbar(aes(ymin = log2(NormProp-sd+1), ymax = log2(NormProp+sd+1)), 
                data = alldf.bas.summary, width = 0.2) +
  scale_y_continuous(limits = c(-0.1,6),breaks = c(0,1,2,3,4,5,6)) +
  ylab("log(Fold change of FOSL1+ cell proportions compared to skin)") +
  geom_hline(yintercept=log2(1+1), linetype="dashed", color = "red") +
  theme_classic()
p1.vu
pdf("CellProportion_allKCs_FOSL1_onlyVU.pdf", useDingbats = FALSE, width = 3, height = 4)
p1.vu
dev.off()

## First calculate the mean value of skin
skin.avg <- alldf %>% 
  #filter(Group %in% c("D0", "D1", "D7", "D30")) %>% 
  filter(Group %in% c("H", "DFU_H", "DFU_NH")) %>% 
  #filter(Group %in% c("NS", "VU")) %>% filter(! Sample %in% c("NS23")) %>% 
  group_by(Sample) %>% mutate(Prop_tot=sum(Prop)) %>% ungroup() %>% 
  group_by(Group) %>% summarize(mean(Prop_tot)) %>% ungroup();skin.avg

alldf.bas <- alldf %>% 
  #filter(Group %in% c("D0", "D1", "D7", "D30")) %>% 
  filter(Group %in% c("H", "DFU_H", "DFU_NH")) %>% 
  #filter(Group %in% c("NS", "VU")) %>% filter(! Sample %in% c("NS23")) %>% 
  group_by(Sample) %>% mutate(Prop_tot=sum(Prop)) %>% ungroup() %>% 
  mutate(NormProp = Prop_tot /skin.avg$`mean(Prop_tot)`[3]) %>% distinct(Group, NormProp)
alldf.bas$Group <- factor(alldf.bas$Group, levels = c("D0", "D1", "D7", "D30"))
alldf.bas$Group <- factor(alldf.bas$Group, levels = c("H", "DFU_H", "DFU_NH"))
#alldf.bas$Group <- factor(alldf.bas$Group, levels = c("NS", "VU"))
alldf.bas.summary <- alldf.bas %>% group_by(Group) %>% 
  summarise(sd = sd(NormProp, na.rm = TRUE), NormProp = mean(NormProp));alldf.bas.summary

p1.dfu <- ggplot(alldf.bas, aes(Group, log2(NormProp+1))) +
  geom_col(data = alldf.bas.summary, fill = NA, color = "black") +
  #scale_y_break(c(30,240)) +
  geom_jitter( position = position_jitter(0.1), color = "black") + 
  geom_errorbar(aes(ymin = log2(NormProp-sd+1), ymax = log2(NormProp+sd+1)), 
                data = alldf.bas.summary, width = 0.2) +
  scale_y_continuous(limits = c(-0.1,6),breaks = c(0,1,2,3,4,5,6)) +
  ylab("log(Fold change of FOSL1+ cell proportions compared to skin)") +
  geom_hline(yintercept=log2(1+1), linetype="dashed", color = "red") +
  theme_classic()
p1.wds
#pdf("wds_bas_mig.pdf", useDingbats = FALSE, width = 2, height = 3)
#pdf("CellProportion_Spi_mig_FOSL1.pdf", useDingbats = FALSE, width = 5, height = 3)
#
pdf("CellProportion_allKCs_FOSL1.pdf", useDingbats = FALSE, width = 5, height = 4)
#p1.wds + p1.dfu + p1.vu + plot_layout(widths = c(4,3,1.5))
p1.wds + p1.dfu + p1.vu + plot_layout(widths = c(2,1.5,1))
dev.off()


 # check the significance of VU
skin.avg <- alldf %>% 
  #filter(Group %in% c("D0", "D1", "D7", "D30")) %>% 
  #filter(Group %in% c("H", "DFU_H", "DFU_NH")) %>% 
  filter(Group %in% c("NS", "VU")) %>% #filter(! Sample %in% c("NS23")) %>% 
  group_by(Sample) %>% mutate(Prop_tot=sum(Prop)) %>% ungroup() %>% 
  group_by(Group) %>% summarize(mean(Prop_tot)) %>% ungroup();skin.avg

wds.bas <- alldf %>% 
  #filter(Group %in% c("D0", "D1", "D7", "D30")) %>% 
  #filter(Group %in% c("H", "DFU_H", "DFU_NH")) %>% 
  filter(Group %in% c("NS", "VU")) %>% filter(! Sample %in% c("NS23")) %>% 
  group_by(Sample) %>% mutate(Prop_tot=sum(Prop)) %>% ungroup() %>% 
  mutate(NormProp = Prop_tot /skin.avg$`mean(Prop_tot)`[1]) %>% distinct(Group, NormProp)
#wds.bas$Group <- factor(wds.bas$Group, levels = c("D0", "D1", "D7", "D30"))
#wds.bas$Group <- factor(wds.bas$Group, levels = c("H", "DFU_H", "DFU_NH"))
wds.bas$Group <- factor(wds.bas$Group, levels = c("NS", "VU"))
# prepare the data for significance test
df.prop <- wds.bas %>% filter(Group %in% c("NS", "VU"))
df.prop$Group <- factor(df.prop$Group, levels = c("NS", "VU"))
test.sum <- aov(NormProp ~ Group, data = df.prop)
summary(test.sum)

t.test(NormProp ~ Group, data = df.prop)
wilcox.test(NormProp ~ Group, data = df.prop)

test.quasi = glm(formula = NormProp ~ Group, data = df.prop, family=quasibinomial)
print(summary(test.quasi))
anova(test.quasi, test = "LRT")$`Pr(>Chi)`[2]




################################################
####---- Expression of interested genes ----####
## Acute wounds
output_path <- "/Users/zhuliu/Desktop/scRNA_STseq/proj_10X_woundhealing/03_results/01-Seurat-PreAnalysis/02_Seurat_BatchCorrection/s2_Seurat_subKeratinocytes/"
kera_sub <- readRDS(paste0(output_path, "allNew_subcluster_keratins_220203.rds"))
# load the annotations from main clusters
mt <- data.table::fread("/Users/zhuliu/Desktop/sc st/plotting scripts/Figure 1/allcombined_wounds_newAnnotation_metadata.txt")
mt_sub <- colnames(kera_sub) %>% as.data.frame() %>% setNames("barcode") %>% 
  left_join(., mt[, c(1, 24, 25, 27)], by=c("barcode"="barcode"))
kera_sub$newCellTypes <- mt_sub$newCellTypes
kera_sub$newMainCellTypes <- mt_sub$newMainCellTypes
kera_sub$doublet_scores <- mt_sub$doublet_scores

# re-annotate the cell clusters
fac_levs <- c("Bas_I", "Bas_prolif", "Bas_mig", 
              "Spi_I", "Spi_II_a", "Spi_II_b",
              "Spi_III", "Spi_mig",
              "Gra_I")
current_cluster_ids <- c(1,5,6,3,0,4,8,2,7) # List of current cluster IDs
tmp <- plyr::mapvalues(x=kera_sub$SCT_snn_res.0.5, from=current_cluster_ids, to=fac_levs)

kera_sub$upCellTypes <- tmp
kera_sub$upCellTypes <- factor(kera_sub$upCellTypes, levels = fac_levs)
DimPlot(kera_sub, group.by = "upCellTypes", label = F) +
  NoAxes() + ggtitle("")
FeaturePlot(kera_sub, features = c("EGFR"), cols = c("grey90", "red"), split.by = "Condition")

# Try to manually calculate the FOSL1 postive cell numbers
smtol <- table(kera_sub$orig.ident) %>% as.data.frame() %>% setNames(c("Sample", "smTot"))
gene_data <- FetchData(kera_sub, vars = c("Condition", "orig.ident", "upCellTypes","FOSL1"), slot = "data")
gene_data <- gene_data %>% rownames_to_column(var = "barcode")
pernt <- gene_data %>% group_by(orig.ident, upCellTypes) %>% 
  summarise(Exp = sum(FOSL1 >0), 
            nonExp = sum(FOSL1 < 0.00001),
            ctTot = sum(FOSL1 > -Inf)) %>% ungroup() %>% 
  left_join(., smtol, by=c("orig.ident" = "Sample"))
data.table::fwrite(pernt, file = "/Users/zhuliu/Desktop/sc st/Versions_SC_ST_Project/V1_Figures_230120/Fig2_3_cellProportion_exp/cellproportion_FOSL1_wounds.txt", sep = "\t")


## DFU
inteData <- readRDS("/Users/zhuliu/Desktop/sc st/compareStudies/Diabetic foot ulcers/step4_DFU_Wound_inteKC.rds")
ct.cols <- c('#807dba','#9e9ac8','#ccebc5',
             '#fe9929','#fec44f','#fee391',
             '#fb8072','#b3de69','#fccde5')

# cell type proportion analysis
fac_levs <- c("Bas_I", "Bas_prolif", "Bas_mig", "Spi_I", "Spi_II_a", "Spi_II_b", "Spi_III", "Spi_mig", "Gra_I", "Other")
DimPlot(inteData, group.by = "newCellType", label = T, split.by = "Condition")
FeaturePlot(inteData, features = c("FOSL1"), cols = c("grey90", "red"), split.by = "Condition")
FeaturePlot(inteData, features = c("EGFR"), cols = c("grey90", "red"), split.by = "Condition")
VlnPlot(inteData, features = c("EGFR"), group.by = "Condition", pt.size = 0)


# Try to manually calculate the FOSL1 postive cell numbers
smtol <- table(inteData$Sample) %>% as.data.frame() %>% setNames(c("Sample", "smTot"))
gene_data <- FetchData(inteData, vars = c("Condition", "Sample", "newCellType","FOSL1"), slot = "data")
gene_data <- gene_data %>% rownames_to_column(var = "barcode")
pernt <- gene_data %>% group_by(Sample, newCellType) %>% 
  summarise(Exp = sum(FOSL1 >0), 
            nonExp = sum(FOSL1 < 0.00001),
            ctTot = sum(FOSL1 > -Inf)) %>% ungroup() %>% 
  left_join(., smtol, by=c("Sample" = "Sample"))
data.table::fwrite(pernt, file = "/Users/zhuliu/Desktop/sc st/Versions_SC_ST_Project/V1_Figures_230120/Fig2_3_cellProportion_exp/cellproportion_FOSL1_dfu.txt", sep = "\t")


## VU
inteData <- readRDS("/Volumes/zhuliu/Groups/Ning Xu Landén/Lihua_Luo/P4_Zhuang/Step2_ComparedVUandWound_SubcellTypes/C1_KC/s4_AfterHarmony.rds")
DimPlot(inteData, group.by = "CellType")
DimPlot(inteData, group.by = "SCT_snn_res.0.5", label = T)

inteData <- subset(inteData, subset = Project == "VU")
DimPlot(inteData, group.by = "SCT_snn_res.0.5", label = T)

FeaturePlot(inteData, features = c("FOSL1"), cols = c("grey90", "red"), split.by = "orig.ident")
FeaturePlot(inteData, features = c("EGFR"), cols = c("grey90", "red"), split.by = "Condition")
VlnPlot(inteData, features = c("EGFR"), group.by = "Condition", pt.size = 0.001)

# Try to manually calculate the FOSL1 postive cell numbers
smtol <- table(inteData$orig.ident) %>% as.data.frame() %>% setNames(c("Sample", "smTot"))
gene_data <- FetchData(inteData, vars = c("Condition", "orig.ident", "SCT_snn_res.0.5","FOSL1"), slot = "data")
gene_data <- gene_data %>% rownames_to_column(var = "barcode")
pernt <- gene_data %>% group_by(orig.ident, SCT_snn_res.0.5) %>% 
  summarise(Exp = sum(FOSL1 >0), 
            nonExp = sum(FOSL1 < 0.00001),
            ctTot = sum(FOSL1 > -Inf)) %>% ungroup() %>% 
  left_join(., smtol, by=c("orig.ident" = "Sample"))
data.table::fwrite(pernt, file = "/Users/zhuliu/Desktop/sc st/Versions_SC_ST_Project/V1_Figures_230120/Fig2_3_cellProportion_exp/cellproportion_FOSL1_vu_addVU6.txt", sep = "\t")

