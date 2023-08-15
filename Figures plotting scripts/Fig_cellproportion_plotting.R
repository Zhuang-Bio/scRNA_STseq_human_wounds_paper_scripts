library(tidyverse)
library(ggplot2)
library(patchwork)
library(ggbreak)
library(ggpubr)
library(Seurat)

setwd("/Users/zhuliu/Desktop/sc st/Versions_SC_ST_Project/V1_Figures_230529/01-Scripts")

################################################
####---- Expression of interested genes ----####
## Acute wounds
inteData <- readRDS("/Users/zhuliu/Desktop/sc st/plotting scripts/Figure 1/allcombined_wounds_newAnnotation.rds")
VlnPlot(inteData, features = c("EGF", "TGFA", "AREG"), group.by = "newCellTypes", split.by = "Condition", ncol = 1)
VlnPlot(inteData, features = c("TGFA", "AREG", "HBEGF", "EREG"), group.by = "CellType", split.by = "Condition", ncol = 1)
VlnPlot(inteData, features = c("EGF", "IFNG"), group.by = "Condition", ncol = 1)

FeaturePlot(inteData, features = c("EGF", "IFNG"), cols = c("grey90", "red"))
FeaturePlot(inteData, features = c("EGF", "EGFR"), cols = c("grey90", "red"), split.by = "Condition")
FeaturePlot(inteData, features = c("TGFA", "AREG", "HBEGF", "EREG"), cols = c("grey90", "red"), split.by = "Condition")
FeaturePlot(inteData, features = c("IFNG", "IFNGR1", "IFNGR2"), cols = c("grey90", "red"), split.by = "Condition")


## DFU Samples
inteData <- readRDS("/Users/zhuliu/Desktop/sc st/compareStudies/Diabetic foot ulcers/step1_diabetic_foot_ulcer_skin_integrated_harmony.rds")
DotPlot(inteData, features = c("IFNG", "IFNGR1", "IFNGR2"), group.by = "Condition")
table(inteData$Donor1)

DimPlot(inteData, group.by = "SCT_snn_res.1", label = T)
# re-annotate the cell clusters
anno.cl <- list()
anno.cl$HE_Fibro = c(12) #MMP3, MMP1, CHI3L2, CHI3L1, TNFAIP6
anno.cl$FB_I_III_DFU = c(2,8,25) #ASPN, POSTN, WISP2, IGFBP6, COMP, APCDD1, PI16, ELN
anno.cl$FB_II_DFU = c(1,18) #APOD, CXCL14, C3, APOE, C2orf40
anno.cl$NKT_Lympho = c(5) #CD69, IL32, CD52, IL7R, PTPRC, NKG7
anno.cl$M1_Macro = c(11) #IL1B,EREG,
anno.cl$M2_Macro = c(9) #C1QA, C1QB, CD14, CD163
anno.cl$B_Lympho = c(23) #CD79A,MS4A1
anno.cl$VasEndo_art = c(16) #"IGFBP3","HEY1","SEMA3G"
anno.cl$VasEndo_ven = c(3) #RGCC
anno.cl$LymphEndo = c(21)
anno.cl$SMC2 = c(15) #MKI67, CENPF
anno.cl$SMC1 = c(0,4,7,13,24)
anno.cl$Melano_Schwann = c(19) #DCT, SOX10, SOX2
anno.cl$Mast = c(20) #TPSAB1, MS4A2
anno.cl$Sweat_Seba = c(17) #DCD, SCGB2A2, MUCL1, KRT19
anno.cl$OtherKera = c(14,22) #KRT17, KRT6B, KRT6A, S100A2, KRT16
anno.cl$BasalKera = c(10) #COL17A1, KRT5
anno.cl$DiffKera = c(6) #KRT1, KRT10, KRT16, KRT6A

trans = rep(names(anno.cl), times = unlist(lapply(anno.cl, length)))
names(trans) = unlist(anno.cl)
inteData$CellType = trans[as.character(inteData$SCT_snn_res.1)]
#DimPlot(inteData, group.by = "SCT_snn_res.1", cols = ct.cols, label = T)
DimPlot(inteData, group.by = "CellType", label = T) + ggtitle("") + NoLegend()

fac_levs <- c("BasalKera", "DiffKera", "OtherKera", "Sweat_Seba", 
              "HE_Fibro", "FB_I_III_DFU", "FB_II_DFU", "SMC1", "SMC2",
              "M1_Macro", "M2_Macro", "B_Lympho", "NKT_Lympho", 
              "LymphEndo", "VasEndo_art", "VasEndo_ven", "Melano_Schwann", "Mast")
ct.cols <- c("#d94701", "#fdae61", "#fd8d3c", "#35978f", 
             "#0570b0", "#33A02C", "#72BF5A", "#B2DF8A", "#3690c0", 
             "#f768a1", "#d4b9da", "#92c5de", "#d1e5f0",
             "#c0a700", "#1a9850", "#fb9a99", "#8d4720", "#fdbe85"
)
names(ct.cols) <- fac_levs

inteData$CellType <- factor(inteData$CellType, levels = fac_levs)
DimPlot(inteData, group.by = "CellType", label = T, cols = ct.cols) + ggtitle("") + NoLegend()
VlnPlot(inteData, features = c("CXCL1", "CXCL5"), group.by = "CellType", split.by = "Condition", ncol = 1)
VlnPlot(inteData, features = c("CXCL1", "CXCL5"), group.by = "Condition", ncol = 1)

VlnPlot(inteData, features = c("CXCR1", "CXCR2"), group.by = "CellType", split.by = "Condition", ncol = 1)
VlnPlot(inteData, features = c("CXCR1", "CXCR2"), group.by = "Condition", ncol = 1)

VlnPlot(inteData, features = c("IL6", "FOSL1"), group.by = "CellType", split.by = "Condition", ncol = 1)
VlnPlot(inteData, features = c("IL6", "FOSL1"), group.by = "Condition", ncol = 1)

FeaturePlot(inteData, features = c("CXCL1", "CXCL5", "CXCR1", "CXCR2", "IL6", "FOSL1"), cols = c("grey90", "red"))
FeaturePlot(inteData, features = c("CXCL1", "CXCL5", "CXCR1"), cols = c("grey90", "red"), split.by = "Condition")
FeaturePlot(inteData, features = c("CXCR2", "IL6", "FOSL1"), cols = c("grey90", "red"), split.by = "Condition")

genes_all <- FetchData(inteData, vars = c("Condition", "Sample", "CellType",
                                          "CXCL1", "CXCL5", "CXCR1", "CXCR2", "IL6", "FOSL1"), slot = "data")
genes_all <- genes_all %>% rownames_to_column(var = "barcode")
data.table::fwrite(genes_all, file = "all_genes_oriexpression_celltype_condition.txt", sep = "\t")

genes_exp <- AverageExpression(inteData, assays = "RNA", features = c("CXCL1", "CXCL5", "CXCR1", "CXCR2", "IL6", "FOSL1"), 
                               group.by = c("CellType", "Condition"))$RNA
genes_exp <- genes_exp %>% t()
genes_exp <- genes_exp %>% as.data.frame() %>% rownames_to_column(var = "Group")
data.table::fwrite(genes_exp, file = "all_genes_avgexpression_celltype_condition.txt", sep = "\t")

genes_exp_condi <- AverageExpression(inteData, assays = "RNA", features = c("CXCL1", "CXCL5", "CXCR1", "CXCR2", "IL6", "FOSL1"),
                                     group.by = c("Condition"))$RNA
genes_exp_condi <- genes_exp_condi %>% as.data.frame() %>% rownames_to_column(var = "Gene")
data.table::fwrite(genes_exp_condi, file = "all_genes_avgexpression_condition.txt", sep = "\t")


## VU Samples
inteData <- readRDS("/Users/zhuliu/Desktop/sc st/Versions_SC_ST_Project/V1_Figures_230120/ComparedStudies/VU_data/D04_annotatedData.rds")
DimPlot(inteData, group.by = "newCellType", label = T)
DotPlot(inteData, features = c("CXCL1", "CXCL5", "CXCR1", "CXCR2", "IL6"), group.by = "Condition") + coord_flip()
FeaturePlot(inteData, features = c("CXCL1"), split.by = "orig.ident")

DimPlot(inteData, group.by = "newCellType", label = T) + ggtitle("") + NoLegend()
VlnPlot(inteData, features = c("CXCL1", "CXCL5"), group.by = "newCellType", split.by = "Condition", ncol = 1)
VlnPlot(inteData, features = c("CXCL1", "CXCL5"), group.by = "Condition", ncol = 1)

VlnPlot(inteData, features = c("CXCR1", "CXCR2"), group.by = "newCellType", split.by = "Condition", ncol = 1)
VlnPlot(inteData, features = c("CXCR1", "CXCR2"), group.by = "Condition", ncol = 1)

VlnPlot(inteData, features = c("IL6", "FOSL1"), group.by = "newCellType", split.by = "Condition", ncol = 1)
VlnPlot(inteData, features = c("IL6", "FOSL1"), group.by = "Condition", ncol = 1)

FeaturePlot(inteData, features = c("CXCL1", "CXCL5", "CXCR1", "CXCR2", "IL6", "FOSL1"), cols = c("grey90", "red"))
FeaturePlot(inteData, features = c("CXCL1", "CXCL5", "CXCR1"), cols = c("grey90", "red"), split.by = "Condition")
FeaturePlot(inteData, features = c("CXCR2", "IL6", "FOSL1"), cols = c("grey90", "red"), split.by = "Condition")

genes_all <- FetchData(inteData, vars = c("Condition", "orig.ident", "newCellType",
                                          "CXCL1", "CXCL5", "CXCR1", "CXCR2", "IL6", "FOSL1"), slot = "data")
genes_all <- genes_all %>% rownames_to_column(var = "barcode")
data.table::fwrite(genes_all, file = "all_genes_oriexpression_celltype_condition.txt", sep = "\t")

genes_exp <- AverageExpression(inteData, assays = "RNA", features = c("CXCL1", "CXCL5", "CXCR1", "CXCR2", "IL6", "FOSL1"), 
                               group.by = c("newCellType", "Condition"))$RNA
genes_exp <- genes_exp %>% t()
genes_exp <- genes_exp %>% as.data.frame() %>% rownames_to_column(var = "Group")
data.table::fwrite(genes_exp, file = "all_genes_avgexpression_celltype_condition.txt", sep = "\t")

genes_exp_condi <- AverageExpression(inteData, assays = "RNA", features = c("CXCL1", "CXCL5", "CXCR1", "CXCR2", "IL6", "FOSL1"),
                                     group.by = c("Condition"))$RNA
genes_exp_condi <- genes_exp_condi %>% as.data.frame() %>% rownames_to_column(var = "Gene")
data.table::fwrite(genes_exp_condi, file = "all_genes_avgexpression_condition.txt", sep = "\t")



#########################################################
####---- Cell proportion analysis of CXCL1/FOSL1 ----####
library(scCustomize)
# acute wound data
pernt <- Percent_Expressing(inteData, features = c("CXCL1", "CXCL5", "FOSL1"), 
                   group_by = "orig.ident") %>% as.matrix() %>% t()
colnames(pernt) <- paste0(colnames(pernt), "_percent")
avg <- AverageExpression(inteData, assays = "RNA", features = c("CXCL1", "CXCL5", "FOSL1"), slot = "data",
                  group.by = c("orig.ident"))$RNA %>% t()
df_wd <- data.frame(avg, pernt) %>% rownames_to_column(var = "Sample") 

# DFU data
pernt <- Percent_Expressing(inteData, features = c("CXCL1", "CXCL5", "FOSL1"), 
                            group_by = "Donor1") %>% as.matrix() %>% t()
colnames(pernt) <- paste0(colnames(pernt), "_percent")
avg <- AverageExpression(inteData, assays = "RNA", features = c("CXCL1", "CXCL5", "FOSL1"), slot = "data",
                         group.by = c("Donor1"))$RNA %>% t()
df_dfu <- data.frame(avg, pernt) %>% rownames_to_column(var = "Sample") 

# VU data (Be careful with the NS23 since the FOSL1 gene is widely expressed [weird])
pernt <- Percent_Expressing(inteData, features = c("CXCL1", "CXCL5", "FOSL1"), 
                            group_by = "orig.ident") %>% as.matrix() %>% t()
colnames(pernt) <- paste0(colnames(pernt), "_percent")
avg <- AverageExpression(inteData, assays = "RNA", features = c("CXCL1", "CXCL5", "FOSL1"), slot = "data",
                         group.by = c("orig.ident"))$RNA %>% t()
df_vu <- data.frame(avg, pernt) %>% rownames_to_column(var = "Sample")

# combine all data
df_gene <- rbind(df_wd, df_dfu, df_vu) %>% 
  mutate(Group = Sample) %>% 
  mutate(Group = gsub("H[0-9]$", "H", Group)) %>% 
  mutate(Group = gsub("^PWH[0-9]{2}", "", Group)) %>% 
  mutate(Group = gsub("NS[0-9]{2}", "NS", Group)) %>% 
  mutate(Group = gsub("VU[0-9]{1}", "VU", Group))

# DFU
df_gene_io <- df_gene %>% 
  select(starts_with(c("CXCL1", "Group"))) %>% 
  #slice(-8,-12) %>% # exclude the outlier of gene expression
  filter(Group %in% c("H", "DFU_H", "DFU_NH")) 
skin.avg <- df_gene_io %>% group_by(Group) %>% summarise(mean(CXCL1_percent));skin.avg

df_gene_io <- df_gene_io %>% mutate(NormExp = CXCL1_percent /skin.avg$`mean(CXCL1_percent)`[3])
df_gene_io$Group <- factor(df_gene_io$Group, levels = c("H", "DFU_H", "DFU_NH"))
wds.bas.summary <- df_gene_io %>% group_by(Group) %>% 
  summarise(sd = sd(NormExp, na.rm = TRUE), NormExp = mean(NormExp))
wds.bas.summary

my_comparisons <- list(c("H", "DFU_H"), c("H", "DFU_NH"), c("DFU_H", "DFU_NH"))
p2.wds <- ggplot(df_gene_io, aes(Group, NormExp)) +
  geom_col(data = wds.bas.summary, fill = NA, color = "black") +
  geom_jitter( position = position_jitter(0.1), color = "black") + 
  geom_errorbar( aes(ymin = NormExp-sd, ymax = NormExp+sd), 
                 data = wds.bas.summary, width = 0.2) +
  #ylab("Fold change of gene expression compared to skin") +
  ylab("Fold change of postive cell proportions compared to skin") +
  #scale_y_continuous(limits = c(0,14),breaks = c(0, 4, 8,12,16)) +
  theme_classic()
p2.wds + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test")
#p1.wds + p2.wds
# add the significance using t.test
#p1 <- p1.wds + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test")
#p2 <- p2.wds + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test")


# VU
df_gene_io <- df_gene %>% 
  select(starts_with(c("CXCL1", "Group"))) %>% 
  filter(Group %in% c("NS", "VU")) 
skin.avg <- df_gene_io %>% group_by(Group) %>% summarise(mean(CXCL1_percent));skin.avg

df_gene_io <- df_gene_io %>% mutate(NormExp = CXCL1_percent /skin.avg$`mean(CXCL1_percent)`[1])
df_gene_io$Group <- factor(df_gene_io$Group, levels = c("NS", "VU"))
wds.bas.summary <- df_gene_io %>% group_by(Group) %>% 
  summarise(sd = sd(NormExp, na.rm = TRUE), NormExp = mean(NormExp))
wds.bas.summary

my_comparisons <- list(c("NS", "VU"))
p2.wds <- ggplot(df_gene_io, aes(Group, NormExp)) +
  geom_col(data = wds.bas.summary, fill = NA, color = "black") +
  geom_jitter( position = position_jitter(0.1), color = "black") + 
  geom_errorbar( aes(ymin = NormExp-sd, ymax = NormExp+sd), 
                 data = wds.bas.summary, width = 0.2) +
  #ylab("Fold change of gene expression compared to skin") +
  ylab("Fold change of postive cell proportions compared to skin") +
  #scale_y_continuous(limits = c(0,14),breaks = c(0, 4, 8,12,16)) +
  theme_classic()
p2.wds + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test")

# wound
df_gene_io <- df_gene %>% 
  select(starts_with(c("CXCL1", "Group"))) %>% 
  filter(Group %in% c("D0", "D1", "D7", "D30")) 
skin.avg <- df_gene_io %>% group_by(Group) %>% summarise(mean(CXCL1_percent));skin.avg

df_gene_io <- df_gene_io %>% mutate(NormExp = CXCL1_percent /skin.avg$`mean(CXCL1_percent)`[1])
df_gene_io$Group <- factor(df_gene_io$Group, levels = c("D0", "D1", "D7", "D30"))
wds.bas.summary <- df_gene_io %>% group_by(Group) %>% 
  summarise(sd = sd(NormExp, na.rm = TRUE), NormExp = mean(NormExp))
wds.bas.summary

my_comparisons <- list(c("NS", "VU"))
p2.wds <- ggplot(df_gene_io, aes(Group, NormExp)) +
  geom_col(data = wds.bas.summary, fill = NA, color = "black") +
  geom_jitter( position = position_jitter(0.1), color = "black") + 
  geom_errorbar( aes(ymin = NormExp-sd, ymax = NormExp+sd), 
                 data = wds.bas.summary, width = 0.2) +
  #ylab("Fold change of gene expression compared to skin") +
  ylab("Fold change of postive cell proportions compared to skin") +
  #scale_y_continuous(limits = c(0,14),breaks = c(0, 4, 8,12,16)) +
  theme_classic()
p2.wds + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test")



# filter out the lowly expressed genes
genes_all_f <- genes_all %>% filter(FOSL1 > 0.1)
table(genes_all_f$newCellType, genes_all_f$Condition)

fosl1exp <- table(genes_all_f$newCellType, genes_all_f$Sample) %>% as.data.frame()

table(inteData$newCellType, inteData$Sample)
table(inteData$Condition)
smTotal <- table(inteData$Sample) %>% as.data.frame()
table(inteData$newCellType, inteData$Sample)

# add the total number into fosl1exp
fosl1_tol <- fosl1exp %>% left_join(., smTotal, by=c("Var2" = "Var1")) %>% 
  setNames(c("CellType", "Sample", "Freq", "Total")) %>% 
  mutate(Group = gsub("H[0-9]$", "H", Sample)) %>% 
  mutate(Group = gsub("^PWH[0-9]{2}", "", Group)) %>% 
  filter(Group %in% c("H", "DFU_H", "DFU_NH")) %>% mutate(Prop = Freq / Total)

# Bas_mig
## First calculate the mean value of skin
skin.avg <- fosl1_tol %>% filter(CellType %in% c("Spi_mig")) %>% group_by(Group) %>% summarize(mean(Prop))
skin.avg
wds.bas <- fosl1_tol %>% filter(CellType %in% c("Spi_mig")) %>% mutate(NormProp = Prop /skin.avg$`mean(Prop)`[3])
wds.bas$Group <- factor(wds.bas$Group, levels = c("H", "DFU_H", "DFU_NH"))
wds.bas.summary <- wds.bas %>% group_by(Group) %>% 
  summarise(sd = sd(NormProp, na.rm = TRUE), NormProp = mean(NormProp))
wds.bas.summary

p2.wds <- ggplot(wds.bas, aes(Group, NormProp)) +
  geom_col(data = wds.bas.summary, fill = NA, color = "black") +
  geom_jitter( position = position_jitter(0.1), color = "black") + 
  geom_errorbar( aes(ymin = NormProp-sd, ymax = NormProp+sd), 
                 data = wds.bas.summary, width = 0.2) +
  #scale_y_continuous(limits = c(0,14),breaks = c(0, 4, 8,12,16)) +
  theme_classic()
p1.wds + p2.wds
# add the significance using t.test
p1 <- p1.wds + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test")
p2 <- p2.wds + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test")
pdf("step4_FOSL1_positiveCells.pdf", useDingbats = FALSE, width = 4, height = 3)
p1 + p2
dev.off()

# add the significance of groups (Pair-wise comparisons)
# prepare the data for significance test
df.prop <- wds.bas %>% filter(Group %in% c("DFU_H", "DFU_NH"))
df.prop$Group <- factor(df.prop$Group, levels = c("DFU_H", "DFU_NH"))

# anova test
test.sum <- aov(NormProp ~ Group, data = df.prop)
summary(test.sum)

#t.test(NormProp ~ Group, data = df.prop)
#wilcox.test(NormProp ~ Group, data = df.prop)

# quasibinomial test
test.quasi = glm(formula = Prop ~ Group, data = df.prop, family=quasibinomial)
print(summary(test.quasi))
anova(test.quasi, test = "LRT")$`Pr(>Chi)`[2]

