library(tidyverse)
library(ggplot2)
library(patchwork)
library(ggbreak)
library(ggpubr)
library(Seurat)
library(scCustomize)

setwd("./Fig2_3_cellProportion_exp")

##################################
####---- Load the ST data ----####
hswound.STcom <- readRDS("/Users/zhuliu/Desktop/sc st/ST plotting scripts/Seurat_STseq_integrate_update.rds")
VlnPlot(hswound.STcom, features = "FOSL1", group.by = "Condition")
AverageExpression(hswound.STcom, assays = "Spatial", features = c("CXCL1"),
                  group.by = "Condition")$Spatial

cxcl1_st <- AverageExpression(hswound.STcom, assays = "Spatial", features = c("CXCL1"),
                              group.by = "Sample_name")$Spatial
cxcl1_st <- cxcl1_st %>% t()
cxcl1_st <- cxcl1_st %>% as.data.frame() %>% 
  rownames_to_column(var = "Sample") %>% setNames(c("Sample", "avgExp"))
data.table::fwrite(cxcl1_st, file = "expression_CXCL1_ST.txt", sep = "\t") 


## DFU data
inteData <- readRDS("/Users/zhuliu/Desktop/sc st/compareStudies/Diabetic foot ulcers/step1_diabetic_foot_ulcer_skin_integrated_harmony.rds")
FeaturePlot(inteData, features = c("CXCL1"))
VlnPlot(inteData, features = "CXCL1", group.by = "Condition")
AverageExpression(inteData, assays = "RNA", features = c("CXCL1"),
                  group.by = "Condition")$RNA

cxcl1_st <- AverageExpression(inteData, assays = "RNA", features = c("CXCL1"),
                              group.by = "Donor1")$RNA
cxcl1_st <- cxcl1_st %>% t()
cxcl1_st <- cxcl1_st %>% as.data.frame() %>% 
  rownames_to_column(var = "Sample") %>% setNames(c("Sample", "avgExp"))
data.table::fwrite(cxcl1_st, file = "expression_CXCL1_DFU.txt", sep = "\t") 

## Acute wounds
inteData <- readRDS("/Users/zhuliu/Desktop/sc st/plotting scripts/Figure 1/allcombined_wounds_newAnnotation.rds")
FeaturePlot(inteData, features = c("CXCL1"))
VlnPlot(inteData, features = "CXCL1", group.by = "Condition")
AverageExpression(inteData, assays = "RNA", features = c("CXCL1"),
                  group.by = "orig.ident")$RNA

cxcl1_st <- AverageExpression(inteData, assays = "RNA", features = c("CXCL1"),
                              group.by = "orig.ident")$RNA
cxcl1_st <- cxcl1_st %>% t()
cxcl1_st <- cxcl1_st %>% as.data.frame() %>% 
  rownames_to_column(var = "Sample") %>% setNames(c("Sample", "avgExp"))
data.table::fwrite(cxcl1_st, file = "expression_CXCL1_AcuteWounds.txt", sep = "\t") 


## VU
#inteData <- readRDS("/Users/zhuliu/Desktop/sc st/Versions_SC_ST_Project/V1_Figures_230120/ComparedStudies/VU_data/D04_annotatedData.rds")
#inteData <- readRDS("/Users/zhuliu/Desktop/sc st/Versions_SC_ST_Project/V1_Figures_230120/ComparedStudies/VU_data/D08_Annotated_VUdata_addVU6.rds")
inteData <- readRDS("/Volumes/zhuliu/Groups/Ning Xu Landén/Zhuang Liu/proj_VU_fromLihua/s02_VU_Inte/D07_Harmony_delOutlierCls_anno.rds")
FeaturePlot(inteData, features = c("CXCL1"))
VlnPlot(inteData, features = "CXCL1", group.by = "Condition")
library(ggpubr)
library(scCustomize)
my_comparisons <- list( c("NS", "VU"))
(VlnPlot_scCustom(seurat_object = inteData, features = "CXCL1", 
                  group.by = "Condition", 
                  plot_median = TRUE,y.max = 6) & NoLegend()) +
  stat_compare_means(comparisons = my_comparisons)


AverageExpression(inteData, assays = "RNA", features = c("CXCL1"),
                  group.by = "Condition")$RNA

cxcl1_st <- AverageExpression(inteData, assays = "RNA", features = c("CXCL1"),
                              group.by = "orig.ident")
cxcl1_st <- cxcl1_st$RNA[,] %>% as.data.frame()
cxcl1_st <- cxcl1_st %>% as.data.frame() %>% 
  rownames_to_column(var = "Sample") %>% setNames(c("Sample", "avgExp"))
data.table::fwrite(cxcl1_st, file = "expression_CXCL1_VU_new.txt", sep = "\t")


## VU (add two more healthy skin)
inteData <- readRDS("/Users/zhuliu/Downloads/D06_AfterHarmony.rds")
FeaturePlot(inteData, features = c("CXCL1"), split.by = "orig.ident")
VlnPlot(inteData, features = "CXCL1", group.by = "Condition")
library(ggpubr)
library(scCustomize)
my_comparisons <- list( c("NS", "VU"))
(VlnPlot_scCustom(seurat_object = inteData, features = "CXCL1", 
                  group.by = "Condition", 
                  plot_median = TRUE,y.max = 6) & NoLegend()) +
  stat_compare_means(comparisons = my_comparisons)

AverageExpression(inteData, assays = "RNA", features = c("CXCL1"),
                  group.by = "Condition")$RNA

cxcl1_st <- AverageExpression(inteData, assays = "RNA", features = c("CXCL1"),
                              group.by = "orig.ident")
cxcl1_st <- cxcl1_st$RNA[,] %>% as.data.frame()
cxcl1_st <- cxcl1_st %>% as.data.frame() %>% 
  rownames_to_column(var = "Sample") %>% setNames(c("Sample", "avgExp"))
data.table::fwrite(cxcl1_st, file = "expression_CXCL1_VU_new_240118.txt", sep = "\t")

cxcl1_st <- cxcl1_st %>% 
  filter(!Sample %in% "VU4") %>% 
  mutate(Sample = gsub("NS[0-9]{2}$", "NS", Sample)) %>% 
  mutate(Sample = gsub("VU[0-9]$", "VU", Sample)) %>% 
  mutate(Sample = gsub("NS[0-9]{2}[A-Z]", "NS", Sample))

my_comparisons <- list( c("VU", "NS") )
vu_plot <- ggboxplot(cxcl1_st, x = "Sample", y = "avgExp",
                     color = "Sample", palette = "npg", add = "jitter") +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") +
  scale_y_continuous(limits = c(0,7),breaks = c(0,1,2,3,4,5,6,7))
vu_plot


## Bulk RNA-seq data
bulk <- readRDS("/Users/zhuliu/Desktop/KI_Derma_zhuliu/05-Archive-Projects/miRNA/a_Manuscript/ShinyApp_miRNA_Profiling/miRNA_Xulab/data/mRNA_exp.rds")
bulk_cxcl1 <- bulk %>% filter(GeneSymbol %in% c("CXCL1")) %>% 
  select(-2:-3) %>% column_to_rownames(var = "GeneSymbol") %>% t() %>% 
  as.data.frame() %>% rownames_to_column(var = "Sample") %>% 
  setNames(c("Sample", "avgExp"))
data.table::fwrite(bulk_cxcl1, file = "expression_CXCL1_bulk.txt", sep = "\t")   


## DFU GEO data download
library(GEOquery)
gse <- getGEO("GSE80178", GSEMatrix = TRUE)
show(gse)
exprSet = exprs(gse[[1]]) %>% as.data.frame() %>%  rownames_to_column(var = "ID")
pdata = pData(gse[[1]]) %>% as.data.frame()
gpl <- gse$GSE80178_series_matrix.txt.gz@featureData@data
colnames(exprSet)[-1] <- c(paste0("DFU", 1:6), paste0("DFU_Skin", 1:3), paste0("Skin", 1:3))
cxcl1_fosl1 <- exprSet %>% filter(ID %in% c("16967794", "16740630")) %>% 
  column_to_rownames(var = "ID") %>% t() %>% as.data.frame() %>% 
  rownames_to_column(var = "Sample") %>% 
  setNames(c("Sample", "FOSL1", "CXCL1"))
data.table::fwrite(cxcl1_fosl1, file = "expression_CXCL1_microarray_FOSL1.txt", sep = "\t")


## DFU scRNA new Data
library(SeuratDisk)
newDFU <- LoadH5Seurat("/Users/zhuliu/Desktop/sc st/Versions_SC_ST_Project/V1_Figures_230120/Fig2_3_cellProportion_exp/DFU_scData_GSE223964/GSE223964_allcell.h5seurat")
DimPlot(newDFU, group.by = "Cell.Type")
FeaturePlot(newDFU, features = c("CXCL1", "FOSL1"), split.by = "condition", 
            cols = c("grey90", "red"))
VlnPlot(newDFU, features = c("CXCL1", "FOSL1"), group.by = "Cell.Type",
        split.by = "condition")
AverageExpression(newDFU, assays = "RNA", features = c("CXCL1"),
                  group.by = "condition")$RNA

cxcl1_st <- AverageExpression(newDFU, assays = "RNA", features = c("CXCL1"),
                              group.by = "sample")$RNA
cxcl1_st <- cxcl1_st %>% t()
cxcl1_st <- cxcl1_st %>% as.data.frame() %>% 
  rownames_to_column(var = "Sample") %>% setNames(c("Sample", "avgExp"))
data.table::fwrite(cxcl1_st, file = "expression_CXCL1_newDFUdata.txt", sep = "\t") 



##################################
##-- Plotting GENE EXPRESSION --##
# Acute wound
df <- data.table::fread("expression_CXCL1_AcuteWounds.txt") %>% 
  mutate(Sample = gsub("^PWH[0-9]{2}", "", Sample))
avg.skin <- df %>% group_by(Sample) %>% 
  summarize(mean(avgExp)) %>% ungroup();avg.skin
df.f <- df %>% mutate(norExp = avgExp/avg.skin$`mean(avgExp)`[1])
df.f$Sample <- factor(df.f$Sample, levels = c("D0", "D1", "D7", "D30"))
alldf.summary <- df.f %>% group_by(Sample) %>% 
  summarise(sd = sd(norExp, na.rm = TRUE), norExp = mean(norExp));alldf.summary

p1 <- ggplot(df.f, aes(Sample, log2(norExp+1))) +
  geom_col(data = alldf.summary, fill = NA, color = "black") +
  geom_jitter( position = position_jitter(0.1), color = "black") + 
  geom_errorbar(aes(ymin = log2(norExp-sd+1), ymax = log2(norExp+sd+1)), 
                data = alldf.summary, width = 0.2) +
  scale_y_continuous(limits = c(-0.1,5),breaks = c(0,1,2,3,4,5)) +
  ylab("Fold change of CXCL1 expression compared to skin") +
  xlab("scRNA acute wound") +
  geom_hline(yintercept=log2(1+1), linetype="dashed", color = "red") +
  theme_classic()
p1
my_comparisons <- list( c("Wound1", "Skin"), c("Wound7", "Skin"), c("Wound30", "Skin") )
#pdf("qRTPCR_CXCL1_exp_sig.pdf", useDingbats = F, width = 4, height = 4)
ggboxplot(df, x = "Sample", y = "NormProp",
          color = "Sample", palette = "npg", add = "jitter") +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")


# Spatial Acute wound
df <- data.table::fread("expression_CXCL1_ST.txt") %>% 
  mutate(Sample = Sample) %>% 
  separate(col = Sample, into = c("Donor1", "Sample"), sep = "_")
#df <- cxcl1_st %>% 
#  mutate(Sample = Sample) %>% 
#  separate(col = Sample, into = c("Donor1", "Sample"), sep = "_")
avg.skin <- df %>% group_by(Sample) %>% 
  summarize(mean(avgExp)) %>% ungroup();avg.skin
df.f <- df %>% mutate(norExp = avgExp/avg.skin$`mean(avgExp)`[1])
df.f$Sample <- factor(df.f$Sample, levels = c("Skin", "Wound1", "Wound7", "Wound30"))
alldf.summary <- df.f %>% group_by(Sample) %>% 
  summarise(sd = sd(norExp, na.rm = TRUE), norExp = mean(norExp));alldf.summary

p2 <- ggplot(df.f, aes(Sample, log2(norExp+1))) +
  geom_col(data = alldf.summary, fill = NA, color = "black") +
  geom_jitter( position = position_jitter(0.1), color = "black") + 
  geom_errorbar(aes(ymin = log2(norExp-sd+1), ymax = log2(norExp+sd+1)), 
                data = alldf.summary, width = 0.2) +
  scale_y_continuous(limits = c(-0.1,5),breaks = c(0,1,2,3,4,5)) +
  ylab("Fold change of CXCL1 expression compared to skin") +
  xlab("STseq acute wound") +
  geom_hline(yintercept=log2(1+1), linetype="dashed", color = "red") +
  theme_classic()
p2
#pdf(file = "FOSL1_exp.pdf", useDingbats = F, width = 2, height = 3)
#dev.off()
skin.avg <- df %>% 
  group_by(Sample) %>% summarize(mean(avgExp));skin.avg
df <- df %>%
  mutate(NormProp = avgExp /skin.avg$`mean(avgExp)`[1])

wds.bas <- df %>% filter(Sample %in% c("Skin", "Wound7")) 
test.sum <- aov(NormProp ~ Sample, data = wds.bas)
summary(test.sum)
t.test(NormProp ~ Sample,  data = wds.bas)
wilcox.test(NormProp ~ Sample, data = wds.bas)
my_comparisons <- list( c("Wound1", "Skin"), c("Wound7", "Skin"), c("Wound30", "Skin") )
#pdf("qRTPCR_CXCL1_exp_sig.pdf", useDingbats = F, width = 4, height = 4)
st_plot <- ggboxplot(df, x = "Sample", y = "log2(NormProp+1)",
          color = "Sample", palette = "npg", add = "jitter") +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") +
  scale_y_continuous(limits = c(0,7),breaks = c(0,1,2,3,4,5,6,7))
st_plot

# bulk RNA-seq
df <- data.table::fread("expression_CXCL1_bulk.txt") %>% 
  mutate(Sample = gsub("[0-9]$", "", Sample)) %>% 
  mutate(Sample = gsub("_$", "", Sample))
avg.skin <- df %>% group_by(Sample) %>% 
  summarize(mean(avgExp)) %>% ungroup();avg.skin
df.f <- df %>% mutate(norExp = avgExp/avg.skin$`mean(avgExp)`[1])
df.f$Sample <- factor(df.f$Sample, levels = c("Skin", "Wound1", "Wound7", "VU"))
alldf.summary <- df.f %>% group_by(Sample) %>% 
  summarise(sd = sd(norExp, na.rm = TRUE), norExp = mean(norExp));alldf.summary

p3 <- ggplot(df.f, aes(Sample, log2(norExp+1))) +
  geom_col(data = alldf.summary, fill = NA, color = "black") +
  geom_jitter( position = position_jitter(0.1), color = "black") + 
  geom_errorbar(aes(ymin = log2(norExp-sd+1), ymax = log2(norExp+sd+1)), 
                data = alldf.summary, width = 0.2) +
  scale_y_continuous(limits = c(-0.1,10),breaks = c(0,1,2,3,4,5,6,7,8,9,10)) +
  ylab("Fold change of CXCL1 expression compared to skin") +
  xlab("BulkSeq acute wound and VU") +
  geom_hline(yintercept=log2(1+1), linetype="dashed", color = "red") +
  theme_classic()
p3
my_comparisons <- list( c("Wound1", "Skin"), c("Wound7", "Skin"), c("VU", "Skin"), c("VU", "Wound1"), c("VU", "Wound7") )
pdf("CXCL1_exp_bulkseq.pdf", useDingbats = F, width = 3, height = 3)
ggboxplot(df, x = "Sample", y = "log2(avgExp+1)", 
          color = "Sample", palette = "npg", add = "jitter") +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") +
  scale_y_continuous(limits = c(0,10),breaks = c(0,2,4,6,8,10))
dev.off()

# VU scRNA-seq (Pay attention to the new Data)
df <- data.table::fread("expression_CXCL1_VU_new_twoMore.txt") %>%  #expression_CXCL1_VU_new.txt
  mutate(Sample = gsub("NS[0-9]{2}$", "NS", Sample)) %>% 
  mutate(Sample = gsub("VU[0-9]$", "VU", Sample)) %>% 
  mutate(Sample = gsub("NS[0-9]{2}[A-Z]", "NS", Sample)) 
avg.skin <- df %>% group_by(Sample) %>% 
  summarize(mean(avgExp)) %>% ungroup();avg.skin
df.f <- df %>% mutate(norExp = avgExp/avg.skin$`mean(avgExp)`[1])
df.f$Sample <- factor(df.f$Sample, levels = c("NS", "VU"))
alldf.summary <- df.f %>% group_by(Sample) %>% 
  summarise(sd = sd(norExp, na.rm = TRUE), norExp = mean(norExp));alldf.summary

p4 <- ggplot(df.f, aes(Sample, log2(norExp+1))) +
  geom_col(data = alldf.summary, fill = NA, color = "black") +
  geom_jitter( position = position_jitter(0.1), color = "black") + 
  geom_errorbar(aes(ymin = log2(norExp-sd+1), ymax = log2(norExp+sd+1)), 
                data = alldf.summary, width = 0.2) +
  scale_y_continuous(limits = c(-0.2,5),breaks = c(0,1,2,3,4,5)) +
  ylab("Fold change of CXCL1 expression compared to skin") +
  xlab("scRNA VU") +
  geom_hline(yintercept=log2(1+1), linetype="dashed", color = "red") +
  theme_classic()
p4
my_comparisons <- list( c("VU", "NS") )
#pdf("qRTPCR_CXCL1_exp_sig.pdf", useDingbats = F, width = 4, height = 4)
vu_plot <- ggboxplot(df, x = "Sample", y = "avgExp",
          color = "Sample", palette = "npg", add = "jitter") +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") +
  scale_y_continuous(limits = c(0,7),breaks = c(0,1,2,3,4,5,6,7))
vu_plot

# DFU scRNA-seq
df <- data.table::fread("expression_CXCL1_DFU.txt") %>% 
  mutate(Sample = gsub("H[0-9]$", "H", Sample)) %>% 
  dplyr::slice(-8,-12)
avg.skin <- df %>% group_by(Sample) %>% 
  summarize(mean(avgExp)) %>% ungroup();avg.skin
df.f <- df %>% mutate(norExp = avgExp/avg.skin$`mean(avgExp)`[3])
df.f$Sample <- factor(df.f$Sample, levels = c("H", "DFU_H", "DFU_NH"))
alldf.summary <- df.f %>% group_by(Sample) %>% 
  summarise(sd = sd(norExp, na.rm = TRUE), norExp = mean(norExp));alldf.summary

p5 <- ggplot(df.f, aes(Sample, log2(norExp+1))) +
  geom_col(data = alldf.summary, fill = NA, color = "black") +
  geom_jitter( position = position_jitter(0.1), color = "black") + 
  geom_errorbar(aes(ymin = log2(norExp-sd+1), ymax = log2(norExp+sd+1)), 
                data = alldf.summary, width = 0.2) +
  scale_y_continuous(limits = c(-0.1,5),breaks = c(0,1,2,3,4,5)) +
  ylab("Fold change of CXCL1 expression compared to skin") +
  xlab("scRNA DFU") +
  geom_hline(yintercept=log2(1+1), linetype="dashed", color = "red") +
  theme_classic()
p5
df$Sample <- factor(df$Sample, levels = c("H", "DFU_H", "DFU_NH"))
my_comparisons <- list( c("DFU_H", "H"), c("DFU_NH", "H"), c("DFU_NH", "DFU_H") )
#pdf("qRTPCR_CXCL1_exp_sig.pdf", useDingbats = F, width = 4, height = 4)
dfu_plot <- ggboxplot(df, x = "Sample", y = "avgExp",
          color = "Sample", palette = "npg", add = "jitter") +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") +
  scale_y_continuous(limits = c(0,7),breaks = c(0,1,2,3,4,5,6,7))
dfu_plot

# newDFU scRNA-seq
df <- data.table::fread("expression_CXCL1_newDFUdata.txt") %>% 
  mutate(Sample = c(rep("DFU",5), rep("nonDFU",3)))
avg.skin <- df %>% group_by(Sample) %>% 
  summarize(mean(avgExp)) %>% ungroup();avg.skin
df.f <- df %>% mutate(norExp = avgExp/avg.skin$`mean(avgExp)`[2])
df.f$Sample <- factor(df.f$Sample, levels = c("nonDFU", "DFU"))
alldf.summary <- df.f %>% group_by(Sample) %>% 
  summarise(sd = sd(norExp, na.rm = TRUE), norExp = mean(norExp));alldf.summary

p6 <- ggplot(df.f, aes(Sample, log2(norExp+1))) +
  geom_col(data = alldf.summary, fill = NA, color = "black") +
  geom_jitter( position = position_jitter(0.1), color = "black") + 
  geom_errorbar(aes(ymin = log2(norExp-sd+1), ymax = log2(norExp+sd+1)), 
                data = alldf.summary, width = 0.2) +
  scale_y_continuous(limits = c(-0.1,5),breaks = c(0,1,2,3,4,5)) +
  ylab("Fold change of CXCL1 expression compared to skin") +
  xlab("scRNA DFU newData") +
  geom_hline(yintercept=log2(1+1), linetype="dashed", color = "red") +
  theme_classic()
p6
my_comparisons <- list( c("DFU", "nonDFU"))
#pdf("qRTPCR_CXCL1_exp_sig.pdf", useDingbats = F, width = 4, height = 4)
ggboxplot(df, x = "Sample", y = "avgExp",
          color = "Sample", palette = "npg", add = "jitter") +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")


# DFU microarray
df <- data.table::fread("expression_CXCL1_microarray_FOSL1.txt") %>% 
  mutate(Sample = gsub("[0-9]$", "", Sample)) %>% 
  dplyr::select(1,3) %>% setNames(c("Sample", "avgExp")) %>% 
  mutate(Sample=gsub("DFU_Skin", "Skin", Sample))
avg.skin <- df %>% group_by(Sample) %>% 
  summarize(mean(avgExp)) %>% ungroup();avg.skin
df.f <- df %>% mutate(norExp = avgExp/avg.skin$`mean(avgExp)`[2])
df.f$Sample <- factor(df.f$Sample, levels = c("Skin", "DFU"))#"DFU_Skin", 
alldf.summary <- df.f %>% group_by(Sample) %>% 
  summarise(sd = sd(norExp, na.rm = TRUE), norExp = mean(norExp));alldf.summary

p7 <- ggplot(df.f, aes(Sample, log2(norExp+1))) +
  geom_col(data = alldf.summary, fill = NA, color = "black") +
  geom_jitter( position = position_jitter(0.1), color = "black") + 
  geom_errorbar(aes(ymin = log2(norExp-sd+1), ymax = log2(norExp+sd+1)), 
                data = alldf.summary, width = 0.2) +
  #scale_y_continuous(limits = c(-0.1,5),breaks = c(0,1,2,3,4,5)) +
  ylab("Fold change of CXCL1 expression compared to skin") +
  xlab("microarray DFU") +
  geom_hline(yintercept=log2(1+1), linetype="dashed", color = "red") +
  theme_classic()
p7
df$Sample <- factor(df$Sample, levels = c("Skin", "DFU")) #"DFU_Skin", 
my_comparisons <- list( c("DFU", "Skin")) #, c("DFU", "DFU_Skin")
#pdf("qRTPCR_CXCL1_exp_sig.pdf", useDingbats = F, width = 4, height = 4)
ggboxplot(df.f, x = "Sample", y = "norExp",
          color = "Sample", palette = "npg", add = "jitter") +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") +
  #scale_y_continuous(limits = c(0,1.5),breaks = c(0,0.3,0.6,0.9,1.2,1.5))
  scale_y_continuous(limits = c(0,4),breaks = c(0,1,2,3,4))
pdf("CXCL1_exp_sc_DFU_microarray2.pdf", useDingbats = FALSE, width = 2, height = 3)
dev.off()

pdf("CXCL1_exp_st_sc_UP_removeVU4.pdf", useDingbats = FALSE, width = 6, height = 3)
st_plot + dfu_plot + vu_plot + plot_layout(widths = c(2,1.5,1))
dev.off()

p1+p2+p3+p4+p5+p6+p7 + plot_layout(widths = c(2,2,2,1,1.5,1,1.5), ncol = 7) +
  plot_annotation("Plotting based on log2(Fold Change+1)")
# reorder the figure
#pdf(file = "CXCL1_EXP.pdf", useDingbats = FALSE, width = 9, height = 4)
p1+p2+p5+p6+p4+p7+p3 + plot_layout(widths = c(2,2,1.5,1,1,1.5,2), ncol = 7) +
  plot_annotation("Plotting based on log2(Fold Change+1)")
#dev.off()
