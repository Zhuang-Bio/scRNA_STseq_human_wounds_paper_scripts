library(tidyverse)
library(ggplot2)
library(patchwork)
#library(ggbreak)
library(ggpubr)

setwd("./")

########################################################
# Figure 6 cell proportions of acute and chronic wounds
# run the cell type one by one
alldf <- data.table::fread("allInte_CellType_Keratinocyte_cellproportion.txt")
alldf <- data.table::fread("allInte_CellType_Fibroblast_cellproportion.txt")
alldf <- data.table::fread("allInte_CellType_Myeloid_cellproportion.txt")
alldf <- data.table::fread("allInte_CellType_Lymphoid_cellproportion.txt")
alldf <- data.table::fread("allInte_CellType_Endothelium_cellproportion.txt")

table(alldf$newCellType)

ct <- unique(alldf$newCellType);ct
# ct <- ct[-c(2,3,9)];ct # Run KC
#ct <- ct[-c(2,4,9)];ct # Run FB
#ct <- ct[-c(5,7,8)];ct # Run Myeloid
ct <- ct[-c(5)];ct # Run Lymphoid


plots <- list()
for (i in seq_along(ct)) {
  # choose the cell types for plotting
  alldf_test <- alldf %>% dplyr::filter(newCellType == ct[i]) 
  
  ##################
  # acute wound
  skin_acute <- alldf_test %>% dplyr::filter(Project == "AcuteWound") %>% 
    group_by(Condition) %>% summarize(mean(Prop));skin_acute
  acw <-  alldf_test %>% dplyr::filter(Project == "AcuteWound") %>%
    mutate(NormProp = Prop /skin_acute$`mean(Prop)`[1])
  acw$Condition <- factor(acw$Condition, levels = c("Skin", "Wound1", "Wound7", "Wound30"))
  acw.summary <- acw %>% group_by(Condition) %>% 
    summarise(sd = sd(NormProp, na.rm = TRUE), NormProp = mean(NormProp))
  acw.summary
  
  #######
  # dfu
  skin_dfu <- alldf_test %>% dplyr::filter(Project == "DFU") %>% 
    group_by(Condition) %>% summarize(mean(Prop));skin_dfu
  dfutest <-  alldf_test %>% dplyr::filter(Project == "DFU") %>%
    mutate(NormProp = Prop /skin_dfu$`mean(Prop)`[3])
  dfutest$Condition <- factor(dfutest$Condition, levels = c("H", "DFU_H", "DFU_NH"))
  dfutest.summary <- dfutest %>% group_by(Condition) %>% 
    summarise(sd = sd(NormProp, na.rm = TRUE), NormProp = mean(NormProp))
  dfutest.summary
  
  #######
  # vu
  skin_vu <- alldf_test %>% dplyr::filter(Project == "VU") %>% 
    group_by(Condition) %>% summarize(mean(Prop));skin_vu
  vutest <-  alldf_test %>% dplyr::filter(Project == "VU") %>%
    mutate(NormProp = Prop /skin_vu$`mean(Prop)`[1])
  vutest$Condition <- factor(vutest$Condition, levels = c("NS", "VU"))
  vutest.summary <- vutest %>% group_by(Condition) %>% 
    summarise(sd = sd(NormProp, na.rm = TRUE), NormProp = mean(NormProp))
  vutest.summary
  
  scalevalue <- floor(max(c(acw$NormProp, dfutest$NormProp, vutest$NormProp)) + 1)
  
  p1.wds <- ggplot(acw, aes(Condition, NormProp, colour=Condition)) +
    geom_col(data = acw.summary, fill = NA, color="black") +
    geom_errorbar( aes(ymin = NormProp-sd, ymax = NormProp+sd), 
                   data = acw.summary, width = 0.3, color="black") +
    geom_jitter( position = position_jitter(0.1), size=1) + 
    scale_color_manual(values = c('#a6cee3','#ff7f00','#af8dc3','#33a02c')) +
    scale_y_continuous(limits = c(0,scalevalue)) +
    theme_classic() +
    theme(legend.position = "none")
  
  p1.dfu <- ggplot(dfutest, aes(Condition, NormProp, colour=Condition)) +
    geom_col(data = dfutest.summary, fill = NA, color="black") + #, color = "black"
    geom_errorbar( aes(ymin = NormProp-sd, ymax = NormProp+sd), 
                   data = dfutest.summary, width = 0.3, color="black") +
    geom_jitter( position = position_jitter(0.1), size=1) + #, color = "black"
    scale_color_manual(values = c('#a6cee3','#a1d76a','#e9a3c9')) +
    scale_y_continuous(limits = c(0,scalevalue)) +
    theme_classic() +
    theme(legend.position = "none")
  
  p1.vu <- ggplot(vutest, aes(Condition, NormProp, colour=Condition)) +
    geom_col(data = vutest.summary, fill = NA, color="black") +
    geom_errorbar( aes(ymin = NormProp-sd, ymax = NormProp+sd), 
                   data = vutest.summary, width = 0.3, color="black") +
    geom_jitter( position = position_jitter(0.1), size=1) + 
    scale_color_manual(values = c('#a6cee3','#fdbf6f')) + 
    scale_y_continuous(limits = c(0,scalevalue)) +
    theme_classic() +
    theme(legend.position = "none")
  
  # combine the plots
  #pdf("Result_CellProportion_allInte_KC_Bas_mig.pdf", useDingbats = FALSE, width = 4, height = 2)
  #pdf("Result_CellProportion_allInte_KC_Spi_mig.pdf", useDingbats = FALSE, width = 4, height = 2)
  plots[[i]] <- (p1.wds + p1.dfu + p1.vu) + plot_layout(widths = c(2,1.5,1)) + labs(caption = ct[i])
  #dev.off()
}

pdf("supple_chronicwounds_cellproportions_Endo.pdf", useDingbats = F, width = 16, height = 6)
wrap_plots(plots, ncol = 4)
dev.off()

#output figures
#plots[[7]] <- mki67kc
# manually add the MKI67 KCs
# then plot together with all KCs

#################################################################
# MKI67 positive cells proportion in all integrated dataset of KC
library(scCustomize)
inteData <- readRDS("allacuteWound_DFU_VU_inteKC_redu.rds")
uplot <- DimPlot(inteData, group.by = "newCellType", label = T) + NoLegend() + NoAxes()

FeaturePlot_scCustom(inteData, features = c("MKI67"), num_columns = 5, split.by = "Condition", combine = T)

FeaturePlot_scCustom(inteData, features = c("MKI67"), num_columns = 10, split.by = "orig.ident", combine = T)


genes_all <- FetchData(inteData, vars = c("Condition", "orig.ident", "Project", "newCellType",
                                          "MKI67"), slot = "data")
genes_all_f <- genes_all %>% filter(MKI67 > 0)

sm.tol <- table(inteData$orig.ident) %>% as.data.frame()

MKI67 <- table(genes_all_f$orig.ident, genes_all_f$newCellType) %>% as.data.frame() %>% 
  setNames(c("orig.ident", "newCellType", "Freq")) %>% 
  left_join(., alldf[,c(1,2,5,8,10)], by=c("orig.ident"="orig.ident", "newCellType"="newCellType")) %>% 
  na.omit() %>% 
  dplyr::mutate(Prop = Freq/newTotal)


allexp=MKI67
## First calculate the mean value of skin
##################
# acute wound
skin_acute <- allexp %>% dplyr::filter(Project == "AcuteWound") %>% 
  group_by(orig.ident) %>% mutate(Prop_tot=sum(Prop)) %>% ungroup() %>% 
  group_by(Condition) %>% summarize(mean(Prop_tot)) %>% ungroup();skin_acute

acw <-  allexp %>% dplyr::filter(Project == "AcuteWound") %>%
  group_by(orig.ident) %>% mutate(Prop_tot=sum(Prop)) %>% ungroup() %>% 
  mutate(NormProp = Prop_tot /skin_acute$`mean(Prop_tot)`[1]) %>% distinct(Condition, NormProp)

acw$Condition <- factor(acw$Condition, levels = c("Skin", "Wound1", "Wound7", "Wound30"))
acw.summary <- acw %>% group_by(Condition) %>% 
  summarise(sd = sd(NormProp, na.rm = TRUE), NormProp = mean(NormProp))
acw.summary

p1.wds <- ggplot(acw, aes(Condition, log2(NormProp+1), colour=Condition)) +
  geom_col(data = acw.summary, fill = NA, color="black") +
  geom_errorbar( aes(ymin = log2(NormProp-sd+1), ymax = log2(NormProp+sd+1)), 
                 data = acw.summary, width = 0.3, color="black") +
  geom_jitter( position = position_jitter(0.1), size=1) + 
  scale_color_manual(values = c('#a6cee3','#ff7f00','#af8dc3','#33a02c')) +
  scale_y_continuous(limits = c(0,3),breaks = c(0,0.5,1,1.5,2,2.5,3)) +
  theme_classic() +
  theme(legend.position = "none")
p1.wds

##################
# DFU
skin_dfu <- allexp %>% dplyr::filter(Project == "DFU") %>% 
  group_by(orig.ident) %>% mutate(Prop_tot=sum(Prop)) %>% ungroup() %>% 
  group_by(Condition) %>% summarize(mean(Prop_tot)) %>% ungroup();skin_dfu

dfutest <-  allexp %>% dplyr::filter(Project == "DFU") %>%
  group_by(orig.ident) %>% mutate(Prop_tot=sum(Prop)) %>% ungroup() %>% 
  mutate(NormProp = Prop_tot /skin_dfu$`mean(Prop_tot)`[3]) %>% distinct(Condition, NormProp)

dfutest$Condition <- factor(dfutest$Condition, levels = c("H", "DFU_H", "DFU_NH"))
dfutest.summary <- dfutest %>% group_by(Condition) %>% 
  summarise(sd = sd(NormProp, na.rm = TRUE), NormProp = mean(NormProp))
dfutest.summary

p1.dfu <- ggplot(dfutest, aes(Condition, log2(NormProp+1), colour=Condition)) +
  geom_col(data = dfutest.summary, fill = NA, color = "black") + #
  geom_errorbar( aes(ymin = log2(NormProp-sd+1), ymax = log2(NormProp+sd+1)), 
                 data = dfutest.summary, width = 0.3, color="black") +
  geom_jitter( position = position_jitter(0.1), size=1) + 
  scale_color_manual(values = c('#a6cee3','#a1d76a','#e9a3c9')) +
  scale_y_continuous(limits = c(0,3),breaks = c(0,0.5,1,1.5,2,2.5,3)) +
  theme_classic() +
  theme(legend.position = "none")
p1.dfu

##################
# VU
skin_vu <- allexp %>% dplyr::filter(Project == "VU") %>% 
  group_by(orig.ident) %>% mutate(Prop_tot=sum(Prop)) %>% ungroup() %>% 
  group_by(Condition) %>% summarize(mean(Prop_tot)) %>% ungroup();skin_vu

vutest <-  allexp %>% dplyr::filter(Project == "VU") %>%
  group_by(orig.ident) %>% mutate(Prop_tot=sum(Prop)) %>% ungroup() %>% 
  mutate(NormProp = Prop_tot /skin_vu$`mean(Prop_tot)`[1]) %>% distinct(Condition, NormProp)

vutest$Condition <- factor(vutest$Condition, levels = c("NS", "VU"))
vutest.summary <- vutest %>% group_by(Condition) %>% 
  summarise(sd = sd(NormProp, na.rm = TRUE), NormProp = mean(NormProp))
vutest.summary

p1.vu <- ggplot(vutest, aes(Condition, log2(NormProp+1), colour=Condition)) +
  geom_col(data = vutest.summary, fill = NA, color = "black") + #
  geom_errorbar( aes(ymin = log2(NormProp-sd+1), ymax = log2(NormProp+sd+1)), 
                 data = vutest.summary, width = 0.3, color="black") +
  geom_jitter( position = position_jitter(0.1), size=1) + 
  scale_color_manual(values = c('#a6cee3','#fdbf6f')) + 
  scale_y_continuous(limits = c(0,3),breaks = c(0,0.5,1,1.5,2,2.5,3)) +
  theme_classic() +
  theme(legend.position = "none")
p1.vu


mki67kc <- p1.wds + p1.dfu + p1.vu + plot_layout(widths = c(2,1.5,1))


#############################################################
# Figure 6i EGF and HGF signaling in acute and chonic wounds
inteData <- readRDS("allacuteWound_DFU_VU_integrated_reduced.rds")
# Make a dotplot show the fold changes of Wounds compared to Skin
egf <- c("TGFA", "AREG", "HBEGF", "EREG", "EGFR", "ERBB2")
hgf <- c("HGF", "MET")

otherLigand <- c(egf, hgf)
exp <- AggregateExpression(inteData, features = otherLigand, group.by = "Condition")$RNA %>% 
  as.data.frame() %>% 
  mutate(Wound1_FC = Wound1/Skin, Wound7_FC= Wound7/Skin, Wound30_FC=Wound30/Skin,
         DFU_H_FC = DFU_H/H, DFU_NH_FC= DFU_NH/H,
         VU_FC = VU/NS)

allegf_exp <- exp %>%
  rownames_to_column(var = "ID") %>% dplyr::select(1,11:16) %>% 
  pivot_longer(!ID, names_to = "Group", values_to = "FoldChange")
allegf_exp$Group <- gsub("_FC", "", allegf_exp$Group)
allegf_exp$Group <- factor(allegf_exp$Group, levels = c("Wound1", "Wound7", "Wound30", "DFU_H", "DFU_NH", "VU"))
allegf_exp$ID <- factor(allegf_exp$ID, levels = rev(otherLigand))

# Do the scale or not
## Setting the color to blue while values are lower than 1 (Down-regulated)
## red colors while values are higher than 1 (Up-regulated) compared to normal skin
ggplot(data = allegf_exp, aes(x=ID, y = Group, 
                              color = FoldChange, size = 10)) + 
  geom_point(size = 6) +
  scale_color_gradient2(limits = c(0,2), low = "blue", mid = "white", high = "red", midpoint = 1, na.value = "red") +
  #scale_size(range = c(0,2), limits = c(0,2))+
  theme_bw() + 
  ylab("") + 
  xlab("") +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1)) +
  coord_flip()
#pdf("Figure6i_EGF_signaling_DFU_VU.pdf", useDingbats = F, width = 4, height = 4)
#dev.off()

