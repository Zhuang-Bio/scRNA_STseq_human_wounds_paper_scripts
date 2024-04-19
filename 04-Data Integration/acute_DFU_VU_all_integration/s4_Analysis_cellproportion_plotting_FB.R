library(tidyverse)
library(ggplot2)
library(patchwork)
library(ggbreak)
library(ggpubr)

setwd("./cellproportion_analysis_results/")

# run the FB
alldf <- data.table::fread("allInte_CellType_Fibroblast_cellproportion.txt")
table(alldf$newCellType)
#FB_ELN_SFRP4 FB_I_POSTN_COL11A1  FB_I_POSTN_COL4A1   FB_I_POSTN_MMP11    FB_I_SFRP4_COMP 
#41                 41                 41                 41                 41 
#FB_II_APOD_ITM2A   FB_II_APOE_CCL19    FB_III_ELN_LEPR          FB_prolif    FB_SFRP1_CRABP1 
#41                 41                 41                 25                 41 

# choose the cell types for plotting
alldf_test <- alldf %>% dplyr::filter(newCellType == "FB_I_POSTN_COL11A1") 
alldf_test <- alldf %>% dplyr::filter(newCellType == "FB_I_POSTN_MMP11") 

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

p1.wds <- ggplot(acw, aes(Condition, NormProp, colour=Condition)) +
  geom_col(data = acw.summary, fill = NA, color="black") +
  geom_errorbar( aes(ymin = NormProp-sd, ymax = NormProp+sd), 
                 data = acw.summary, width = 0.3, color="black") +
  geom_jitter( position = position_jitter(0.1), size=1) + 
  scale_color_manual(values = c('#a6cee3','#ff7f00','#af8dc3','#33a02c')) +
  #scale_y_continuous(limits = c(0,10),breaks = c(0, 2,4,6,8,10)) +
  #scale_y_continuous(limits = c(0,5),breaks = c(0, 1,2,3,4,5,6)) +
  theme_classic() +
  theme(legend.position = "none")
p1.wds

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

p1.dfu <- ggplot(dfutest, aes(Condition, NormProp, colour=Condition)) +
  geom_col(data = dfutest.summary, fill = NA, color="black") + #, color = "black"
  geom_errorbar( aes(ymin = NormProp-sd, ymax = NormProp+sd), 
                 data = dfutest.summary, width = 0.3, color="black") +
  geom_jitter( position = position_jitter(0.1), size=1) + #, color = "black"
  scale_color_manual(values = c('#a6cee3','#a1d76a','#e9a3c9')) +
  #scale_y_continuous(limits = c(0,10),breaks = c(0, 2,4,6,8,10)) +
  #scale_y_continuous(limits = c(0,5),breaks = c(0, 1,2,3,4,5,6)) +
  theme_classic() +
  theme(legend.position = "none")
p1.dfu

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

p1.vu <- ggplot(vutest, aes(Condition, NormProp, colour=Condition)) +
  geom_col(data = vutest.summary, fill = NA, color="black") +
  geom_errorbar( aes(ymin = NormProp-sd, ymax = NormProp+sd), 
                 data = vutest.summary, width = 0.3, color="black") +
  geom_jitter( position = position_jitter(0.1), size=1) + 
  scale_color_manual(values = c('#a6cee3','#fdbf6f')) + 
  #scale_y_continuous(limits = c(0,10),breaks = c(0, 2,4,6,8,10)) +
  #scale_y_continuous(limits = c(0,5),breaks = c(0, 1,2,3,4,5,6)) +
  theme_classic() +
  theme(legend.position = "none")
p1.vu

# combine the plots
#pdf("Result_CellProportion_allInte_FB_POSTN_COL11A1.pdf", useDingbats = FALSE, width = 4, height = 2)
#pdf("Result_CellProportion_allInte_FB_POSTN_MMP11.pdf", useDingbats = FALSE, width = 4, height = 2)
p1.wds + p1.dfu + p1.vu + plot_layout(widths = c(2,1.5,1))
dev.off()


##########################################################
# MKI67 positive cells proportion in all integrated dataset
library(scCustomize)
alldf <- data.table::fread("allInte_CellType_Fibroblast_cellproportion.txt")

inteData <- readRDS("../02-fibroblast/allacuteWound_DFU_VU_inteFB_redu.rds")
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
  #scale_y_continuous(limits = c(0,7),breaks = c(0,1,2,3,4,5,6,7)) +
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
  #scale_y_continuous(limits = c(0,7),breaks = c(0,1,2,3,4,5,6,7)) +
  theme_classic() +
  theme(legend.position = "none")
p1.dfu

##################
# VU
skin_vu <- allexp %>% dplyr::filter(Project == "VU") %>% 
  dplyr::filter(! orig.ident %in% c("NS23")) %>% # weird expression of MKI67 in this sample
  group_by(orig.ident) %>% mutate(Prop_tot=sum(Prop)) %>% ungroup() %>% 
  group_by(Condition) %>% summarize(mean(Prop_tot)) %>% ungroup();skin_vu

vutest <-  allexp %>% dplyr::filter(Project == "VU") %>%
  dplyr::filter(! orig.ident %in% c("NS23")) %>% # weird expression of MKI67 in this sample
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
  #scale_y_continuous(limits = c(0,7),breaks = c(0,1,2,3,4,5,6,7)) +
  theme_classic() +
  theme(legend.position = "none")
p1.vu

#pdf("Result_CellProportion_allInte_FB_MKI67positive.pdf", useDingbats = FALSE, width = 4, height = 2)
p1.wds + p1.dfu + p1.vu + plot_layout(widths = c(2,1.5,1))
dev.off()


################################################
# Run the statistics for cell proportion analysis
# load the data
cellproportion_statis <- function(data_cp=NULL, project = "AcuteWound"){
  step3 <- data_cp %>% dplyr::filter(Project == project) %>% mutate(norProp = 100*Prop)
  
  ct.stat <- unique(step3$Project) %>% as.character();print(ct.stat)
  
  stat_sig <- list()
  for (i in seq_along(ct.stat)){
    ct_sel <- ct.stat[i]
    ct_por <- step3 %>% 
      #dplyr::filter(newCellType == ct_sel) %>% 
      arrange(Condition)
    
    stats_all <- function(GroupInfo = NULL, paired = TRUE){
      ct_por_gr <- ct_por %>% dplyr::filter(Condition %in% GroupInfo)
      # paired samples Wilcoxon test
      if(paired){
        #res <- wilcox.test(norProp ~ Condition, data = ct_por_gr, paired = TRUE)
        res <- t.test(norProp ~ Condition, data = ct_por_gr)
      }else{
        #res <- wilcox.test(norProp ~ Condition, data = ct_por_gr, paired = FALSE)
        res <- t.test(norProp ~ Condition, data = ct_por_gr)
      }
      print(res)
      wilcoxon_test = res$p.value
      # quasibinomial test
      test.quasi = glm(formula = Prop ~ Condition, data = ct_por_gr, family=quasibinomial)
      print(summary(test.quasi))
      quasibinomial=anova(test.quasi, test = "LRT")$`Pr(>Chi)`[2]
      return(c(wilcoxon_test, quasibinomial))
    }
    
    if(project == "AcuteWound"){
      res_all_1=stats_all(GroupInfo = c("Skin", "Wound1"))
      res_all_2=stats_all(GroupInfo = c("Skin", "Wound7"))
      res_all_3=stats_all(GroupInfo = c("Skin", "Wound30"))
      res_all_4=stats_all(GroupInfo = c("Wound1", "Wound7"))
      res_all_5=stats_all(GroupInfo = c("Wound1", "Wound30"))
      res_all_6=stats_all(GroupInfo = c("Wound7", "Wound30"))
      stat_sig[[ct_sel]] <- append(stat_sig[[ct_sel]], res_all_1)
      stat_sig[[ct_sel]] <- append(stat_sig[[ct_sel]], res_all_2)
      stat_sig[[ct_sel]] <- append(stat_sig[[ct_sel]], res_all_3)
      stat_sig[[ct_sel]] <- append(stat_sig[[ct_sel]], res_all_4)
      stat_sig[[ct_sel]] <- append(stat_sig[[ct_sel]], res_all_5)
      stat_sig[[ct_sel]] <- append(stat_sig[[ct_sel]], res_all_6)
    }
    if(project == "DFU"){
      res_all_1=stats_all(GroupInfo = c("H", "DFU_H"), paired = FALSE)
      res_all_2=stats_all(GroupInfo = c("H", "DFU_NH"), paired = FALSE)
      res_all_3=stats_all(GroupInfo = c("DFU_NH", "DFU_H"), paired = FALSE)
      stat_sig[[ct_sel]] <- append(stat_sig[[ct_sel]], res_all_1)
      stat_sig[[ct_sel]] <- append(stat_sig[[ct_sel]], res_all_2)
      stat_sig[[ct_sel]] <- append(stat_sig[[ct_sel]], res_all_3)
    }
    if(project == "VU"){
      res_all_1=stats_all(GroupInfo = c("NS", "VU"), paired = FALSE)
      stat_sig[[ct_sel]] <- append(stat_sig[[ct_sel]], res_all_1)
    }
  }
  # output the results
  if(project == "AcuteWound"){
    stat_sig_all <- do.call("rbind", stat_sig) %>% as.data.frame() %>% 
      rownames_to_column(var = "newCellType") %>% 
      setNames(c("newCellType", 
                 "Wilco_Skin_Wound1", "QuasiBino_Skin_Wound1", "Wilco_Skin_Wound7", "QuasiBino_Skin_Wound7", 
                 "Wilco_Skin_Wound30", "QuasiBino_Skin_Wound30", "Wilco_Wound1_Wound7", "QuasiBino_Wound1_Wound7",
                 "Wilco_Wound1_Wound30", "QuasiBino_Wound1_Wound30", "Wilco_Wound7_Wound30", "QuasiBino_Wound7_Wound30"))
  }
  if(project == "DFU"){
    stat_sig_all <- do.call("rbind", stat_sig) %>% as.data.frame() %>% 
      rownames_to_column(var = "newCellType") %>% 
      setNames(c("newCellType", "Wilco_H_DFU_H", "QuasiBino_H_DFU_H", 
                 "Wilco_H_DFU_NH", "QuasiBino_H_DFU_NH", "Wilco_DFU_H_DFU_NH", "QuasiBino_DFU_H_DFU_NH"))
  }
  if(project == "VU"){
    stat_sig_all <- do.call("rbind", stat_sig) %>% as.data.frame() %>% 
      rownames_to_column(var = "newCellType") %>% 
      setNames(c("newCellType",  "Wilco_NS_VU", "QuasiBino_NS_VU"))
  }
  return(stat_sig_all)
}

# run the functions
## Keratinocyte
ct_acute <- cellproportion_statis(data_cp = allexp, project = "AcuteWound")
ct_dfu <- cellproportion_statis(data_cp = allexp, project = "DFU")
ct_vu <- cellproportion_statis(data_cp = allexp, project = "VU")

exp_sig <- rbind(t(ct_acute), t(ct_dfu), t(ct_vu)) %>% as.data.frame() %>% rownames_to_column(var = "Group")
data.table::fwrite(exp_sig, "Result_CellProportion_allInte_FB_MKI67positive_sig.txt", sep = "\t")


AverageExpression(inteData, features = c("MKI67"), assays = "RNA", group.by = "Condition")


# FB-prolif log2FoldChange
alldf_test <- alldf %>% dplyr::filter(newCellType == "FB_prolif") 

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

p1.wds <- ggplot(acw, aes(Condition, log2(NormProp+1), colour=Condition)) +
  geom_col(data = acw.summary, fill = NA, color="black") +
  geom_errorbar( aes(ymin = log2(NormProp-sd+1), ymax = log2(NormProp+sd+1)), 
                 data = acw.summary, width = 0.3, color="black") +
  geom_jitter( position = position_jitter(0.1), size=1) + 
  scale_color_manual(values = c('#a6cee3','#ff7f00','#af8dc3','#33a02c')) +
  #scale_y_continuous(limits = c(0,10),breaks = c(0, 2,4,6,8,10)) +
  scale_y_continuous(limits = c(0,6),breaks = c(0, 1,2,3,4,5,6)) +
  theme_classic() +
  theme(legend.position = "none")
p1.wds

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

p1.dfu <- ggplot(dfutest, aes(Condition, log2(NormProp+1), colour=Condition)) +
  geom_col(data = dfutest.summary, fill = NA, color = "black") + #
  geom_errorbar( aes(ymin = log2(NormProp-sd+1), ymax = log2(NormProp+sd+1)), 
                 data = dfutest.summary, width = 0.3, color="black") +
  geom_jitter( position = position_jitter(0.1), size=1) + 
  scale_color_manual(values = c('#a6cee3','#a1d76a','#e9a3c9')) +
  #scale_y_continuous(limits = c(0,10),breaks = c(0, 2,4,6,8,10)) +
  #scale_y_continuous(limits = c(0,6),breaks = c(0, 1,2,3,4,5,6)) +
  theme_classic() +
  theme(legend.position = "none")
p1.dfu

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

p1.vu <- ggplot(vutest, aes(Condition, log2(NormProp+1), colour=Condition)) +
  geom_col(data = vutest.summary, fill = NA, color = "black") + #
  geom_errorbar( aes(ymin = log2(NormProp-sd+1), ymax = log2(NormProp+sd+1)), 
                 data = vutest.summary, width = 0.3, color="black") +
  geom_jitter( position = position_jitter(0.1), size=1) + 
  scale_color_manual(values = c('#a6cee3','#fdbf6f')) + 
  #scale_y_continuous(limits = c(0,10),breaks = c(0, 2,4,6,8,10)) +
  scale_y_continuous(limits = c(0,6),breaks = c(0, 1,2,3,4,5,6)) +
  theme_classic() +
  theme(legend.position = "none")
p1.vu

# combine the plots
#pdf("Result_CellProportion_allInte_FB_prolif.pdf", useDingbats = FALSE, width = 4, height = 2)
p1.wds + p1.dfu + p1.vu + plot_layout(widths = c(2,1.5,1))
dev.off()

