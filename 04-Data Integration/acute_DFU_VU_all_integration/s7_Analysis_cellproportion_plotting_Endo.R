library(tidyverse)
library(ggplot2)
library(patchwork)
library(ggbreak)
library(ggpubr)

setwd("./cellproportion_analysis_results/")

# run the Endothelium
alldf <- data.table::fread("allInte_CellType_Endothelium_cellproportion.txt")
table(alldf$newCellType)
#LE    Pericytes          SMC VE_arteriole VE_capillary   VE_venule1   VE_venule2 
#41           41           41           41           41           41           41 

# choose the cell types for plotting
alldf_test <- alldf %>% dplyr::filter(newCellType == "VE_capillary")

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
  scale_y_continuous(limits = c(0,4),breaks = c(0, 1,2,3,4,5,6)) +
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
  scale_y_continuous(limits = c(0,4),breaks = c(0, 1,2,3,4,5,6)) +
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
  scale_y_continuous(limits = c(0,4),breaks = c(0, 1,2,3,4,5,6)) +
  theme_classic() +
  theme(legend.position = "none")
p1.vu

# combine the plots
#pdf("Result_CellProportion_allInte_Endo_VEcapillary.pdf", useDingbats = FALSE, width = 4, height = 2)
p1.wds + p1.dfu + p1.vu + plot_layout(widths = c(2,1.5,1))
dev.off()

