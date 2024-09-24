############################
# This script is updated from 
# "/Users/zhuliu/Desktop/sc st/Versions_SC_ST_Project/V2_Figures_231017/01-Scripts/Figure_cell2location"
# Author: zhuang.liu@ki.se 
############################
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(directlabels)
library(patchwork)


###############################
# load the original NMF results
## Change the niche names to the same as in the manuscript
nmf_n15 <- readxl::read_xlsx("microenvironment_NMF_n15.xlsx", sheet = 1)

df0 <- nmf_n15 %>% select(2:3, 4) %>% group_by(Sample, Condition) %>% 
  summarise(sum(mean_nUMI_factorsfact_0)) %>% ungroup() %>% 
  setNames(c("Sample", "Condition", "fact")) %>% 
  dplyr::mutate(Factors = "Fact0")
df1 <- nmf_n15 %>% select(2:3, 5) %>% group_by(Sample, Condition) %>% 
  summarise(sum(mean_nUMI_factorsfact_1)) %>% ungroup() %>% 
  setNames(c("Sample", "Condition", "fact")) %>% 
  dplyr::mutate(Factors = "Fact1")
df2 <- nmf_n15 %>% select(2:3, 6) %>% group_by(Sample, Condition) %>% 
  summarise(sum(mean_nUMI_factorsfact_2)) %>% ungroup() %>% 
  setNames(c("Sample", "Condition", "fact")) %>% 
  dplyr::mutate(Factors = "Fact2")
df3 <- nmf_n15 %>% select(2:3, 7) %>% group_by(Sample, Condition) %>% 
  summarise(sum(mean_nUMI_factorsfact_3)) %>% ungroup() %>% 
  setNames(c("Sample", "Condition", "fact")) %>% 
  dplyr::mutate(Factors = "Fact3")
df4 <- nmf_n15 %>% select(2:3, 8) %>% group_by(Sample, Condition) %>% 
  summarise(sum(mean_nUMI_factorsfact_4)) %>% ungroup() %>% 
  setNames(c("Sample", "Condition", "fact")) %>% 
  dplyr::mutate(Factors = "Fact4")
df5 <- nmf_n15 %>% select(2:3, 9) %>% group_by(Sample, Condition) %>% 
  summarise(sum(mean_nUMI_factorsfact_5)) %>% ungroup() %>% 
  setNames(c("Sample", "Condition", "fact")) %>% 
  dplyr::mutate(Factors = "Fact5")
df6 <- nmf_n15 %>% select(2:3, 10) %>% group_by(Sample, Condition) %>% 
  summarise(sum(mean_nUMI_factorsfact_6)) %>% ungroup() %>% 
  setNames(c("Sample", "Condition", "fact")) %>% 
  dplyr::mutate(Factors = "Fact6")
df7 <- nmf_n15 %>% select(2:3, 11) %>% group_by(Sample, Condition) %>% 
  summarise(sum(mean_nUMI_factorsfact_7)) %>% ungroup() %>% 
  setNames(c("Sample", "Condition", "fact")) %>% 
  dplyr::mutate(Factors = "Fact7")
df8 <- nmf_n15 %>% select(2:3, 12) %>% group_by(Sample, Condition) %>% 
  summarise(sum(mean_nUMI_factorsfact_8)) %>% ungroup() %>% 
  setNames(c("Sample", "Condition", "fact")) %>% 
  dplyr::mutate(Factors = "Fact8")
df9 <- nmf_n15 %>% select(2:3, 13) %>% group_by(Sample, Condition) %>% 
  summarise(sum(mean_nUMI_factorsfact_9)) %>% ungroup() %>% 
  setNames(c("Sample", "Condition", "fact")) %>% 
  dplyr::mutate(Factors = "Fact9")
df10 <- nmf_n15 %>% select(2:3, 14) %>% group_by(Sample, Condition) %>% 
  summarise(sum(mean_nUMI_factorsfact_10)) %>% ungroup() %>% 
  setNames(c("Sample", "Condition", "fact")) %>% 
  dplyr::mutate(Factors = "Fact10")
df11 <- nmf_n15 %>% select(2:3, 15) %>% group_by(Sample, Condition) %>% 
  summarise(sum(mean_nUMI_factorsfact_11)) %>% ungroup() %>% 
  setNames(c("Sample", "Condition", "fact")) %>% 
  dplyr::mutate(Factors = "Fact11")
df12 <- nmf_n15 %>% select(2:3, 16) %>% group_by(Sample, Condition) %>% 
  summarise(sum(mean_nUMI_factorsfact_12)) %>% ungroup() %>% 
  setNames(c("Sample", "Condition", "fact")) %>% 
  dplyr::mutate(Factors = "Fact12")
df13 <- nmf_n15 %>% select(2:3, 17) %>% group_by(Sample, Condition) %>% 
  summarise(sum(mean_nUMI_factorsfact_13)) %>% ungroup() %>% 
  setNames(c("Sample", "Condition", "fact")) %>% 
  dplyr::mutate(Factors = "Fact13")
df14 <- nmf_n15 %>% select(2:3, 18) %>% group_by(Sample, Condition) %>% 
  summarise(sum(mean_nUMI_factorsfact_14)) %>% ungroup() %>% 
  setNames(c("Sample", "Condition", "fact")) %>% 
  dplyr::mutate(Factors = "Fact14")

alldf <- rbind(df0, df1, df2, df3, df4, df5, df6, df7, df8, df9, df10,
               df11, df12, df13, df14)

newFactor = paste0("Fact",c(1,2,3,4,5,6,
                            7,8,9,10,11,
                            12,13,14,15))
rawFactor = paste0("Fact",c(0,1,2,3,4,5,
                            6,14,12,13,8,
                            11,9,10,7))

alldf$newFactors <- plyr::mapvalues(alldf$Factors, from = rawFactor, to = newFactor)


###################
# Plotting
## Prepare the data and plot all niches at one panel
alldf_re <- alldf %>% select(-4) %>% pivot_wider(names_from = newFactors, values_from = fact)
ggboxplot(df0, x = "Condition", y = "fact",
          add = "jitter")
alldf_re$Condition <- factor(alldf_re$Condition, levels = c("Skin", "Wound1",
                                                            "Wound7", "Wound30"))
ggboxplot(alldf_re, x = "Condition",
          y = c("Fact1","Fact2","Fact3", "Fact4","Fact5","Fact6","Fact7",
                "Fact8","Fact9", "Fact10","Fact11","Fact12","Fact13","Fact14", "Fact15"),
          combine = TRUE,
          ylab = "Total number",
          color = "Condition", palette = "jco",
          add = "jitter"             # column containing point labels
)

## Prepare the data and combine the plots of niches
alldf$Condition <- factor(alldf$Condition, levels = c("Skin", "Wound1", "Wound7", "Wound30"))

niches_plot <- list()
for (i in seq_along(newFactor)) {
  niche_io <- newFactor[i]
  alldf_tmp <- alldf %>% dplyr::filter(newFactors %in% niche_io)
  tmp_pl <- ggplot(data = alldf_tmp, aes(x=Condition, y=fact, group=Sample, color=Sample)) +
              geom_line(linewidth=1) + geom_point(size=1) +
              xlab("") + ylab("")+
              theme_bw() +
              theme(
                aspect.ratio = 1,
                legend.position = 'right',
                panel.grid.major = element_line(colour = "grey90",size = 0.2,linetype = 1),
                axis.text.x = element_text(colour = "black", size = 12, angle = 30),
                axis.text.y = element_text(colour = "black", size = 12)) +
              facet_wrap(~newFactors) +
              theme(strip.background = element_rect(colour = "black", fill = "white"))
  
  niches_plot[[i]] <- tmp_pl
}

#pdf("spatial_all_microNiches_1.pdf", useDingbats = FALSE, width = 16, height = 16)
# plot all the figures
wrap_plots(niches_plot, ncol = 5, nrow = 3, guides = "collect", axis_titles = "collect_y") 
#ylab(title="Total UMI counts for each microarray environment") 
#dev.off()

###################################################################
# Re-draw the spatial feature plots of niche7, niche8, and niche13
# load the niche data
library(Seurat)
library(tidyverse)
library(scCustomize)
library(RColorBrewer)
library(viridis)

nmf_n15 <- readxl::read_xlsx("microenvironment_NMF_n15.xlsx", sheet = 1)
colnames(nmf_n15)
newFactor = paste0("Niche",c(1,2,3,4,5,6,
                            7,8,9,10,11,
                            12,13,14,15))
rawFactor = paste0("mean_nUMI_factorsfact_",c(0,1,2,3,4,5,
                            6,14,12,13,8,
                            11,9,10,7))
# rename the factor, keep it same as in Figure 1f
tmp_match <- plyr::mapvalues(colnames(nmf_n15)[4:18], from = rawFactor, to = newFactor)
colnames(nmf_n15)[4:18] <- tmp_match

# load the spatial data
st_seu <- readRDS("/Users/zhuliu/Desktop/sc st/ST plotting scripts/Seurat_STseq_integrate_update.rds")
mt <- st_seu@meta.data %>% rownames_to_column(var = "barcode")
mt <- mt %>% left_join(., nmf_n15[, -c(2:3)], by=c("barcode" = "SpotIndex" ))
identical(mt$barcode, nmf_n15$SpotIndex)

# modified the metadata of spatial seurat object
mt_f <- mt %>% column_to_rownames(var = "barcode")
st_seu@meta.data <- mt_f


########################
# plot interested niches
img.names <- unique(st_seu@meta.data$Sample_name) %>% as.character();img.names
# choose which donor
plot_donors <- img.names[grep("Donor4_", img.names)];plot_donors

ratio_list <- list()
for (i in seq_along(plot_donors)) {
  coord <- GetTissueCoordinates(object = st_seu, image = plot_donors[i])
  # calculate the aspect ratio of rows to columns
  myratio <- (max(coord$imagerow) - min(coord$imagerow)) / (max(coord$imagecol) - min(coord$imagecol))
  ratio_list[[i]] <- myratio
}
ratio_list

colours = magma(10)
plot_list <- list()
for (i in seq_along(plot_donors)) {
  plot_list[[i]] <- SpatialFeaturePlot(st_seu, 
                                       features = "Niche1",
                                       images = plot_donors[i], 
                                       stroke = 0, combine = T, 
                                       pt.size.factor = 1.2, 
                                       #pt.size.factor = ratio_list[[i]], 
                                       alpha = c(1, 1), 
                                       max.cutoff = "q98",
                                       image.alpha = 0, crop = FALSE) +
    #theme(aspect.ratio = ratio_list[[i]]) + 
    #scale_fill_gradientn(limits=c(0, 1), colours=magma(10), na.value = "#FCFDBFFF") # niche13
    #scale_fill_gradientn(limits=c(0, 5), colours=magma(10), na.value = "#FCFDBFFF") # niche8
    #scale_fill_gradientn(limits=c(0, 3), colours=magma(10), na.value = "#FCFDBFFF") # niche7
    scale_fill_gradientn(limits=c(0, 10), colours=magma(10), na.value = "#FCFDBFFF") # niche10
}
patchwork::wrap_plots(plot_list, ncol = 2)

#pdf("ST_Donor4_genePlot_Niche10_defaultColor.pdf", useDingbats = F, width = 8, height = 8)
pdf("ST_Donor4_genePlot_Niche1_custColor.pdf", useDingbats = F, width = 18, height = 18)
patchwork::wrap_plots(plot_list, ncol = 4)
dev.off()

