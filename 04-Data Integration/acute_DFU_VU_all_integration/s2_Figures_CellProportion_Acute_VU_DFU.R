#######################
# This script was used for cell proportion analysis of integrated acute wound, VU, and DFU cell clusters
# And the integrated clusters were further annotated according to annotation of sole acute wound scRNA-seq

## PS: For figures of cell proportion in paper is only based on acute wound data (Slightly different 
## to the analysis here, won't change our conclusion)
### For the normalization of cell numbers of each cell type, we normalized the cell numbers by the total number 
### of epidermal and dermal parts of acute wounds since we put 1:1 ratio of epidermis and dermis when constructing the library
### Regarding the VU and DFU, we normalized the cell numbers by the total number of each sample

# If you have any questions, please contact me: zhuang.liu@ki.se 

library(ggplot2)
library(patchwork)
library(tidyverse)


#########################
# cell types and colors
#########################
kc.ct <- c("Bas_I", "Bas_prolif", "Bas_mig", 
           "Spi_I", "Spi_II_a", "Spi_II_b",
            "Spi_III", "Spi_mig", "Gra_I")
kc.cols <- c('#807dba','#9e9ac8','#ccebc5',
                '#fe9929','#fec44f','#fee391',
                '#fb8072','#b3de69','#fccde5')
names(kc.cols) <- kc.ct

fb.ct <- c("FB_I_POSTN_COL11A1", "FB_I_POSTN_MMP11", 
           "FB_I_POSTN_COL4A1", "FB_I_SFRP4_COMP", 
           "FB_II_APOD_ITM2A", "FB_II_APOE_CCL19",
           "FB_III_ELN_LEPR", 
           "FB_prolif",
           "FB_SFRP1_CRABP1", "FB_ELN_SFRP4")
fb.cols <- c("#57a9e2", "#31c783", "#c6cc85", "#c072fd", "#e36767", 
             "#aa6d60", "#ff9a41", "#eca1d5", "#31d7e8", "#b0b0b0")
names(fb.cols) <- fb.ct

mye.ct <- c("Mac_inf", "Mac1","Mac2", "Mac3", 
            "pDC", "cDC1", "cDC2", "DC3", 
            "LC", "Apoptotic", "Cycling")
mye.cols <- c('#f1b6da','#df65b0', '#41ab5d', '#addd8e',
             "#810f7c", "#807dba", "#9970ab", "#b2abd2", 
             "#bf812d", "#fee0b6", "#f1a340"
)
names(mye.cols) <- mye.ct

lym.ct <- c("Treg", "Th", "ILCs", "Tc", 
              "ILC1_NK", "NK", "Ttol", 
              "Plasma", "Bcell")
lym.cols <- c('#807dba','#9970ab','#810f7c', '#41ae76', 
             '#74add1', '#8dd3c7','#76d72f',
             '#f768a1', '#fdb462'
)
names(lym.cols) <- lym.ct

endo.ct <- c("SMC", "Pericytes", "LE",
              "VE_arteriole", "VE_capillary",
              "VE_venule1", "VE_venule2")
endo.cols <- c('#12a22d', '#ce8b1a', '#fe8Bac', 
             '#d93860', '#3ba4db', 
             '#807dba', '#d4b9da')
names(endo.cols) <- endo.ct

# combine all the cell types and colors
ct.cols <- list(Keratinocyte = list(kc.ct, kc.cols),
                Fibroblast =list(fb.ct, fb.cols),
                Myeloid= list(mye.ct, mye.cols),
                Lymphoid= list(lym.ct, lym.cols),
                Endothelium = list(endo.ct, endo.cols)
                )


###########################
# Cell proportion analysis
# load the total number of epidermis and dermis for each sample
sepa_sm <- data.table::fread("./cellproportion_analysis_data/allInte_Dermis_Epi_perSample_Totcells.txt")
# load the cell number for each cell type
all_ct <- data.table::fread("./cellproportion_analysis_data/allInte_celltype_cellnumber.txt")

all_ct$newCellType <- gsub("SFRP1_CRABP1", "FB_SFRP1_CRABP1", all_ct$newCellType)
all_ct$newCellType <- gsub("ELN_SFRP4", "FB_ELN_SFRP4", all_ct$newCellType)

# get the part and group information
all_ct_part <- all_ct %>% select(5,7) %>% distinct()
groupinfo <- all_ct %>% select(2,3,4) %>% distinct()

unique(all_ct$newCellType)
celltypes <- unique(all_ct$mainCellType)

#i=1
for (i in seq_along(celltypes)) {
  ct <- celltypes[i]
  # extract the main cell type (do it one by one)
  all_ct_sub <- all_ct %>% filter(mainCellType == ct)
  
  # step 1. Samples divided into cell types and their numbers
  sm.tol <- table(all_ct_sub$orig.ident) %>% as.data.frame() %>% setNames(c("orig.ident", "TotalNumber"))
  
  step1 <- table(all_ct_sub$orig.ident, all_ct_sub$newCellType) %>% as.data.frame() %>% setNames(c("orig.ident", "newCellType", "Freq")) %>% 
    left_join(., all_ct_part, by=c("newCellType" = "newCellType")) %>% 
    left_join(., groupinfo[,c(1,3)],by=c("orig.ident"="orig.ident"))
  step1$newCellType <- as.character(step1$newCellType) 
  #step1 <- step1 %>% filter(Sample != "NS63D") # remove NS63D 
  
  # step 2. calculate the proportion of each cell type for each individual
  step2 <- step1 %>% left_join(., sepa_sm, by=c("orig.ident"="orig.ident", "sepa"="sepa")) %>% 
    left_join(., sm.tol, by=c("orig.ident"="orig.ident")) %>% 
    distinct() %>%
    mutate(newTotal=ifelse(.$Project == "AcuteWound", .$Total, .$TotalNumber)) %>% # For acute wound, normalize the cell number by separating dermis and epidermis
    #mutate(newTotal=TotalNumber) %>% # set the acute wound same as DFU and VU
    mutate(Prop=Freq/newTotal) %>% 
    left_join(., groupinfo[,c(1,2)], by=c("orig.ident" = "orig.ident")) 
  
  # step 3. calculate the total normalized proportions of each cell type per condition
  df.group <- step2 %>% 
    group_by(Condition, newCellType) %>% summarise(Freq=sum(Prop)) %>% 
    ungroup() %>% group_by(Condition) %>% 
    mutate(Freq_new = Freq/sum(Freq), lbl = scales::percent(Freq_new)) %>% ungroup()
  
  df.group$Condition <- factor(df.group$Condition, levels = c("Skin", "Wound1", "Wound7", "Wound30",
                                                              "H", "DFU_H", "DFU_NH", "NS", "VU"))
  df.group$newCellType <- factor(df.group$newCellType, levels = ct.cols[[ct]][[1]])
  
  p_cp_5 <- ggplot(df.group, aes(x = Condition, y = Freq_new, fill=newCellType)) +
    geom_bar(position = "fill", stat = "identity") +
    scale_fill_manual(values = ct.cols[[ct]][[2]]) + 
    xlab('') +
    scale_y_continuous(breaks = seq(0, 1, .2), 
                       expand = c(0, 0.01),
                       labels = c("0%", "20%", "40%", "60%", "80%", "100%"),
                       name = 'Percentage') +
    geom_text(aes(label = lbl), 
              size = 2, 
              position = position_stack(vjust = 0.5)) +
    theme_bw() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1),
      axis.text.y = element_text(size = 12, color = "black"),
      axis.title.y = element_text(size = 12, color = "black"),
      legend.position = "right",
      legend.text = element_text(size = 12, color = "black")
    )
  print(p_cp_5)
  pdf(file = paste0("allInte_CellType_", ct, "_cellproportion.pdf"), useDingbats = FALSE, width = 5, height = 5)
  print(p_cp_5)
  dev.off()
  
  # export the cell proportion
  data.table::fwrite(step2, file = paste0("allInte_CellType_", ct, "_cellproportion.txt"), sep="\t")
}


################################################
# Run the statistics for cell proportion analysis
# load the data
cellproportion_statis <- function(data_cp=NULL, project = "AcuteWound"){
  step3 <- data.table::fread(file = data_cp) %>% dplyr::filter(Project == project) %>% mutate(norProp = 100*Prop) #%>% dplyr::filter(!newCellType %in% c("FB_prolif"))
  
  ct.stat <- unique(step3$newCellType) %>% as.character();print(ct.stat)
  
  stat_sig <- list()
  for (i in seq_along(ct.stat)){
    ct_sel <- ct.stat[i]
    ct_por <- step3 %>% dplyr::filter(newCellType == ct_sel) %>% arrange(Condition)
    
    stats_all <- function(GroupInfo = NULL, paired = TRUE){
      ct_por_gr <- ct_por %>% dplyr::filter(Condition %in% GroupInfo)
      
      # check the normality of the variables
      print(paste0("Shapiro.test on ", GroupInfo[1]))
      ct_por %>% dplyr::filter(Condition %in% GroupInfo[1]) %>% pull(norProp) %>% shapiro.test() %>% print()
      print(paste0("Shapiro.test on ", GroupInfo[2]))
      ct_por %>% dplyr::filter(Condition %in% GroupInfo[2]) %>% pull(norProp) %>% shapiro.test() %>% print()
      
      # paired samples Wilcoxon test
      if(paired){
        #res <- wilcox.test(norProp ~ Condition, data = ct_por_gr, paired = FALSE, exact=FALSE, correct= TRUE) # try the unpaired test
        res <- t.test(norProp ~ Condition, data = ct_por_gr)
      }else{
        #res <- wilcox.test(norProp ~ Condition, data = ct_por_gr, paired = FALSE, exact=FALSE, correct= TRUE)
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
ct_acute <- cellproportion_statis(data_cp = "allInte_CellType_Keratinocyte_cellproportion.txt", 
                                  project = "AcuteWound")

tmp_df <- data.table::fread("allInte_CellType_Keratinocyte_cellproportion.txt")
tmp_df %>% filter(newCellType == "Bas_I") %>% filter(Condition == "Skin") %>% pull(Prop) %>% shapiro.test()

ct_acute4 <- cellproportion_statis(data_cp = "allInte_CellType_Keratinocyte_cellproportion.txt", 
                                  project = "AcuteWound")
data.table::fwrite(ct_acute, "allInte_CellType_Keratinocyte_cellproportion_sig_Acute.txt", sep="\t")

ct_dfu <- cellproportion_statis(data_cp = "allInte_CellType_Keratinocyte_cellproportion.txt", 
                                  project = "DFU")
data.table::fwrite(ct_dfu, "allInte_CellType_Keratinocyte_cellproportion_sig_DFU.txt", sep="\t")

ct_vu <- cellproportion_statis(data_cp = "allInte_CellType_Keratinocyte_cellproportion.txt", 
                                  project = "VU")
data.table::fwrite(ct_vu, "allInte_CellType_Keratinocyte_cellproportion_sig_VU.txt", sep="\t")

## Fibroblast
ct_acute <- cellproportion_statis(data_cp = "allInte_CellType_Fibroblast_cellproportion.txt", 
                                  project = "AcuteWound")
data.table::fwrite(ct_acute, "allInte_CellType_Fibroblast_cellproportion_sig_Acute.txt", sep="\t")

ct_dfu <- cellproportion_statis(data_cp = "allInte_CellType_Fibroblast_cellproportion.txt", 
                                project = "DFU")
data.table::fwrite(ct_dfu, "allInte_CellType_Fibroblast_cellproportion_sig_DFU.txt", sep="\t")

ct_vu <- cellproportion_statis(data_cp = "allInte_CellType_Fibroblast_cellproportion.txt", 
                               project = "VU")
data.table::fwrite(ct_vu, "allInte_CellType_Fibroblast_cellproportion_sig_VU.txt", sep="\t")

## Myeloid
ct_acute <- cellproportion_statis(data_cp = "allInte_CellType_Myeloid_cellproportion.txt", 
                                  project = "AcuteWound")
data.table::fwrite(ct_acute, "allInte_CellType_Myeloid_cellproportion_sig_Acute.txt", sep="\t")

ct_dfu <- cellproportion_statis(data_cp = "allInte_CellType_Myeloid_cellproportion.txt", 
                                project = "DFU")
data.table::fwrite(ct_dfu, "allInte_CellType_Myeloid_cellproportion_sig_DFU.txt", sep="\t")

ct_vu <- cellproportion_statis(data_cp = "allInte_CellType_Myeloid_cellproportion.txt", 
                               project = "VU")
data.table::fwrite(ct_vu, "allInte_CellType_Myeloid_cellproportion_sig_VU.txt", sep="\t")

## Lymphoid
ct_acute <- cellproportion_statis(data_cp = "allInte_CellType_Lymphoid_cellproportion.txt", 
                                  project = "AcuteWound")
data.table::fwrite(ct_acute, "allInte_CellType_Lymphoid_cellproportion_sig_Acute.txt", sep="\t")

ct_dfu <- cellproportion_statis(data_cp = "allInte_CellType_Lymphoid_cellproportion.txt", 
                                project = "DFU")
data.table::fwrite(ct_dfu, "allInte_CellType_Lymphoid_cellproportion_sig_DFU.txt", sep="\t")

ct_vu <- cellproportion_statis(data_cp = "allInte_CellType_Lymphoid_cellproportion.txt", 
                               project = "VU")
data.table::fwrite(ct_vu, "allInte_CellType_Lymphoid_cellproportion_sig_VU.txt", sep="\t")

## Endothelium
ct_acute <- cellproportion_statis(data_cp = "allInte_CellType_Endothelium_cellproportion.txt", 
                                  project = "AcuteWound")
data.table::fwrite(ct_acute, "allInte_CellType_Endothelium_cellproportion_sig_Acute.txt", sep="\t")

ct_dfu <- cellproportion_statis(data_cp = "allInte_CellType_Endothelium_cellproportion.txt", 
                                project = "DFU")
data.table::fwrite(ct_dfu, "allInte_CellType_Endothelium_cellproportion_sig_DFU.txt", sep="\t")

ct_vu <- cellproportion_statis(data_cp = "allInte_CellType_Endothelium_cellproportion.txt", 
                               project = "VU")
data.table::fwrite(ct_vu, "allInte_CellType_Endothelium_cellproportion_sig_VU.txt", sep="\t")


## organize the files
df.name <- dir(path = "cellproportion_analysis_results/", pattern = ".txt")
df.path <- paste0("cellproportion_analysis_results/", df.name)

df.list <- lapply(df.path, FUN = function(x) data.table::fread(x))
df.name.mv <- gsub("(^allInte_CellType_)|(.txt)|(_cellproportion)", "", df.name)
names(df.list) <- df.name.mv


