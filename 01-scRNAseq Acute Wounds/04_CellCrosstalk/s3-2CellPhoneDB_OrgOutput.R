setwd("/Users/zhuliu/Desktop/scRNA_STseq/proj_10X_woundhealing/03_results/01-Seurat-PreAnalysis/02_Seurat_BatchCorrection/s1_Seurat_allSample_harmony/s8_Cell_Cell_L_R_analysis_upCellTypes/CellPhoneDB_221107/")

library(tidyverse)

rm(list = ls())
output <- dir(path = "countsData_Wound30/", pattern = ".txt$")
output.name <- gsub(".txt", "", output)
output.path <- paste0(getwd(), "/countsData_Wound30/", output)

output.files <- lapply(output.path, FUN = function(x){
  read.table(x, header=T, stringsAsFactors = F, sep="\t", comment.char = '', check.names=F)
})  
names(output.files) <- output.name
rm(output, output.name, output.path)

sig_int <- output.files$significant_means[, -c(1:12)]

int_num <- list()
for (i in seq_along(colnames(sig_int))) {
  int_num[[colnames(sig_int)[i]]] <- table(sig_int[, i]>0) %>% sum() 
}
inte_num_com <- data.frame(CellType = names(int_num), Number = unlist(int_num)) %>% separate(col = "CellType", into = c("CellType1", "CellType2"), sep = "\\|")

cts <- list()
for (i in 1:nrow(inte_num_com)) {
  cts[[i]] <- c(inte_num_com[i,]$CellType1, inte_num_com[i,]$CellType2)
}

# sort the data
cts_sort <- lapply(cts, sort)
cts_sort <- do.call("rbind", cts_sort) %>% as.data.frame()
cts_sort$newCellType <- paste0(cts_sort$V1, "|", cts_sort$V2)
table(cts_sort$newCellType)

inte_num_com$newCellType <- cts_sort$newCellType

inte_num_com_f <- inte_num_com %>% group_by(newCellType) %>% summarise(Total = sum(Number)) %>% ungroup() %>% 
  separate(col = newCellType, into = c("CellType1", "CellType2"), sep = "\\|") %>%
  pivot_wider(names_from = CellType2, values_from = Total) %>% column_to_rownames(var = "CellType1")


library(viridis)
pheatmap::pheatmap(inte_num_com_f, cluster_rows = F, cluster_cols = F, 
                   #color = magma(20), breaks = seq(0,400, by=20), 
                   na_col="white", border_color ="white")





# load the clean metadata
mt <- data.table::fread("/Users/zhuliu/Desktop/sc st/plotting scripts/newsubCluster_annotation/Subclusters.txt")
unique(mt$CellType) %>% sort()
  

#########################
# Dotplot or heatmap
source("plot_functions.R")

lr_pairs <- c("CXCR1_CXCL5",  "CXCR2_CXCL5",  "CXCL1_CXCR1",  "CXCL1_CXCR2")

ct_io <- c("M1|Bas_mig", "M1|Spi_mig", "M2|Bas_mig", "M2|Spi_mig", "Mac_inf|Bas_mig", "Mac_inf|Spi_mig", 
           "Mac_mig|Bas_mig", "Mac_mig|Spi_mig", "Bas_mig|M1", "Spi_mig|M1", "Bas_mig|M2", "Spi_mig|M2")

dot_plot(selected_rows = lr_pairs, selected_columns = ct_io, 
         filename = "CXCL1_5_LR.pdf", means_path = "countsData_Wound30/means.txt", pvalues_path = "countsData_Wound30/pvalues.txt")


#####################
# reorganize the data
org_lrs <- function(selected_rows = NULL, selected_columns = NULL,
                    means_path = './means.txt', pvalues_path = './pvalues.txt'
){
  all_pval = read.table(pvalues_path, header=T, stringsAsFactors = F, sep="\t", comment.char = '', check.names=F)
  all_means = read.table(means_path, header=T, stringsAsFactors = F, sep="\t", comment.char = '', check.names=F)
  
  intr_pairs = all_pval$interacting_pair
  all_pval = all_pval[,-c(1:11)]
  all_means = all_means[,-c(1:11)]
  
  if(is.null(selected_rows)){
    selected_rows = intr_pairs
  }
  
  if(is.null(selected_columns)){
    selected_columns = colnames(all_pval)
  }
  
  sel_pval = all_pval[match(selected_rows, intr_pairs), selected_columns]
  sel_means = all_means[match(selected_rows, intr_pairs), selected_columns]
  
  df_names = expand.grid(selected_rows, selected_columns)
  pval = unlist(sel_pval)
  pval[pval==0] = 0.0009
  plot.data = cbind(df_names, pval)
  pr = unlist(as.data.frame(sel_means))
  #pr[pr==0] = 1
  plot.data = cbind(plot.data, pr)
  colnames(plot.data) = c('pair', 'clusters', 'pvalue', 'mean')
  return(plot.data)
}

skin_data <- org_lrs(selected_rows = lr_pairs, selected_columns = ct_io, 
                     means_path = "./countsData_Skin/means.txt", pvalues_path = "./countsData_Skin/pvalues.txt")
wound1_data <- org_lrs(selected_rows = lr_pairs, selected_columns = ct_io, 
                     means_path = "./countsData_Wound1/means.txt", pvalues_path = "./countsData_Wound1/pvalues.txt")
wound7_data <- org_lrs(selected_rows = lr_pairs, selected_columns = ct_io, 
                     means_path = "./countsData_Wound7/means.txt", pvalues_path = "./countsData_Wound7/pvalues.txt")
wound30_data <- org_lrs(selected_rows = lr_pairs, selected_columns = ct_io, 
                     means_path = "./countsData_Wound30/means.txt", pvalues_path = "./countsData_Wound30/pvalues.txt")


ggplot(skin_data,aes(x=clusters, y=pair)) +
  geom_point(aes(size=-log10(pvalue), color=mean)) +
  scale_color_gradientn('Interaction Score', colors=magma(10), limits=c(0,1)) +
  #scale_size(range = c(0,3)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text=element_text(size=14, colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(size=12, colour = "black"),
        axis.title=element_blank(),
        panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))

