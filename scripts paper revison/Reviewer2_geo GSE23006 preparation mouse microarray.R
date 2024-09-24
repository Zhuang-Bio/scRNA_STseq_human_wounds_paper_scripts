# Function to download the GEO data
downGSE <- function(GSE_ID = NULL, destdir = ".") {
  
  require(GEOquery)
  eSet <- getGEO(GSE_ID, destdir = destdir, getGPL = T) #getGPL:Download the annotation file
  
  exprSet = exprs(eSet[[1]])
  pdata = pData(eSet[[1]])
  
  write.csv(exprSet, paste0(GSE_ID, "_exprSet.csv"))
  write.csv(pdata, paste0(GSE_ID, "_metadata.csv"))
  return(eSet)
}
downGSE(GSE_ID = "GSE23006", destdir = ".")

GSE23006_gpl <- getGEO("GSE23006", AnnotGPL=TRUE)
GSE23006_gpl <- GSE23006_gpl$GSE23006_series_matrix.txt.gz@featureData@data

####----Step 2 Process microarray data----####
######--GSE23006--######
library(tidyverse)
library(GEOquery)

GSE23006_exp <- read_csv("GSE23006_exprSet.csv")
colnames(GSE23006_exp)[1] <- "probeID"
dim(GSE23006_exp);table(is.na(GSE23006_exp))
GSE23006_met <- read_csv("GSE23006_metadata.csv")
GSE23006_gpl <- getGEO("GSE23006", AnnotGPL=TRUE)
GSE23006_gpl <- GSE23006_gpl$GSE23006_series_matrix.txt.gz@featureData@data

GSE23006_exp_anno <- GSE23006_exp %>% 
  left_join(.,GSE23006_gpl[,c(1,3)],by=c("probeID"="ID")) %>% 
  dplyr::select(50,2:49) %>% distinct(`Gene symbol`, .keep_all = TRUE)
colnames(GSE23006_exp_anno)[1] <- "Gene"

GSE23006_exp_anno <- GSE23006_exp_anno %>% arrange(Gene) %>% slice(-1)



# Extract the metadata
GroupInfo <- GSE23006_met[,c(3,45,46)] %>% setNames(c("GEO", "Condition", "Tissue")) %>% 
  filter(Tissue == "skin")
length(unique(GroupInfo$GEO))
GroupInfo$Condition <- gsub(" ", "", GroupInfo$Condition)
GroupInfo$Sample <- paste0(GroupInfo$GEO,"_",GroupInfo$Condition)
GroupInfo <- GroupInfo %>% select(4, everything())

microarrayData <- GSE23006_exp_anno[, GroupInfo$GEO]
microarrayData <- cbind(GSE23006_exp_anno$Gene, microarrayData)

colnames(microarrayData)[1] <- "GeneID"

identical(colnames(microarrayData)[-1], GroupInfo$GEO)
colnames(microarrayData)[-1] <- GroupInfo$Sample

data.table::fwrite(microarrayData, "GSE23006_mouse_skin_healing_exp.txt", sep = "\t")
data.table::fwrite(GroupInfo, "GSE23006_mouse_skin_healing_meta.txt", sep = "\t")


microarrayData <- data.table::fread("GSE23006_mouse_skin_healing_exp.txt")
GroupInfo <- data.table::fread("GSE23006_mouse_skin_healing_meta.txt")
# load the human and mouse genes
geneanno <- data.table::fread("Private_human_mouse_GeneSymbol_conversion.txt")

microarrayData_f <- microarrayData %>% left_join(., geneanno, by=c("GeneID" = "ms_gene")) %>% 
  select(26,2:25) %>% distinct(hs_gene, .keep_all = TRUE) %>% slice(-1)

colnames(microarrayData_f)[1] <- "GeneID"
data.table::fwrite(microarrayData_f, "GSE23006_mouse_skin_healing_exp.txt", sep = "\t")
