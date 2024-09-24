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
downGSE(GSE_ID = "GSE209609", destdir = ".")

GSE209609_gpl <- getGEO("GSE209609", AnnotGPL=TRUE)
GSE209609_gpl <- GSE209609_gpl$GSE209609_series_matrix.txt.gz@featureData@data

####----Step 2 Process microarray data----####
######--GSE209609--######
library(tidyverse)
library(GEOquery)

GSE209609_exp <- read_csv("GSE209609_exprSet.csv")
colnames(GSE209609_exp)[1] <- "probeID"
dim(GSE209609_exp);table(is.na(GSE209609_exp))
GSE209609_met <- read_csv("GSE209609_metadata.csv")
GSE209609_gpl <- getGEO("GSE209609", AnnotGPL=TRUE)
GSE209609_gpl <- GSE209609_gpl$GSE209609_series_matrix.txt.gz@featureData@data

GSE209609_exp_anno <- GSE209609_exp %>% 
  left_join(.,GSE209609_gpl[,c(1,2)],by=c("probeID"="ID"))

# annotate the ensembl genes
library(biomaRt)
# Create a character vector of Ensembl IDs		
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset(mart = mart, dataset = "hsapiens_gene_ensembl")
# extract the genes
genedata <- getBM(mart=mart, attributes=c("ensembl_gene_id", "external_gene_name",
                                          "gene_biotype", "description"), uniqueRows=T, useCache=FALSE)
genedata <- genedata %>% distinct(ensembl_gene_id, .keep_all = TRUE) #remove duplicates of ensembl ID

GSE209609_exp_anno <- GSE209609_exp %>% 
  left_join(.,GSE209609_gpl[,c(1,2)],by=c("probeID"="ID")) %>% 
  left_join(., genedata[,1:2], by=c("ORF" = "ensembl_gene_id")) %>% 
  dplyr::select(99,2:97) %>% distinct(external_gene_name, .keep_all = TRUE) %>% 
  arrange(external_gene_name) %>% slice(-1)

# Extract the metadata
GroupInfo <- GSE209609_met[,c(2,3,40,41)] %>% setNames(c("Sample", "GEO", "Condition", "Tissue")) %>% 
  filter(Tissue == "skin")
length(unique(GroupInfo$Sample))
GroupInfo$Condition <- gsub(" ", "", GroupInfo$Condition)
GroupInfo$Condition <- gsub("0hours", "Hour0", GroupInfo$Condition)
GroupInfo$Condition <- gsub("6hours", "Hour6", GroupInfo$Condition)
GroupInfo$Sample <- gsub("Skin_", "s", GroupInfo$Sample)

microarrayData <- GSE209609_exp_anno[, GroupInfo$GEO]
microarrayData <- cbind(GSE209609_exp_anno$external_gene_name, microarrayData)

colnames(microarrayData)[1] <- "GeneID"

identical(colnames(microarrayData)[-1], GroupInfo$GEO)
colnames(microarrayData)[-1] <- GroupInfo$Sample

data.table::fwrite(microarrayData, "GSE209609_palate_skin_healing_exp.txt", sep = "\t")
data.table::fwrite(GroupInfo, "GSE209609_palate_skin_healing_meta.txt", sep = "\t")

