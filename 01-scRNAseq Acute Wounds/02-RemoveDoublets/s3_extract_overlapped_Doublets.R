library(tidyverse)
library(magrittr)
require(VennDiagram)
require(grDevices)
require(venn)

setwd("/Users/zhuliu/Desktop/scRNA_STseq/proj_10X_woundhealing/03_results/01-Seurat-PreAnalysis/01_Seurat_qc_Doublets")

######################################
####---Step 1. load the results---####
######################################
#read the DoubletFinder results
samples <- list.files("./DoubletFinder/", pattern = ".txt$")
samples.path <- paste0("./DoubletFinder/", samples)
sm_files <- lapply(samples.path, FUN = function(x){
  data.table::fread(x)
})
names(sm_files) <- gsub("_.*", "", samples) 
for (i in seq_along(sm_files)) {
  sm_files[[i]] <-  sm_files[[i]] %>% 
    mutate(barcodename = paste0(names(sm_files[i]), "_", barcode)) %>% 
    mutate(Samplename = names(sm_files[i])) %>% select(-2:-7)
}
rm(samples, samples.path, i)

#read the pyScrublet results
samples.py <- list.files("./pyScrublet/", pattern = ".txt$")
samples.py.path <- paste0("./pyScrublet/", samples.py)
sm_pyfiles <- lapply(samples.py.path, FUN = function(x){
  data.table::fread(x)
})
names(sm_pyfiles) <- gsub("_.*", "", samples.py) 
for (i in seq_along(sm_pyfiles)) {
 sm_pyfiles[[i]] <-  sm_pyfiles[[i]] %>% 
   mutate(barcodename = paste0(names(sm_pyfiles[i]), "_", barcode)) %>% 
   mutate(Samplename = names(sm_pyfiles[i])) %>% 
   mutate(predicted_doublets = ifelse(predicted_doublets == "TRUE", "Doublet", "Singlet"))
}
rm(samples.py, samples.py.path, i)
sm_doublets <- do.call("rbind", sm_pyfiles)

##############################################################
####---Step 2. Check the distribution of doublet scores---####
##############################################################
####PS: Be careful with three variables of trait, doublet_scores thresholds and optnames in DoubletFinder results
print(names(sm_files)[1]); print(names(sm_pyfiles)[1])
trait="PWH26D0"   ##replaced with different sample nams
# density plot of all data
d <- density(sm_pyfiles[[trait]]$doublet_scores)
plot(d)
polygon(d, col="red", border="blue")

sm_pyfiles[[trait]] %>% 
  group_by(predicted_doublets) %>% 
  summarise(Mean =mean(doublet_scores),
            Median = median(doublet_scores),
            Min = min(doublet_scores),
            Max = max(doublet_scores),
            n = n())

#sm_pyfiles[[trait]] %>% 
#  group_by(NewDoublets) %>% #with new thresholds
#  summarise(Mean =mean(doublet_scores), 
#            Median = median(doublet_scores),
#            Min = min(doublet_scores),
#            Max = max(doublet_scores),
#            n = n())

#reset the doublet thresholds for each sample
# PWH26D0 = 0.05 # new output  0.1
# PWH26D1 = 0.1 # new output  0.15
# PWH26D30 = 0.05 # new output  0.1
# PWH26D7 = 0.1
# PWH27D0 = 0.1
# PWH27D1 = 0.15 # new output  0.2
# PWH27D30 = 0.1 # new output  0.15
# PWH27D7 = 0.1 # new output  0.15
# PWH28D0 = 0.1 # new output  0.15
# PWH28D1 = 0.1 # new output  0.15
# PWH28D30 = 0.1 # new output  0.15
# PWH28D7 = 0.1

trait="PWH26D0"
# density plot of two types (singlet, doublet)
ggplot(sm_pyfiles[[trait]], aes(x=doublet_scores, color=predicted_doublets)) +
  geom_density() +
  geom_histogram(aes(y=..density..), alpha=0.5, 
                 position="identity", binwidth = 0.05)

sm_pyfiles[[trait]] %<>% 
  mutate(NewDoublets = ifelse(doublet_scores > 0.15, "Doublet", "Singlet")) ##Change the threshold values to define the doublets 
table(sm_pyfiles[[trait]]$NewDoublets); table(sm_pyfiles[[trait]]$predicted_doublets)
data.table::fwrite(sm_pyfiles[[trait]], file = paste0("Scrublet_", trait, ".txt"), sep = "\t")
#Check the interactions between pyScrublet and DoubletFinder
#ven plot, first remember to source the function
optname=colnames(sm_files[[trait]])[grep("DF\\.classification", colnames(sm_files[[trait]]))]
#try all the index [1...n] of optname and choose the most overlapped one
dev.off(); vennPlotFunction(DFopt = optname[2], trait = trait) %>% grid.draw() 

df_res <- sm_files[[trait]][, c(1,6,7,8:9)] #optname[3] = 6:7, optname[2] = 4:5, optname[1] = 2:3
data.table::fwrite(df_res, file = paste0("DoubletFinder_", trait, ".txt"), sep = "\t")

#write the overlapped results
doublets.df <- which(sm_files[[trait]][[optname[2]]] == "Doublet")
doublet.pyS <- which(sm_pyfiles[[trait]][["NewDoublets"]] == "Doublet")
tmps.doublet <- intersect(sm_files[[trait]]$barcode[doublets.df], 
                          sm_pyfiles[[trait]]$barcode[doublet.pyS]) %>% 
                as.data.frame() %>% 
                rename(Doublets = ".") %>% 
                mutate(sampleID = trait,
                       mapID = paste0(trait, "_", Doublets))
data.table::fwrite(tmps.doublet, file = paste0("Shared_doublets_", trait, ".txt"), sep = "\t")

vennPlotFunction <- function(DFopt = NULL, trait = trait){
  tmp1 <- which(sm_files[[trait]][[DFopt]] == "Singlet")
  tmp2 <- which(sm_files[[trait]][[DFopt]] == "Doublet")
  tmp3 <- which(sm_pyfiles[[trait]][["NewDoublets"]] == "Singlet")
  tmp4 <- which(sm_pyfiles[[trait]][["NewDoublets"]] == "Doublet")
  
  temp.venn <- venn.diagram(list(DF.singlet = sm_files[[trait]]$barcode[tmp1],
                                 DF.doublet = sm_files[[trait]]$barcode[tmp2],
                                 pySc.singlet = sm_pyfiles[[trait]]$barcode[tmp3],
                                 pySc.doublet = sm_pyfiles[[trait]]$barcode[tmp4]),
                            fill = c("dodgerblue", "seagreen", "darkorange1", "orchid3"), 
                            alpha = c(0.5,0.5,0.5,0.5), cex=1.2, 
                            cat.fontface = 2, filename = NULL, scaled = T, euler.d = TRUE)
  return(temp.venn)
}


###Output Doublets from either pyScrublet or DoubletFinder
#Union doublets from either Scrublet or DoubletFinder
#Doublets <- unique(c(sm_doublets_f$barcodename, sm_doublets_df_f$barcodename)) %>% as.data.frame() %>% rename("barcodename" = ".") %>% mutate(newDoublet = "Doublet")
Doublets <- data.table::fread("Doublets_either_Scrublet_DoubletFinder.txt")
sm_doublets_f <- sm_doublets %>% left_join(., Doublets, by=c("barcodename" = "barcodename")) %>% 
  select(4,2,6) %>% 
  mutate(newDoublet = ifelse(is.na(newDoublet), "Singlet", newDoublet)) %>% 
  rename('predicted_doublets' = 'newDoublet')
data.table::fwrite(sm_doublets_f, file = "Doublets_either_Scrublet_DoubletFinder.txt", sep = "\t")


