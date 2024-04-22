library(Seurat)
library(SeuratObject)
library(tidyverse)
library(magrittr)
library(ggplot2)
library(scCustomize)
library(patchwork)
library(RColorBrewer)
library(viridis)
library(ComplexHeatmap)
library(circlize)

setwd("")

rm(list = ls())

kera_sub <- readRDS("allNew_subcluster_keratins_220203.rds")

# re-annotate the cell clusters
fac_levs <- c("Bas_I", "Bas_prolif", "Bas_mig", 
              "Spi_I", "Spi_II_a", "Spi_II_b",
              "Spi_III", "Spi_mig",
              "Gra_I")
current_cluster_ids <- c(1,5,6,3,0,4,8,2,7) # List of current cluster IDs
tmp <- plyr::mapvalues(x=kera_sub$SCT_snn_res.0.5, from=current_cluster_ids, to=fac_levs)

kera_sub$upCellTypes <- tmp
kera_sub$upCellTypes <- factor(kera_sub$upCellTypes, levels = fac_levs)

ct.cols <- c('#807dba','#9e9ac8','#ccebc5',
               '#fe9929','#fec44f','#fee391',
               '#fb8072','#b3de69','#fccde5')
names(ct.cols) <- fac_levs

##################################################
####---- Fig 2 UMAP plot of keratinocytes ----####
##################################################
p.umap.kera <- DimPlot(kera_sub, group.by = "upCellTypes", label = F, cols = ct.cols) +
  NoAxes() + ggtitle("")
p.umap.kera

pdf("Fig2_UMAP_kera.pdf", useDingbats = F, width = 5, height = 3)
p.umap.kera
dev.off()

library(scCustomize)
FeaturePlot_scCustom(kera_sub, features = c("MKI67"), split.by = "Condition") & NoAxes()
#pdf("suppleFig3_UMAP_kera_MKI67.pdf", useDingbats = F, width = 15, height = 3)
#dev.off()

# Check the expression of KRT10, KRT14, and KRT16
VlnPlot_scCustom(seurat_object = kera_sub, features = "KRT10", 
                 group.by = "upCellTypes", 
                 plot_median = TRUE,y.max = 8) & NoLegend()
VlnPlot_scCustom(seurat_object = kera_sub, features = "KRT14", 
                 group.by = "upCellTypes", 
                 plot_median = TRUE,y.max = 8) & NoLegend()
VlnPlot_scCustom(seurat_object = kera_sub, features = "KRT16", 
                 group.by = "upCellTypes", 
                 plot_median = TRUE,y.max = 8) & NoLegend()

###############################################
####---- Fig 2 dotplot of marker genes ----####
###############################################
# load the marker genes from the FindAllMarker step
markers <- data.table::fread("allNew_subcluster_keratins_220203_mgCellTypes.txt")
table(markers$cluster)
top200markers <- markers %>% group_by(cluster) %>% top_n(n = 200, wt = avg_log2FC) %>% arrange(cluster)
data.table::fwrite(top200markers, "Markers_top200_subKC.txt", sep = "\t")

# calculate the scale data
DefaultAssay(kera_sub) <- "RNA"
kera_sub <- NormalizeData(kera_sub, verbose = TRUE,  normalization.method = "LogNormalize", scale.factor = 10000)
kera_sub <- ScaleData(object = kera_sub, features = rownames(kera_sub))

top20 <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) %>% 
  mutate(pct.dif = pct.1 - pct.2) %>% 
  arrange(cluster, desc(pct.dif), desc(avg_log2FC)) %>% ungroup()
top20 %>% filter(cluster == "Gra_I") %>% pull(gene)
FeaturePlot(kera_sub, features = c("KRT31"))

# manually select the top 5-8 representative markers according to the top20 marker genes
top_repre_markers <- c("KRT15", "KRT31", "COL17A1", "ASS1", "POSTN", 
                       "NUSAP1", "STMN1", "TOP2A", "MKI67", "CENPF",
                       "MMP3", "MMP1", "AREG", "TNFRSF12A", "FGFBP1", "S100A2", "S100A10",
                       "ITPRIPL2", "NRARP", "MT1X", "MT1G", "APOE",
                       "CHP2", "KCNK7", "DEFB1", 
                       "THEM5", "DSC1", "KRTDAP",
                       "TNFAIP2", "IRF1", "TNFSF10", "NFKBIA", "CCL27",
                       "KRT6A", "KRT6B", "KRT6C", "KRT17", "KRT16", "S100A7", "S100A8", "S100A9",
                       "SLURP2", "KLK7", "CNFN", "FLG", "LORICRIN")

#----Dotplot plotting----#
kera_sub$upCellTypes <- factor(kera_sub$upCellTypes, levels = rev(fac_levs))
plot_marker <- DotPlot(kera_sub, features = top_repre_markers, 
                       group.by = "upCellTypes", cols = c("white", "#cb181d"), 
                       dot.scale = 4, col.min = 0, dot.min = 0.1
) + theme(axis.text.x = element_text(angle = 90, hjust = 1),
          panel.border = element_rect(colour = "black", fill=NA)) + labs(x="", y="")
#pdf("Fig2_Dotplot_markers.pdf", useDingbats = F, width = 10, height = 3)
plot_marker
#dev.off()


################################################
####---- Fig 2 cell proportion analysis ----####
################################################
# step 1. 12 individuals divided into cell types and their numbers
step1 <- table(kera_sub$orig.ident, kera_sub$upCellTypes) %>% as.data.frame()
step1$sepa <- "Epidemis"

# step 2. read the cell numbers of dermis and epidermis from the main clustering analysis
## since we use a 1:1 ratio epidermis:dermis cells for scRNA-seq (more unbiased way)
step2 <- data.table::fread("all_epidermis_dermis_perIndi_220613.txt")

# step 3. calculate the proportion of each cell type for each individual
step3 <- step1 %>% left_join(., step2, by=c("sepa"="sepa", "Var1"="Var1")) %>% distinct() %>% 
  mutate(Prop=Freq/a_sum)

# step 4. calculate the total normalized proportions of each cell type per condition
df.group <- step3 %>% mutate(Sample = gsub("PWH..", "", Var1)) %>% group_by(Sample, Var2) %>% summarise(Freq=sum(Prop)) %>% 
  ungroup() %>% group_by(Sample) %>% mutate(Freq_new = Freq/sum(Freq), lbl = scales::percent(Freq_new)) %>% ungroup() %>% 
  rename("Cluster" = "Var2")

df.group$Sample <- gsub("D0", "Skin", df.group$Sample)
df.group$Sample <- gsub("D1", "Wound1", df.group$Sample)
df.group$Sample <- gsub("D7", "Wound7", df.group$Sample)
df.group$Sample <- gsub("D30", "Wound30", df.group$Sample)
df.group$Sample <- factor(df.group$Sample, levels = c("Skin", "Wound1", "Wound7", "Wound30"))
df.group$Cluster <- factor(df.group$Cluster, levels = fac_levs)

p_cp_5 <- ggplot(df.group, aes(x = Sample, y = Freq_new, fill = Cluster)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = ct.cols) + 
  xlab('') +
  scale_y_continuous(breaks = seq(0, 1, .2), 
                     expand = c(0, 0.01),
                     labels = c("0%", "20%", "40%", "60%", "80%", "100%"),
                     name = 'Percentage') +
  geom_text(aes(label = lbl), 
            size = 3, 
            position = position_stack(vjust = 0.5)) +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    axis.text.x = element_text(size = 14, color = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.y = element_text(size = 16, color = "black"),
    legend.position = "right",
    legend.text = element_text(size = 12, color = "black")
  ) 
p_cp_5

#pdf("Fig2_cellproportion_kera.pdf", useDingbats = F, width = 5, height = 6)
p_cp_5
#dev.off()

#########################################################
# run the statistics of cell proportions compared to Skin
df.group.out <- step3 %>% mutate(Group = gsub("PWH..", "", Var1)) 
data.table::fwrite(df.group.out, file = "Fig2/KCcellproportion_normalized.txt", sep = "\t")

df.group.out <- data.table::fread("Fig2/KCcellproportion_normalized.txt")
# make a function to run all the test using the
# quasibinomial and normal Wilcoxon test
ct <- unique(df.group.out$Var2) %>% as.character();ct

stat_sig <- list()
for (i in seq_along(ct)){
  ct_sel <- ct[i]
  ct_por <- df.group.out %>% dplyr::filter(Var2 == ct_sel) %>% arrange(Group)
  
  stats_all <- function(GroupInfo = NULL){
    ct_por_gr <- ct_por %>% dplyr::filter(Group %in% GroupInfo)
    # check the normality of the variables
    print(paste0("Shapiro.test on ", GroupInfo[1]))
    ct_por %>% dplyr::filter(Group %in% GroupInfo[1]) %>% pull(Prop) %>% shapiro.test() %>% print()
    print(paste0("Shapiro.test on ", GroupInfo[2]))
    ct_por %>% dplyr::filter(Group %in% GroupInfo[2]) %>% pull(Prop) %>% shapiro.test() %>% print()
    
    # paired samples Wilcoxon test
    if(TRUE){
      res <- t.test(Prop ~ Group, data = ct_por_gr)
    }else{
      res <- t.test(Prop ~ Group, data = ct_por_gr)
    }
    print(res)
    wilcoxon_paired = res$p.value
    # quasibinomial test
    test.quasi = glm(formula = Prop ~ Group, data = ct_por_gr, family=quasibinomial)
    print(summary(test.quasi))
    quasibinomial=anova(test.quasi, test = "LRT")$`Pr(>Chi)`[2]
    return(c(wilcoxon_paired, quasibinomial))
  }
  # Test of D1 versus D0
  res_all_1=stats_all(GroupInfo = c("D0", "D1"))
  res_all_2=stats_all(GroupInfo = c("D0", "D7"))
  res_all_3=stats_all(GroupInfo = c("D0", "D30"))
  stat_sig[[ct_sel]] <- append(stat_sig[[ct_sel]], res_all_1)
  stat_sig[[ct_sel]] <- append(stat_sig[[ct_sel]], res_all_2)
  stat_sig[[ct_sel]] <- append(stat_sig[[ct_sel]], res_all_3)
}
# combine all the data
stat_sig_all <- do.call("rbind", stat_sig) %>% as.data.frame() %>% 
  rownames_to_column(var = "CellType") %>% 
  setNames(c("CellType", "WilcoPair_D1D0", "QuasiBino_D1D0", "WilcoPair_D7D0", "QuasiBino_D7D0", "WilcoPair_D30D0", "QuasiBino_D30D0"))
data.table::fwrite(stat_sig_all, "Fig2_KCcellproportion_normalized_addSig.txt", sep = "\t")


##########################################################
####---- Fig 2 GO analysis of interested clusters ----####
##########################################################
# load the custom GO functions
source("../Functions/fig2Functions.R")
# Select the top200 marker genes 
top200 <- markers %>% group_by(cluster) %>% top_n(n = 200, wt = avg_log2FC) %>% 
  mutate(pct.dif = pct.1 - pct.2) %>% 
  arrange(cluster, desc(pct.dif), desc(avg_log2FC)) %>% ungroup()
#data.table::fwrite(top200, file = "Fig 2 tmp files/Top200_marker_genes_cluster.txt", sep = "\t")

bas_gos <- GO_BP(marksDF = top200, i="Bas_III") #Bas_III = Bas_mig
bas_gos <- bas_gos %>% dplyr::select(-3,-4) %>% arrange(qvalue, desc(Count))
spi_gos <- GO_BP(marksDF = top200, i="Spi_V") #Spi_V = Spi_mig
spi_gos <- spi_gos %>% dplyr::select(-3,-4) %>% arrange(qvalue, desc(Count))
spi3_gos <- GO_BP(marksDF = top200, i="Spi_IV") #Spi_IV = Spi_III
spi3_gos <- spi3_gos %>% dplyr::select(-3,-4) %>% arrange(qvalue, desc(Count))
#data.table::fwrite(bas_gos, file = "Fig 2 tmp files/GO_BPs_Bas_mig.txt", sep = "\t")

# read the interesting GOs manually selected from excel file
bas_mig1 <- readxl::read_xlsx("Fig 2 tmp files/GO_BPs.xlsx", sheet = 1)
spi_mig1 <- readxl::read_xlsx("Fig 2 tmp files/GO_BPs.xlsx", sheet = 2)
spi_III1 <- readxl::read_xlsx("Fig 2 tmp files/GO_BPs.xlsx", sheet = 3)

GOs_io <- readxl::read_xlsx("Fig 2 tmp files/GO_BPs.xlsx", sheet = 4) %>% arrange(orderid) %>% drop_na()
bas_mig1_f <- bas_mig1[bas_mig1$ID %in% GOs_io$ID, ] %>% mutate(Celltype = "Bas_mig")
spi_mig1_f <- spi_mig1[spi_mig1$ID %in% GOs_io$ID, ] %>% mutate(Celltype = "Spi_mig")
#spi_III1_f <- spi_III1[spi_III1$ID %in% GOs_io$ID, ] %>% mutate(Celltype = "Spi_III")

gos_comb <- rbind(bas_mig1_f, spi_mig1_f)#, spi_III1_f)
#data.table::fwrite(gos_comb, file = "Fig 2 tmp files/GO_combined.txt", sep = "\t")

# read the new order and reset the factor leves
gos_comb$Celltype <- factor(gos_comb$Celltype, levels = c("Bas_mig", "Spi_mig"))
gos_comb$Description <- factor(gos_comb$Description, levels = rev(GOs_io$Description))

goplot <- ggplot(data = gos_comb, aes(x = Celltype, y = Description, 
                            color = Celltype, size = -log10(qvalue))) + 
  geom_point() +
  scale_color_manual(values = c("#ccebc5", "#b3de69", "#fb8072")) +
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  scale_size(range = c(0,6), breaks = c(0, 2, 4, 6)) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
goplot

#pdf("Fig2_GOplot_migClusters.pdf", useDingbats = F, width = 5, height = 4)
goplot
#dev.off()


###############################################################
####---- Fig 2 DE analysis across different conditions ----####
###############################################################
# load the custom the DE analysis function
source("/custom functions/fig2Functions.R")
celltypes <- unique(kera_sub$upCellTypes) %>% as.character()
kera_sub@active.ident <- kera_sub$upCellTypes

# perform the DE analysis across conditions for each cell type
for (i in seq_along(celltypes)) {
  cellident <- celltypes[i]
  grps <- unique(kera_sub$Condition) %>% as.character()
  for (j in seq_along(grps)) {
    cond <- grps[j]
    degs <- DE.genes(objdata = kera_sub, group1 = cond, 
                     groupby = "Condition", subsetIdent = cellident)
    data.table::fwrite(degs, file = paste0("DE_", cellident, "_", cond, ".txt"), sep = "\t")
  }
}
#DE_wound1_skin <- DE.genes(objdata = kera_sub, group1 = "Wound1", group2 = "Skin",
#                           groupby = "Condition", subsetIdent = "Bas_mig")

# reload the DE results again
de.file <- list.files("DEGs_kera_subcluster/original DE results", pattern = ".txt")
de.file.path <- file.path("DEGs_kera_subcluster/original DE results", de.file)

de.files <- lapply(de.file.path, FUN = function(x) data.table::fread(x))
de.name <- gsub("(DE_)|(.txt)", "", de.file)
names(de.files) <- de.name

for (i in seq_along(de.name)) {
  de.files[[i]] <- de.files[[i]] %>% mutate(Group = de.name[i]) %>% mutate(diff.pct = pct.1 - pct.2) %>% dplyr::filter(abs(avg_log2FC) > 0.5)
  # filter the diff.pct to keep the same direction to the comparison
  tmp1 <- de.files[[i]] %>% dplyr::filter(Type == "Up") %>% dplyr::filter(diff.pct > 0)
  tmp2 <- de.files[[i]] %>% dplyr::filter(Type == "Down") %>% dplyr::filter(diff.pct < 0)
  de.files[[i]] <- rbind(tmp1, tmp2)
  #data.table::fwrite(de.files[[i]], file = paste0("DEGs_kera_subcluster/up organized DE results/", de.name[i], ".txt"), sep = "\t")
}

# run the GO analysis
i= 1:12 ; 13:27 ; 28:31
for (i in 32:36) {
  tmp_go <- GO_BP_genes(genelist = de.files[[de.name[i]]] %>% dplyr::filter(Type == "Up") %>% dplyr::slice(1:200) %>% pull(gene), i.name = paste0(de.name[i], "_up"))
  if(is.null(tmp_go)){
    print("No GOs")
  }else{
    data.table::fwrite(tmp_go, file = paste0(de.name[i], "_GO_BP_up.txt"), sep = "\t")
  }
  
  tmp_go_down <- GO_BP_genes(genelist = de.files[[de.name[i]]] %>% dplyr::filter(Type == "Down") %>% dplyr::slice(1:200) %>% pull(gene), i.name = paste0(de.name[i], "_down"))
  if(is.null(tmp_go_down)){
    print("No GOs")
  }else{
    data.table::fwrite(tmp_go_down, file = paste0(de.name[i], "_GO_BP_down.txt"), sep = "\t")
  }
  #rm(tmp_go, tmp_go_down)
}

###############################################################
# draw the dotplot of top DE genes across different conditions
top15genes <- c(skin_basmig = de.files[[5]] %>% dplyr::slice(1:15) %>% pull(gene), 
                w1_basmig = de.files[[6]] %>% dplyr::slice(1:15) %>% pull(gene),
                w7_basmig = de.files[[8]] %>% dplyr::slice(1:15) %>% pull(gene),
                w30_basmig = de.files[[7]] %>% dplyr::slice(1:15) %>% pull(gene)) %>% unique()

top15genes <- c(skin_spimig = de.files[[33]] %>% dplyr::slice(1:15) %>% pull(gene), 
                w1_spimig = de.files[[34]] %>% dplyr::slice(1:15) %>% pull(gene),
                w7_spimig = de.files[[36]] %>% dplyr::slice(1:15) %>% pull(gene),
                w30_spimig = de.files[[35]] %>% dplyr::slice(1:15) %>% pull(gene)) %>% unique()


plot_DEs <- DotPlot(object = subset(kera_sub, subset = upCellTypes == "Bas_mig"), #change the cell type
                    features = rev(top15genes), cols = c("white", "#cb181d"), 
                    dot.scale = 3, col.min = -1, col.max = 1, 
                    group.by = "Condition") + coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10),
        panel.border = element_rect(colour = "black", fill=NA)) + labs(x="", y="")

#pdf("Fig2_DEgenes_conditions_Bas_mig.pdf", useDingbats = F, width = 5, height = 10)
plot_DEs
#dev.off()


#####################################
# Draw the expression of DEGs (Based on the GO terms of Figure 2f)
library(scRNAtoolVis)
Idents(kera_sub) <- kera_sub$upCellTypes

# Bas_mig DEGs
## Day 1
mig_gene <- c("MMP1", "MMP3", "ITGA3", "HBEGF", "ANXA3", "PRSS3", "S100A2", "FGFBP1", "ANXA1", "TMSB4X") 
immune <- c("SERPINB1", "SERPINB3","SERPINB4", "CAP1", "S100A11", "S100A8", "S100A9")
#metallothionein <- c("MT1G", "MT1H","MT2A")
## Day 7
ecm_day7 <- c("LAMC2","LAMA3", "LAMB3", "COL17A1")
cell_adhesion <- c("ITGAV", "ITGB6", "ITGB4", "GJB2", "GJB6")

higenes <- c(mig_gene, immune, ecm_day7, cell_adhesion) #metallothionein,
#pdf("DEG_Bas_mig_dotplot.pdf", useDingbats = F, width = 10, height = 6)
jjDotPlot(object = subset(kera_sub, subset = upCellTypes == "Bas_mig"), #subset(inteData, idents = c("Bas_mig")
          gene = higenes,
          xtree = FALSE, ytree = FALSE,
          id = 'Condition',
          cluster.order = rev(c("Skin","Wound1","Wound7","Wound30")),
          rescale = T,
          #dot.col = c('blue','white','red'),
          #midpoint = 0,
          rescale.min = 0,
          rescale.max = 1) + scale_size(breaks = c(0, 25, 50, 75, 100), range = c(1,4))
#dev.off()

# Spi_mig DEGs
## Day 1
mig_gene <- c("KRT16", "HBEGF", "ANXA3", "PRSS3", "S100A2", "FGFBP1", "MYH9", "PLK2","GIPC1", "RAB25", "GFUS") 
immune <- c("SERPINB1", "SERPINB3","SERPINB4", "CTSC", "S100A9", "S100A7", "S100A8")
#metallothionein <- c("MT1G", "MT1H","MT2A")
## Day 7
differ <- c("SPRR1B", "IVL", "DSG3", "DSC2", "CDH3", "SFN")
metabolic <- c("SULT2B1", "ELOVL4", "UPP1", "PYGL", "VCP", "ATP5F1D", "NME1", "PGAM1", "PGK1", "EIF6", "TPI1", "PKM", "ENO1", "RAN", "GAPDH", "LDHA", "ATP5MC3")

higenes <- c(mig_gene, immune, differ) #metallothionein,
#pdf("DEG_Spi_mig_dotplot.pdf", useDingbats = F, width = 10, height = 6)
jjDotPlot(object = subset(kera_sub, subset = upCellTypes == "Spi_mig"), #subset(inteData, idents = c("Bas_mig")
          gene = higenes,
          xtree = FALSE, ytree = FALSE,
          id = 'Condition',
          cluster.order = rev(c("Skin","Wound1","Wound7","Wound30")),
          rescale = T,
          #dot.col = c('blue','white','red'),
          #midpoint = 0,
          rescale.min = 0,
          rescale.max = 1) + scale_size(breaks = c(0, 25, 50, 75, 100), range = c(1,4))
#dev.off()

