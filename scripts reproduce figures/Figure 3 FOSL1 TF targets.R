library(tidyverse)
library(magrittr)
library(SCENIC)
library(SCopeLoomR)
library(AUCell)

setwd("")

############################
# Figure 3
## Draw the plot of Regulon Z-scores for wound clusters
regulons <- readxl::read_xlsx("/pySCENIC/z_hswoundWound_cellType_RSSs_zscores.xlsx", sheet = 1) %>% 
  select(-1) %>% group_by(CellTypes) %>% arrange(CellTypes, desc(Z)) %>% 
  ungroup() %>% mutate(Rank = rep(1:384, 9))

#regulons$regulon <- gsub("_","", regulons$regulon)
zscore_plot <- function(Cluster = NULL, title=NULL){
  require(ggrepel)
  inputDF <- regulons %>% filter(CellTypes == Cluster)
  ggplot(inputDF, aes(x=Rank, y= Z)) +
    geom_point(colour = "#3b5cc4", size = 0.8) +
    geom_point(data = inputDF %>% dplyr::slice(1:15), mapping = aes(x = Rank, y=Z), colour = "red", size = 1.2) +
    ylab('Regulon specificity Z-scores') +
    ggtitle(paste0(title)) +
    geom_hline(yintercept=0, linetype="dashed", color = "black") +
    #scale_y_continuous(limits = c(-3,3), breaks = c(-3,-2,-1,0,1,2,3)) +
    scale_x_continuous(limits = c(0, 400), breaks = c(0, 100, 200, 300)) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(size= 12),
      axis.title.x = element_text(size= 12),
      axis.text.y = element_text(size= 12),
      axis.title.y = element_text(size= 12),
      axis.title=element_text(size=12), #face = "bold"
      plot.title = element_text(size = 12)
    ) +
    geom_text_repel(data = inputDF %>% dplyr::slice(1:5), 
                    mapping = aes(x = Rank, y=Z, label = regulon),
                    color = "red",
                    nudge_x      = 100,
                    direction    = "y",
                    hjust        = 0,
                    segment.size = 0.2,
                    size = 2,
                    segment.color = "red"
    )
}

p1 <- zscore_plot(Cluster = "Bas_III", title = "Bas_mig") #Bas_mig
p2 <- zscore_plot(Cluster = "Spi_V", title = "Spi_mig") #Spi_mig

#pdf("Fig3_TF_Regulons_mig_bas_spi.pdf", useDingbats = F, width = 5, height = 3)
require(patchwork)
p1+p2
#dev.off()

# draw all the regulon TFs
z.plots <- lapply(unique(regulons$CellTypes), function(x) zscore_plot(x))

#####################################
# Plotting the FOSL1 regulon activity
regulon_mt <- data.table::fread("pySCENIC/z_hswoundWound_cellType_metadata_addRegulons.csv")
identical(colnames(kera_sub), regulon_mt$V1)
grep("FOSL1", colnames(regulon_mt))
reac <- regulon_mt[, 96]
kera_sub$Regulon_FOSL1 <- reac$`Regulon_FOSL1_(+)`

#fosl1 <- data.table::fread("FOSL1_regulonValues.txt")
hswound.krt.com <- readRDS("p_try_KRT6A+_wound_cells.rds")
mt <- hswound.krt.com@meta.data %>% rownames_to_column(var = "barcode") %>% 
  left_join(., fosl1[,c(1,29,30)])
hswound.krt.com$Regulon_FOSL1 <- mt$Regulon_FOSL1
hswound.krt.com$Regulon_FOSL1_scale <- mt$Regulon_FOSL1_scale
DimPlot(hswound.krt.com)
#pdf("FOSL1_regulon_umap_onlyWoundCells.pdf", useDingbats = FALSE, width = 6, height = 4)
FeaturePlot(
  hswound.krt.com,
  features=c("Regulon_FOSL1"), 
  #split.by = "Condition",
  combine=T, cols = c("#850aff","#400080","yellow"), pt.size = 0.2, max.cutoff = "q95"#, pt.size = 1.2
) + NoAxes() + ggtitle("Regulon_FOSL1(+)")
#dev.off()


###########################
# FOSL1 regulated targets
scenicLoomPath='hswound_SCENIC_AUC.loom' #here loom file contains multi-annotations from Step 2.2
loom <- open_loom(scenicLoomPath, mode="r") #only read mode

# Read regulons from loom file: 
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")

# extract the regulons and its activity scores
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom, column.attr.name="RegulonsAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]

# load the motifs of each TF from GRN
enrich_motif_TFs <- read_csv("hswound_reg.csv") #modified the header
enrich_motif_TFs[1,];enrich_motif_TFs[2,]
colnames(enrich_motif_TFs) <- c("TF", "MotifID", "AUC", "NES", "MotifSimilarityQvalue",
                                "OrthologousIdentity", "Annotation", "Context",
                                "TargetGenes", "RankAtMax")
enrich_motif_TFs <- enrich_motif_TFs[-1:-2, ]

# Extract the TF FOSL1
fosl1 <- enrich_motif_TFs[enrich_motif_TFs$TF %in% c("FOSL1", "ASCL2"),] %>% 
  mutate(AUC = parse_number(AUC)) %>% 
  mutate(RankAtMax = parse_number(RankAtMax)) #%>% 
  #dplyr::filter(AUC < 0.05 & Annotation == "gene is directly annotated")

grep("KRT6A", fosl1$TargetGenes) #[1]  1  2  3  5  6  7  8  9 10 11 12 13

fosl1$TargetGenes[1]

# load the marker genes of each cluster
mg <- data.table::fread("Top200_marker_genes_cluster.txt")
mg_spimig <- mg %>% dplyr::filter(cluster == "Spi_V") %>% arrange(desc(avg_log2FC)) %>% slice(1:50)
mg_basmig <- mg %>% dplyr::filter(cluster == "Bas_III") %>% arrange(desc(avg_log2FC)) %>% slice(1:50)

spi <- intersect(mg_spimig$gene, unlist(regulons['FOSL1(+)'])) %>% sort() %>% as.data.frame()
data.table::fwrite(spi, "FOSL1_target_genes_Spi_mig.txt")

bas <- intersect(mg_basmig$gene, unlist(regulons['FOSL1(+)'])) %>% sort() %>% as.data.frame()
data.table::fwrite(bas, "FOSL1_target_genes_Bas_mig.txt")


# GO plot of FOSL1 regulated targets
library(tidyverse)
library(ggplot2)
go_df <- read.csv("GO_hsapiens_20-12-2022_15-39-34__intersections.csv")
go_bp <- go_df %>% slice(1:3,5,6)
go_bp$term_name <- factor(go_bp$term_name, levels = rev(go_bp$term_name))
p3 <- ggplot(data = go_bp, aes(x=term_name, y = -log10(adjusted_p_value))) + 
  geom_bar(stat = "identity") +
  theme_bw() + 
  ylab("") + 
  xlab("") +
  coord_flip()
p3
#pdf("Fig3_FOSL1_targets_GOs.pdf", useDingbats = FALSE, width = 5, height = 2)
p3
#dev.off()

#########################################
# KC subclustering supplementary figures
## Wound and non-wound KCs classification
kera_sub <- readRDS("allNew_subcluster_keratins_220203.rds")

# load the classification of wound and healthy cells based on the expression of gene KRT6A
mt <- read.csv("p_try_KRT6A_wound_cells_GaussianSelection_metadata.csv")

identical(colnames(kera_sub), mt$barcode) #if TRUE, add the classification
kera_sub$KRT6A_class <- mt$KRT6A_class
#pdf("P_try_KRT6A_classification_gene_expression.pdf", height = 2.5, width = 12, useDingbats = F)
(DimPlot(object = kera_sub, 
         reduction = "umap", 
         label = F, group.by = "KRT6A_class", label.size = 5, split.by = "Condition") + NoAxes() + ggtitle(""))
#dev.off()

# set the color to keep the same as Figure 2a
cl.colors <- c('#807dba','#9e9ac8','#ccebc5',
               '#fe9929','#fec44f','#fee391',
               '#fb8072','#b3de69','#fccde5')
fac_levs <- c("Bas_I", "Bas_prolif", "Bas_mig", 
              "Spi_I", "Spi_II_a", "Spi_II_b",
              "Spi_III", "Spi_mig",
              "Gra_I")

# Wound cells UMAP (clusters)
wounds <- readRDS("p_try_KRT6A+_wound_cells.rds")
mt.io <- colnames(wounds) %>% as.data.frame() %>% rename("barcode" = ".") %>% 
  left_join(., mt[, c(1,28:31)], by=c("barcode" = "barcode")) %>% column_to_rownames(var = "barcode")
identical(colnames(wounds), rownames(mt.io))

wounds$doublet_scores <- mt.io$doublet_scores
wounds$upCellTypes <- mt.io$upCellTypes
wounds$upCellTypes <- factor(wounds$upCellTypes, levels = fac_levs)
#pdf("Fig3_woundcell.pdf", useDingbats = F, width = 16, height = 4)
DimPlot(wounds, group.by = "upCellTypes", cols = cl.colors) +
  DimPlot(wounds, group.by = "SCT_snn_res.0.8") +
  DimPlot(wounds, group.by = "SCT_snn_res.0.5")
#dev.off()

# Non-wound cells UMAP (clusters)
healthy <- readRDS("p_try_KRT6A-dim_healthy_cells.rds")
mt.io <- colnames(healthy) %>% as.data.frame() %>% rename("barcode" = ".") %>% 
  left_join(., mt[, c(1,28:31)], by=c("barcode" = "barcode")) %>% column_to_rownames(var = "barcode")
identical(colnames(healthy), rownames(mt.io))

healthy$doublet_scores <- mt.io$doublet_scores
healthy$upCellTypes <- mt.io$upCellTypes
healthy$upCellTypes <- factor(healthy$upCellTypes, levels = fac_levs)

#pdf("Fig3_healthycell.pdf", useDingbats = F, width = 16, height = 4)
DimPlot(healthy, group.by = "upCellTypes", cols = cl.colors) +
  DimPlot(healthy, group.by = "SCT_snn_res.0.8") +
  DimPlot(healthy, group.by = "SCT_snn_res.0.5")
#dev.off()

########################################
# Fig3e-h monocle3 pseudotime analysis
## please see script "subclustering/s3-3Wound_nonWound_cells_monocle3_CellOracle.rmd"

# Fig3i-j CellOracle in silico TF analysis
## please see python notebooks

###############################################
# FOSL1 expression in ST and scRNA-seq datasets
## ST-seq (only three epidermal clusters)
hswound.STcom <- readRDS("Seurat_STseq_integrate_update.rds")
DimPlot(hswound.STcom, label = T, label.size = 10) + NoAxes()
# re-annotate the cell clusters
fac_levs <- c("Basal_epi", "Suprabasal_epi", "Wound_edges",
              "Hair_follicle", "Papillary_dermis", 
              "Sebaceous_gland", "Sweat_gland", "Sweatgland_FB",
              "FB1_rich", "FB2_rich", 
              "Dermis1", "Dermis2","Smooth_muscle", 
              "Immune1_rich",  "Immune_endo", "Lymphatic_endo", "Mast_cell")

current_cluster_ids <- c(10, 8, 5,
                         16, 9, 13, 14, 4,
                         7, 3,
                         1, 15, 11, 
                         6, 2, 12, 17) # List of current cluster IDs
tmp <- plyr::mapvalues(x=hswound.STcom$SCT_snn_res.0.5, from=current_cluster_ids, to=fac_levs)
hswound.STcom$AnnoType <- tmp
hswound.STcom$AnnoType <- factor(hswound.STcom$AnnoType, levels = fac_levs)
DimPlot(hswound.STcom, label = T, group.by = "AnnoType") + NoAxes()

# check the expression of FOSL1 in epidermal clusters (excluding Hair Follicle cluster)
hswound.STcom <- subset(hswound.STcom, subset = AnnoType %in% c("Basal_epi","Suprabasal_epi", "Wound_edges"))
VlnPlot(hswound.STcom, features = "FOSL1", group.by = "Condition")

fosl1_st <- AverageExpression(hswound.STcom, assays = "Spatial", features = c("FOSL1"),
                              group.by = "Sample_name")$Spatial
fosl1_st <- fosl1_st %>% t()
fosl1_st <- fosl1_st %>% as.data.frame() %>% 
  rownames_to_column(var = "Sample") %>% setNames(c("Sample", "avgExp"))

df <- fosl1_st %>% 
  mutate(Sample = Sample) %>% 
  separate(col = Sample, into = c("Donor1", "Sample"), sep = "_")
avg.skin <- df %>% group_by(Sample) %>% 
  summarize(mean(avgExp)) %>% ungroup();avg.skin
df.f <- df %>% mutate(norExp = avgExp/avg.skin$`mean(avgExp)`[1]) %>% 
  mutate(lognorExp = log2(norExp+1))
df.f$Sample <- factor(df.f$Sample, levels = c("Skin", "Wound1", "Wound7", "Wound30"))
alldf.summary <- df.f %>% group_by(Sample) %>% 
  summarise(sd = sd(norExp, na.rm = TRUE), norExp = mean(norExp));alldf.summary

# plotting using ggpubr
require(ggpubr)
df.f$Sample <- factor(df.f$Sample, levels = c("Skin", "Wound1", "Wound7", "Wound30"))
my_comparisons <- list(c("Skin", "Wound1"), c("Skin", "Wound7"), c("Skin", "Wound30"))

pdf("FOSL1_exp_in_STdata.pdf", useDingbats = FALSE, width = 3, height = 4)
ggboxplot(df.f, "Sample", "lognorExp", color = "Sample",
          add = "jitter", ylab = "log2(fold change) of FOSL1 expression compared to skin in epidermal clusters", title = "",
          #ylim=c(0,7),
          xlab = "" #palette = c("#d95f02", "#e7298a"), 
) +
  stat_compare_means(comparisons = my_comparisons)
dev.off()

