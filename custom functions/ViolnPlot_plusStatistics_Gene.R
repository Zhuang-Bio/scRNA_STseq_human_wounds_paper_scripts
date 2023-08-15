hswound.sub <- readRDS("../s2_Seurat_allSample_subclustering/allNew_subcluster_keratins_220203.rds")
FeaturePlot(hswound.sub, features = "WAKMAR2")

hswound.sub.gra <- subset(x = hswound.sub, CellTypes == "Gra_I")
VlnPlot(hswound.sub.gra, features = "WAKMAR1", pt.size = 0.0001, group.by = "Condition") + NoLegend()

wakmar2_conditions <- subset(x = hswound.sub.gra, features = "WAKMAR1")
wakmar2 <- FetchData(object = wakmar2_conditions, vars = c("Condition", "WAKMAR1"), slot = "data")
FetchData(object = wakmar2_conditions, vars = c("Condition", "WAKMAR1"), slot = "data") %>% 
  group_by(Condition) %>%  
  summarise(avg = mean(`WAKMAR1`, na.rm = TRUE)) %>% ungroup()


library(ggpubr)
my_comparisons <- list(c("Wound1", "Skin"),c("Wound30", "Wound7"), c("Wound7", "Wound1"),
                       c("Wound30", "Wound1"), c("Wound30", "Skin"), c("Wound7", "Skin"))
ggviolin(wakmar2, x = "Condition", y = "WAKMAR1", fill = "Condition",
         palette = c('#e41a1c','#377eb8','#4daf4a','#984ea3'),
         add = "boxplot", add.params = list(fill = "white")) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") +
  theme(legend.position = "none")
