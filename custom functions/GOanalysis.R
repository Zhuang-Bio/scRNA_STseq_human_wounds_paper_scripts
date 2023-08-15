GO_BP <- function(data=NULL, regulation=NULL){
  suppressMessages(require(clusterProfiler))  #loading the package if it is not loaded
  suppressMessages(require(org.Hs.eg.db)) 
  
  theme_cus <- theme(
    #panel.grid.major.x = element_blank(),
    #panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(size = 0.5, color = "black"),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"))
  
  if(is.null(regulation)){
    print("Using all the DE genes for GO analysis")
    Enriched_gene <- as.character(unique(data$entrezgene_id))
    ego_BP_up <- enrichGO(gene = Enriched_gene,
                          OrgDb         = org.Hs.eg.db,
                          ont           = "BP",
                          keyType = "ENTREZID",
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.05,
                          qvalueCutoff  = 0.2,
                          minGSSize = 10,
                          maxGSSize = 500,
                          readable = T)
    ego_BP_up_f <- as.data.frame(ego_BP_up)
    data.table::fwrite(ego_BP_up_f, file = paste0(condition, celltypename, "_GOterms.txt"), sep="\t")
    if (nrow(ego_BP_up_f) > 0) {
      dotplot(ego_BP_up, showCategory=10, label_format = 80) + ggtitle("Top 10 GO-BPs") +
        theme_cus
    } else{
      print("No enriched terms")
    }
      
  }else {
    print(paste0("Using all the ", regulation, " DE genes for GO analysis"))
    tmpdata <- data %>% dplyr::filter(Type == regulation)
    Enriched_gene <- as.character(unique(tmpdata$entrezgene_id))
    ego_BP_up <- enrichGO(gene = Enriched_gene,
                          OrgDb         = org.Hs.eg.db,
                          ont           = "BP",
                          keyType = "ENTREZID",
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.05,
                          qvalueCutoff  = 1,
                          minGSSize = 10,
                          maxGSSize = 500,
                          readable = T)
    ego_BP_up_f <- as.data.frame(ego_BP_up)
    data.table::fwrite(ego_BP_up_f, file = paste0(condition, celltypename, "_GOterms.txt"), sep="\t")
    if (nrow(ego_BP_up_f) > 0) {
      dotplot(ego_BP_up, showCategory=10, label_format = 80) + ggtitle("Top 10 GO-BPs") +
        theme_cus
    } else{
      print("No enriched terms")
    }
  }
}


KEGG <- function(data=NULL, regulation=NULL){
  suppressMessages(require(clusterProfiler))  #loading the package if it is not loaded
  suppressMessages(require(org.Hs.eg.db)) 
  
  theme_cus <- theme(
    #panel.grid.major.x = element_blank(),
    #panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(size = 0.5, color = "black"),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"))

  if(is.null(regulation)){
    print("Using all the DE genes for KEGG analysis")
    Enriched_gene <- as.character(unique(data$entrezgene_id))
    kk_up <- enrichKEGG(gene         = Enriched_gene,
                        organism     = 'hsa',
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff  = 1,
                        minGSSize = 10,
                        maxGSSize = 500)
    kk_up_f <- as.data.frame(kk_up)
    if (nrow(kk_up_f) > 0) {
      dotplot(kk_up, showCategory=5, label_format = 80) + ggtitle("Top 10 KEGG") +
        theme_cus
    } else{
      print("No enriched terms")
    }
    
  }else {
    print(paste0("Using all the ", regulation, " DE genes for KEGG analysis"))
    tmpdata <- data %>% dplyr::filter(Type == regulation)
    Enriched_gene <- as.character(unique(tmpdata$entrezgene_id))
    kk_up <- enrichKEGG(gene         = Enriched_gene,
                        organism     = 'hsa',
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff  = 1,
                        minGSSize = 10,
                        maxGSSize = 500)
    kk_up_f <- as.data.frame(kk_up)
    if (nrow(kk_up_f) > 0) {
      dotplot(kk_up, showCategory=5, label_format = 80) + ggtitle("Top 10 KEGG") +
        theme_cus
    } else{
      print("No enriched terms")
    }
  }  
}
