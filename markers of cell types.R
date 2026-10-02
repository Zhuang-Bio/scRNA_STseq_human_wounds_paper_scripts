# Dotplot of marker genes
makers <- c("KRT15", "COL17A1", "CENPK", "MKI67", "MMP3", "GJB2",
            "KRT1", "DMKN", "KRT6B", "KRT6C", "SLURP1", "FLG",
            "SOX9", "ADGRL3","DACH1", "GPC5", "ESRRG", "PNLIPRP3",
            "DCT", "TYRP1", # MEL
            "COL1A1", "COL23A1", "PI16", "CCN5", "C3", "CCL19",
            "NEGR1", "LAMA2", "DIAPH3", "INHBA", "ADAM12", "FAP",
            "SOX2", "NRXN1", "MYOCD", "MYH11", "CCDC102B", "GUCY1A2",
            "CCL21", "MMRN1", "PODXL", "RASGRF2", 
            "SELE", "ACKR1", "CCL14", "VWF", "KDR", "FGF12",
            "ADGRF5", "RGS5", # Periva_stroma (this cluster close to pericyte and FBs)
            "EREG", "IL1B", "C1QA", "RNASE1", "MMP9", "ITGAX", # Mac3 is different to ATAC-seq Mac3
            "IDO1", "CLEC9A", "CCR7", "FCER1A", "LAMP3", "IL15",
            "CD207", "FCGBP", "TPSAB1", "TPSB2",
            "IL7R", "CAMK4", "GZMA", "GZMK","GNLY", "XCL1",
            "IGKC", "JCHAIN")

ct.names <- c("Bas", "Bas_prolif", "Bas_mig", "Spi", 
              "Spi_mig", "Gra", "HF_epi", "HF_bulgeCP",
              "Sweat_sebaceous","Melanocyte",
              "F1_SupPapi", "F2_UniReti", "F3_FRClike",
              "F4_5_HFnerve", "F6_InflaProlif", "F7_Myofib",
              "Schwann",
              "SMC", "Pericyte", "LE", "VE_art", 
              "VE_ven1", "VE_ven2", "VE_capi", 
              "Periva_stroma",
              "Mac_Inf", "Mac_antiInf", "Mac3",
              "cDC1", "cDC2", "DC3",
              "LC", "Mast",
              "T_help", "NKcell", "T_cytotoxic",
              "Plasma_B"
)


marker_list <- list(
  Bas = c("POSTN", "COL17A1", "KRT15", "LAMB4", "ASS1"),
  Bas_prolif = c("BRIP1", "DIAPH3", "CENPK", "BRCA2", "MKI67"),
  Bas_mig = c("TNC", "MMP3", "MMP1", "AREG", "MIR31HG"),
  Spi_I = c("KRT1", "KCNK7", "KRTDAP", "CHP2", "IL34", "IL33", "NRARP"),
  Spi_mig = c("KRT16", "KRT6B", "KRT6C", "S100A7", "SERPINB4"),
  Gra = c("KRT78", "KLK7", "SLURP1", "SLURP2", "FLG"),
  HF_ORS = c("SUSD4", "IGFL4", "PTN", "CIDEA", "CARD18", "SOX9"),
  HF_basal = c("ADGRL3", "DACH1", "MOXD1", "RUNX1", "TNC"),
  HF_LGR5 = c("GPC6", "GPC5", "CTNND2", "LGR5", "LHX2"),
  Sebaceous = c("MGST1", "PPARG", "ACSBG1", "CYP4F8", "GLDC"),
  SG_I = c("ESRRG", "PNLIPRP3", "TG", "CDH11", "PPARGC1A"),
  SG_II = c("DCD", "PDE4B", "PLCB4", "SCGB1B2P", "SCGB2A2")
)


fb_marker_list <- list(
  F1_SupPapi1 = c("COL18A1", "COL23A1", "NCKAP5", "LEPR", "TRABD2B"),
  F2_UniReti = c("PI16", "CCN5", "TNXB", "CD34", "SCARA5"),
  F3_PeriVasc = c("C3", "APOD", "FGF7", "IL6", "TNFSF14"),
  F3_FRClike = c("APOE", "CCL19", "TMEM132C", "CH25H", "TCIM"),
  F4_HFasso = c("COCH", "TNN", "PCDH15", "NRG3", "COL11A1"),
  F5_NerveLike = c("A2M","TM4SF1","MCTP1","EBF2","ITGA6"),
  F6_InflaProlif = c("DIAPH3", "RRM2", "TOP2A", "ASPM", "MKI67"),
  F7_preMyoFB = c("SFRP4", "BCAT1", "ALDH1A3", "GPC6", "RUNX2"),
  F7_ScarMyoFB = c("ADAM12", "POSTN", "ASPN", "PLPP4", "NRG1", "MMP11"),
  #F7_MyoFB = c("MDK", "FABP5", "PLPP4", "TAGLN"),
  F8_lipidFB = c("ABCA10", "ABCA6", "ABCA8", "ABCA9", "HPSE2"),  
  F9_Cilium = c("SNX31", "SAXO1", "SLC2A14", "ERI2", "IL19", "C9")
)


marker.list <- list(
  Mac1_IL1B       = c("CD68","CD163","THBS1","EREG","SERPINB2","IL1A","IL1B"),
  Mac1_SPP1       = c("AQP9","SEMA3C","APOE","SPP1","PLD1"),
  Mac2_COLEC12    = c("COLEC12","ME1","DAB2","FGFR1"),
  Mac2_F13A1      = c("MAFB","F13A1","RHOB","EGR1","SLC40A1"),
  Mac_FOLR2       = c("IL10","MMP19","MMP9","FOLR2","FGL2"),
  Mac_stress      = c("JUN","ATF3","HSPB1","HSPA1B","DNAJB1"),
  cDC1            = c("CLEC9A","FRY","NEGR1","CADM1","IDO1"),
  cDC2            = c("IL1R2","FCER1A","CD1C","CLEC10A"),
  DC_LAMP3        = c("LAMP3","ENOX1","IL15","CCL22","CD200"),
  DC_cycling      = c("CDC45","MKI67","BIRC5","DIAPH3","RRM2"),
  PlasmacytoidDC  = c("TCF4","IL3RA","COL24A1","GLT1D1","P2RY6","PTPRS"),
  LC              = c("CD207","CDH20","FCGBP","CLDN1","COL21A1"),
  Mast            = c("TPSB2","TPSAB1","CPA3","IL1RL1","MS4A2"),
  CD4T_helper     = c("CD40LG","TRAT1","BCL11B","CAMK4","CD6"),
  Treg            = c("IKZF2","CTLA4","TIGIT","CCDC141","FOXP3"),
  T_ECMhigh       = c("COL1A2","COL1A1","SPARC","LUM"),
  CD8T_cytotoxic  = c("GZMK","GZMA","CD8A","CD8B","CCL5"),
  NKcell          = c("NKG7","GNLY","KLRC1","KLRD1","SH2D1B"),
  ILCs            = c("XCL1","XCL2","PDE7B","DPF3","TNFSF4"),
  Bcell           = c("MS4A1","BANK1","EBF1","CD79A","PAX5"),
  Plasma          = c("IGHG1","IGHA1","IGHG3","MZB1","JCHAIN")
)


marker.list <- list(
  SMC          = c("MYH11","RGS6","NTRK2","PIP5K1B","FRAS1"),
  Pericyte     = c("COL6A3","NCKAP5","CCDC102B","GUCY1A2","THY1"),
  Vascu_stroma = c("DCN","TNFAIP6","THBS2","TWIST2","MFAP5"),
  LE           = c("CCL21","PROX1","PKHD1L1","SEMA3D","PIEZO2"),
  VE_art       = c("NEBL","PCSK5","IGFBP3","HEY1","SRGN"),
  VE_ven1      = c("SELE","ACKR1","COL15A1","HLA-DQA1","CSF3"),
  VE_ven2      = c("VWF","CCL14","ID1","THSD7A","MYRIP"),
  VE_capi      = c("RGCC","ARHGAP18","LXN","NOX4","BTNL9"),
  VE_cyclic    = c("DIAPH3","TOP2A","KNL1","RRM2","MKI67")
)

