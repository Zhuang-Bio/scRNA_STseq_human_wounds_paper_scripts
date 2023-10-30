library(shiny) 
library(shinyhelper) 
library(data.table) 
library(Matrix) 
library(DT) 
library(magrittr) 
library(ggplot2) 
library(ggrepel) 
library(hdf5r) 
library(ggdendro) 
library(gridExtra) 
sc1conf = readRDS("sc1conf.rds")
sc1def  = readRDS("sc1def.rds")
sc1gene = readRDS("sc1gene.rds")
sc1meta = readRDS("sc1meta.rds")



### Useful stuff 
# Colour palette 
cList = list(c("grey85","#FFF7EC","#FEE8C8","#FDD49E","#FDBB84", 
               "#FC8D59","#EF6548","#D7301F","#B30000","#7F0000"), 
             c("#4575B4","#74ADD1","#ABD9E9","#E0F3F8","#FFFFBF", 
               "#FEE090","#FDAE61","#F46D43","#D73027")[c(1,1:9,9)], 
             c("#FDE725","#AADC32","#5DC863","#27AD81","#21908C", 
               "#2C728E","#3B528B","#472D7B","#440154")) 
names(cList) = c("White-Red", "Blue-Yellow-Red", "Yellow-Green-Purple") 
 
# Panel sizes 
pList = c("400px", "600px", "800px") 
names(pList) = c("Small", "Medium", "Large") 
pList2 = c("500px", "700px", "900px") 
names(pList2) = c("Small", "Medium", "Large") 
pList3 = c("600px", "800px", "1000px") 
names(pList3) = c("Small", "Medium", "Large") 
sList = c(18,24,30) 
names(sList) = c("Small", "Medium", "Large") 
lList = c(5,6,7) 
names(lList) = c("Small", "Medium", "Large") 
 
# Function to extract legend 
g_legend <- function(a.gplot){  
  tmp <- ggplot_gtable(ggplot_build(a.gplot))  
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")  
  legend <- tmp$grobs[[leg]]  
  legend 
}  
 
# Plot theme 
sctheme <- function(base_size = 24, XYval = TRUE, Xang = 0, XjusH = 0.5){ 
  oupTheme = theme( 
    text =             element_text(size = base_size, family = "Helvetica"), 
    panel.background = element_rect(fill = "white", colour = NA), 
    axis.line =   element_line(colour = "black"), 
    axis.ticks =  element_line(colour = "black", size = base_size / 20), 
    axis.title =  element_text(face = "bold"), 
    axis.text =   element_text(size = base_size), 
    axis.text.x = element_text(angle = Xang, hjust = XjusH), 
    legend.position = "bottom", 
    legend.key =      element_rect(colour = NA, fill = NA) 
  ) 
  if(!XYval){ 
    oupTheme = oupTheme + theme( 
      axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
      axis.text.y = element_blank(), axis.ticks.y = element_blank()) 
  } 
  return(oupTheme) 
} 
 
### Common plotting functions 
# Plot cell information on dimred 
scDRcell <- function(inpConf, inpMeta, inpdrX, inpdrY, inp1, inpsub1, inpsub2, 
                     inpsiz, inpcol, inpord, inpfsz, inpasp, inptxt, inplab){ 
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]} 
  # Prepare ggData 
  ggData = inpMeta[, c(inpConf[UI == inpdrX]$ID, inpConf[UI == inpdrY]$ID, 
                       inpConf[UI == inp1]$ID, inpConf[UI == inpsub1]$ID),  
                   with = FALSE] 
  colnames(ggData) = c("X", "Y", "val", "sub") 
  rat = (max(ggData$X) - min(ggData$X)) / (max(ggData$Y) - min(ggData$Y)) 
  bgCells = FALSE 
  if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){ 
    bgCells = TRUE 
    ggData2 = ggData[!sub %in% inpsub2] 
    ggData = ggData[sub %in% inpsub2] 
  } 
  if(inpord == "Max-1st"){ 
    ggData = ggData[order(val)] 
  } else if(inpord == "Min-1st"){ 
    ggData = ggData[order(-val)] 
  } else if(inpord == "Random"){ 
    ggData = ggData[sample(nrow(ggData))] 
  } 
  
  # Do factoring if required 
  if(!is.na(inpConf[UI == inp1]$fCL)){ 
    ggCol = strsplit(inpConf[UI == inp1]$fCL, "\\|")[[1]] 
    names(ggCol) = levels(ggData$val) 
    ggLvl = levels(ggData$val)[levels(ggData$val) %in% unique(ggData$val)] 
    ggData$val = factor(ggData$val, levels = ggLvl) 
    ggCol = ggCol[ggLvl] 
  } 
 
  # Actual ggplot 
  ggOut = ggplot(ggData, aes(X, Y, color = val)) 
  if(bgCells){ 
    ggOut = ggOut + 
      geom_point(data = ggData2, color = "snow2", size = inpsiz, shape = 16) 
  } 
  ggOut = ggOut + 
    geom_point(size = inpsiz, shape = 16) + xlab(inpdrX) + ylab(inpdrY) + 
    sctheme(base_size = sList[inpfsz], XYval = inptxt) 
  if(is.na(inpConf[UI == inp1]$fCL)){ 
    ggOut = ggOut + scale_color_gradientn("", colours = cList[[inpcol]]) + 
      guides(color = guide_colorbar(barwidth = 15)) 
  } else { 
    sListX = min(nchar(paste0(levels(ggData$val), collapse = "")), 200) 
    sListX = 0.75 * (sList - (1.5 * floor(sListX/50))) 
    ggOut = ggOut + scale_color_manual("", values = ggCol) + 
      guides(color = guide_legend(override.aes = list(size = 5),  
                                  nrow = inpConf[UI == inp1]$fRow)) + 
      theme(legend.text = element_text(size = sListX[inpfsz])) 
    if(inplab){ 
      ggData3 = ggData[, .(X = mean(X), Y = mean(Y)), by = "val"] 
      lListX = min(nchar(paste0(ggData3$val, collapse = "")), 200) 
      lListX = lList - (0.25 * floor(lListX/50)) 
      ggOut = ggOut + 
        geom_text_repel(data = ggData3, aes(X, Y, label = val), 
                        color = "grey10", bg.color = "grey95", bg.r = 0.15, 
                        size = lListX[inpfsz], seed = 42) 
    } 
  } 
  if(inpasp == "Square") { 
    ggOut = ggOut + coord_fixed(ratio = rat) 
  } else if(inpasp == "Fixed") { 
    ggOut = ggOut + coord_fixed() 
  } 
  return(ggOut) 
} 
 
scDRnum <- function(inpConf, inpMeta, inp1, inp2, inpsub1, inpsub2, 
                    inpH5, inpGene, inpsplt){ 
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]} 
  # Prepare ggData 
  ggData = inpMeta[, c(inpConf[UI == inp1]$ID, inpConf[UI == inpsub1]$ID), 
                   with = FALSE] 
  colnames(ggData) = c("group", "sub") 
  h5file <- H5File$new(inpH5, mode = "r") 
  h5data <- h5file[["grp"]][["data"]] 
  ggData$val2 = h5data$read(args = list(inpGene[inp2], quote(expr=))) 
  ggData[val2 < 0]$val2 = 0 
  h5file$close_all() 
  if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){ 
    ggData = ggData[sub %in% inpsub2] 
  } 
  
  # Split inp1 if necessary 
  if(is.na(inpConf[UI == inp1]$fCL)){ 
    if(inpsplt == "Quartile"){nBk = 4} 
    if(inpsplt == "Decile"){nBk = 10} 
    ggData$group = cut(ggData$group, breaks = nBk) 
  } 
  
  # Actual data.table 
  ggData$express = FALSE 
  ggData[val2 > 0]$express = TRUE 
  ggData1 = ggData[express == TRUE, .(nExpress = .N), by = "group"] 
  ggData = ggData[, .(nCells = .N), by = "group"] 
  ggData = ggData1[ggData, on = "group"] 
  ggData = ggData[, c("group", "nCells", "nExpress"), with = FALSE] 
  ggData[is.na(nExpress)]$nExpress = 0 
  ggData$pctExpress = 100 * ggData$nExpress / ggData$nCells 
  ggData = ggData[order(group)] 
  colnames(ggData)[3] = paste0(colnames(ggData)[3], "_", inp2) 
  return(ggData) 
} 
# Plot gene expression on dimred 
scDRgene <- function(inpConf, inpMeta, inpdrX, inpdrY, inp1, inpsub1, inpsub2, 
                     inpH5, inpGene, 
                     inpsiz, inpcol, inpord, inpfsz, inpasp, inptxt){ 
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]} 
  # Prepare ggData 
  ggData = inpMeta[, c(inpConf[UI == inpdrX]$ID, inpConf[UI == inpdrY]$ID, 
                       inpConf[UI == inpsub1]$ID),  
                   with = FALSE] 
  colnames(ggData) = c("X", "Y", "sub") 
  rat = (max(ggData$X) - min(ggData$X)) / (max(ggData$Y) - min(ggData$Y)) 
  
  h5file <- H5File$new(inpH5, mode = "r") 
  h5data <- h5file[["grp"]][["data"]] 
  ggData$val = h5data$read(args = list(inpGene[inp1], quote(expr=))) 
  ggData[val < 0]$val = 0 
  h5file$close_all() 
  bgCells = FALSE 
  if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){ 
    bgCells = TRUE 
    ggData2 = ggData[!sub %in% inpsub2] 
    ggData = ggData[sub %in% inpsub2] 
  } 
  if(inpord == "Max-1st"){ 
    ggData = ggData[order(val)] 
  } else if(inpord == "Min-1st"){ 
    ggData = ggData[order(-val)] 
  } else if(inpord == "Random"){ 
    ggData = ggData[sample(nrow(ggData))] 
  } 
   
  # Actual ggplot 
  ggOut = ggplot(ggData, aes(X, Y, color = val)) 
  if(bgCells){ 
    ggOut = ggOut + 
      geom_point(data = ggData2, color = "snow2", size = inpsiz, shape = 16) 
  } 
  ggOut = ggOut + 
    geom_point(size = inpsiz, shape = 16) + xlab(inpdrX) + ylab(inpdrY) + 
    sctheme(base_size = sList[inpfsz], XYval = inptxt) +  
    scale_color_gradientn(inp1, colours = cList[[inpcol]]) + 
      guides(color = guide_colorbar(barwidth = 15)) 
  if(inpasp == "Square") { 
    ggOut = ggOut + coord_fixed(ratio = rat) 
  } else if(inpasp == "Fixed") { 
    ggOut = ggOut + coord_fixed() 
  } 
  return(ggOut) 
} 
 
# Plot gene coexpression on dimred 
bilinear <- function(x,y,xy,Q11,Q21,Q12,Q22){ 
  oup = (xy-x)*(xy-y)*Q11 + x*(xy-y)*Q21 + (xy-x)*y*Q12 + x*y*Q22 
  oup = oup / (xy*xy) 
  return(oup) 
} 
scDRcoex <- function(inpConf, inpMeta, inpdrX, inpdrY, inp1, inp2, 
                     inpsub1, inpsub2, inpH5, inpGene, 
                     inpsiz, inpcol, inpord, inpfsz, inpasp, inptxt){ 
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]} 
  # Prepare ggData 
  ggData = inpMeta[, c(inpConf[UI == inpdrX]$ID, inpConf[UI == inpdrY]$ID, 
                       inpConf[UI == inpsub1]$ID),  
                   with = FALSE] 
  colnames(ggData) = c("X", "Y", "sub") 
  rat = (max(ggData$X) - min(ggData$X)) / (max(ggData$Y) - min(ggData$Y)) 
  
  h5file <- H5File$new(inpH5, mode = "r") 
  h5data <- h5file[["grp"]][["data"]] 
  ggData$val1 = h5data$read(args = list(inpGene[inp1], quote(expr=))) 
  ggData[val1 < 0]$val1 = 0 
  ggData$val2 = h5data$read(args = list(inpGene[inp2], quote(expr=))) 
  ggData[val2 < 0]$val2 = 0 
  h5file$close_all() 
  bgCells = FALSE 
  if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){ 
    bgCells = TRUE 
    ggData2 = ggData[!sub %in% inpsub2] 
    ggData = ggData[sub %in% inpsub2] 
  } 
  
  # Generate coex color palette 
  cInp = strsplit(inpcol, "; ")[[1]] 
  if(cInp[1] == "Red (Gene1)"){ 
    c10 = c(255,0,0) 
  } else if(cInp[1] == "Orange (Gene1)"){ 
    c10 = c(255,140,0) 
  } else { 
    c10 = c(0,255,0) 
  } 
  if(cInp[2] == "Green (Gene2)"){ 
    c01 = c(0,255,0) 
  } else { 
    c01 = c(0,0,255) 
  } 
  c00 = c(217,217,217) ; c11 = c10 + c01 
  nGrid = 16; nPad = 2; nTot = nGrid + nPad * 2 
  gg = data.table(v1 = rep(0:nTot,nTot+1), v2 = sort(rep(0:nTot,nTot+1))) 
  gg$vv1 = gg$v1 - nPad ; gg[vv1 < 0]$vv1 = 0; gg[vv1 > nGrid]$vv1 = nGrid 
  gg$vv2 = gg$v2 - nPad ; gg[vv2 < 0]$vv2 = 0; gg[vv2 > nGrid]$vv2 = nGrid 
  gg$cR = bilinear(gg$vv1, gg$vv2, nGrid, c00[1], c10[1], c01[1], c11[1]) 
  gg$cG = bilinear(gg$vv1, gg$vv2, nGrid, c00[2], c10[2], c01[2], c11[2]) 
  gg$cB = bilinear(gg$vv1, gg$vv2, nGrid, c00[3], c10[3], c01[3], c11[3]) 
  gg$cMix = rgb(gg$cR, gg$cG, gg$cB, maxColorValue = 255) 
  gg = gg[, c("v1", "v2", "cMix")] 
  
  # Map colours 
  ggData$v1 = round(nTot * ggData$val1 / max(ggData$val1)) 
  ggData$v2 = round(nTot * ggData$val2 / max(ggData$val2)) 
  ggData$v0 = ggData$v1 + ggData$v2 
  ggData = gg[ggData, on = c("v1", "v2")] 
  if(inpord == "Max-1st"){ 
    ggData = ggData[order(v0)] 
  } else if(inpord == "Min-1st"){ 
    ggData = ggData[order(-v0)] 
  } else if(inpord == "Random"){ 
    ggData = ggData[sample(nrow(ggData))] 
  } 
  
  # Actual ggplot 
  ggOut = ggplot(ggData, aes(X, Y)) 
  if(bgCells){ 
    ggOut = ggOut + 
      geom_point(data = ggData2, color = "snow2", size = inpsiz, shape = 16) 
  } 
  ggOut = ggOut + 
    geom_point(size = inpsiz, shape = 16, color = ggData$cMix) + 
    xlab(inpdrX) + ylab(inpdrY) + 
    sctheme(base_size = sList[inpfsz], XYval = inptxt) + 
    scale_color_gradientn(inp1, colours = cList[[1]]) + 
    guides(color = guide_colorbar(barwidth = 15)) 
  if(inpasp == "Square") { 
    ggOut = ggOut + coord_fixed(ratio = rat) 
  } else if(inpasp == "Fixed") { 
    ggOut = ggOut + coord_fixed() 
  } 
  return(ggOut) 
} 
 
scDRcoexLeg <- function(inp1, inp2, inpcol, inpfsz){ 
  # Generate coex color palette 
  cInp = strsplit(inpcol, "; ")[[1]] 
  if(cInp[1] == "Red (Gene1)"){ 
    c10 = c(255,0,0) 
  } else if(cInp[1] == "Orange (Gene1)"){ 
    c10 = c(255,140,0) 
  } else { 
    c10 = c(0,255,0) 
  } 
  if(cInp[2] == "Green (Gene2)"){ 
    c01 = c(0,255,0) 
  } else { 
    c01 = c(0,0,255) 
  } 
  c00 = c(217,217,217) ; c11 = c10 + c01 
  nGrid = 16; nPad = 2; nTot = nGrid + nPad * 2 
  gg = data.table(v1 = rep(0:nTot,nTot+1), v2 = sort(rep(0:nTot,nTot+1))) 
  gg$vv1 = gg$v1 - nPad ; gg[vv1 < 0]$vv1 = 0; gg[vv1 > nGrid]$vv1 = nGrid 
  gg$vv2 = gg$v2 - nPad ; gg[vv2 < 0]$vv2 = 0; gg[vv2 > nGrid]$vv2 = nGrid 
  gg$cR = bilinear(gg$vv1, gg$vv2, nGrid, c00[1], c10[1], c01[1], c11[1]) 
  gg$cG = bilinear(gg$vv1, gg$vv2, nGrid, c00[2], c10[2], c01[2], c11[2]) 
  gg$cB = bilinear(gg$vv1, gg$vv2, nGrid, c00[3], c10[3], c01[3], c11[3]) 
  gg$cMix = rgb(gg$cR, gg$cG, gg$cB, maxColorValue = 255) 
  gg = gg[, c("v1", "v2", "cMix")] 
  
  # Actual ggplot 
  ggOut = ggplot(gg, aes(v1, v2)) + 
    geom_tile(fill = gg$cMix) + 
    xlab(inp1) + ylab(inp2) + coord_fixed(ratio = 1) + 
    scale_x_continuous(breaks = c(0, nTot), label = c("low", "high")) + 
    scale_y_continuous(breaks = c(0, nTot), label = c("low", "high")) + 
    sctheme(base_size = sList[inpfsz], XYval = TRUE) 
  return(ggOut) 
} 
 
scDRcoexNum <- function(inpConf, inpMeta, inp1, inp2, 
                        inpsub1, inpsub2, inpH5, inpGene){ 
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]} 
  # Prepare ggData 
  ggData = inpMeta[, c(inpConf[UI == inpsub1]$ID), with = FALSE] 
  colnames(ggData) = c("sub") 
  h5file <- H5File$new(inpH5, mode = "r") 
  h5data <- h5file[["grp"]][["data"]] 
  ggData$val1 = h5data$read(args = list(inpGene[inp1], quote(expr=))) 
  ggData[val1 < 0]$val1 = 0 
  ggData$val2 = h5data$read(args = list(inpGene[inp2], quote(expr=))) 
  ggData[val2 < 0]$val2 = 0 
  h5file$close_all() 
  if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){ 
    ggData = ggData[sub %in% inpsub2] 
  } 
  
  # Actual data.table 
  ggData$express = "none" 
  ggData[val1 > 0]$express = inp1 
  ggData[val2 > 0]$express = inp2 
  ggData[val1 > 0 & val2 > 0]$express = "both" 
  ggData$express = factor(ggData$express, levels = unique(c("both", inp1, inp2, "none"))) 
  ggData = ggData[, .(nCells = .N), by = "express"] 
  ggData$percent = 100 * ggData$nCells / sum(ggData$nCells) 
  ggData = ggData[order(express)] 
  colnames(ggData)[1] = "expression > 0" 
  return(ggData) 
} 
 
# Plot violin / boxplot 
scVioBox <- function(inpConf, inpMeta, inp1, inp2, 
                     inpsub1, inpsub2, inpH5, inpGene, 
                     inptyp, inppts, inpsiz, inpfsz){ 
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]} 
  # Prepare ggData 
  ggData = inpMeta[, c(inpConf[UI == inp1]$ID, inpConf[UI == inpsub1]$ID), 
                   with = FALSE] 
  colnames(ggData) = c("X", "sub") 
  
  # Load in either cell meta or gene expr
  if(inp2 %in% inpConf$UI){ 
    ggData$val = inpMeta[[inpConf[UI == inp2]$ID]] 
  } else { 
    h5file <- H5File$new(inpH5, mode = "r") 
    h5data <- h5file[["grp"]][["data"]] 
    ggData$val = h5data$read(args = list(inpGene[inp2], quote(expr=))) 
    ggData[val < 0]$val = 0 
    set.seed(42) 
    tmpNoise = rnorm(length(ggData$val)) * diff(range(ggData$val)) / 1000 
    ggData$val = ggData$val + tmpNoise 
    h5file$close_all() 
  } 
  if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){ 
    ggData = ggData[sub %in% inpsub2] 
  } 
  
  # Do factoring 
  ggCol = strsplit(inpConf[UI == inp1]$fCL, "\\|")[[1]] 
  names(ggCol) = levels(ggData$X) 
  ggLvl = levels(ggData$X)[levels(ggData$X) %in% unique(ggData$X)] 
  ggData$X = factor(ggData$X, levels = ggLvl) 
  ggCol = ggCol[ggLvl] 
  
  # Actual ggplot 
  if(inptyp == "violin"){ 
    ggOut = ggplot(ggData, aes(X, val, fill = X)) + geom_violin(scale = "width") 
  } else { 
    ggOut = ggplot(ggData, aes(X, val, fill = X)) + geom_boxplot() 
  } 
  if(inppts){ 
    ggOut = ggOut + geom_jitter(size = inpsiz, shape = 16) 
  } 
  ggOut = ggOut + xlab(inp1) + ylab(inp2) + 
    sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1) +  
    scale_fill_manual("", values = ggCol) +
    theme(legend.position = "none")
  return(ggOut) 
} 
 
# Plot proportion plot 
scProp <- function(inpConf, inpMeta, inp1, inp2, inpsub1, inpsub2, 
                   inptyp, inpflp, inpfsz){ 
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]} 
  # Prepare ggData 
  ggData = inpMeta[, c(inpConf[UI == inp1]$ID, inpConf[UI == inp2]$ID, 
                       inpConf[UI == inpsub1]$ID),  
                   with = FALSE] 
  colnames(ggData) = c("X", "grp", "sub") 
  if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){ 
    ggData = ggData[sub %in% inpsub2] 
  } 
  ggData = ggData[, .(nCells = .N), by = c("X", "grp")] 
  ggData = ggData[, {tot = sum(nCells) 
                      .SD[,.(pctCells = 100 * sum(nCells) / tot, 
                             nCells = nCells), by = "grp"]}, by = "X"] 
  
  # Do factoring 
  ggCol = strsplit(inpConf[UI == inp2]$fCL, "\\|")[[1]] 
  names(ggCol) = levels(ggData$grp) 
  ggLvl = levels(ggData$grp)[levels(ggData$grp) %in% unique(ggData$grp)] 
  ggData$grp = factor(ggData$grp, levels = ggLvl) 
  ggCol = ggCol[ggLvl] 
  
  # Actual ggplot 
  if(inptyp == "Proportion"){ 
    ggOut = ggplot(ggData, aes(X, pctCells, fill = grp)) + 
      geom_col() + ylab("Cell Proportion (%)") 
  } else { 
    ggOut = ggplot(ggData, aes(X, nCells, fill = grp)) + 
      geom_col() + ylab("Number of Cells") 
  } 
  if(inpflp){ 
    ggOut = ggOut + coord_flip() 
  } 
  ggOut = ggOut + xlab(inp1) + 
    sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1) +  
    scale_fill_manual("", values = ggCol) + 
    theme(legend.position = "right") 
  return(ggOut) 
} 
 
# Get gene list 
scGeneList <- function(inp, inpGene){ 
  geneList = data.table(gene = unique(trimws(strsplit(inp, ",|;|
")[[1]])), 
                        present = TRUE) 
  geneList[!gene %in% names(inpGene)]$present = FALSE 
  return(geneList) 
} 
 
# Plot gene expression bubbleplot / heatmap 
scBubbHeat <- function(inpConf, inpMeta, inp, inpGrp, inpPlt, 
                       inpsub1, inpsub2, inpH5, inpGene, inpScl, inpRow, inpCol, 
                       inpcols, inpfsz, save = FALSE){ 
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]} 
  # Identify genes that are in our dataset 
  geneList = scGeneList(inp, inpGene) 
  geneList = geneList[present == TRUE] 
  shiny::validate(need(nrow(geneList) <= 50, "More than 50 genes to plot! Please reduce the gene list!")) 
  shiny::validate(need(nrow(geneList) > 1, "Please input at least 2 genes to plot!")) 
   
  # Prepare ggData 
  h5file <- H5File$new(inpH5, mode = "r") 
  h5data <- h5file[["grp"]][["data"]] 
  ggData = data.table() 
  for(iGene in geneList$gene){ 
    tmp = inpMeta[, c("SampleID", inpConf[UI == inpsub1]$ID), with = FALSE] 
    colnames(tmp) = c("SampleID", "sub") 
    tmp$grpBy = inpMeta[[inpConf[UI == inpGrp]$ID]] 
    tmp$geneName = iGene 
    tmp$val = h5data$read(args = list(inpGene[iGene], quote(expr=))) 
    ggData = rbindlist(list(ggData, tmp)) 
  } 
  h5file$close_all() 
  if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){ 
    ggData = ggData[sub %in% inpsub2] 
  } 
  shiny::validate(need(uniqueN(ggData$grpBy) > 1, "Only 1 group present, unable to plot!")) 
   
  # Aggregate 
  ggData$val = expm1(ggData$val) 
  ggData = ggData[, .(val = mean(val), prop = sum(val>0) / length(SampleID)), 
                  by = c("geneName", "grpBy")] 
  ggData$val = log1p(ggData$val) 
   
  # Scale if required 
  colRange = range(ggData$val) 
  if(inpScl){ 
    ggData[, val:= scale(val), keyby = "geneName"] 
    colRange = c(-max(abs(range(ggData$val))), max(abs(range(ggData$val)))) 
  } 
   
  # hclust row/col if necessary 
  ggMat = dcast.data.table(ggData, geneName~grpBy, value.var = "val") 
  tmp = ggMat$geneName 
  ggMat = as.matrix(ggMat[, -1]) 
  rownames(ggMat) = tmp 
  if(inpRow){ 
    hcRow = dendro_data(as.dendrogram(hclust(dist(ggMat)))) 
    ggRow = ggplot() + coord_flip() + 
      geom_segment(data = hcRow$segments, aes(x=x,y=y,xend=xend,yend=yend)) + 
      scale_y_continuous(breaks = rep(0, uniqueN(ggData$grpBy)), 
                         labels = unique(ggData$grpBy), expand = c(0, 0)) + 
      scale_x_continuous(breaks = seq_along(hcRow$labels$label), 
                         labels = hcRow$labels$label, expand = c(0, 0.5)) + 
      sctheme(base_size = sList[inpfsz]) + 
      theme(axis.title = element_blank(), axis.line = element_blank(), 
            axis.ticks = element_blank(), axis.text.y = element_blank(), 
            axis.text.x = element_text(color="white", angle = 45, hjust = 1)) 
    ggData$geneName = factor(ggData$geneName, levels = hcRow$labels$label) 
  } else { 
    ggData$geneName = factor(ggData$geneName, levels = rev(geneList$gene)) 
  } 
  if(inpCol){ 
    hcCol = dendro_data(as.dendrogram(hclust(dist(t(ggMat))))) 
    ggCol = ggplot() + 
      geom_segment(data = hcCol$segments, aes(x=x,y=y,xend=xend,yend=yend)) + 
      scale_x_continuous(breaks = seq_along(hcCol$labels$label), 
                         labels = hcCol$labels$label, expand = c(0.05, 0)) + 
      scale_y_continuous(breaks = rep(0, uniqueN(ggData$geneName)), 
                         labels = unique(ggData$geneName), expand=c(0,0)) + 
      sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1) + 
      theme(axis.title = element_blank(), axis.line = element_blank(), 
            axis.ticks = element_blank(), axis.text.x = element_blank(), 
            axis.text.y = element_text(color = "white")) 
    ggData$grpBy = factor(ggData$grpBy, levels = hcCol$labels$label) 
  } 
   
  # Actual plot according to plottype 
  if(inpPlt == "Bubbleplot"){ 
    # Bubbleplot 
    ggOut = ggplot(ggData, aes(grpBy, geneName, color = val, size = prop)) + 
      geom_point() +  
      sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1) +  
      scale_x_discrete(expand = c(0.05, 0)) +  
      scale_y_discrete(expand = c(0, 0.5)) + 
      scale_size_continuous("proportion", range = c(0, 8), 
                            limits = c(0, 1), breaks = c(0.00,0.25,0.50,0.75,1.00)) + 
      scale_color_gradientn("expression", limits = colRange, colours = cList[[inpcols]]) + 
      guides(color = guide_colorbar(barwidth = 15)) + 
      theme(axis.title = element_blank(), legend.box = "vertical") 
  } else { 
    # Heatmap 
    ggOut = ggplot(ggData, aes(grpBy, geneName, fill = val)) + 
      geom_tile() +  
      sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1) + 
      scale_x_discrete(expand = c(0.05, 0)) +  
      scale_y_discrete(expand = c(0, 0.5)) + 
      scale_fill_gradientn("expression", limits = colRange, colours = cList[[inpcols]]) + 
      guides(fill = guide_colorbar(barwidth = 15)) + 
      theme(axis.title = element_blank()) 
  } 
     
  # Final tidy 
  ggLeg = g_legend(ggOut) 
  ggOut = ggOut + theme(legend.position = "none") 
  if(!save){ 
    if(inpRow & inpCol){ggOut =  
      grid.arrange(ggOut, ggLeg, ggCol, ggRow, widths = c(7,1), heights = c(1,7,2),  
                   layout_matrix = rbind(c(3,NA),c(1,4),c(2,NA)))  
    } else if(inpRow){ggOut =  
      grid.arrange(ggOut, ggLeg, ggRow, widths = c(7,1), heights = c(7,2),  
                   layout_matrix = rbind(c(1,3),c(2,NA)))  
    } else if(inpCol){ggOut =  
      grid.arrange(ggOut, ggLeg, ggCol, heights = c(1,7,2),  
                   layout_matrix = rbind(c(3),c(1),c(2)))  
    } else {ggOut =  
      grid.arrange(ggOut, ggLeg, heights = c(7,2),  
                   layout_matrix = rbind(c(1),c(2)))  
    }  
  } else { 
    if(inpRow & inpCol){ggOut =  
      arrangeGrob(ggOut, ggLeg, ggCol, ggRow, widths = c(7,1), heights = c(1,7,2),  
                  layout_matrix = rbind(c(3,NA),c(1,4),c(2,NA)))  
    } else if(inpRow){ggOut =  
      arrangeGrob(ggOut, ggLeg, ggRow, widths = c(7,1), heights = c(7,2),  
                  layout_matrix = rbind(c(1,3),c(2,NA)))  
    } else if(inpCol){ggOut =  
      arrangeGrob(ggOut, ggLeg, ggCol, heights = c(1,7,2),  
                  layout_matrix = rbind(c(3),c(1),c(2)))  
    } else {ggOut =  
      arrangeGrob(ggOut, ggLeg, heights = c(7,2),  
                  layout_matrix = rbind(c(1),c(2)))  
    }  
  } 
  return(ggOut) 
} 
 
 
 
 
 
### Start server code 
shinyServer(function(input, output, session) { 
  ### For all tags and Server-side selectize 
  observe_helpers() 
 optCrt="{ option_create: function(data,escape) {return('<div class=\"create\"><strong>' + '</strong></div>');} }" 
  updateSelectizeInput(session, "sc1a1inp2", choices = names(sc1gene), server = TRUE, 
                       selected = sc1def$gene1, options = list( 
                         maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt))) 
  updateSelectizeInput(session, "sc1a3inp1", choices = names(sc1gene), server = TRUE, 
                       selected = sc1def$gene1, options = list( 
                         maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt))) 
  updateSelectizeInput(session, "sc1a3inp2", choices = names(sc1gene), server = TRUE, 
                       selected = sc1def$gene2, options = list( 
                         maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt))) 
  updateSelectizeInput(session, "sc1b2inp1", choices = names(sc1gene), server = TRUE, 
                       selected = sc1def$gene1, options = list( 
                         maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt))) 
  updateSelectizeInput(session, "sc1b2inp2", choices = names(sc1gene), server = TRUE, 
                       selected = sc1def$gene2, options = list( 
                         maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt))) 
  updateSelectizeInput(session, "sc1c1inp2", server = TRUE, 
                       choices = c(names(sc1gene)), 
                       selected = "FOSL1", options = list( 
                         maxOptions = 10, 
                         create = TRUE, persist = TRUE, render = I(optCrt))) 
 
 
  ### Plots for tab a1 
  output$sc1a1sub1.ui <- renderUI({ 
    sub = strsplit(sc1conf[UI == input$sc1a1sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("sc1a1sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$sc1a1sub1non, { 
    sub = strsplit(sc1conf[UI == input$sc1a1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "sc1a1sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$sc1a1sub1all, { 
    sub = strsplit(sc1conf[UI == input$sc1a1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "sc1a1sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$sc1a1oup1 <- renderPlot({ 
    scDRcell(sc1conf, sc1meta, input$sc1a1drX, input$sc1a1drY, input$sc1a1inp1,  
             input$sc1a1sub1, input$sc1a1sub2, 
             input$sc1a1siz, input$sc1a1col1, input$sc1a1ord1, 
             input$sc1a1fsz, input$sc1a1asp, input$sc1a1txt, input$sc1a1lab1) 
  }) 
  output$sc1a1oup1.ui <- renderUI({ 
    plotOutput("sc1a1oup1", height = pList[input$sc1a1psz]) 
  }) 
  output$sc1a1oup1.pdf <- downloadHandler( 
    filename = function() { paste0("sc1",input$sc1a1drX,"_",input$sc1a1drY,"_",  
                                   input$sc1a1inp1,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$sc1a1oup1.h, width = input$sc1a1oup1.w, useDingbats = FALSE, 
      plot = scDRcell(sc1conf, sc1meta, input$sc1a1drX, input$sc1a1drY, input$sc1a1inp1,   
                      input$sc1a1sub1, input$sc1a1sub2, 
                      input$sc1a1siz, input$sc1a1col1, input$sc1a1ord1,  
                      input$sc1a1fsz, input$sc1a1asp, input$sc1a1txt, input$sc1a1lab1) ) 
  }) 
  output$sc1a1oup1.png <- downloadHandler( 
    filename = function() { paste0("sc1",input$sc1a1drX,"_",input$sc1a1drY,"_",  
                                   input$sc1a1inp1,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$sc1a1oup1.h, width = input$sc1a1oup1.w, 
      plot = scDRcell(sc1conf, sc1meta, input$sc1a1drX, input$sc1a1drY, input$sc1a1inp1,   
                      input$sc1a1sub1, input$sc1a1sub2, 
                      input$sc1a1siz, input$sc1a1col1, input$sc1a1ord1,  
                      input$sc1a1fsz, input$sc1a1asp, input$sc1a1txt, input$sc1a1lab1) ) 
  }) 
  output$sc1a1.dt <- renderDataTable({ 
    ggData = scDRnum(sc1conf, sc1meta, input$sc1a1inp1, input$sc1a1inp2, 
                     input$sc1a1sub1, input$sc1a1sub2, 
                     "sc1gexpr.h5", sc1gene, input$sc1a1splt) 
    datatable(ggData, rownames = FALSE, extensions = "Buttons", 
              options = list(pageLength = -1, dom = "tB", buttons = c("copy", "csv", "excel"))) %>% 
      formatRound(columns = c("pctExpress"), digits = 2) 
  }) 
  
  output$sc1a1oup2 <- renderPlot({ 
    scDRgene(sc1conf, sc1meta, input$sc1a1drX, input$sc1a1drY, input$sc1a1inp2,  
             input$sc1a1sub1, input$sc1a1sub2, 
             "sc1gexpr.h5", sc1gene, 
             input$sc1a1siz, input$sc1a1col2, input$sc1a1ord2, 
             input$sc1a1fsz, input$sc1a1asp, input$sc1a1txt) 
  }) 
  output$sc1a1oup2.ui <- renderUI({ 
    plotOutput("sc1a1oup2", height = pList[input$sc1a1psz]) 
  }) 
  output$sc1a1oup2.pdf <- downloadHandler( 
    filename = function() { paste0("sc1",input$sc1a1drX,"_",input$sc1a1drY,"_",  
                                   input$sc1a1inp2,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$sc1a1oup2.h, width = input$sc1a1oup2.w, useDingbats = FALSE, 
      plot = scDRgene(sc1conf, sc1meta, input$sc1a1drX, input$sc1a1drY, input$sc1a1inp2,  
                      input$sc1a1sub1, input$sc1a1sub2, 
                      "sc1gexpr.h5", sc1gene, 
                      input$sc1a1siz, input$sc1a1col2, input$sc1a1ord2, 
                      input$sc1a1fsz, input$sc1a1asp, input$sc1a1txt) ) 
  }) 
  output$sc1a1oup2.png <- downloadHandler( 
    filename = function() { paste0("sc1",input$sc1a1drX,"_",input$sc1a1drY,"_",  
                                   input$sc1a1inp2,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$sc1a1oup2.h, width = input$sc1a1oup2.w, 
      plot = scDRgene(sc1conf, sc1meta, input$sc1a1drX, input$sc1a1drY, input$sc1a1inp2,  
                      input$sc1a1sub1, input$sc1a1sub2, 
                      "sc1gexpr.h5", sc1gene, 
                      input$sc1a1siz, input$sc1a1col2, input$sc1a1ord2, 
                      input$sc1a1fsz, input$sc1a1asp, input$sc1a1txt) ) 
  }) 
   
   
  ### Plots for tab a2 
  output$sc1a2sub1.ui <- renderUI({ 
    sub = strsplit(sc1conf[UI == input$sc1a2sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("sc1a2sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$sc1a2sub1non, { 
    sub = strsplit(sc1conf[UI == input$sc1a2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "sc1a2sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$sc1a2sub1all, { 
    sub = strsplit(sc1conf[UI == input$sc1a2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "sc1a2sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$sc1a2oup1 <- renderPlot({ 
    scDRcell(sc1conf, sc1meta, input$sc1a2drX, input$sc1a2drY, input$sc1a2inp1,  
             input$sc1a2sub1, input$sc1a2sub2, 
             input$sc1a2siz, input$sc1a2col1, input$sc1a2ord1, 
             input$sc1a2fsz, input$sc1a2asp, input$sc1a2txt, input$sc1a2lab1) 
  }) 
  output$sc1a2oup1.ui <- renderUI({ 
    plotOutput("sc1a2oup1", height = pList[input$sc1a2psz]) 
  }) 
  output$sc1a2oup1.pdf <- downloadHandler( 
    filename = function() { paste0("sc1",input$sc1a2drX,"_",input$sc1a2drY,"_",  
                                   input$sc1a2inp1,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$sc1a2oup1.h, width = input$sc1a2oup1.w, useDingbats = FALSE, 
      plot = scDRcell(sc1conf, sc1meta, input$sc1a2drX, input$sc1a2drY, input$sc1a2inp1,   
                      input$sc1a2sub1, input$sc1a2sub2, 
                      input$sc1a2siz, input$sc1a2col1, input$sc1a2ord1,  
                      input$sc1a2fsz, input$sc1a2asp, input$sc1a2txt, input$sc1a2lab1) ) 
  }) 
  output$sc1a2oup1.png <- downloadHandler( 
    filename = function() { paste0("sc1",input$sc1a2drX,"_",input$sc1a2drY,"_",  
                                   input$sc1a2inp1,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$sc1a2oup1.h, width = input$sc1a2oup1.w, 
      plot = scDRcell(sc1conf, sc1meta, input$sc1a2drX, input$sc1a2drY, input$sc1a2inp1,   
                      input$sc1a2sub1, input$sc1a2sub2, 
                      input$sc1a2siz, input$sc1a2col1, input$sc1a2ord1,  
                      input$sc1a2fsz, input$sc1a2asp, input$sc1a2txt, input$sc1a2lab1) ) 
  }) 
   
  output$sc1a2oup2 <- renderPlot({ 
    scDRcell(sc1conf, sc1meta, input$sc1a2drX, input$sc1a2drY, input$sc1a2inp2,  
             input$sc1a2sub1, input$sc1a2sub2, 
             input$sc1a2siz, input$sc1a2col2, input$sc1a2ord2, 
             input$sc1a2fsz, input$sc1a2asp, input$sc1a2txt, input$sc1a2lab2) 
  }) 
  output$sc1a2oup2.ui <- renderUI({ 
    plotOutput("sc1a2oup2", height = pList[input$sc1a2psz]) 
  }) 
  output$sc1a2oup2.pdf <- downloadHandler( 
    filename = function() { paste0("sc1",input$sc1a2drX,"_",input$sc1a2drY,"_",  
                                   input$sc1a2inp2,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$sc1a2oup2.h, width = input$sc1a2oup2.w, useDingbats = FALSE, 
      plot = scDRcell(sc1conf, sc1meta, input$sc1a2drX, input$sc1a2drY, input$sc1a2inp2,   
                      input$sc1a2sub1, input$sc1a2sub2, 
                      input$sc1a2siz, input$sc1a2col2, input$sc1a2ord2,  
                      input$sc1a2fsz, input$sc1a2asp, input$sc1a2txt, input$sc1a2lab2) ) 
  }) 
  output$sc1a2oup2.png <- downloadHandler( 
    filename = function() { paste0("sc1",input$sc1a2drX,"_",input$sc1a2drY,"_",  
                                   input$sc1a2inp2,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$sc1a2oup2.h, width = input$sc1a2oup2.w, 
      plot = scDRcell(sc1conf, sc1meta, input$sc1a2drX, input$sc1a2drY, input$sc1a2inp2,   
                      input$sc1a2sub1, input$sc1a2sub2, 
                      input$sc1a2siz, input$sc1a2col2, input$sc1a2ord2,  
                      input$sc1a2fsz, input$sc1a2asp, input$sc1a2txt, input$sc1a2lab2) ) 
  }) 
   
   
  ### Plots for tab a3 
  output$sc1a3sub1.ui <- renderUI({ 
    sub = strsplit(sc1conf[UI == input$sc1a3sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("sc1a3sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$sc1a3sub1non, { 
    sub = strsplit(sc1conf[UI == input$sc1a3sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "sc1a3sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$sc1a3sub1all, { 
    sub = strsplit(sc1conf[UI == input$sc1a3sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "sc1a3sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$sc1a3oup1 <- renderPlot({ 
    scDRgene(sc1conf, sc1meta, input$sc1a3drX, input$sc1a3drY, input$sc1a3inp1,  
             input$sc1a3sub1, input$sc1a3sub2, 
             "sc1gexpr.h5", sc1gene, 
             input$sc1a3siz, input$sc1a3col1, input$sc1a3ord1, 
             input$sc1a3fsz, input$sc1a3asp, input$sc1a3txt) 
  }) 
  output$sc1a3oup1.ui <- renderUI({ 
    plotOutput("sc1a3oup1", height = pList[input$sc1a3psz]) 
  }) 
  output$sc1a3oup1.pdf <- downloadHandler( 
    filename = function() { paste0("sc1",input$sc1a3drX,"_",input$sc1a3drY,"_",  
                                   input$sc1a3inp1,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$sc1a3oup1.h, width = input$sc1a3oup1.w, useDingbats = FALSE, 
      plot = scDRgene(sc1conf, sc1meta, input$sc1a3drX, input$sc1a3drY, input$sc1a3inp1,  
                      input$sc1a3sub1, input$sc1a3sub2, 
                      "sc1gexpr.h5", sc1gene, 
                      input$sc1a3siz, input$sc1a3col1, input$sc1a3ord1, 
                      input$sc1a3fsz, input$sc1a3asp, input$sc1a3txt) ) 
  }) 
  output$sc1a3oup1.png <- downloadHandler( 
    filename = function() { paste0("sc1",input$sc1a3drX,"_",input$sc1a3drY,"_",  
                                   input$sc1a3inp1,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$sc1a3oup1.h, width = input$sc1a3oup1.w, 
      plot = scDRgene(sc1conf, sc1meta, input$sc1a3drX, input$sc1a3drY, input$sc1a3inp1,  
                      input$sc1a3sub1, input$sc1a3sub2, 
                      "sc1gexpr.h5", sc1gene, 
                      input$sc1a3siz, input$sc1a3col1, input$sc1a3ord1, 
                      input$sc1a3fsz, input$sc1a3asp, input$sc1a3txt) ) 
  }) 
   
  output$sc1a3oup2 <- renderPlot({ 
    scDRgene(sc1conf, sc1meta, input$sc1a3drX, input$sc1a3drY, input$sc1a3inp2,  
             input$sc1a3sub1, input$sc1a3sub2, 
             "sc1gexpr.h5", sc1gene, 
             input$sc1a3siz, input$sc1a3col2, input$sc1a3ord2, 
             input$sc1a3fsz, input$sc1a3asp, input$sc1a3txt) 
  }) 
  output$sc1a3oup2.ui <- renderUI({ 
    plotOutput("sc1a3oup2", height = pList[input$sc1a3psz]) 
  }) 
  output$sc1a3oup2.pdf <- downloadHandler( 
    filename = function() { paste0("sc1",input$sc1a3drX,"_",input$sc1a3drY,"_",  
                                   input$sc1a3inp2,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$sc1a3oup2.h, width = input$sc1a3oup2.w, useDingbats = FALSE, 
      plot = scDRgene(sc1conf, sc1meta, input$sc1a3drX, input$sc1a3drY, input$sc1a3inp2,  
                      input$sc1a3sub1, input$sc1a3sub2, 
                      "sc1gexpr.h5", sc1gene, 
                      input$sc1a3siz, input$sc1a3col2, input$sc1a3ord2, 
                      input$sc1a3fsz, input$sc1a3asp, input$sc1a3txt) ) 
  }) 
  output$sc1a3oup2.png <- downloadHandler( 
    filename = function() { paste0("sc1",input$sc1a3drX,"_",input$sc1a3drY,"_",  
                                   input$sc1a3inp2,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$sc1a3oup2.h, width = input$sc1a3oup2.w, 
      plot = scDRgene(sc1conf, sc1meta, input$sc1a3drX, input$sc1a3drY, input$sc1a3inp2,  
                      input$sc1a3sub1, input$sc1a3sub2, 
                      "sc1gexpr.h5", sc1gene, 
                      input$sc1a3siz, input$sc1a3col2, input$sc1a3ord2, 
                      input$sc1a3fsz, input$sc1a3asp, input$sc1a3txt) ) 
  }) 
     
   
  ### Plots for tab b2 
  output$sc1b2sub1.ui <- renderUI({ 
    sub = strsplit(sc1conf[UI == input$sc1b2sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("sc1b2sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$sc1b2sub1non, { 
    sub = strsplit(sc1conf[UI == input$sc1b2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "sc1b2sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$sc1b2sub1all, { 
    sub = strsplit(sc1conf[UI == input$sc1b2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "sc1b2sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$sc1b2oup1 <- renderPlot({ 
    scDRcoex(sc1conf, sc1meta, input$sc1b2drX, input$sc1b2drY,   
             input$sc1b2inp1, input$sc1b2inp2, input$sc1b2sub1, input$sc1b2sub2, 
             "sc1gexpr.h5", sc1gene, 
             input$sc1b2siz, input$sc1b2col1, input$sc1b2ord1, 
             input$sc1b2fsz, input$sc1b2asp, input$sc1b2txt) 
  }) 
  output$sc1b2oup1.ui <- renderUI({ 
    plotOutput("sc1b2oup1", height = pList2[input$sc1b2psz]) 
  }) 
  output$sc1b2oup1.pdf <- downloadHandler( 
    filename = function() { paste0("sc1",input$sc1b2drX,"_",input$sc1b2drY,"_",  
                                    input$sc1b2inp1,"_",input$sc1b2inp2,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$sc1b2oup1.h, width = input$sc1b2oup1.w, useDingbats = FALSE, 
      plot = scDRcoex(sc1conf, sc1meta, input$sc1b2drX, input$sc1b2drY,  
                      input$sc1b2inp1, input$sc1b2inp2, input$sc1b2sub1, input$sc1b2sub2, 
                      "sc1gexpr.h5", sc1gene, 
                      input$sc1b2siz, input$sc1b2col1, input$sc1b2ord1, 
                      input$sc1b2fsz, input$sc1b2asp, input$sc1b2txt) ) 
  }) 
  output$sc1b2oup1.png <- downloadHandler( 
    filename = function() { paste0("sc1",input$sc1b2drX,"_",input$sc1b2drY,"_",  
                                    input$sc1b2inp1,"_",input$sc1b2inp2,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$sc1b2oup1.h, width = input$sc1b2oup1.w, 
      plot = scDRcoex(sc1conf, sc1meta, input$sc1b2drX, input$sc1b2drY,  
                      input$sc1b2inp1, input$sc1b2inp2, input$sc1b2sub1, input$sc1b2sub2, 
                      "sc1gexpr.h5", sc1gene, 
                      input$sc1b2siz, input$sc1b2col1, input$sc1b2ord1, 
                      input$sc1b2fsz, input$sc1b2asp, input$sc1b2txt) ) 
  }) 
  output$sc1b2oup2 <- renderPlot({ 
    scDRcoexLeg(input$sc1b2inp1, input$sc1b2inp2, input$sc1b2col1, input$sc1b2fsz) 
  }) 
  output$sc1b2oup2.ui <- renderUI({ 
    plotOutput("sc1b2oup2", height = "300px") 
  }) 
  output$sc1b2oup2.pdf <- downloadHandler( 
    filename = function() { paste0("sc1",input$sc1b2drX,"_",input$sc1b2drY,"_",  
                                    input$sc1b2inp1,"_",input$sc1b2inp2,"_leg.pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = 3, width = 4, useDingbats = FALSE, 
      plot = scDRcoexLeg(input$sc1b2inp1, input$sc1b2inp2, input$sc1b2col1, input$sc1b2fsz) ) 
  }) 
  output$sc1b2oup2.png <- downloadHandler( 
    filename = function() { paste0("sc1",input$sc1b2drX,"_",input$sc1b2drY,"_",  
                                    input$sc1b2inp1,"_",input$sc1b2inp2,"_leg.png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = 3, width = 4, 
      plot = scDRcoexLeg(input$sc1b2inp1, input$sc1b2inp2, input$sc1b2col1, input$sc1b2fsz) ) 
  }) 
  output$sc1b2.dt <- renderDataTable({ 
    ggData = scDRcoexNum(sc1conf, sc1meta, input$sc1b2inp1, input$sc1b2inp2, 
                         input$sc1b2sub1, input$sc1b2sub2, "sc1gexpr.h5", sc1gene) 
    datatable(ggData, rownames = FALSE, extensions = "Buttons", 
              options = list(pageLength = -1, dom = "tB", buttons = c("copy", "csv", "excel"))) %>% 
      formatRound(columns = c("percent"), digits = 2) 
  }) 
     
   
  ### Plots for tab c1 
  output$sc1c1sub1.ui <- renderUI({ 
    sub = strsplit(sc1conf[UI == input$sc1c1sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("sc1c1sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$sc1c1sub1non, { 
    sub = strsplit(sc1conf[UI == input$sc1c1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "sc1c1sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$sc1c1sub1all, { 
    sub = strsplit(sc1conf[UI == input$sc1c1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "sc1c1sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$sc1c1oup <- renderPlot({ 
    scVioBox(sc1conf, sc1meta, input$sc1c1inp1, input$sc1c1inp2, 
             input$sc1c1sub1, input$sc1c1sub2, 
             "sc1gexpr.h5", sc1gene, input$sc1c1typ, input$sc1c1pts, 
             input$sc1c1siz, input$sc1c1fsz) 
  }) 
  output$sc1c1oup.ui <- renderUI({ 
    plotOutput("sc1c1oup", height = pList2[input$sc1c1psz]) 
  }) 
  output$sc1c1oup.pdf <- downloadHandler( 
    filename = function() { paste0("sc1",input$sc1c1typ,"_",input$sc1c1inp1,"_",  
                                   input$sc1c1inp2,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$sc1c1oup.h, width = input$sc1c1oup.w, useDingbats = FALSE, 
      plot = scVioBox(sc1conf, sc1meta, input$sc1c1inp1, input$sc1c1inp2, 
                      input$sc1c1sub1, input$sc1c1sub2, 
                      "sc1gexpr.h5", sc1gene, input$sc1c1typ, input$sc1c1pts, 
                      input$sc1c1siz, input$sc1c1fsz) ) 
  }) 
  output$sc1c1oup.png <- downloadHandler( 
    filename = function() { paste0("sc1",input$sc1c1typ,"_",input$sc1c1inp1,"_",  
                                   input$sc1c1inp2,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$sc1c1oup.h, width = input$sc1c1oup.w, 
      plot = scVioBox(sc1conf, sc1meta, input$sc1c1inp1, input$sc1c1inp2, 
                      input$sc1c1sub1, input$sc1c1sub2, 
                      "sc1gexpr.h5", sc1gene, input$sc1c1typ, input$sc1c1pts, 
                      input$sc1c1siz, input$sc1c1fsz) ) 
  }) 
     
   
### Plots for tab c2 
  output$sc1c2sub1.ui <- renderUI({ 
    sub = strsplit(sc1conf[UI == input$sc1c2sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("sc1c2sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$sc1c2sub1non, { 
    sub = strsplit(sc1conf[UI == input$sc1c2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "sc1c2sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$sc1c2sub1all, { 
    sub = strsplit(sc1conf[UI == input$sc1c2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "sc1c2sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
output$sc1c2oup <- renderPlot({ 
  scProp(sc1conf, sc1meta, input$sc1c2inp1, input$sc1c2inp2,  
         input$sc1c2sub1, input$sc1c2sub2, 
         input$sc1c2typ, input$sc1c2flp, input$sc1c2fsz) 
}) 
output$sc1c2oup.ui <- renderUI({ 
  plotOutput("sc1c2oup", height = pList2[input$sc1c2psz]) 
}) 
output$sc1c2oup.pdf <- downloadHandler( 
  filename = function() { paste0("sc1",input$sc1c2typ,"_",input$sc1c2inp1,"_",  
                                 input$sc1c2inp2,".pdf") }, 
  content = function(file) { ggsave( 
    file, device = "pdf", height = input$sc1c2oup.h, width = input$sc1c2oup.w, useDingbats = FALSE, 
    plot = scProp(sc1conf, sc1meta, input$sc1c2inp1, input$sc1c2inp2,  
                  input$sc1c2sub1, input$sc1c2sub2, 
                  input$sc1c2typ, input$sc1c2flp, input$sc1c2fsz) ) 
  }) 
output$sc1c2oup.png <- downloadHandler( 
  filename = function() { paste0("sc1",input$sc1c2typ,"_",input$sc1c2inp1,"_",  
                                 input$sc1c2inp2,".png") }, 
  content = function(file) { ggsave( 
    file, device = "png", height = input$sc1c2oup.h, width = input$sc1c2oup.w, 
    plot = scProp(sc1conf, sc1meta, input$sc1c2inp1, input$sc1c2inp2,  
                  input$sc1c2sub1, input$sc1c2sub2, 
                  input$sc1c2typ, input$sc1c2flp, input$sc1c2fsz) ) 
  }) 
     
   
  ### Plots for tab d1 
  output$sc1d1sub1.ui <- renderUI({ 
    sub = strsplit(sc1conf[UI == input$sc1d1sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("sc1d1sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$sc1d1sub1non, { 
    sub = strsplit(sc1conf[UI == input$sc1d1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "sc1d1sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$sc1d1sub1all, { 
    sub = strsplit(sc1conf[UI == input$sc1d1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "sc1d1sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$sc1d1oupTxt <- renderUI({ 
    geneList = scGeneList(input$sc1d1inp, sc1gene) 
    if(nrow(geneList) > 50){ 
      HTML("More than 50 input genes! Please reduce the gene list!") 
    } else { 
      oup = paste0(nrow(geneList[present == TRUE]), " genes OK and will be plotted") 
      if(nrow(geneList[present == FALSE]) > 0){ 
        oup = paste0(oup, "<br/>", 
                     nrow(geneList[present == FALSE]), " genes not found (", 
                     paste0(geneList[present == FALSE]$gene, collapse = ", "), ")") 
      } 
      HTML(oup) 
    } 
  }) 
  output$sc1d1oup <- renderPlot({ 
    scBubbHeat(sc1conf, sc1meta, input$sc1d1inp, input$sc1d1grp, input$sc1d1plt, 
               input$sc1d1sub1, input$sc1d1sub2, "sc1gexpr.h5", sc1gene, 
               input$sc1d1scl, input$sc1d1row, input$sc1d1col, 
               input$sc1d1cols, input$sc1d1fsz) 
  }) 
  output$sc1d1oup.ui <- renderUI({ 
    plotOutput("sc1d1oup", height = pList3[input$sc1d1psz]) 
  }) 
  output$sc1d1oup.pdf <- downloadHandler( 
    filename = function() { paste0("sc1",input$sc1d1plt,"_",input$sc1d1grp,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$sc1d1oup.h, width = input$sc1d1oup.w, 
      plot = scBubbHeat(sc1conf, sc1meta, input$sc1d1inp, input$sc1d1grp, input$sc1d1plt, 
                        input$sc1d1sub1, input$sc1d1sub2, "sc1gexpr.h5", sc1gene, 
                        input$sc1d1scl, input$sc1d1row, input$sc1d1col, 
                        input$sc1d1cols, input$sc1d1fsz, save = TRUE) ) 
  }) 
  output$sc1d1oup.png <- downloadHandler( 
    filename = function() { paste0("sc1",input$sc1d1plt,"_",input$sc1d1grp,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$sc1d1oup.h, width = input$sc1d1oup.w, 
      plot = scBubbHeat(sc1conf, sc1meta, input$sc1d1inp, input$sc1d1grp, input$sc1d1plt, 
                        input$sc1d1sub1, input$sc1d1sub2, "sc1gexpr.h5", sc1gene, 
                        input$sc1d1scl, input$sc1d1row, input$sc1d1col, 
                        input$sc1d1cols, input$sc1d1fsz, save = TRUE) ) 
  }) 
   
   
   
}) 
 
 
 
 