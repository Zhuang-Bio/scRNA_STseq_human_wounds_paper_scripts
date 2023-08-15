#Functions
convertDPTperc = function(.sce){
  .sce = .sce[,order(.sce$dpt_pseudotime)]
  .sce$dpt_order = 1:ncol(.sce)
  .sce$dpt_order_perc = .sce$dpt_order / max(.sce$dpt_order)
  return(.sce)
}

testOneGene = function(.gene, .sce){
  dat = data.frame(exp=as.vector(assay(.sce[.gene,], "logcounts")), time=.sce$dpt_order_perc)
  dat = dat[!is.na(dat$exp),]
  #fit = gam(exp ~ lo(time), data=dat)
  fit = gam(exp ~ time, data=dat)
  summ = summary(fit)
  coef = unname(coef(fit)[2])
  p = summ[4][[1]][1,5]
  m_sq = summ[4][[1]][1,3]
  df = data.frame(geneSymbol=.gene, mean_sq=m_sq, coef=coef, pval=p)
  return(df)
}

testProcess = function(.sce, .tfs){
  # correaltion test
  gam.test = ldply(rownames(.sce), testOneGene, .sce, .parallel=T)
  colnames(gam.test) = c("geneSymbol","mean_sq", "coef", "pval")
  gam.test$geneSymbol = as.character(gam.test$geneSymbol)
  rownames(gam.test) = gam.test$geneSymbol
  gam.test = gam.test[order(gam.test$pval,decreasing=F),]
  gam.test$adj.p = p.adjust(gam.test$pval, method="BH")
  gam.test$is.TF = ifelse(as.character(gam.test$geneSymbol) %in% .tfs, T, F)
  return(gam.test)
}
