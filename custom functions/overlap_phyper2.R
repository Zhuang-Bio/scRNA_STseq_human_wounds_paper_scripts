overlap_phyper2<-function(L1,L2,bg=length(unique(c(unlist(L1),unlist(L2)))), plot=FALSE, title="Overlap", remove.diag = FALSE, with.total = TRUE, to.file = NULL, silent=FALSE){
  # takes two list with indices and creates a matrix with pvalue for overlap between elements
  # phyper test uses all entries in L as background if bg is not specified.
  nL1<-length(L1)
  nL2<-length(L2)
  M<-mat.or.vec(nL1,nL2)
  P<-mat.or.vec(nL1,nL2)
  P[,]<-1
  for (i in 1:nL1){
    for (j in  1:nL2){
      M[i,j]<-length(intersect(L1[[i]],L2[[j]]))
      if (M[i,j] ==0) {
        P[i,j] <- 1
      }else{
        P[i,j]<-phyper(M[i,j]-1,length(L1[[i]]),bg-length(L1[[i]]),length(L2[[j]]), lower.tail = FALSE)
      }
    }
  }
  colnames(P)<-names(L2)
  rownames(P)<-names(L1)
  colnames(M)<-names(L2)
  rownames(M)<-names(L1)
  
  P[M==0]<-1
  # still may have values that are zero
  pseudo = min(P[P>0])*0.1
  
  if(remove.diag){
    diag(P) = NA
  }
  
  # ad column/row with total genes per list.
  # lower square is total unique genes in the 2 lists.
  if (with.total){
    nTotal1 = unlist(lapply(L1,length))
    nTotal2 = c(unlist(lapply(L2,length)),length(unique(c(unlist(L1),unlist(L2)))))
    M = cbind(M,nTotal1)
    M = rbind(M,nTotal2)
    P = cbind(P, rep(NA,nrow(P)))
    P = rbind(P, rep(NA,ncol(P)))     
    
  }
  suppressMessages(require(pheatmap))
  if (is.null(to.file)){
    pl = pheatmap(-log10(P+pseudo), cluster_rows = F, cluster_cols = F, display_numbers = M, fontsize = 20, main = title, silent = silent)
  }else{
    pheatmap(-log10(P+pseudo), cluster_rows = F, cluster_cols = F, display_numbers = M, fontsize = 20, main = title, filename = to.file, silent = TRUE)
  }
  return(plot=pl) #return(list(P=P,M=M, plot=pl)) #old results
}


pred_phyper<-function(L1,L2,bg=length(unique(c(unlist(L1),unlist(L2)))), cutoff=0.01, nM=3, with.cl.names=FALSE){
  # takes two list with gene names and creates a matrix with pvalue for overlap 
  # phyper test uses all entries in L as background if bg is not specified.
  # will predict for each group in L1, the one with best overlap to L2
  # predict "Unass" if phyper pvalue < cutoff and if no gene overlap is >=nM
  # with.cl.names: add in cluster names in the new annotation names
  
  nL1<-length(L1)
  nL2<-length(L2)
  M<-mat.or.vec(nL1,nL2)
  P<-mat.or.vec(nL1,nL2)
  P[,]<-1
  for (i in 1:nL1){
    for (j in  1:nL2){
      M[i,j]<-length(intersect(L1[[i]],L2[[j]]))
      if (M[i,j] < nM) {
        P[i,j] <- 1    
      }else{
        P[i,j]<-phyper(M[i,j]-1,length(L1[[i]]),bg-length(L1[[i]]),length(L2[[j]]), lower.tail = FALSE)
      }
    }
  }
  
  # find lowest pvalue
  lowest = apply(P,1,which.min)
  assign = names(L2)[lowest]
  assign[apply(P,1,min)>cutoff]="Unass"
  if (with.cl.names){
    assign = paste(names(L1),assign,sep=":")
  }
  names(assign)=names(L1)  
  
  return(assign)
}
