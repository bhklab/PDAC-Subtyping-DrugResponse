


pam50_subtyping <- function(dataset, pam50_centroid ){
  
  #rownames(dataset)=dataset[1:nrow(dataset), 1]
  
  
  print(paste("Number of genes mapped out of 50 pam50 genes ",length(which(rownames(dataset) %in%  rownames(pam50_centroid))), sep=" " ))
  print("Genes not mapped")
  print(paste( setdiff(rownames(pam50_centroid), rownames(dataset) ), sep=" "))
  
  #not_mapped =  setdiff(xx[2:nrow(xx), 1],pcsi[1:nrow(pcsi), 1] )

  dataset_pam50=dataset[which(rownames(dataset) %in% rownames(pam50_centroid)),]
  dataset_pam50_scaled <-scale(t(data.matrix(dataset_pam50)))
  
  #clusters = cutree(hclust(as.dist(1-cor(t(dataset_pam50_scaled ), method="spearman"))), k=2) # get 5 clusters
  
  
  
  clusters = ConsensusClusterPlus(t(dataset_pam50_scaled),maxK=4,reps=100,pItem=1,pFeature=1,clusterAlg="hc",distance="spearman",seed=1262118388.71279,corUse="pairwise.complete.obs")[[2]][["consensusClass"]]
     clusters1 =  clusters
  
  
  
  dataset_centroid=sapply(unique(clusters), clust.centroid, dataset_pam50_scaled , clusters)
  
  pam50_cent = pam50_centroid[which(rownames(pam50_centroid) %in% rownames(dataset_centroid)),]
  
  dataset_centroid = dataset_centroid[sort(rownames(dataset_centroid)),]
  pam50_cent = pam50_cent[sort(rownames(pam50_cent)),]
  
  print(rownames(dataset_centroid) ==   rownames(pam50_cent))
  
  ### cluster 1 in pam50 is "BASAL" and cluster 2 is "CLASSICAL"  
  
 xxx= which(is.na(rowSums(dataset_centroid)))
  if(length(xxx) >0){
    dataset_centroid = dataset_centroid[-xxx,]
    pam50_cent =pam50_cent[-xxx,]
    }
  
  cc= cor(dataset_centroid, pam50_cent)
  
  ss= sapply(1:dim(cc)[1], function(x) which(cc[x,] == max(cc[x,])) )
  corlist = sapply(1:dim(cc)[1], function(x) max(cc[x,])) 
  clusters = mgsub( clusters, unique(clusters), names(ss))
  
  zz=cbind(names(ss), corlist)
  
  ret=list(dataset_centroid, names(clusters1) ,clusters,zz)
  
  return(ret)
}















