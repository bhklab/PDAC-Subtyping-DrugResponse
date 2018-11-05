bailey_subtyping <- function(dataset, bailey_centroid ){
  
  
  #rownames(dataset)=dataset[1:nrow(dataset), 1]
  print(paste("Number of genes mapped out of 613 bailey genes ",length(which(rownames(dataset) %in%  rownames(bailey_centroid))), sep=" " ))
  print("Genes not mapped")
  print(paste( setdiff(rownames(bailey_centroid), rownames(dataset) ), sep=" "))
  
  #not_mapped =  setdiff(xx[2:nrow(xx), 1],pcsi[1:nrow(pcsi), 1] )
  
  dataset_bailey=dataset[which(rownames(dataset) %in% rownames(bailey_centroid)),]
  dataset_bailey_scaled <-scale(t(data.matrix(dataset_bailey)))
  
  #clusters = cutree(hclust(as.dist(1-cor(t(dataset_bailey_scaled ), method="spearman"))), k=2) # get 5 clusters
  
  
  clusters = ConsensusClusterPlus(t(dataset_bailey_scaled),maxK=4,reps=100,pItem=0.8,pFeature=1,clusterAlg="hc",distance="spearman",seed=1262118388.71279,corUse="pairwise.complete.obs")[[3]][["consensusClass"]]
  clusters1 = clusters 
  if(length(which(is.na( dataset_bailey_scaled )))> 0){
    
    dataset_bailey_scaled= as.data.frame(dataset_bailey_scaled)
dataset_bailey_scaled[ is.na(dataset_bailey_scaled)] <- 0
  }
  
  dataset_centroid=sapply(unique(clusters), clust.centroid, dataset_bailey_scaled , clusters)
  
  bailey_cent = bailey_centroid[which(rownames(bailey_centroid) %in% rownames(dataset_centroid)),]
  dataset_centroid = dataset_centroid[sort(rownames(dataset_centroid)),]
  bailey_cent = bailey_cent[sort(rownames(bailey_cent)),]
  
  print(rownames(dataset_centroid) ==   rownames(bailey_cent))
  
  
  xxx= which(is.na(rowSums(dataset_centroid)))
  if(length(xxx) >0){
    dataset_centroid = dataset_centroid[-xxx,]
    bailey_cent =bailey_cent[-xxx,]
  }
  
  ### cluster 1 in bailey is "BASAL" and cluster 2 is "CLASSICAL"  
  cc= cor(dataset_centroid,bailey_cent)
  ss= sapply(1:dim(cc)[1], function(x) which(cc[x,] == max(cc[x,])) )
  corlist = sapply(1:dim(cc)[1], function(x) max(cc[x,])) 
  clusters = mgsub( clusters, unique(clusters), names(ss))
  
  zz=cbind(names(ss), corlist)
  
  ret=list(dataset_centroid, names(clusters1) ,clusters,zz)
  
  return(ret)
  
  
}















