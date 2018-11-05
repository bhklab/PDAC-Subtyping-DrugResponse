


collisson_subtyping <- function(dataset, collisson_centroid ){
  
  #rownames(dataset)=dataset[1:nrow(dataset), 1]
  
  
  print(paste("Number of genes mapped out of 63 collisson genes ",length(which(rownames(dataset) %in%  rownames(collisson_centroid))), sep=" " ))
  print("Genes not mapped")
  print(paste( setdiff(rownames(collisson_centroid), rownames(dataset) ), sep=" "))
  
  #not_mapped =  setdiff(xx[2:nrow(xx), 1],pcsi[1:nrow(pcsi), 1] )
  
  dataset_collisson=dataset[which(rownames(dataset) %in% rownames(collisson_centroid)),]
  dataset_collisson_scaled <-scale(t(data.matrix(dataset_collisson)))
  
  #clusters = cutree(hclust(as.dist(1-cor(t(dataset_collisson_scaled ), method="spearman"))), k=2) # get 5 clusters
  
  #dataset_collisson_scaled[is.na(dataset_collisson_scaled)] <- 0
  
  clusters = ConsensusClusterPlus(t(dataset_collisson_scaled),maxK=4,reps=100,pItem=0.8,pFeature=1,clusterAlg="hc",distance="spearman",seed=1262118388.71279,corUse="pairwise.complete.obs")[[3]][["consensusClass"]]
 
   clusters1 =  clusters
  
  dataset_centroid=sapply(unique(clusters), clust.centroid, dataset_collisson_scaled , clusters)
  
  collisson_cent = collisson_centroid[which(rownames(collisson_centroid) %in% rownames(dataset_centroid)),]
  
  dataset_centroid = dataset_centroid[sort(rownames(dataset_centroid)),]
  collisson_cent = collisson_cent[sort(rownames(collisson_cent)),]
  
  print(rownames(dataset_centroid) ==   rownames(collisson_cent))
  
  
  ### cluster 1 in bailey is "BASAL" and cluster 2 is "CLASSICAL"  
  
  xxx= which(is.na(rowSums(dataset_centroid)))
  if(length(xxx) >0){
    dataset_centroid = dataset_centroid[-xxx,]
    collisson_cent =collisson_cent[-xxx,]
    }
  
  cc= cor(dataset_centroid,collisson_cent)
  
  ss= sapply(1:dim(cc)[1], function(x) which(cc[x,] == max(cc[x,])) )
  corlist = sapply(1:dim(cc)[1], function(x) max(cc[x,])) 
  clusters = mgsub( clusters, unique(clusters), names(ss))
  
  zz=cbind(names(ss), corlist)
  
  ret=list(dataset_centroid, names(clusters1) ,clusters,zz)
  
  return(ret)
  
  
  
  
}















