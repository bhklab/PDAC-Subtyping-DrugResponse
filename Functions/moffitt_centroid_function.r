


moffitt_subtyping <- function(dataset, moff_centroid ){
  
 #rownames(dataset)=dataset[1:nrow(dataset), 1]


print(paste("Number of genes mapped out of 50 Moffitt genes ",length(which(rownames(dataset) %in%  rownames(moff_centroid))), sep=" " ))
print("Genes not mapped")
print(paste( setdiff(rownames(moff_centroid), rownames(dataset) ), sep=" "))

#not_mapped =  setdiff(xx[2:nrow(xx), 1],pcsi[1:nrow(pcsi), 1] )
  ### ATAD4 === PRR15L  ### LOC400573 not found ### ANXA8L2 ===ANXA8L1 ### CTSL2 == CTSV
  
dataset_moff=dataset[which(rownames(dataset) %in% rownames(moff_centroid)),]
dataset_moff_scaled <-scale(t(data.matrix(dataset_moff)))

#clusters = cutree(hclust(as.dist(1-cor(t(dataset_moff_scaled ), method="spearman"))), k=2) # get 5 clusters



clusters = ConsensusClusterPlus(t(dataset_moff_scaled),maxK=6,reps=100,pItem=0.8,pFeature=1,clusterAlg="hc",distance="spearman",seed=1262118388.71279,corUse="pairwise.complete.obs")[[2]][["consensusClass"]]

   clusters1 =  clusters


dataset_centroid=sapply(unique(clusters), clust.centroid, dataset_moff_scaled , clusters)
  
moff_cent = moff_centroid[which(rownames(moff_centroid) %in% rownames(dataset_centroid)),]
  
dataset_centroid = dataset_centroid[sort(rownames(dataset_centroid)),]
moff_cent = moff_cent[sort(rownames(moff_cent)),]

  print(rownames(dataset_centroid) ==   rownames(moff_cent))

### cluster 1 in moffitt is "BASAL" and cluster 2 is "CLASSICAL"  

  xxx= which(is.na(rowSums(dataset_centroid)))
  if(length(xxx) >0){
    dataset_centroid = dataset_centroid[-xxx,]
    moff_cent =moff_cent[-xxx,]
    }
  
  cc= cor(dataset_centroid,moff_cent)
  
  ss= sapply(1:dim(cc)[1], function(x) which(cc[x,] == max(cc[x,])) )
  corlist = sapply(1:dim(cc)[1], function(x) max(cc[x,])) 
  clusters = mgsub( clusters, unique(clusters), names(ss))
  
  zz=cbind(names(ss), corlist)
  
  ret=list(dataset_centroid, names(clusters1) ,clusters,zz)
  
  return(ret)
  


}















