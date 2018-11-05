clust.centroid = function(i, dat, clusters) {
  ind = (clusters == i)
  
  if(length( which(ind == TRUE))==1){
    dat[ind,]
  }else{
    colMeans(dat[ind,])
  }}