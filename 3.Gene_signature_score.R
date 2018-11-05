
load("/Users/vandanasandhu/Desktop/Subtyping_PDACs/RS_New_Analysis/Datasets/common_genes_cohorts_new.RData")
load("/Users/vandanasandhu/Desktop/Subtyping_PDACs/meta_clusters.RData")

signature_score <- function(gene_list, List_name){
  my_comparisons <- list( c("Basal","Classical"), c("Classical", "Exocrine"), c("Exocrine", "Basal"))
  
  gene_list_score=list()
  sample_names=list()
  
  for(i in 1:(length(cohorts1))){
    gene_list_score[[i]] = colMeans( t(scale(t(data.matrix(cohorts1[[i]][gene_list,])))), na.rm = TRUE)
    sample_names[[i]] = names(colMeans(cohorts1[[i]][gene_list,], na.rm = TRUE))
  }
  
  zz= data.frame(gene_list_score=unlist(gene_list_score), class=dd$meta_class)
  REM = which(zz$class %in% c(NA))
  zz$class = mgsub(c("1","2","3"),  c("Basal","Exocrine","Classical"), zz$class)
  
  zz=zz[-REM,]
  #zz= zz[-which(zz$class =="4"),]
  
  p <- ggboxplot(zz, x = "class", y = "gene_list_score",
                 color = "class", palette = "jco",ylab=List_name, 
                 add = "jitter", xlab = "", legend="none")
  
 p= p + stat_compare_means(comparisons = my_comparisons)
  
  
  
  return(p)
}

#1.##################### Cell cycle proliferation score 
pp= read.table("/Users/vandanasandhu/Desktop/Subtyping_PDACs/Datasets/CCP_datasets.txt", sep="\t")

signature_score(as.character(pp[,1]), "Proliferation signature") 

#2. ##################### Double Strand break repair
pp= read.table("/Users/vandanasandhu/Desktop/Subtyping_PDACs/Datasets/dsbr_signature.txt", sep="\t")
signature_score(as.character(pp[,1]), " Double Strand break repair") 

########################### Prolifertaion score
pp=c("CCNB1","UBE2C","BIRC5","CDC20","PTTG1","RRM2","MKI67","TYMS","CEP55","KNTC2","CDCA1")
signature_score(as.character(pp), "Proliferation signature") 

########################### Hypoxia
pp= read.table("/Users/vandanasandhu/Desktop/Subtyping_PDACs/Datasets/HYPOXIA_buffa.txt", sep="\t")
signature_score(as.character(pp[,1]), "Hypoxia signature") 

#3. ##################### iMMUNE CYTOSOLIC ACTIVVITY

pp=c("GZMA","PRF1")
signature_score(as.character(pp), "Immune signature") 

