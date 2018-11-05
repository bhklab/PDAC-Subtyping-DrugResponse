library(randomForest)
library(limma)

load("/Users/vandanasandhu/Desktop/Subtyping_PDACs/RS_New_Analysis/Datasets/gg.RData")
load("/Users/vandanasandhu/Desktop/Subtyping_PDACs/RS_New_Analysis/Datasets/binary_rf_model.RData")
load("/Users/vandanasandhu/Desktop/Subtyping_PDACs/RS_New_Analysis/Datasets/RF_response.RData")
load("/Users/vandanasandhu/Desktop/Subtyping_PDACs/RS_New_Analysis/Datasets/output1_RF.RData")
source("/Users/vandanasandhu/Desktop/Subtyping_PDACs/RS_New_Analysis/Functions/meta_subtypes_celllines.R")
source("/Users/vandanasandhu/Desktop/Subtyping_PDACs/Scripts/mgsub_function.R")
load("/Users/vandanasandhu/Desktop/Subtyping_PDACs/Datasets/pharmacogx_pdac_celllines1.RData")

##################################################################################################################
#### CELLLINES

dataset_format <- function(dataset){
  
  dataset1=avereps(dataset)
  dataset1=data.frame(dataset1)
  dataset2= sapply(dataset1, function(x) as.numeric(as.character(x)))
 rownames(dataset2)=rownames(dataset1)
  return(dataset2)
}

gcsi = dataset_format(panc_celllines1$GCSI)
gdsc = dataset_format(panc_celllines1$GDSC)
ccle = dataset_format(panc_celllines1$CCLE)
ctrpv2= dataset_format(panc_celllines1$CTRPV2)
#rownames(gcsi)[9189] ="MIA2"  ### Aliase for CTAGE5 present in metagene
#rownames(gcsi)[4363] ="APOC2"  ### Aliase for APC2 present in metagene



gcsi_subtypes = meta_subtype_cellines(gcsi) 
ccle_subtypes = meta_subtype_cellines(ccle) 
gdsc_subtypes = meta_subtype_cellines(gdsc) 
ctrpv2_subtypes = meta_subtype_cellines(ctrpv2) 

#### COMPARE CELLLINES SUBTYPING CALLS



aa=merge(gcsi_subtypes, ccle_subtypes, by.x="id", by.y="id")
length(which(as.character(aa[,2]) == as.character(aa[,6])))/dim(aa)[1]
aa[,1][which(as.character(aa[,2]) != as.character(aa[,6]))]
aa[which(as.character(aa[,2]) != as.character(aa[,6])),]

aa=merge( ccle_subtypes, ctrpv2_subtypes, by.x="id", by.y="id")
length(which(as.character(aa[,2]) == as.character(aa[,6])))/dim(aa)[1]
aa[,1][which(as.character(aa[,2]) != as.character(aa[,6]))]
aa[which(as.character(aa[,2]) != as.character(aa[,6])),]


aa=merge(gcsi_subtypes,  ctrpv2_subtypes, by.x="id", by.y="id")
length(which(as.character(aa[,2]) == as.character(aa[,6])))/dim(aa)[1]
aa[,1][which(as.character(aa[,2]) != as.character(aa[,6]))]
aa[which(as.character(aa[,2]) != as.character(aa[,6])),]

####### COMPASS SAMPLES 
load("/Users/vandanasandhu/Desktop/RS_Remodelling_TSP_project/COMPASS/log_transformed_compass.RData")
com_sur= read.csv("/Users/vandanasandhu/Desktop/RS_Remodelling_TSP_project/COMPASS/compass_survival", sep="\t", header = TRUE)

aa=meta_subtype_cellines(compass_mat) 
table(aa$subtypes)


####### PDXE clustering
data <- readRDS("/Users/vandanasandhu/Downloads/PDXE_PDAC_RNASeq.Rda")
library(Biobase)
dd=exprs(data)


pdxe_subtypes = meta_subtype_cellines(dataset_format(dd)) 
table(pdxe_subtypes$subtypes)


######################  PanCuRx Organoids and Xenografts #####################################################

xx= read.table("/Users/vandanasandhu/Desktop/Subtyping_PDACs/Datasets/PanCuRx_rna_O_X.txt", sep="\t", header = TRUE)

mat=xx[, 4:ncol(xx)]
log_mat= as.matrix(log2(mat+1))
rownames(log_mat) = xx[,2]
avg_mat=avereps(log_mat)
classes= meta_subtype_cellines(avg_mat) 
