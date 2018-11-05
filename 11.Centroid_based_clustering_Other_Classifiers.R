load("/Users/vandanasandhu/Desktop/Subtyping_PDACs/Datasets/pharmacogx_pdac_celllines1.RData")
load("/Users/vandanasandhu/Desktop/subtyping_using_other_classifier/DATA/bailey_centroid.RData")
load("/Users/vandanasandhu/Desktop/subtyping_using_other_classifier/DATA/pam50_centroid.RData")
load("/Users/vandanasandhu/Desktop/subtyping_using_other_classifier/DATA/collisson_centroid.RData")
load("/Users/vandanasandhu/Desktop/subtyping_using_other_classifier/DATA/moffitt_centroid.RData")

source("/Users/vandanasandhu/Desktop/subtyping_using_other_classifier/CODES/Functions/collisson_centroid_function.r")
source("/Users/vandanasandhu/Desktop/subtyping_using_other_classifier/CODES/Functions/bailey_centroid_function.R")
source("/Users/vandanasandhu/Desktop/subtyping_using_other_classifier/CODES/Functions/moffitt_centroid_function.r")
source("/Users/vandanasandhu/Desktop/subtyping_using_other_classifier/CODES/Functions/pam50_centroid_function.r")
source("/Users/vandanasandhu/Desktop/Subtyping_PDACs/RS_New_Analysis/Functions/meta_subtypes_celllines.R")
source("/Users/vandanasandhu/Desktop/subtyping_using_other_classifier/CODES/Functions/clust.centroid_function.R")

load("/Users/vandanasandhu/Desktop/Subtyping_PDACs/RS_New_Analysis/Datasets/binary_rf_model.RData")
load("/Users/vandanasandhu/Desktop/Subtyping_PDACs/RS_New_Analysis/Datasets/output1_RF.RData")
load("/Users/vandanasandhu/Desktop/Subtyping_PDACs/RS_New_Analysis/Datasets/RF_response.RData")
load("/Users/vandanasandhu/Desktop/Subtyping_PDACs/RS_New_Analysis/Datasets/gg.RData")


####### Pharmacogx celllines 
dataset_format <- function(dataset){
  
  dataset1=avereps(dataset)
  dataset1=data.frame(dataset1)
  dataset2= sapply(dataset1, function(x) as.numeric(as.character(x)))
  rownames(dataset2)=rownames(dataset1)
  return(dataset2)
}

##CELL-LINES
gcsi = dataset_format(panc_celllines1$GCSI)
gdsc = dataset_format(panc_celllines1$GDSC)
ccle = dataset_format(panc_celllines1$CCLE)
ctrpv2= dataset_format(panc_celllines1$CTRPV2)

##COMPASS
xx=read.table("/Users/vandanasandhu/Documents/Projects/Average_based_analysis_VALIDATION_COHORT/Datasets_Average_Based/compass_data.txt", sep="\t", header=TRUE)
compass=xx
rownames(compass)= compass[,1]
compass=compass[,-1]

##PDXE
data <- readRDS("/Users/vandanasandhu/Downloads/PDXE_PDAC_RNASeq.Rda")
library(Biobase)
PDXE=exprs(data)

### OICR PDXE
xx= read.table("/Users/vandanasandhu/Desktop/Subtyping_PDACs/Datasets/PanCuRx_rna_O_X.txt", sep="\t", header = TRUE)
mat=xx[, 4:ncol(xx)]
log_mat= as.matrix(log2(mat+1))
rownames(log_mat) = xx[,2]

OICR_PDXE =avereps(log_mat)


##### PDACs

load("/Users/vandanasandhu/Desktop/Subtyping_PDACs/RS_New_Analysis/Datasets/cohorts_new.RData")


#######################################################################################
#################################### Clustering using all the methods

dataset<-ctrpv2

dataset_collisson = collisson_subtyping(dataset,collisson_centroid)
dataset_bailey = bailey_subtyping(dataset,bailey_centroid)
dataset_moffitt = moffitt_subtyping(dataset,moff_centroid)
dataset_pam50 = pam50_subtyping(dataset,pam50_centroid)
dataset_metasubtypes = mgsub( meta_subtype_cellines(dataset)$subtypes  , c(1,2,3), c("Meta_Basal","Meta_Exocrine","Meta_Classical"))

################################# Making the heatmap 

Moffitt = mgsub(dataset_moffitt[[3]], c("Classical","Basal"), c("Moffitt_Classical", "Moffitt_Basal") )
Collisson = mgsub(dataset_collisson[[3]], c("Classical PDA","QM-PDA","Exocrine-like PDA"), c("Collisson_classical","Collisson-QM-PDA","Collisson-Exocrine-like PDA") )
PAM50 = dataset_pam50[[3]]
Bailey = mgsub( dataset_bailey[[3]], c("Pancreatic Progenitor","Squamous","ADEX","Immunogenic"), c("Bailey_PP", "Bailey_Squamous","Bailey_ADEX","Bailey_Immune") )
Metaclusters = dataset_metasubtypes

dat<- rbind(Metaclusters, Bailey,  Moffitt,Collisson, PAM50)
dat1<- data.frame(dat)
colnames(dat1) = colnames(dataset)
rownames(dat1) = c("1.Metaclusters", "2.Bailey","3.Moffitt", "4.Collisson", "5.PAM50" )

dat1$id= rownames(dat1)

colr= c(brewer.pal(n=12, "Paired"),"black","grey")
dat3 <- melt(dat1, id.var = 'id')
ggplot(dat3, aes(variable, id)) + geom_tile(aes(fill = value), colour = "white") + scale_fill_manual(values=colr)


################################# Confusion matrix
cc= cbind(Metaclusters, Bailey,  Moffitt,Collisson, PAM50)
zz=expand.grid(x =1:dim(cc)[2], y = 1:dim(cc)[2])
a1= sapply(1:dim(zz)[1], function(x) as.numeric(summary(assocstats(table(cc[,zz[x,1]], cc[,zz[x,2]])))[[2]][5]))
a2= sapply(1:dim(zz)[1], function(x) as.numeric(summary(assocstats(table(cc[,zz[x,1]], cc[,zz[x,2]])))[[1]]$p.value))
  
##MCC
#a1= sapply(1:dim(zz)[1], function(x)  mcc(as.factor(cc[,zz[x,1]]), as.factor( cc[,zz[x,2]]),nperm = 1000, setseed = 12345, nthread = 1)$estimate)
#a2= sapply(1:dim(zz)[1], function(x)  mcc(as.factor(cc[,zz[x,1]]), as.factor( cc[,zz[x,2]]),nperm = 1000, setseed = 12345, nthread = 1)$p.value)

cramer_matrix=matrix(a1, ncol=5, nrow=5)
p_value=matrix(a2, ncol=5, nrow=5)

new <- matrix(NA, nrow = 5, ncol = 5)
new[upper.tri(new)] <- cramer_matrix[upper.tri(cramer_matrix)]
new[lower.tri(new)] <- p_value[upper.tri(p_value)]

rownames(new) = rownames(dat1)
colnames(new)= rownames(new) 
new

###
colnames(cramer_matrix) = rownames(dat1)
rownames(cramer_matrix) = rownames(dat1)
cramer_matrix= round(cramer_matrix,2)
colfunc <- colorRampPalette(c("black", "white", "red"))
heatmap(cramer_matrix,margin=c(5,8),cexRow=1,cexCol =0.8)

# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}

# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

library(reshape2)

cormat <- reorder_cormat(cramer_matrix)
upper_tri <- get_upper_tri(cormat)
# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)

library(ggplot2)
ggheatmap <- ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Cramer's V index") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()

gg_heatmap =ggheatmap + 
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),legend.justification = c(1, 0),
    legend.position = c(0.5, 0.8),
    legend.direction = "horizontal") +
  guides(fill = guide_colorbar(barwidth = 10, barheight = 3,
                               title.position = "top", title.hjust = 0.5))



####

dd.row <- as.dendrogram(hclust(dist(cramer_matrix)))
row.ord <- order.dendrogram(dd.row)

dd.col <- as.dendrogram(hclust(dist(t(cramer_matrix))))
col.ord <- order.dendrogram(dd.col)

ddata_x <- dendro_data(dd.row)
ddata_y <- dendro_data(dd.col)

### Set up a blank theme
theme_none <- theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  axis.title.x = element_text(colour=NA),
  axis.title.y = element_blank(),
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.line = element_blank()
  #axis.ticks.length = element_blank()
)
# Dendrogram 2
p3 <- ggplot(segment(ddata_y)) + 
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) + 
  coord_flip() + theme_none

grid.newpage()
print(gg_heatmap , vp=viewport(0.8, 0.8, x=0.4, y=0.4))
print(p3, vp=viewport(0.1, 0.65, x=0.8, y=0.4))

