library(randomForest)
library(vcdExtra)
library(switchBox)
library(impute)

a=load("/Users/vandanasandhu/Desktop/Subtyping_PDACs/RS_New_Analysis/Datasets/common_genes_cohorts_new.RData")
meta_genes=load("/Users/vandanasandhu/Desktop/Subtyping_PDACs/meta_genes.RData")

pcsi=data.matrix(cohorts1$pcsi[as.character(meta_genes),])
tcga=data.matrix(cohorts1$tcga[as.character(meta_genes),])
icgc=data.matrix(cohorts1$icgc_seq[as.character(meta_genes),])
kirby=data.matrix(cohorts1$kirby[as.character(meta_genes),])
ouh=data.matrix(cohorts1$ouh[as.character(meta_genes),])
winter=data.matrix(cohorts1$winter[as.character(meta_genes),])
collisson=data.matrix(cohorts1$collisson[as.character(meta_genes),])
zhang=data.matrix(cohorts1$zhang[as.character(meta_genes),])
chen=data.matrix(cohorts1$chen[as.character(meta_genes),])
unc=data.matrix(cohorts1$unc[as.character(meta_genes),])
icgc_arr=data.matrix(cohorts1$icgc_arr[as.character(meta_genes),])
balagurunathan=data.matrix(cohorts1$balagurunathan[as.character(meta_genes),])
pei=data.matrix(cohorts1$pei[as.character(meta_genes),])
grutzmann=data.matrix(cohorts1$grutzmann[as.character(meta_genes),])
badea=data.matrix(cohorts1$badea[as.character(meta_genes),])

hamidi=data.matrix(cohorts1$hamidi[as.character(meta_genes),])
haider=data.matrix(cohorts1$haider[as.character(meta_genes),])
bauer=data.matrix(cohorts1$bauer[as.character(meta_genes),])
yang=data.matrix(cohorts1$yang[as.character(meta_genes),])
#van_den_broeck=data.matrix(cohorts1$van_den_broeck[as.character(meta_genes),])
lunardi=data.matrix(cohorts1$lunardi[as.character(meta_genes),])
janky=data.matrix(cohorts1$janky[as.character(meta_genes),])
hamidi=impute.knn(hamidi)[[1]]
lunardi=impute.knn(lunardi)[[1]]


#clusters = list(PCSI= z_PCSI, TCGA=z_TCGA, ICGC_seq=z_ICGC, 
#                Kirby = z_Kirby, OUH= z_OUH, Winter= z_winter,
#               Collisson=z_Collisson,Zhang=z_zhang,Chen=z_chen,
#                UNC=z_unc,ICGC_arr=z_icgc_arr,Balagurunathan=z_balagurunathan, 
#                Grutzmann=z_grutzmann,Badea= z_badea,Pei=z_pei)


##########################################################################
##########################################################################
##########################################################################
##########################################################################
##########################################################################
##########################################################################

clusters1=clusters
classes=list()
sample_names=list()
cohorts=list()
for(i in 1:(length(clusters1))){
  
  classes[[i]]=clusters1[[i]]$meta_classes
  sample_names[[i]]=names(clusters1[[i]]$meta_classes)
  cohorts[[i]]=rep(names(clusters1)[i],length(names(clusters1[[i]]$meta_classes)))
}

meta_Cluster= data.frame(sample = unlist(sample_names), meta_class= unlist(classes), cohorts=unlist(cohorts))

##########################################################################
##########################################################################
##########################################################################
##########################################################################
##########################################################################
##########################################################################

#########


com_matrix= cbind( pcsi, tcga, icgc, kirby, ouh, winter, collisson, zhang, chen, unc, icgc_arr, balagurunathan, grutzmann, badea , pei, hamidi, yang, lunardi, janky, bauer, haider)
length(which(colnames(com_matrix) == meta_Cluster$sample))
#com_matrix= NaRV.omit(com_matrix)
#meta_Cluster$meta_class
#########################################################################


dat <- com_matrix
Groups <- meta_Cluster$meta_class
y <- as.numeric(as.character(Groups))




pheno_grp <- function(y, com){
  
  z=setdiff(y,com)
  class_com = unlist(sapply(z, function(x) which(Groups ==x)))
  na_index= which(y %in% NA)
  y1=as.numeric(y)
  y1[c(class_com, na_index)] =0
  return(as.factor(y1))  
}

n=3
#dat[is.na(dat)] <- 0
dat[is.na(dat)] <- 0

model=lapply(1:n, function(x) SWAP.KTSP.Train(dat, pheno_grp(y, x) , k=50))

gene1 = unlist(sapply(1:n, function(x) model[[x]]$TSPs[,1]))
gene2 = unlist(sapply(1:n, function(x) model[[x]]$TSPs[,2]))
genepairs <- as.character(paste(gene1,gene2,sep=">"))

feat <- c(gene1, gene2)
xx= dat[feat,]


ll=lapply(1: dim(xx)[2], function(x)  ifelse(xx[gene1,x] > xx[gene2,x], 1,0))
output <- matrix(unlist(ll), nrow = ncol(xx), byrow = TRUE)

rownames(output)=paste('Patient',1:ncol(xx),sep='') 
colnames(output)=genepairs

y <- as.factor(y)
y_NA = which(y %in% NA)
y1= y[-y_NA]
output1= output[-y_NA,]
output1=t(output1)


cores=detectCores()
cl <- makeCluster(cores[1]-1, outfile="") #not to overload your computer
registerDoParallel(cl)

pred=list()
pred= foreach(i = 1: dim(output1)[2], .packages = c( "randomForest","foreach")) %dopar% {         
  
  output2 <- t(output1)
  print(i)
  train <- output2[-i,]
  test <- output2[i,]
  ytrain <- y1[-i]
  model <- randomForest(train, ytrain)
  pred[[i]] <- as.character(predict(model, test))
  
}

binary_rf_model <- randomForest( t(output1), y1)
gg=list(gene1=gene1, gene2= gene2, genepairs =genepairs)

save(binary_rf_model, file="/Users/vandanasandhu/Desktop/Subtyping_PDACs/RS_New_Analysis/Datasets/binary_rf_model.RData")
save(output1, file="/Users/vandanasandhu/Desktop/Subtyping_PDACs/RS_New_Analysis/Datasets/output1_RF.RData")
save(y1, file="/Users/vandanasandhu/Desktop/Subtyping_PDACs/RS_New_Analysis/Datasets/RF_response.RData")
save(gg, file="/Users/vandanasandhu/Desktop/Subtyping_PDACs/RS_New_Analysis/Datasets/gg.RData")


pred1 <- do.call("c", pred)
table(pred1, as.character(y1))
a <- data.frame(predictions = pred1,group = y1)
a1 <- xtabs(~a[,1] + a[,2], data = a)
summary(assocstats(a1))
length(which(pred1 == as.character(y1)))/length(pred1)

##########################################################################
##########################################################################
############################ LEAVE ONE COHORT OUT CROSS VALIDATION ##############################################
##########################################################################
##########################################################################

com_matrix= cbind( pcsi, tcga, icgc, kirby, ouh, winter, collisson, zhang, chen, unc, icgc_arr, balagurunathan, grutzmann, badea , pei, hamidi, yang, lunardi, janky, bauer, haider)
#com_matrix= NaRV.omit(com_matrix)
length(which(colnames(com_matrix) == meta_Cluster$sample))

meta_Cluster$meta_class
uu=unique(meta_Cluster$cohorts)
pred1=list()
for(i in 1:length(unique(meta_Cluster$cohorts)) ){
  print(uu[i])
  z1 = which(meta_Cluster$cohorts == uu[i])
  
  dat <- com_matrix[,-z1]
  Groups <- meta_Cluster$meta_class[-z1]
  y <- as.numeric(as.character(Groups))
  n=3
  dat[is.na(dat)] <- 0
  model=lapply(1:n, function(x) SWAP.KTSP.Train(dat, pheno_grp(y, x) ,k=50))
  
  gene1 = unlist(sapply(1:n, function(x) model[[x]]$TSPs[,1]))
  gene2 = unlist(sapply(1:n, function(x) model[[x]]$TSPs[,2]))
  genepairs <- as.character(paste(gene1,gene2,sep=">"))
  
  feat <- c(gene1, gene2)
  xx= dat[feat,]
  
  ll=lapply(1: dim(xx)[2], function(x)  ifelse(xx[gene1,x] > xx[gene2,x], 1,0))
  output <- matrix(unlist(ll), nrow = ncol(xx), byrow = TRUE)
  
  rownames(output)=paste('Patient',1:ncol(xx),sep='') 
  colnames(output)=genepairs
  
  y <- as.factor(y)
  y_NA = which(y %in% NA)
  y1= y[-y_NA]
  output1= output[-y_NA,]
  output1=output1
  
  model <- randomForest(output1, y1)
  test= com_matrix[,z1]
  
  ll1 = lapply(1: dim(test)[2], function(x)  ifelse(test[gene1,x] > test[gene2,x], 1,0))
  output <- matrix(unlist(ll1), nrow = ncol(test), byrow = TRUE)
  rownames(output)=colnames(com_matrix)[z1] 
  colnames(output)=genepairs
  
  pred1[[i]] <- as.character(predict(model, output,na.action =  na.exclude))
}

pred2 <- do.call("c", pred1)
table(pred2, as.character(meta_Cluster$meta_class))
a <- data.frame(predictions = pred2,group = as.character(meta_Cluster$meta_class))
a1 <- xtabs(~a[,1] + a[,2], data = a)
summary(assocstats(a1))
length(which(pred2 == as.character(meta_Cluster$meta_class)))/length(pred2)


##########################################################################

mm1 = meta_Cluster[-y_NA,]
mm1$tsp_classes = pred1

arr = mm1[which(mm1$cohorts == "ICGC_arr"),]
rownames(arr)= arr$sample
seq = mm1[which(mm1$cohorts == "ICGC_seq"),]
rownames(seq)= seq$sample


################ ICGC-array and ICGC-seq

common_seq = which(rownames(seq) %in% rownames(arr))
common_arr = which(rownames(arr) %in% rownames(seq))

sort(rownames(seq)[common_seq]) == sort(rownames(arr)[common_arr])


arr_classes= arr[sort(rownames(seq)[common_seq]),]
seq_classes= seq[sort(rownames(arr)[common_arr]),]

table(arr_classes$tsp_classes, seq_classes$tsp_classes)
length(which(arr_classes$tsp_classes == seq_classes$tsp_classes))

######### Survival 
load("/Users/vandanasandhu/Desktop/RS_Remodelling_TSP_project/Github/RData/PDAC_Expression_dataset.RData")
yang_survival=read.csv("/Users/vandanasandhu/Desktop/Subtyping_PDACs/RS_New_Analysis/Datasets/yang_survival.txt", sep="\t", header=T)

samples=c(rownames(rs_coh$PCSI_new), rownames(rs_coh$TCGA), rownames(rs_coh$Kirby),
          rownames(rs_coh$ICGC_arr),rownames(rs_coh$ICGC_seq), rownames(rs_coh$UNC),
          rownames(rs_coh$Chen), rownames(rs_coh$Collisson), rownames(rs_coh$Zhang),
          rownames(rs_coh$OUH), rownames(rs_coh$Winter), as.character(yang_survival$ID))


os= c(as.numeric(as.character(rs_coh$PCSI_new$OS)), as.numeric(as.character(rs_coh$TCGA$OS)), 
      as.numeric(as.character(rs_coh$Kirby$OS)),as.numeric(as.character(rs_coh$ICGC_arr$OS)),
      as.numeric(as.character(rs_coh$ICGC_seq$OS)), as.numeric(as.character(rs_coh$UNC$OS)),
      as.numeric(as.character(rs_coh$Chen$OS)), as.numeric(as.character(rs_coh$Collisson$OS)), 
      as.numeric(as.character(rs_coh$Zhang$OS)), as.numeric(as.character(rs_coh$OUH$OS)), 
      as.numeric(as.character(rs_coh$Winter$OS)), as.numeric(as.character(yang_survival$OS)))


os_status= c(as.numeric(as.character(rs_coh$PCSI_new$OS_Status)), as.numeric(as.character(rs_coh$TCGA$OS_Status)), 
             as.numeric(as.character(rs_coh$Kirby$OS_Status)),as.numeric(as.character(rs_coh$ICGC_arr$OS_Status)),
             as.numeric(as.character(rs_coh$ICGC_seq$OS_Status)), as.numeric(as.character(rs_coh$UNC$OS_Status)),
             as.numeric(as.character(rs_coh$Chen$OS_Status)), as.numeric(as.character(rs_coh$Collisson$OS_Status)), 
             as.numeric(as.character(rs_coh$Zhang$OS_Status)), as.numeric(as.character(rs_coh$OUH$OS_Status)), 
             as.numeric(as.character(rs_coh$Winter$OS_Status)),as.numeric(as.character(yang_survival$OS_Status)))



survival= data.frame(samples=samples, os=os, os_status=os_status)


##########################################################################
##########################################################################

duplicate_id = which(duplicated(mm1$sample) == TRUE)
mm2=mm1[-duplicate_id,]


merged_mat= merge(survival, mm2, by.x="samples", by.y="sample")


fit <- survfit(Surv(merged_mat$os, merged_mat$os_status == 1) ~ merged_mat$tsp_classes + strata(cohorts), data=merged_mat)
ggsurvplot(fit, data =merged_mat, risk.table = TRUE, pval = TRUE,ggtheme = theme_minimal())


su=coxph(Surv(merged_mat$os, merged_mat$os_status == 1) ~ merged_mat$tsp_classes+ strata(cohorts), data=merged_mat)
summary(su)

##########################################################################
##########################################################################
unique_cohorts=unique(merged_mat$cohorts)
par(mfrow=c(3,4))

for(i in 1:length(unique(merged_mat$cohorts))){
  zz=merged_mat[which(merged_mat$cohorts == unique_cohorts[i] ),]
  
  fit <- survfit(Surv(zz$os, zz$os_status == 1) ~ zz$tsp_classes, data=zz)
  
  sd=survdiff(Surv(zz$os, zz$os_status == 1) ~  zz$tsp_classes, data=zz)
  p.val <- 1 - pchisq(sd$chisq, length(sd$n) - 1)
  plot(fit,col=c("tomato","green","cyan","purple","yellow","black"), lwd=2, main=unique_cohorts[i], p.val=TRUE)
  legend("bottomleft",paste("P =", round(p.val,2), sep = ""), bty = "n") 
}

##########################################################################
##########################################################################

