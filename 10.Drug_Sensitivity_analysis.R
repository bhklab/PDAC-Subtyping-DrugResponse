library(survcomp)
library(gridExtra)

load("/Users/vandanasandhu/Desktop/Subtyping_PDACs/Datasets/pharmacogx_pdac_celllines.RData")
load("/Users/vandanasandhu/Desktop/Subtyping_PDACs/Datasets/celllines_subtyping_RF.RData")
load("/Users/vandanasandhu/Desktop/Subtyping_PDACs/Datasets/sensitivity_celllines_auc_recomputed.RData")
source("/Users/vandanasandhu/Desktop/Subtyping_PDACs/Scripts/mgsub_function.R")



length(which(rownames(sensitivity_cellines$GCSI) %in% rownames(sensitivity_cellines$CCLE)))
drugs= intersect(intersect(rownames(sensitivity_cellines$GCSI),rownames(sensitivity_cellines$CCLE)),rownames(sensitivity_cellines$CTRPv2))

drugs1=drugs
########### GCSI drugs, subtypes, 17009 genes

rownames(celline_subtypes$gcsi_subtypes) = colnames(panc_celllines1$GCSI)
A1=mgsub("[-]",".", colnames(sensitivity_cellines$GCSI))
A2=mgsub(" ",".",A1)
colnames(sensitivity_cellines$GCSI)=A2

length(which(rownames(celline_subtypes$gcsi_subtypes) %in% colnames(sensitivity_cellines$GCSI)))
celllines= rownames(celline_subtypes$gcsi_subtypes)[which(rownames(celline_subtypes$gcsi_subtypes)  %in% colnames(sensitivity_cellines$GCSI))]

gcsi_subtypes = celline_subtypes$gcsi_subtypes[celllines,]
gcsi_ic50 = sensitivity_cellines$GCSI[drugs, celllines]
a=mgsub(c("1","3"), c("Basal","Classical"),  gcsi_subtypes[,2])
gcsi_subtypes$meta_class= a


############ CCLE  drugs, subtypes, 17009 genes

rownames(celline_subtypes$ccle_subtypes) = colnames(panc_celllines1$CCLE)
length(which(rownames(celline_subtypes$ccle_subtypes) %in% colnames(sensitivity_cellines$CCLE)))


celllines= rownames(celline_subtypes$ccle_subtypes)[which(rownames(celline_subtypes$ccle_subtypes) %in% colnames(sensitivity_cellines$CCLE))]


ccle_subtypes = celline_subtypes$ccle_subtypes[celllines,]
ccle_ic50 = sensitivity_cellines$CCLE[drugs, celllines]
#ctrpv2_ic50= -1 * log10(ctrpv2_ic50)


############ CTRPV2  drugs, subtypes, 17009 genes

rownames(celline_subtypes$ctrpv2_subtypes) = colnames(panc_celllines1$CTRPV2)
A1=mgsub("[-]",".", colnames(sensitivity_cellines$CTRPv2))
A2=mgsub(" ",".",A1)
colnames(sensitivity_cellines$CTRPv2)=A2


length(which(rownames(celline_subtypes$ctrpv2_subtypes) %in% colnames(sensitivity_cellines$CTRPv2)))
celllines= rownames(celline_subtypes$ctrpv2_subtypes)[which(rownames(celline_subtypes$ctrpv2_subtypes) %in% colnames(sensitivity_cellines$CTRPv2))]


ctrpv2_subtypes = celline_subtypes$ctrpv2_subtypes[celllines,]
ctrpv2_ic50 = sensitivity_cellines$CTRPv2[drugs, celllines]
#ctrpv2_ic50= -1 * log10(ctrpv2_ic50)

p1=sapply(1:length(drugs), function(x)wilcox.test(ccle_ic50[x,]~ ccle_subtypes$subtypes)$p.value)
p2=sapply(1:length(drugs), function(x)wilcox.test(ctrpv2_ic50[x,]~ ctrpv2_subtypes$subtypes)$p.value)
p3=sapply(1:length(drugs), function(x)wilcox.test(gcsi_ic50[x,]~ gcsi_subtypes$subtypes)$p.value)


#meta_p=sapply(1:length(p1), function(x) combine.test(p=c(p1[x], p2[x],p3[x]), weight=c(length(colnames(ccle_ic50)),length( colnames(ctrpv2_ic50)), length( colnames(gcsi_ic50))),method="z.transform"))
#data.frame(meta_p, drugs)




name=function(x){
  ss=ifelse(x==1,"Basal","Classical")
return(ss)  
}

abc1 = sapply(1:length(drugs), function(x) concordance.index(as.numeric(ccle_ic50[x,]), cl= ifelse(ccle_subtypes$subtypes == 1, 1, 0),na.rm = TRUE, method="noether"))
abc2 = sapply(1:length(drugs), function(x) concordance.index(as.numeric(ctrpv2_ic50[x,]), cl= ifelse(ctrpv2_subtypes$subtypes == 1, 1, 0),na.rm = TRUE, method="noether"))
abc3 = sapply(1:length(drugs), function(x) concordance.index(as.numeric(gcsi_ic50[x,]), cl= ifelse(gcsi_subtypes$subtypes == 1, 1, 0),na.rm = TRUE, method="noether"))

meta_ci=vector()
meta_se=vector()
meta_pvalue=vector()

for(i in 1:length(rownames(ccle_ic50))){
mat= data.frame(cellines= colnames(ccle_ic50), drug=ccle_ic50[i,], subtype= name(ccle_subtypes$subtypes) )
mat1= mat[ !is.na(mat$drug),]
p1=ggboxplot(mat1 , x = "subtype", y = "drug", fill = "subtype", add = "jitter",ylab = "AAC", xlab = "Subtype", title= paste(rownames(ccle_ic50)[i]," " ,"CI: ", round(abc1[,i][[1]],2), " (P: ",  round(abc1[,i][[5]],2),")", sep="") ,  legend = "none" )+
 stat_compare_means(method = "wilcox.test")+  scale_fill_brewer(palette="Set2")

mat= data.frame(cellines= colnames(ctrpv2_ic50), drug=ctrpv2_ic50[i,], subtype= name(ctrpv2_subtypes$subtypes ))
mat1= mat[ !is.na(mat$drug),]
p2=ggboxplot(mat1 , x = "subtype", y = "drug", fill = "subtype", add = "jitter",ylab = "AAC", xlab = "Subtype",title= paste("CI: ", round(abc2[,i][[1]],2), " (P: ",  round(abc2[,i][[5]],2),")", sep="") ,  legend = "none" )+
   stat_compare_means(method = "wilcox.test")+  scale_fill_brewer(palette="Set2")

mat= data.frame(cellines= colnames(gcsi_ic50), drug=gcsi_ic50[i,], subtype= name(gcsi_subtypes$subtypes) )
mat1= mat[ !is.na(mat$drug),]
p3=ggboxplot(mat1 , x = "subtype", y = "drug", fill = "subtype", add = "jitter",ylab = "AAC", xlab = "Subtype", title= paste("CI: ", round(abc3[,i][[1]],2), " (P: ",  round(abc3[,i][[5]],2),")", sep="") ,  legend = "none" )+
stat_compare_means(method = "wilcox.test")+
  scale_fill_brewer(palette="Set2")
grid.arrange(p1,p2,p3, nrow=1, ncol=3)

meta_ci[i]=combine.est(x= c( abc1[,i][[1]],  abc2[,i][[1]], abc3[,i][[1]]), x.se=  c( abc1[,i][[2]],  abc2[,i][[2]], abc3[,i][[2]]), hetero = TRUE, na.rm = TRUE)$estimate
meta_se[i]=combine.est(x= c( abc1[,i][[1]],  abc2[,i][[1]], abc3[,i][[1]]), x.se=  c( abc1[,i][[2]],  abc2[,i][[2]], abc3[,i][[2]]), hetero = TRUE, na.rm = TRUE)$se

#meta_pvalue[i]= combine.test(c( abc1[,i][[5]],  abc2[,i][[5]], abc3[,i][[5]]), c(length(colnames(ccle_ic50)),length( colnames(ctrpv2_ic50)), length( colnames(gcsi_ic50))),method="z.transform")
meta_pvalue[i] <- pnorm((meta_ci[i] - 0.5)/meta_se[i], lower.tail = meta_ci[i] < 0.5) * 2

}

data.frame(drugs, meta_ci, meta_pvalue)



##############################################################################################################################################################################################

#### CCLE vs. GCSI #### #### CCLE vs. GCSI #### #### CCLE vs. GCSI #### #### CCLE vs. GCSI #### #### CCLE vs. GCSI #### #### CCLE vs. GCSI #### 
##################################################################################################################################################
##################################################################################################################################################
######################################################################################################################################################


length(which(rownames(sensitivity_cellines$GCSI) %in% rownames(sensitivity_cellines$CCLE)))
drugs=intersect(rownames(sensitivity_cellines$GCSI),rownames(sensitivity_cellines$CCLE))
drugs2=drugs

########### GCSI drugs, subtypes, 17009 genes

rownames(celline_subtypes$gcsi_subtypes) = colnames(panc_celllines1$GCSI)
A1=mgsub("[-]",".", colnames(sensitivity_cellines$GCSI))
A2=mgsub(" ",".",A1)
colnames(sensitivity_cellines$GCSI)=A2

length(which(rownames(celline_subtypes$gcsi_subtypes) %in% colnames(sensitivity_cellines$GCSI)))
celllines= rownames(celline_subtypes$gcsi_subtypes)[which(rownames(celline_subtypes$gcsi_subtypes)  %in% colnames(sensitivity_cellines$GCSI))]

gcsi_subtypes = celline_subtypes$gcsi_subtypes[celllines,]
gcsi_ic50 = sensitivity_cellines$GCSI[drugs, celllines]
a=mgsub(c("1","3"), c("Basal","Classical"),  gcsi_subtypes[,2])
gcsi_subtypes$meta_class= a


############ CCLE  drugs, subtypes, 17009 genes

rownames(celline_subtypes$ccle_subtypes) = colnames(panc_celllines1$CCLE)
length(which(rownames(celline_subtypes$ccle_subtypes) %in% colnames(sensitivity_cellines$CCLE)))


celllines= rownames(celline_subtypes$ccle_subtypes)[which(rownames(celline_subtypes$ccle_subtypes) %in% colnames(sensitivity_cellines$CCLE))]


ccle_subtypes = celline_subtypes$ccle_subtypes[celllines,]
ccle_ic50 = sensitivity_cellines$CCLE[drugs, celllines]
#ctrpv2_ic50= -1 * log10(ctrpv2_ic50)


p1=sapply(1:length(drugs), function(x)wilcox.test(ccle_ic50[x,]~ ccle_subtypes$subtypes)$p.value)
p3=sapply(1:length(drugs), function(x)wilcox.test(gcsi_ic50[x,]~ gcsi_subtypes$subtypes)$p.value)


abc1 = sapply(1:length(drugs), function(x) concordance.index(as.numeric(ccle_ic50[x,]), cl= ifelse(ccle_subtypes$subtypes == 1, 1, 0),na.rm = TRUE, method="noether"))
abc3 = sapply(1:length(drugs), function(x) concordance.index(as.numeric(gcsi_ic50[x,]), cl= ifelse(gcsi_subtypes$subtypes == 1, 1, 0),na.rm = TRUE, method="noether"))

meta_ci=vector()
meta_pvalue=vector()
meta_se=vector()

for(i in 1:length(rownames(ccle_ic50))){
  mat= data.frame(cellines= colnames(ccle_ic50), drug=ccle_ic50[i,], subtype= name(ccle_subtypes$subtypes) )
  mat1= mat[ !is.na(mat$drug),]
  p1=ggboxplot(mat1 , x = "subtype", y = "drug", fill = "subtype", add = "jitter",ylab = "AAC", xlab = "Subtype", title= paste(rownames(ccle_ic50)[i]," " ,"CI: ", round(abc1[,i][[1]],2), " (P: ",  round(abc1[,i][[5]],2),")", sep="") ,  legend = "none" )+
    stat_compare_means(method = "wilcox.test")+  scale_fill_brewer(palette="Set2")
  
  
  mat= data.frame(cellines= colnames(gcsi_ic50), drug=gcsi_ic50[i,], subtype= name(gcsi_subtypes$subtypes) )
  mat1= mat[ !is.na(mat$drug),]
  p3=ggboxplot(mat1 , x = "subtype", y = "drug", fill = "subtype", add = "jitter",ylab = "AAC", xlab = "Subtype", title= paste("CI: ", round(abc3[,i][[1]],2), " (P: ",  round(abc3[,i][[5]],2),")", sep="") ,  legend = "none" )+
    stat_compare_means(method = "wilcox.test")+
    scale_fill_brewer(palette="Set2")
  grid.arrange(p1,p3, nrow=1, ncol=2)
  
  meta_ci[i]=combine.est(x= c( abc1[,i][[1]],  abc3[,i][[1]]), x.se=  c( abc1[,i][[2]],   abc3[,i][[2]]), hetero = TRUE, na.rm = TRUE)$estimate
  #meta_pvalue[i]= combine.test(c( abc1[,i][[5]],  abc3[,i][[5]]), c(length(colnames(ccle_ic50)), length( colnames(gcsi_ic50))),method="z.transform")
  meta_se[i]=combine.est(x= c( abc1[,i][[1]],  abc3[,i][[1]]), x.se=  c( abc1[,i][[2]],   abc3[,i][[2]]), hetero = TRUE, na.rm = TRUE)$se
   #meta_pvalue[i]= combine.test(c( abc1[,i][[5]],  abc2[,i][[5]], abc3[,i][[5]]), c(length(colnames(ccle_ic50)),length( colnames(ctrpv2_ic50)), length( colnames(gcsi_ic50))),method="z.transform")
  meta_pvalue[i] <- pnorm((meta_ci[i] - 0.5)/meta_se[i], lower.tail = meta_ci[i] < 0.5) * 2
  
}

data.frame(drugs, meta_ci, meta_pvalue)

##############################################################################################################################################################################################

#### CCLE vs. ctrpv2 #### #### CCLE vs.ctrpv2 #### #### CCLE vs. ctrpv2 #### #### CCLE vs. ctrpv2 #### #### CCLE vs. ctrpv2#### #### CCLE vs.ctrpv2 #### 
##################################################################################################################################################
##################################################################################################################################################
######################################################################################################################################################



length(which(rownames(sensitivity_cellines$CTRPv2) %in% rownames(sensitivity_cellines$CCLE)))
drugs= intersect(rownames(sensitivity_cellines$CTRPv2),rownames(sensitivity_cellines$CCLE))
drugs3=drugs

############ CCLE  drugs, subtypes, 17009 genes

rownames(celline_subtypes$ccle_subtypes) = colnames(panc_celllines1$CCLE)
length(which(rownames(celline_subtypes$ccle_subtypes) %in% colnames(sensitivity_cellines$CCLE)))


celllines= rownames(celline_subtypes$ccle_subtypes)[which(rownames(celline_subtypes$ccle_subtypes) %in% colnames(sensitivity_cellines$CCLE))]


ccle_subtypes = celline_subtypes$ccle_subtypes[celllines,]
ccle_ic50 = sensitivity_cellines$CCLE[drugs, celllines]
#ctrpv2_ic50= -1 * log10(ctrpv2_ic50)


############ CTRPV2  drugs, subtypes, 17009 genes

rownames(celline_subtypes$ctrpv2_subtypes) = colnames(panc_celllines1$CTRPV2)
A1=mgsub("[-]",".", colnames(sensitivity_cellines$CTRPv2))
A2=mgsub(" ",".",A1)
colnames(sensitivity_cellines$CTRPv2)=A2


length(which(rownames(celline_subtypes$ctrpv2_subtypes) %in% colnames(sensitivity_cellines$CTRPv2)))
celllines= rownames(celline_subtypes$ctrpv2_subtypes)[which(rownames(celline_subtypes$ctrpv2_subtypes) %in% colnames(sensitivity_cellines$CTRPv2))]


ctrpv2_subtypes = celline_subtypes$ctrpv2_subtypes[celllines,]
ctrpv2_ic50 = sensitivity_cellines$CTRPv2[drugs, celllines]
#ctrpv2_ic50= -1 * log10(ctrpv2_ic50)

p1=sapply(1:length(drugs), function(x)wilcox.test(ccle_ic50[x,]~ ccle_subtypes$subtypes)$p.value)
p2=sapply(1:length(drugs), function(x)wilcox.test(ctrpv2_ic50[x,]~ ctrpv2_subtypes$subtypes)$p.value)


meta_p=sapply(1:length(p1), function(x) combine.test(p=c(p1[x], p2[x]), weight=c(length(colnames(ccle_ic50)),length( colnames(ctrpv2_ic50))),method="z.transform"))
data.frame(meta_p, drugs)


abc1 = sapply(1:length(drugs), function(x) concordance.index(as.numeric(ccle_ic50[x,]), cl= ifelse(ccle_subtypes$subtypes == 1, 1, 0),na.rm = TRUE, method="noether"))
abc2 = sapply(1:length(drugs), function(x) concordance.index(as.numeric(ctrpv2_ic50[x,]), cl= ifelse(ctrpv2_subtypes$subtypes == 1, 1, 0),na.rm = TRUE, method="noether"))

meta_ci=vector()
meta_pvalue=vector()
meta_se=vector()

for(i in 1:length(rownames(ccle_ic50))){
  mat= data.frame(cellines= colnames(ccle_ic50), drug=ccle_ic50[i,], subtype= name(ccle_subtypes$subtypes) )
  mat1= mat[ !is.na(mat$drug),]
  p1=ggboxplot(mat1 , x = "subtype", y = "drug", fill = "subtype", add = "jitter",ylab = "AAC", xlab = "Subtype", title= paste(rownames(ccle_ic50)[i]," " ,"CI: ", round(abc1[,i][[1]],2), " (P: ",  round(abc1[,i][[5]],2),")", sep="") ,  legend = "none" )+
    stat_compare_means(method = "wilcox.test")+  scale_fill_brewer(palette="Set2")
  
  mat= data.frame(cellines= colnames(ctrpv2_ic50), drug=ctrpv2_ic50[i,], subtype= name(ctrpv2_subtypes$subtypes ))
  mat1= mat[ !is.na(mat$drug),]
  p2=ggboxplot(mat1 , x = "subtype", y = "drug", fill = "subtype", add = "jitter",ylab = "AAC", xlab = "Subtype",title= paste("CI: ", round(abc2[,i][[1]],2), " (P: ",  round(abc2[,i][[5]],2),")", sep="") ,  legend = "none" )+
    stat_compare_means(method = "wilcox.test")+  scale_fill_brewer(palette="Set2")
  
  grid.arrange(p1,p2, nrow=1, ncol=2)
  
  meta_ci[i]=combine.est(x= c( abc1[,i][[1]],  abc2[,i][[1]]), x.se=  c( abc1[,i][[2]],  abc2[,i][[2]]), hetero = TRUE, na.rm = TRUE)$estimate
  #meta_pvalue[i]= combine.test(c( abc1[,i][[5]],  abc2[,i][[5]]), c(length(colnames(ccle_ic50)),length( colnames(ctrpv2_ic50))),method="z.transform")
  meta_se[i]=combine.est(x= c( abc1[,i][[1]],  abc2[,i][[1]]), x.se=  c( abc1[,i][[2]],  abc2[,i][[2]]), hetero = TRUE, na.rm = TRUE)$se
  #meta_pvalue[i]= combine.test(c( abc1[,i][[5]],  abc2[,i][[5]], abc3[,i][[5]]), c(length(colnames(ccle_ic50)),length( colnames(ctrpv2_ic50)), length( colnames(gcsi_ic50))),method="z.transform")
  meta_pvalue[i] <- pnorm((meta_ci[i] - 0.5)/meta_se[i], lower.tail = meta_ci[i] < 0.5) * 2
  
}

data.frame(drugs, meta_ci, meta_pvalue)




##############################################################################################################################################################################################

#### GCSI vs. ctrpv2 #### #### GCSI vs.ctrpv2 #### #### GCSI vs. ctrpv2 #### #### GCSI vs. ctrpv2 #### #### GCSI vs. ctrpv2#### #### GCSI vs.ctrpv2 #### 
##################################################################################################################################################
##################################################################################################################################################
######################################################################################################################################################




length(which(rownames(sensitivity_cellines$GCSI) %in% rownames(sensitivity_cellines$CTRPv2)))
drugs= intersect(rownames(sensitivity_cellines$GCSI),rownames(sensitivity_cellines$CTRPv2))

drugs4=drugs
########### GCSI drugs, subtypes, 17009 genes

rownames(celline_subtypes$gcsi_subtypes) = colnames(panc_celllines1$GCSI)
A1=mgsub("[-]",".", colnames(sensitivity_cellines$GCSI))
A2=mgsub(" ",".",A1)
colnames(sensitivity_cellines$GCSI)=A2

length(which(rownames(celline_subtypes$gcsi_subtypes) %in% colnames(sensitivity_cellines$GCSI)))
celllines= rownames(celline_subtypes$gcsi_subtypes)[which(rownames(celline_subtypes$gcsi_subtypes)  %in% colnames(sensitivity_cellines$GCSI))]

gcsi_subtypes = celline_subtypes$gcsi_subtypes[celllines,]
gcsi_ic50 = sensitivity_cellines$GCSI[drugs, celllines]
a=mgsub(c("1","3"), c("Basal","Classical"),  gcsi_subtypes[,2])
gcsi_subtypes$meta_class= a


############ CTRPV2  drugs, subtypes, 17009 genes


rownames(celline_subtypes$ctrpv2_subtypes) = colnames(panc_celllines1$CTRPV2)
A1=mgsub("[-]",".", colnames(sensitivity_cellines$CTRPv2))
A2=mgsub(" ",".",A1)
colnames(sensitivity_cellines$CTRPv2)=A2


length(which(rownames(celline_subtypes$ctrpv2_subtypes) %in% colnames(sensitivity_cellines$CTRPv2)))
celllines= rownames(celline_subtypes$ctrpv2_subtypes)[which(rownames(celline_subtypes$ctrpv2_subtypes) %in% colnames(sensitivity_cellines$CTRPv2))]


ctrpv2_subtypes = celline_subtypes$ctrpv2_subtypes[celllines,]
ctrpv2_ic50 = sensitivity_cellines$CTRPv2[drugs, celllines]
#ctrpv2_ic50= -1 * log10(ctrpv2_ic50)

p2=sapply(1:length(drugs), function(x)wilcox.test(ctrpv2_ic50[x,]~  ctrpv2_subtypes$subtypes)$p.value)
p3=sapply(1:length(drugs), function(x)wilcox.test(gcsi_ic50[x,]~  gcsi_subtypes$subtypes)$p.value)



abc2 = sapply(1:length(drugs), function(x) concordance.index(as.numeric(ctrpv2_ic50[x,]), cl= ifelse(ctrpv2_subtypes$subtypes == 1, 1, 0),na.rm = TRUE, method="noether"))
abc3 = sapply(1:length(drugs), function(x) concordance.index(as.numeric(gcsi_ic50[x,]), cl= ifelse(gcsi_subtypes$subtypes == 1, 1, 0),na.rm = TRUE, method="noether"))

meta_ci=vector()
meta_pvalue=vector()
meta_se=vector()
for(i in 1:length(rownames(gcsi_ic50))){

  mat= data.frame(cellines= colnames(ctrpv2_ic50), drug=ctrpv2_ic50[i,], subtype= name(ctrpv2_subtypes$subtypes ))
  mat1= mat[ !is.na(mat$drug),]
  p2=ggboxplot(mat1 , x = "subtype", y = "drug", fill = "subtype", add = "jitter",ylab = "AAC", xlab = "Subtype",title= paste(rownames(gcsi_ic50)[i]," " ,"CI: ", round(abc2[,i][[1]],2), " (P: ",  round(abc2[,i][[5]],2),")", sep="") ,  legend = "none" )+
    stat_compare_means(method = "wilcox.test")+  scale_fill_brewer(palette="Set2")
  
  mat= data.frame(cellines= colnames(gcsi_ic50), drug=gcsi_ic50[i,], subtype= name(gcsi_subtypes$subtypes) )
  mat1= mat[ !is.na(mat$drug),]
  p3=ggboxplot(mat1 , x = "subtype", y = "drug", fill = "subtype", add = "jitter",ylab = "AAC", xlab = "Subtype", title= paste("CI: ", round(abc3[,i][[1]],2), " (P: ",  round(abc3[,i][[5]],2),")", sep="") ,  legend = "none" )+
    stat_compare_means(method = "wilcox.test")+
    scale_fill_brewer(palette="Set2")
  grid.arrange(p2,p3, nrow=1, ncol=2)
  
  meta_ci[i]=combine.est(x= c(  abc2[,i][[1]], abc3[,i][[1]]), x.se=  c(   abc2[,i][[2]], abc3[,i][[2]]), hetero = TRUE, na.rm = TRUE)$estimate
  #meta_pvalue[i]= combine.test(c(  abc2[,i][[5]], abc3[,i][[5]]), c(length( colnames(ctrpv2_ic50)), length( colnames(gcsi_ic50))),method="z.transform")
  
  meta_se[i]=combine.est(x= c( abc2[,i][[1]],  abc3[,i][[1]]), x.se=  c( abc2[,i][[2]],  abc3[,i][[2]]), hetero = TRUE, na.rm = TRUE)$se
  #meta_pvalue[i]= combine.test(c( abc1[,i][[5]],  abc2[,i][[5]], abc3[,i][[5]]), c(length(colnames(ccle_ic50)),length( colnames(ctrpv2_ic50)), length( colnames(gcsi_ic50))),method="z.transform")
  meta_pvalue[i] <- pnorm((meta_ci[i] - 0.5)/meta_se[i], lower.tail = meta_ci[i] < 0.5) * 2
  
}

data.frame(drugs, meta_ci, meta_pvalue)






################### ALL CTRPV2 DRUGS
dd=unique(c(drugs1, drugs2, drugs3, drugs4))

dd1=setdiff(rownames(sensitivity_cellines$CTRPv2), dd)

rownames(celline_subtypes$ctrpv2_subtypes) = colnames(panc_celllines1$CTRPV2)
A1=mgsub("[-]",".", rownames(celline_subtypes$ctrpv2_subtypes))
A2=mgsub(" ",".",A1)
rownames(celline_subtypes$ctrpv2_subtypes)=A2


length(which(rownames(celline_subtypes$ctrpv2_subtypes) %in% colnames(sensitivity_cellines$CTRPv2)))
celllines= rownames(celline_subtypes$ctrpv2_subtypes)[which(rownames(celline_subtypes$ctrpv2_subtypes) %in% colnames(sensitivity_cellines$CTRPv2))]

ctrpv2_subtypes = celline_subtypes$ctrpv2_subtypes[celllines,]
ctrpv2_ic50 = sensitivity_cellines$CTRPv2[dd1, celllines]
ctrpv2_drugs_pvalue= vector()

for( i in 1: dim(ctrpv2_ic50 )[1]){
  
  if(length(which(ctrpv2_ic50[i,][ctrpv2_subtypes$subtypes==3] %in% NA)) <= 5){
  #ctrpv2_drugs_pvalue[i] = wilcox.test(ctrpv2_ic50[i,]~ ctrpv2_subtypes$subtypes)$p.value
  ctrpv2_drugs_pvalue[i]= concordance.index(as.numeric(ctrpv2_ic50[ i,]), cl= ifelse(ctrpv2_subtypes$subtypes == 1, 1, 0),na.rm = TRUE, method="noether")$p.value
  

  }
  else{
    ctrpv2_drugs_pvalue[i] = NA
  }
}

adjusted_p= p.adjust(ctrpv2_drugs_pvalue, method = 'fdr')
DD=data.frame(rownames(ctrpv2_ic50),ctrpv2_drugs_pvalue, adjusted_p)
DD1= DD[which(DD$ctrpv2_drugs_pvalue <0.05),]
Z=which(rownames(ctrpv2_ic50) %in% DD1$rownames.ctrpv2_ic50.)
for(y in 1: length( Z)){
xxx= concordance.index(as.numeric(ctrpv2_ic50[ Z[y],]), cl= ifelse(ctrpv2_subtypes$subtypes == 1, 1, 0),na.rm = TRUE, method="noether")
mat= data.frame(cellines= colnames(ctrpv2_ic50), drug=ctrpv2_ic50[Z[y],], subtype= name(ctrpv2_subtypes$subtypes) )
mat1= mat[ !is.na(mat$drug),]
p=ggboxplot(mat1 , x = "subtype", y = "drug", fill = "subtype", add = "jitter",ylab = "AAC", xlab = "Subtype", title= paste(DD1$rownames.ctrpv2_ic50.[y]," CI: ", round(xxx$c.index,2), " (P: ",  round(xxx$p.value,2),")", sep="") ,  legend = "none" )+
  stat_compare_means(method = "wilcox.test")+
  scale_fill_brewer(palette="Set2")
plot(p)
}

### GCSI
dd=unique(c(drugs1, drugs2, drugs3, drugs4))

dd1=setdiff(rownames(sensitivity_cellines$GCSI), dd)


rownames(celline_subtypes$gcsi_subtypes) = colnames(panc_celllines1$GCSI)
A1=mgsub("[-]",".", colnames(sensitivity_cellines$GCSI))
A2=mgsub(" ",".",A1)
colnames(sensitivity_cellines$GCSI)=A2

length(which(rownames(celline_subtypes$gcsi_subtypes) %in% colnames(sensitivity_cellines$GCSI)))
celllines= rownames(celline_subtypes$gcsi_subtypes)[which(rownames(celline_subtypes$gcsi_subtypes)  %in% colnames(sensitivity_cellines$GCSI))]

gcsi_subtypes = celline_subtypes$gcsi_subtypes[celllines,]
gcsi_ic50 = sensitivity_cellines$GCSI[dd1, celllines]

gcsi_drugs_pvalue= vector()

for( i in 1: dim(gcsi_ic50 )[1]){
  
  if(length(which(gcsi_ic50[i,][gcsi_subtypes$subtypes==3] %in% NA)) <= 6){
    gcsi_drugs_pvalue[i] = wilcox.test(gcsi_ic50[i,]~ gcsi_subtypes$subtypes)$p.value
  }
  else{
    gcsi_drugs_pvalue[i] = NA
  }
}

adjusted_p= p.adjust(gcsi_drugs_pvalue, method = 'fdr')
DD=data.frame(rownames(gcsi_ic50),gcsi_drugs_pvalue, adjusted_p)

DD[which(DD$gcsi_drugs_pvalue <0.05),]

xxx= concordance.index(as.numeric(gcsi_ic50["Doxorubicin",]), cl= ifelse(gcsi_subtypes$subtypes == 1, 1, 0),na.rm = TRUE, method="noether")
mat= data.frame(cellines= colnames(gcsi_ic50), drug=gcsi_ic50["Doxorubicin",], subtype= name(gcsi_subtypes$subtypes) )
mat1= mat[ !is.na(mat$drug),]
ggboxplot(mat1 , x = "subtype", y = "drug", fill = "subtype", add = "jitter",ylab = "AAC", xlab = "Subtype", title= paste("Doxorubicin CI: ", round(xxx$c.index,2), " (P: ",  round(xxx$p.value,2),")", sep="") ,  legend = "none" )+
  stat_compare_means(method = "wilcox.test")+
  scale_fill_brewer(palette="Set2")


### CCLE
dd=unique(c(drugs1, drugs2, drugs3, drugs4))

dd1=setdiff(rownames(sensitivity_cellines$CCLE), dd)


#rownames(celline_subtypes$ccle_subtypes) = colnames(panc_celllines1$CCLE)
#A1=mgsub("[-]",".", rownames(celline_subtypes$ccle_subtypes) )
#A2=mgsub(" ",".",A1)
#rownames(celline_subtypes$ccle_subtypes) =A2
rownames(celline_subtypes$ccle_subtypes) = colnames(panc_celllines1$CCLE)
length(which(rownames(celline_subtypes$ccle_subtypes) %in% colnames(sensitivity_cellines$CCLE)))
celllines= rownames(celline_subtypes$ccle_subtypes)[which(rownames(celline_subtypes$ccle_subtypes)  %in% colnames(sensitivity_cellines$CCLE))]

ccle_subtypes = celline_subtypes$ccle_subtypes[celllines,]
ccle_ic50 = sensitivity_cellines$CCLE[dd1, celllines]

ccle_drugs_pvalue= vector()

for( i in 1: dim(ccle_ic50 )[1]){
  
  if(length(which(ccle_ic50[i,][ccle_subtypes$subtypes==3] %in% NA)) <= 6){
    ccle_drugs_pvalue[i] = wilcox.test(ccle_ic50[i,]~ ccle_subtypes$subtypes)$p.value
  }
  else{
    ccle_drugs_pvalue[i] = NA
  }
}

adjusted_p= p.adjust(ccle_drugs_pvalue, method = 'fdr')
DD=data.frame(rownames(ccle_ic50),ccle_drugs_pvalue, adjusted_p)

DD[which(DD$ccle_drugs_pvalue <0.05),]




xxx= concordance.index(as.numeric(ccle_ic50["AEW541",]), cl= ifelse(ccle_subtypes$subtypes == 1, 1, 0),na.rm = TRUE, method="noether")
mat= data.frame(cellines= colnames(ccle_ic50), drug=ccle_ic50["AEW541",], subtype= name(ccle_subtypes$subtypes) )
mat1= mat[ !is.na(mat$drug),]
ggboxplot(mat1 , x = "subtype", y = "drug", fill = "subtype", add = "jitter",ylab = "AAC", xlab = "Subtype", title= paste("AEW541 CI: ", round(xxx$c.index,2), " (P: ",  round(xxx$p.value,2),")", sep="") ,  legend = "none" )+
  stat_compare_means(method = "wilcox.test")+
  scale_fill_brewer(palette="Set2")



##################################################################################################

unique_drugs= unique(c(rownames(ccle_ic50),rownames(ctrpv2_ic50), rownames(gcsi_ic50)))
c3 <- cbind(ifelse(unique_drugs %in% rownames(ccle_ic50), TRUE, FALSE), ifelse(unique_drugs %in% rownames(ctrpv2_ic50), TRUE, FALSE),ifelse(unique_drugs %in%  rownames(gcsi_ic50), TRUE, FALSE))
a <- vennCounts(c3)

vennDiagram(a, names = c("CCLE","CTRPv2","GCSI"), counts.col = "black", circle.col=brewer.pal(3,"Dark2"))



unique_drugs= unique(c(rownames(ccle_ic50),rownames(ctrpv2_ic50), rownames(gcsi_ic50)))
c3 <- cbind(ifelse(unique_drugs %in% rownames(ccle_ic50), TRUE, FALSE), ifelse(unique_drugs %in% rownames(ctrpv2_ic50), TRUE, FALSE),ifelse(unique_drugs %in%  rownames(gcsi_ic50), TRUE, FALSE))
a <- vennCounts(c3)

vennDiagram(a, names = c("CCLE","CTRPv2","GCSI"), counts.col = "black", circle.col=brewer.pal(3,"Dark2"))


##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################
######################  PDXE drug comparison #####################################################

################################################################################################################################
#load("/Users/vandanasandhu/Desktop/Subtyping_PDACs/Datasets/PDXe_classes.RData")
load("/Users/vandanasandhu/Desktop/Subtyping_PDACs/Datasets/PDXe_drug_response.RData")

pdxe_drugs = unique(meta$drug)
pdxe_subtypes$subtypes=mgsub(c("1","3"), c("Basal","Classical"),  pdxe_subtypes$subtypes)


for(i in 1: length(pdxe_drugs)){
  pdxe_mm = meta[which(meta$drug == pdxe_drugs[i]),]
  pdxe_mm$patient.id=mgsub("-",".",pdxe_mm$patient.id)
  
  mm = merge(pdxe_mm, pdxe_subtypes, by.x="patient.id", by.y="id")
  
  p=ggboxplot(mm , x = "subtypes", y = "AAC", fill = "subtypes", add = "jitter",ylab = "AAC", xlab = "Subtype", title= mm$drug[1])+ stat_compare_means(method = "kruskal.test")+theme_bw(base_family = 'Helvetica')
  plot(p)
}


pdxe_mm = meta[which(meta$drug == "BKM120 + binimetinib"),]
mm = merge(pdxe_mm, pdxe_subtypes, by.x="patient.id", by.y="id")


ggboxplot(mm ,  x = "subtypes", y = "AAC", fill = "subtypes", add = "jitter",ylab = "AAC", xlab = "Subtype", title= mm$drug[1])+ stat_compare_means(method = "kruskal.test")
mm$color =mgsub(c("Basal","Classical"), c("gold", "blue"),  mm$subtypes )
mm= mm[order(mm$AAC),]
barplot(mm$AAC, col=mm$color,  ylab= "AAC", border = "white")
legend("topleft", c("Basal", "Classical","P=0.01"), fill   = c("gold", "blue","white"), bty = "n", border = FALSE)
########



length(which(rownames(sensitivity_cellines$CTRPv2) %in% pdxe_drugs))




pdxe_mm = meta[which(meta$drug == "abraxane"),]
mm = merge(pdxe_mm, pdxe_subtypes, by.x="patient.id", by.y="id")

mm$AAC=mm$AAC/100
ggboxplot(mm ,  x = "subtypes", y = "AAC", fill = "subtypes", color=c("orange","pink"),add = "jitter",ylab = "AAC", xlab = "Subtype", title= mm$drug[1])+ stat_compare_means(method = "kruskal.test")
mm$color =mgsub(c("Basal","Classical"), c("gold", "blue"),  mm$subtypes )
mm= mm[order(mm$AAC),]
barplot(mm$AAC, col=mm$color,  ylab= "AAC", border = "white")
legend("topleft", c("Basal", "Classical","P=0.01"), fill   = c("gold", "blue","white"), bty = "n", border = FALSE)
########




