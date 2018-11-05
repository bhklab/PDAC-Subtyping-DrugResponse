library(piano)
library(GSEA)

## META EFFECT SIZE ########################################################################################

a=load("/Users/vandanasandhu/Desktop/Subtyping_PDACs/RS_New_Analysis/Datasets/common_genes_cohorts_new.RData")

names(clusters)
names(cohorts1)

meta_effect=vector('list', 22)
sd=vector('list', 22)



for(j in 1:length(clusters)){
  for (i in 1:3){
    z=sort(unique(clusters[[j]]$meta_classes))
    k=1:4
    if(i %in% z){
      s= setdiff(unique(clusters[[j]]$meta_classes),i)
      
      response=mapvalues(clusters[[j]]$meta_classes,s, rep(0,length(s)))
      
      meta_effect[[j]][[i]] = sapply(1:dim(cohorts1[[j]])[1], function(x) cohen.d(as.numeric(cohorts1[[j]][x,]),response,hedges.correction=TRUE)$estimate)
      sd[[j]][[i]]= sapply(1:dim(cohorts1[[j]])[1], function(x) sd(as.numeric(cohorts1[[j]][x,which(response!=0)])))
      
    }
    else{
      meta_effect[[j]][[i]] = NA
      sd[[j]][[i]]= NA
    }
  }
  print(names(cohorts1)[j])
}

n=10331  #no. of genes
mean_effect_size=vector('list', 4)

for(j in 1: 3){
  for(i in 1:n){
    
    val=sapply(1:22, function(x) meta_effect[[x]][[j]][i])
    wt=sapply(1:22, function(x) 1/sd[[x]][[j]][i])
    
    
    mean_effect_size[[j]][[i]] = weighted.mean(unlist(val),unlist(wt), na.rm = TRUE)
  } 
}


mm=matrix(unlist(mean_effect_size), ncol=3)
rownames(mm)=rownames(cohorts1$pcsi)
colnames(mm)=c("Metacluster1","Metacluster2","Metacluster3")
mm=data.frame(mm)

save(mm,file="/Users/vandanasandhu/Desktop/Subtyping_PDACs/RS_New_Analysis/Datasets/meta_effectsize.RData")  


#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
load("/Users/vandanasandhu/Desktop/Subtyping_PDACs/RS_New_Analysis/Datasets/meta_effectsize.RData")
metacluster1_genes=abs(mm$Metacluster1) >0.5
metacluster2_genes=abs(mm$Metacluster2) >0.5
metacluster3_genes=abs(mm$Metacluster3) >0.5
#metacluster4_genes=abs(mm$Metacluster4) >0.5


c3 <- cbind(metacluster1_genes, metacluster2_genes, metacluster3_genes)
a <- vennCounts(c3)

vennDiagram(a, names = c("Metacluster1","Metacluster2","Metacluster3"), counts.col = "black",
            circle.col=c("pink","green","gold"))



### Druggable genes
druggable_genes=read.table("/Users/vandanasandhu/Desktop/Subtyping_PDACs/Datasets/druggable_genes.txt")
drug_gene_int= read.csv("/Users/vandanasandhu/Desktop/Subtyping_PDACs/Datasets/drug_gene_interaction.txt",  header=TRUE, sep="\t")

drug_names = unique(drug_gene_int$drug[which(drug_gene_int$drug %in% colnames(ctrpv2_ic50))])
dd=  drug_gene_int[which(drug_gene_int$drug %in% colnames(ctrpv2_ic50)),]

metacluster1_genes=mm[order(abs(mm$Metacluster1), decreasing = TRUE)[1:1000],]
metacluster2_genes=mm[order(abs(mm$Metacluster2), decreasing = TRUE)[1:1000],]
metacluster3_genes=mm[order(abs(mm$Metacluster3), decreasing = TRUE)[1:1000],]
#metacluster4_genes=mm[order(abs(mm$Metacluster4), decreasing = TRUE)[1:100],]

mm1=metacluster1_genes[which(rownames(metacluster1_genes) %in% druggable_genes[,1]),]
mm2=metacluster2_genes[which(rownames(metacluster2_genes) %in% druggable_genes[,1]),]
mm3=metacluster3_genes[which(rownames(metacluster3_genes) %in% druggable_genes[,1]),]
#mm4=metacluster4_genes[which(rownames(metacluster4_genes) %in% druggable_genes[,1]),]

mm1=mm1[order(abs(mm1$Metacluster1), decreasing = TRUE),]
mm2=mm2[order(abs(mm2$Metacluster2), decreasing = TRUE),]
mm3=mm3[order(abs(mm3$Metacluster3), decreasing = TRUE),]
#mm4=mm4[order(abs(mm4$Metacluster4), decreasing = TRUE),]

##################### Gene enrichment


reference_genes= rownames(cohorts1$pcsi)

hallmark= loadGSC("/Users/vandanasandhu/Desktop/RS_Remodelling_TSP_project/Gene-sets/h.all.v6.1.symbols.gmt")
results_hallmark=runGSAhyper(rownames(metacluster1_genes),gsc= hallmark,  adjMethod="fdr",  universe= reference_genes)
results_hallmark$resTab[which(results_hallmark$resTab[,2]<0.05),]
clust1=qnorm(1 - (results_hallmark$pvalues/2)) 

results_hallmark=runGSAhyper(rownames(metacluster2_genes),gsc= hallmark,  adjMethod="fdr",  universe= reference_genes)
results_hallmark$resTab[which(results_hallmark$resTab[,2]<0.05),]
clust2=qnorm(1 - (results_hallmark$pvalues/2)) 

results_hallmark=runGSAhyper(rownames(metacluster3_genes),gsc= hallmark,  adjMethod="fdr",  universe= reference_genes)
results_hallmark$resTab[which(results_hallmark$resTab[,1]<0.05),]
clust3=qnorm(1 - (results_hallmark$pvalues/2)) 

#results_hallmark=runGSAhyper(rownames(metacluster4_genes),gsc= hallmark,  adjMethod="fdr",  universe= reference_genes)
#results_hallmark$resTab[which(results_hallmark$resTab[,2]<0.05),]
#clust4=qnorm(1 - (results_hallmark$pvalues/2)) 


zz=data.frame(clust1=clust1, clust2=clust2, clust3=clust3)
zz[,1][is.infinite(zz[,1])] = NA
rownames(zz)=sub("HALLMARK_","",rownames(zz))

colfunc <- colorRampPalette(c("gold", "white", "black"))
heatmap(data.matrix(zz),Rowv=NA, Colv=NA, col = cm.colors(256), scale="column", margins=c(5,10))

heatmap(data.matrix(zz),Rowv=NA, Colv=NA, col=colfunc(15), scale="column", margins=c(5,10))

heatmap.2(data.matrix(zz),Rowv=TRUE, Colv=NA, col=colfunc(15), scale="column", margins=c(5,10), trace="none",dendrogram="row")

####### Conanical 
genesets=loadGSC("/Users/vandanasandhu/Desktop/RS_Remodelling_TSP_project/Gene-sets/c2.cp.v6.1.symbols.gmt", type="auto")

results=runGSAhyper(rownames(metacluster1_genes),gsc= genesets,  adjMethod="fdr", universe= reference_genes)
results$resTab[which(results$p.adj < 0.05),]
clust1=qnorm(1 - (results$pvalues/2)) 
pval1=results$pvalues

results=runGSAhyper(rownames(metacluster2_genes),gsc= genesets,  adjMethod="fdr", universe= reference_genes)
results$resTab[which(results$p.adj < 0.05),]
clust2=qnorm(1 - (results$pvalues/2)) 
pval2=results$pvalues

results=runGSAhyper(rownames(metacluster3_genes),gsc= genesets,  adjMethod="fdr", universe= reference_genes)
results$resTab[which(results$pvalues < 0.01),]
clust3=qnorm(1 - (results$pvalues/2)) 
pval3=results$pvalues
# 
# results=runGSAhyper(rownames(metacluster4_genes),gsc= genesets,  adjMethod="fdr", universe= reference_genes)
# results$resTab[which(results$p.adj < 0.05),]
# clust4=qnorm(1 - (results$pvalues/2)) 
# pval4=results$pvalues


zz=data.frame(clust1=clust1, clust2=clust2, clust3=clust3)
pp=data.frame(clust1=pval1, clust2=pval2, clust3=pval3)

count=0
ss=vector()
for(i in 1:dim(pp)[1]){
  if( length(which(pp[i,1:3]>0.05))==3){
    ss[count]=i
    count=count+1
  }
}

pp=pp[-ss,]
zz=zz[-ss,]

zz[,1][is.infinite(zz[,1])] = NA
#rownames(zz)=mapvalues(rownames(zz), c("KEGG_","REACTOME_","PID_","BIOCARTA_"),c("","","",""))
colnames(zz)= c("Basal","Exocrine","Classical")

colfunc <- colorRampPalette(c("gold", "white", "black"))
heatmap(data.matrix(zz),Rowv=NA, Colv=NA, col = cm.colors(256), scale="column", margins=c(5,10))

heatmap(data.matrix(zz),Rowv=TRUE, Colv=NA, col=colfunc(15), scale="column", margins=c(5,10))
pdf("/Users/vandanasandhu/Desktop/Analysis_results/conanical_pathways.pdf", width = 20, height = 100)
heatmap.2(data.matrix(zz),Rowv=TRUE, Colv=NA, col=colfunc(15), scale="column", margins=c(5,10),lhei = c(0.05, 0.20), trace="none",dendrogram="row")


heatmap.2(data.matrix(zz),Rowv=TRUE, Colv=NA, col=colfunc(15), scale="column", margins=c(5,10), trace="none",dendrogram="row")

######IPA

path_Set=load("/Users/vandanasandhu/Desktop/os_predictor_project/IPA_pathways.RData")
path_SetS=loadGSC(path_set)
results=runGSAhyper(rownames(metacluster1_genes),gsc=path_SetS,  adjMethod="fdr", universe= reference_genes)
clust1=qnorm(1 - (results$pvalues/2)) 
pval1=results$pvalues

results=runGSAhyper(rownames(metacluster2_genes),gsc=path_SetS,  adjMethod="fdr", universe= reference_genes)
clust2=qnorm(1 - (results$pvalues/2)) 
pval2=results$pvalues

results=runGSAhyper(rownames(metacluster3_genes),gsc=path_SetS,  adjMethod="fdr", universe= reference_genes)
clust3=qnorm(1 - (results$pvalues/2)) 
pval3=results$pvalues

results=runGSAhyper(rownames(metacluster4_genes),gsc=path_SetS,  adjMethod="fdr", universe= reference_genes)
clust4=qnorm(1 - (results$pvalues/2)) 
pval4=results$pvalues

zz=data.frame(clust1=clust1, clust2=clust2, clust3=clust3, clust4=clust4)
pp=data.frame(clust1=pval1, clust2=pval2, clust3=pval3, clust4=pval4)

count=0
ss=vector()
for(i in 1:dim(pp)[1]){
  if( length(which(pp[i,1:4]>0.05))==4){
    ss[count]=i
    count=count+1
  }
}

pp=pp[-ss,]
zz=zz[-ss,]

zz[,1][is.infinite(zz[,1])] = NA
#rownames(zz)=mapvalues(rownames(zz), c("KEGG_","REACTOME_","PID_","BIOCARTA_"),c("","","",""))

colfunc <- colorRampPalette(c("gold", "white", "black"))
heatmap(data.matrix(zz),Rowv=NA, Colv=NA, col = cm.colors(256), scale="column", margins=c(5,10))

heatmap(data.matrix(zz),Rowv=TRUE, Colv=NA, col=colfunc(15), scale="column", margins=c(5,10))
pdf("/Users/vandanasandhu/Desktop/Subtyping_PDACs/Figures/ipa_pathways.pdf", width = 20, height = 100)
heatmap.2(data.matrix(zz),Rowv=TRUE, Colv=NA, col=colfunc(15), scale="column", margins=c(5,10),lhei = c(0.05, 0.20), trace="none",dendrogram="row")


heatmap.2(data.matrix(zz),Rowv=TRUE, Colv=NA, col=colfunc(15), scale="column", margins=c(5,10), trace="none",dendrogram="row")



##### GO Molecular function
genesets=loadGSC("/Users/vandanasandhu/Desktop/RS_Remodelling_TSP_project/Gene-sets/c5.mf.v6.1.symbols.gmt", type="auto")
results=runGSAhyper(rownames(metacluster1_genes),gsc= genesets,  adjMethod="fdr", universe= reference_genes)
results$resTab[which(results$pvalues < 0.05),]
clust1=qnorm(1 - (results$pvalues/2)) 
pval1=results$pvalues



results=runGSAhyper(rownames(metacluster2_genes),gsc= genesets,  adjMethod="fdr", universe= reference_genes)
results$resTab[which(results$pvalues < 0.05),]
clust2=qnorm(1 - (results$pvalues/2)) 
pval2=results$pvalues



results=runGSAhyper(rownames(metacluster3_genes),gsc= genesets,  adjMethod="fdr", universe= reference_genes)
results$resTab[which(results$pvalues < 0.05),]

clust3=qnorm(1 - (results$pvalues/2)) 
pval3=results$pvalues



zz=data.frame(clust1=clust1, clust2=clust2, clust3=clust3)
pp=data.frame(clust1=pval1, clust2=pval2, clust3=pval3)

count=0
ss=vector()
for(i in 1:dim(pp)[1]){
  if( length(which(pp[i,1:3]>0.05))==3){
    ss[count]=i
    count=count+1
  }
}

pp=pp[-ss,]
zz=zz[-ss,]

zz[,1][is.infinite(zz[,1])] = NA
#rownames(zz)=mapvalues(rownames(zz), c("KEGG_","REACTOME_","PID_","BIOCARTA_"),c("","","",""))
colnames(zz)= c("Basal","Exocrine","Classical")
colfunc <- colorRampPalette(c("gold", "white", "black"))
heatmap(data.matrix(zz),Rowv=NA, Colv=NA, col = cm.colors(256), scale="column", margins=c(5,10))

heatmap(data.matrix(zz),Rowv=TRUE, Colv=NA, col=colfunc(15), scale="column", margins=c(5,10), col)
pdf("/Users/vandanasandhu/Desktop/Subtyping_PDACs/RS_New_Analysis/GO_MF.pdf", width = 20, height = 100)
heatmap.2(data.matrix(zz),Rowv=TRUE, Colv=NA, col=colfunc(15), scale="column", margins=c(5,10),lhei = c(0.05, 0.20), trace="none",dendrogram="row")


heatmap.2(data.matrix(zz),Rowv=TRUE, Colv=NA, col=colfunc(15), scale="column", margins=c(5,10), trace="none",dendrogram="row")




##### GO Cellular component
genesets=loadGSC("/Users/vandanasandhu/Desktop/RS_Remodelling_TSP_project/Gene-sets/c5.cc.v6.1.symbols.gmt", type="auto")
results=runGSAhyper(rownames(metacluster1_genes),gsc= genesets,  adjMethod="fdr", universe= reference_genes)
results$resTab[which(results$pvalues < 0.05),]


results=runGSAhyper(rownames(metacluster2_genes),gsc= genesets,  adjMethod="fdr", universe= reference_genes)
results$resTab[which(results$pvalues < 0.05),]

results=runGSAhyper(rownames(metacluster3_genes),gsc= genesets,  adjMethod="fdr", universe= reference_genes)
results$resTab[which(results$pvalues < 0.05),]











