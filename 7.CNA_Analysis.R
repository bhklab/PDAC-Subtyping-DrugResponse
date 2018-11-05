
a=load("/Users/vandanasandhu/Desktop/Subtyping_PDACs/RS_New_Analysis/CNA_OUH_TCGA/Oslo_TCGA_score.RData")
SNPpos <- read.table("/Users/vandanasandhu/Desktop/Subtyping_PDACs/RS_New_Analysis/CNA_OUH_TCGA/SNPpos.txt",header=T,sep="\t",row.names=1, stringsAsFactors=F)
ss=read.table("/Users/vandanasandhu/Desktop/Subtyping_PDACs/RS_New_Analysis/CNA_OUH_TCGA/rs-snp-annot.txt",sep="\t",header=T)

score_file = output_score_sc
rownames(score_file) = score_file$ProbeID
score_file = score_file[,2:ncol(score_file)]

resp_in= data.frame(SampleID= names(clusters$OUH$meta_classes), Resp= as.numeric(clusters$OUH$meta_classes))
names(clusters$OUH$meta_classes)= mgsub( "X","",names(clusters$OUH$meta_classes))

resp_in=  resp_in[which(resp_in$SampleID %in%  colnames(score_file)),]

p1.sampleID = resp_in[which(resp_in$Resp==1), 1]
p2.sampleID = resp_in[which(resp_in$Resp==2), 1]
p3.sampleID = resp_in[which(resp_in$Resp==3), 1]
#p4.sampleID = resp_in[which(resp_in$Resp==4), 1]
#p5.sampleID = resp_in[which(resp_in$Resp==5), 1]

countNonzero <- function(data) {
  countUp = apply(data, 1, function(r) { sum(r == 1) } )
  countDown = apply(data, 1, function(r) { sum(r == -1) } )
  res = data.frame(ProbeID = rownames(data), up=countUp, down=countDown)
  res$up = 100.0 * (res$up / ncol(data))
  res$down = 100.0 * (res$down / ncol(data))
  res = merge(res, SNPpos, by.x="ProbeID", by.y="row.names")
  return(res)
}



p1.count = countNonzero(score_file[, p1.sampleID])
p2.count = countNonzero(score_file[, p2.sampleID])
p3.count = countNonzero(score_file[, p3.sampleID])


p1.count$what=rep("Basal", nrow(p1.count))
p2.count$what=rep("Exocrine", nrow(p2.count))
p3.count$what=rep("Classical", nrow(p3.count))


ptotal = rbind(p1.count, p2.count,p3.count)
ptotal$Chr = factor(as.character(ptotal$Chr), levels=c(as.character(1:22), "X"))

write.table(ptotal,"ptotal_all_gp.txt")

library(ggplot2)

png(filename="/Users/vandanasandhu/Desktop/OUH.png", width=4096, height=1024)
ggplot(ptotal) + 
  geom_segment(color="red", aes(x=Position, xend=Position, y=0, yend = up)) +
  geom_segment(color="green", aes(x=Position, xend=Position, y=0, yend = -down)) +
  facet_grid(what~Chr , scales="free_x", space="free_x") +
  theme_bw() +
  theme(axis.text.x = element_blank() ) +
  ggtitle("PDAC subtypes") + xlab("Genomic position") + ylab("Frequency (%)")
#dev.off()

####################################
resp_in= data.frame(SampleID= names(clusters$TCGA$meta_classes), Resp= as.numeric(clusters$TCGA$meta_classes))

resp_in=  resp_in[which(resp_in$SampleID %in%  colnames(score_file)),]

p1.sampleID = resp_in[which(resp_in$Resp==1), 1]
p2.sampleID = resp_in[which(resp_in$Resp==2), 1]
p3.sampleID = resp_in[which(resp_in$Resp==3), 1]
#p4.sampleID = resp_in[which(resp_in$Resp==4), 1]
#p5.sampleID = resp_in[which(resp_in$Resp==5), 1]

countNonzero <- function(data) {
  countUp = apply(data, 1, function(r) { sum(r == 1) } )
  countDown = apply(data, 1, function(r) { sum(r == -1) } )
  res = data.frame(ProbeID = rownames(data), up=countUp, down=countDown)
  res$up = 100.0 * (res$up / ncol(data))
  res$down = 100.0 * (res$down / ncol(data))
  res = merge(res, SNPpos, by.x="ProbeID", by.y="row.names")
  return(res)
}



p1.count = countNonzero(score_file[, p1.sampleID])
p2.count = countNonzero(score_file[, p2.sampleID])
p3.count = countNonzero(score_file[, p3.sampleID])


p1.count$what=rep("Basal", nrow(p1.count))
p2.count$what=rep("Exocrine", nrow(p2.count))
p3.count$what=rep("Classical", nrow(p3.count))


ptotal = rbind(p1.count, p2.count,p3.count)
ptotal$Chr = factor(as.character(ptotal$Chr), levels=c(as.character(1:22), "X"))

write.table(ptotal,"ptotal_all_gp.txt")

library(ggplot2)

png(filename="/Users/vandanasandhu/Desktop/TCGA.png", width=4096, height=1024)
ggplot(ptotal) + 
  geom_segment(color="red", aes(x=Position, xend=Position, y=0, yend = up)) +
  geom_segment(color="green", aes(x=Position, xend=Position, y=0, yend = -down)) +
  facet_grid(what~Chr , scales="free_x", space="free_x") +
  theme_bw() +
  theme(axis.text.x = element_blank() ) +
  ggtitle("PDAC subtypes") + xlab("Genomic position") + ylab("Frequency (%)")
#dev.off()
###################

rr= cbind(ptotal[1:901042,], ptotal[901043:1802084,], ptotal[1802085:2703126,])
cc=rr

p=vector()

for( i in 1:nrow(cc)){
  x=c(cc[i,2],cc[i,8],cc[i,14],cc[i,3],cc[i,9],cc[i,15])
  if(sum(x) != 0){
  x=matrix(x, nrow=3,ncol=2)
  p[i]=chisq.test(x)$p.value
  }
  else{
    p[i]=NA
  }
  
} 
adj=p.adjust(p)
zz=cbind(cc,p,adj)
colnames(zz)[1]="id"
zz=zz[which(adj <0.05),]



mm=merge(zz,ss, by="id")
write.table(mm,"/Users/vandanasandhu/Desktop/1.txt", sep="\t")

