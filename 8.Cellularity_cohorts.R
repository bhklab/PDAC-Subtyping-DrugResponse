
pcsi_cellularity = read.table("/Users/vandanasandhu/Desktop/Subtyping_PDACs/RS_New_Analysis/Datasets/purity/pcsi_cellularity.txt", header=TRUE, sep="\t")
pp = data.frame(id= names(clusters$PCSI$meta_classes), classes= clusters$PCSI$meta_classes)
mm = merge(pcsi_cellularity, pp, by.x="Tumor", by.y = "id")
mm$cohorts= rep("PCSI", length(mm$classes))
colnames(mm) = c("Id","Tumor percentage", "Subtypes", "Cohorts")


icgc_cellularity = read.table("/Users/vandanasandhu/Desktop/Subtyping_PDACs/RS_New_Analysis/Datasets/purity/icgc_seq_cellularity.txt", header=TRUE, sep="\t")
pp = data.frame(id= names(clusters$ICGC_seq$meta_classes), classes= clusters$ICGC_seq$meta_classes)
mm1 = merge(icgc_cellularity, pp, by.x="Icgc_id", by.y = "id")
mm1$cohorts= rep("ICGC", length(mm1$classes))
mm1=mm1[,c(-2,-3,-4,-6)]
colnames(mm1) = c("Id","Tumor percentage", "Subtypes", "Cohorts")



tcga_cellularity = read.table("/Users/vandanasandhu/Desktop/Subtyping_PDACs/RS_New_Analysis/Datasets/purity/tcga_cellularity.txt", header=TRUE, sep="\t")
pp = data.frame(id= names(clusters$TCGA$meta_classes), classes= clusters$TCGA$meta_classes)
mm2 = merge(tcga_cellularity, pp, by.x="Sample_ID", by.y = "id")
mm2$cohorts= rep("TCGA", length(mm2$classes))
colnames(mm2) = c("Id","Tumor percentage", "Subtypes", "Cohorts")


ouh_cellularity = read.table("/Users/vandanasandhu/Desktop/Subtyping_PDACs/RS_New_Analysis/Datasets/ouh_cellularity.txt", header=TRUE, sep="\t")
pp = data.frame(id= names(clusters$OUH$meta_classes), classes= clusters$OUH$meta_classes)
mm3 = merge(ouh_cellularity, pp, by.x="id", by.y = "id")
mm3$cohorts= rep("OUH", length(mm3$classes))
mm3$tumor_percentage= mm3$tumor_percentage/100
colnames(mm3) = c("Id","Tumor percentage", "Subtypes", "Cohorts")

haider_cellularity = read.table("/Users/vandanasandhu/Desktop/Subtyping_PDACs/RS_New_Analysis/Datasets/purity/haider_cellularity.txt", header=TRUE, sep="\t")
pp = data.frame(id= names(clusters$haider$meta_classes), classes= clusters$haider$meta_classes)
mm4 = merge(haider_cellularity, pp, by.x="id", by.y = "id")
mm4$cellularity =mm4$cellularity/100
mm4$cohorts= rep("Haider", length(mm4$classes))
colnames(mm4) = c("Id","Tumor percentage", "Subtypes", "Cohorts")


mat= rbind(mm,mm1,mm2,mm3,mm4)
colnames(mat)= c("Id","Cellularity", "Subtypes", "Cohorts")
mat$Subtypes = mgsub(c("1","2","3"), c("Basal","Exocrine","Classical"), mat$Subtypes)


ggplot(mat, aes(Subtypes, Cellularity, fill=Cohorts, group=Cohorts )) + 
  geom_bar(stat='identity', position='dodge') +   scale_fill_brewer(palette = "Paired")

ggplot(mat, aes(Cohorts, Cellularity, fill=Subtypes, group=Subtypes )) + 
  geom_bar(stat='identity', position='dodge') +   scale_fill_brewer(palette = "Pastel1")


ggplot(data = mat, aes(Cohorts, Cellularity)) + geom_boxplot(aes(fill=Subtypes), width=0.5)+
  scale_fill_brewer(palette = "Pastel1") + ylab("Cellularity")+ 
  guides(fill=guide_legend(title="Cellularity"))+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

############################################################################
############################################################################
##### If you want to plot single cohorts for cellularity comparison
############################################################################

pcsi_cellularity = read.table("/Users/vandanasandhu/Desktop/Subtyping_PDACs/RS_New_Analysis/Datasets/purity/pcsi_cellularity.txt", header=TRUE, sep="\t")

pp = data.frame(id= names(clusters$PCSI$meta_classes), classes= clusters$PCSI$meta_classes)
mm = merge(pcsi_cellularity, pp, by.x="Tumor", by.y = "id")


kruskal.test(mm$Cellularity, mm$classes)

#mm= mm[-which(mm$classes == c(4)),]

mm$classes=mgsub(c("1","2","3"), c("Basal","Exocrine","Classical"), mm$classes)
my_comparisons <- list( c("Exocrine","Classical"), c("Classical","Basal"), c("Exocrine","Basal"))
###



p=ggboxplot(mm , x = "classes", y = "Cellularity", color= "classes",palette = "jco", add = "jitter",ylab = "Cellularity", xlab = "Subtype", title= "PCSI", label=NULL) 
p+stat_compare_means(comparisons = my_comparisons)


mm$color =mgsub(c(1,2), c("gold", "pink"), mm$classes )
mm=mm[order(mm$Cellularity),]
barplot(mm$Cellularity, col=mm$color, ylab= "cellularity", border = "white")


##################################################

icgc_cellularity = read.table("/Users/vandanasandhu/Desktop/Subtyping_PDACs/RS_New_Analysis/Datasets/purity/icgc_seq_cellularity.txt", header=TRUE, sep="\t")

pp = data.frame(id= names(clusters$ICGC_seq$meta_classes), classes= clusters$ICGC_seq$meta_classes)
mm = merge(icgc_cellularity, pp, by.x="Icgc_id", by.y = "id")


kruskal.test(mm$TumorPurity, mm$classes)


mm$classes=mgsub(c("1","2","3"),c("Basal","Exocrine","Classical"), mm$classes)
my_comparisons <- list( c("Exocrine","Classical"), c("Classical","Basal"), c("Exocrine","Basal"))
###


p=ggboxplot(mm , x = "classes", y = "TumorPurity", color= "classes",palette = "jco", add = "jitter",ylab = "Cellularity", xlab = "Subtype", title= "ICGC-Seq", label=NULL) 
p+stat_compare_means(comparisons = my_comparisons)

mm$color =mgsub(c(1,2,3), c("gold", "pink","cyan"), mm$classes )
mm=mm[order(mm$TumorPurity),]
barplot(mm$TumorPurity, col=mm$color, ylab= "cellularity", border = "white")



##################################################

tcga_cellularity = read.table("/Users/vandanasandhu/Desktop/Subtyping_PDACs/RS_New_Analysis/Datasets/purity/tcga_cellularity.txt", header=TRUE, sep="\t")

pp = data.frame(id= names(clusters$TCGA$meta_classes), classes= clusters$TCGA$meta_classes)
mm = merge(tcga_cellularity, pp, by.x="Sample_ID", by.y = "id")

#mm= mm[-which(mm$classes == c(4)),]
mm$classes=mgsub(c("1","2","3"), c("Basal","Exocrine","Classical"), mm$classes)

kruskal.test(mm$ABSOLUTE_Purity, mm$classes)

my_comparisons <- list( c("Exocrine","Classical"), c("Classical","Basal"), c("Exocrine","Basal"))
###



p=ggboxplot(mm , x = "classes", y = "ABSOLUTE_Purity", color= "classes",palette = "jco", add = "jitter",ylab = "Cellularity", xlab = "Subtype", title= "TCGA", label=NULL) 
p+stat_compare_means(comparisons = my_comparisons)


mm$color =mgsub(c(1,2,3), c("gold", "pink","cyan"), mm$classes )
mm=mm[order(mm$ABSOLUTE_Purity),]
barplot(mm$ABSOLUTE_Purity, col=mm$color, ylab= "cellularity", border = "white")

##################################################

haider_cellularity = read.table("/Users/vandanasandhu/Desktop/Subtyping_PDACs/RS_New_Analysis/Datasets/purity/haider_cellularity.txt", header=TRUE, sep="\t")

pp = data.frame(id= names(clusters$haider$meta_classes), classes= clusters$haider$meta_classes)
mm = merge(haider_cellularity, pp, by.x="id", by.y = "id")
mm$classes=mgsub(c("1","2","3"), c("Basal","Exocrine","Classical"), mm$classes)


kruskal.test(mm$cellularity, mm$classes)

my_comparisons <- list(  c("Exocrine","Basal"))
###

p=ggboxplot(mm , x = "classes", y = "cellularity", color= "classes",palette = "jco", add = "jitter",ylab = "Cellularity", xlab = "Subtype", title= "haider", label=NULL) 
p+stat_compare_means(comparisons = my_comparisons)


mm$color =mgsub(c(1,2,3), c("gold", "pink","cyan"), mm$classes )
mm=mm[order(mm$cellularity),]
barplot(mm$cellularity, col=mm$color, ylab= "cellularity", border = "white")


################

ouh_cellularity = read.table("/Users/vandanasandhu/Desktop/Subtyping_PDACs/RS_New_Analysis/Datasets/ouh_cellularity.txt", header=TRUE, sep="\t")

pp = data.frame(id= names(clusters$OUH$meta_classes), classes= clusters$OUH$meta_classes)
mm = merge(ouh_cellularity, pp, by.x="id", by.y = "id")
mm$classes=mgsub(c("1","2","3"), c("Basal","Exocrine","Classical"), mm$classes)


kruskal.test(mm$tumor_percentage, mm$classes)

my_comparisons <- list( c("Exocrine","Classical"), c("Classical","Basal"), c("Exocrine","Basal"))
###

p=ggboxplot(mm , x = "classes", y = "tumor_percentage", color= "classes",palette = "jco", add = "jitter",ylab = "tumor_percentage", xlab = "Subtype", title= "OUH", label=NULL) 
p+stat_compare_means(comparisons = my_comparisons)


mm$color =mgsub(c(1,2,3), c("gold", "pink","cyan"), mm$classes )
mm=mm[order(mm$cellularity),]
barplot(mm$cellularity, col=mm$color, ylab= "tumor_percentage", border = "white")

