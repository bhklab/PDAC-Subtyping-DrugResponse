
load("/Users/vandanasandhu/Desktop/Subtyping_PDACs/Datasets/pharmacogx_pdac_celllines1.RData")

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

gcsi_subtypes = meta_subtype_cellines(gcsi)
ctrpv2_subtypes =meta_subtype_cellines(ctrpv2)
gdsc_subtypes = meta_subtype_cellines(gdsc)
ccle_subtypes = meta_subtype_cellines(ccle)


aa=merge(gcsi_subtypes, ccle_subtypes, by.x="id", by.y="id")
length(which(as.character(aa[,2]) == as.character(aa[,6])))/dim(aa)[1]

aa=merge( ccle_subtypes, ctrpv2_subtypes, by.x="id", by.y="id")
length(which(as.character(aa[,2]) == as.character(aa[,6])))/dim(aa)[1]


aa=merge(gcsi_subtypes,  ctrpv2_subtypes, by.x="id", by.y="id")
length(which(as.character(aa[,2]) == as.character(aa[,6])))/dim(aa)[1]


aa=merge(gdsc_subtypes, ccle_subtypes, by.x="id", by.y="id")
 length(which(as.character(aa[,2]) == as.character(aa[,6])))/dim(aa)[1]
# aa[,1][which(as.character(aa[,2]) != as.character(aa[,3]))]
# 
aa=merge(gcsi_subtypes, gdsc_subtypes, by.x="id", by.y="id")
length(which(as.character(aa[,2]) == as.character(aa[,6])))/dim(aa)[1]
# aa[,1][which(as.character(aa[,2]) != as.character(aa[,3]))]
# 
aa=merge(ctrpv2_subtypes, gdsc_subtypes, by.x="id", by.y="id")
length(which(as.character(aa[,2]) == as.character(aa[,6])))/dim(aa)[1]


seq=collisson_subtyping(cohorts$icgc_seq, collisson_centroid)
arr=collisson_subtyping(cohorts$icgc_arr, collisson_centroid)
seq= data.frame(id=seq[[2]], subtype=seq[[3]])
arr= data.frame(id=arr[[2]], subtype=arr[[3]])


aa=merge(seq, arr, by.x="id", by.y="id")
length(which(as.character(aa[,2]) == as.character(aa[,3])))
length(which(as.character(aa[,2]) == as.character(aa[,3])))/dim(aa)[1]


############################# 2. Bailey

gcsi_subtypes = bailey_subtyping(gcsi,bailey_centroid)
gdsc_subtypes = bailey_subtyping(gdsc,bailey_centroid)
ccle_subtypes = bailey_subtyping(ccle,bailey_centroid)
ctrpv2_subtypes = bailey_subtyping(ctrpv2,bailey_centroid)

gcsi_subtypes= data.frame(id=gcsi_subtypes[[2]], subtype=gcsi_subtypes[[3]])
gdsc_subtypes= data.frame(id=gdsc_subtypes[[2]], subtype=gdsc_subtypes[[3]])
ccle_subtypes= data.frame(id=ccle_subtypes[[2]], subtype=ccle_subtypes[[3]])
ctrpv2_subtypes= data.frame(id=ctrpv2_subtypes[[2]], subtype=ctrpv2_subtypes[[3]])


aa=merge(gcsi_subtypes, ccle_subtypes, by.x="id", by.y="id")
length(which(as.character(aa[,2]) == as.character(aa[,3])))/dim(aa)[1]

aa=merge( ccle_subtypes, ctrpv2_subtypes, by.x="id", by.y="id")
length(which(as.character(aa[,2]) == as.character(aa[,3])))/dim(aa)[1]


aa=merge(gcsi_subtypes,  ctrpv2_subtypes, by.x="id", by.y="id")
length(which(as.character(aa[,2]) == as.character(aa[,3])))/dim(aa)[1]


aa=merge(gdsc_subtypes, ccle_subtypes, by.x="id", by.y="id")
length(which(as.character(aa[,2]) == as.character(aa[,3])))/dim(aa)[1]
# aa[,1][which(as.character(aa[,2]) != as.character(aa[,3]))]
# 
aa=merge(gcsi_subtypes, gdsc_subtypes, by.x="id", by.y="id")
length(which(as.character(aa[,2]) == as.character(aa[,3])))/dim(aa)[1]
# aa[,1][which(as.character(aa[,2]) != as.character(aa[,3]))]
# 
aa=merge(ctrpv2_subtypes, gdsc_subtypes, by.x="id", by.y="id")
length(which(as.character(aa[,2]) == as.character(aa[,3])))/dim(aa)[1]

seq=meta_subtype_cellines(cohorts$icgc_seq)
arr=meta_subtype_cellines(cohorts$icgc_arr)

aa=merge(seq, arr, by.x="id", by.y="id")
length(which(as.character(aa[,2]) == as.character(aa[,6])))
length(which(as.character(aa[,2]) == as.character(aa[,6])))/dim(aa)[1]


############################# 3. Collisson

gcsi_subtypes = collisson_subtyping(gcsi,collisson_centroid)
gdsc_subtypes = collisson_subtyping(gdsc,collisson_centroid)
ccle_subtypes = collisson_subtyping(ccle,collisson_centroid)
ctrpv2_subtypes = collisson_subtyping(ctrpv2,collisson_centroid)

gcsi_subtypes= data.frame(id=gcsi_subtypes[[2]], subtype=gcsi_subtypes[[3]])
gdsc_subtypes= data.frame(id=gdsc_subtypes[[2]], subtype=gdsc_subtypes[[3]])
ccle_subtypes= data.frame(id=ccle_subtypes[[2]], subtype=ccle_subtypes[[3]])
ctrpv2_subtypes= data.frame(id=ctrpv2_subtypes[[2]], subtype=ctrpv2_subtypes[[3]])


aa=merge(gcsi_subtypes, ccle_subtypes, by.x="id", by.y="id")
length(which(as.character(aa[,2]) == as.character(aa[,3])))/dim(aa)[1]

aa=merge( ccle_subtypes, ctrpv2_subtypes, by.x="id", by.y="id")
length(which(as.character(aa[,2]) == as.character(aa[,3])))/dim(aa)[1]


aa=merge(gcsi_subtypes,  ctrpv2_subtypes, by.x="id", by.y="id")
length(which(as.character(aa[,2]) == as.character(aa[,3])))/dim(aa)[1]


aa=merge(gdsc_subtypes, ccle_subtypes, by.x="id", by.y="id")
length(which(as.character(aa[,2]) == as.character(aa[,3])))/dim(aa)[1]
# aa[,1][which(as.character(aa[,2]) != as.character(aa[,3]))]
# 
aa=merge(gcsi_subtypes, gdsc_subtypes, by.x="id", by.y="id")
length(which(as.character(aa[,2]) == as.character(aa[,3])))/dim(aa)[1]
# aa[,1][which(as.character(aa[,2]) != as.character(aa[,3]))]
# 
aa=merge(ctrpv2_subtypes, gdsc_subtypes, by.x="id", by.y="id")
length(which(as.character(aa[,2]) == as.character(aa[,3])))/dim(aa)[1]

seq=collisson_subtyping(cohorts$icgc_seq, collisson_centroid)
arr=collisson_subtyping(cohorts$icgc_arr, collisson_centroid)
seq= data.frame(id=seq[[2]], subtype=seq[[3]])
arr= data.frame(id=arr[[2]], subtype=arr[[3]])


aa=merge(seq, arr, by.x="id", by.y="id")
length(which(as.character(aa[,2]) == as.character(aa[,3])))
length(which(as.character(aa[,2]) == as.character(aa[,3])))/dim(aa)[1]

############################# 4.Moffitt

gcsi_subtypes = moffitt_subtyping(gcsi,moff_centroid)
gdsc_subtypes = moffitt_subtyping(gdsc,moff_centroid)
ccle_subtypes = moffitt_subtyping(ccle,moff_centroid)
ctrpv2_subtypes = moffitt_subtyping(ctrpv2,moff_centroid)

gcsi_subtypes= data.frame(id=gcsi_subtypes[[2]], subtype=gcsi_subtypes[[3]])
gdsc_subtypes= data.frame(id=gdsc_subtypes[[2]], subtype=gdsc_subtypes[[3]])
ccle_subtypes= data.frame(id=ccle_subtypes[[2]], subtype=ccle_subtypes[[3]])
ctrpv2_subtypes= data.frame(id=ctrpv2_subtypes[[2]], subtype=ctrpv2_subtypes[[3]])


aa=merge(gcsi_subtypes, ccle_subtypes, by.x="id", by.y="id")
length(which(as.character(aa[,2]) == as.character(aa[,3])))/dim(aa)[1]

aa=merge( ccle_subtypes, ctrpv2_subtypes, by.x="id", by.y="id")
length(which(as.character(aa[,2]) == as.character(aa[,3])))/dim(aa)[1]


aa=merge(gcsi_subtypes,  ctrpv2_subtypes, by.x="id", by.y="id")
length(which(as.character(aa[,2]) == as.character(aa[,3])))/dim(aa)[1]


aa=merge(gdsc_subtypes, ccle_subtypes, by.x="id", by.y="id")
length(which(as.character(aa[,2]) == as.character(aa[,3])))/dim(aa)[1]
# aa[,1][which(as.character(aa[,2]) != as.character(aa[,3]))]
# 
aa=merge(gcsi_subtypes, gdsc_subtypes, by.x="id", by.y="id")
length(which(as.character(aa[,2]) == as.character(aa[,3])))/dim(aa)[1]
# aa[,1][which(as.character(aa[,2]) != as.character(aa[,3]))]
# 
aa=merge(ctrpv2_subtypes, gdsc_subtypes, by.x="id", by.y="id")
length(which(as.character(aa[,2]) == as.character(aa[,3])))/dim(aa)[1]

#########
seq=moffitt_subtyping(cohorts$icgc_seq, moff_centroid)
arr=moffitt_subtyping(cohorts$icgc_arr, moff_centroid)
seq= data.frame(id=seq[[2]], subtype=seq[[3]])
arr= data.frame(id=arr[[2]], subtype=arr[[3]])


aa=merge(seq, arr, by.x="id", by.y="id")
length(which(as.character(aa[,2]) == as.character(aa[,3])))
length(which(as.character(aa[,2]) == as.character(aa[,3])))/dim(aa)[1]

#############################5. PAM50

gcsi_subtypes = pam50_subtyping(gcsi,pam50_centroid)
gdsc_subtypes = pam50_subtyping(gdsc,pam50_centroid)
ccle_subtypes = pam50_subtyping(ccle,pam50_centroid)
ctrpv2_subtypes = pam50_subtyping(ctrpv2,pam50_centroid)

gcsi_subtypes= data.frame(id=gcsi_subtypes[[2]], subtype=gcsi_subtypes[[3]])
gdsc_subtypes= data.frame(id=gdsc_subtypes[[2]], subtype=gdsc_subtypes[[3]])
ccle_subtypes= data.frame(id=ccle_subtypes[[2]], subtype=ccle_subtypes[[3]])
ctrpv2_subtypes= data.frame(id=ctrpv2_subtypes[[2]], subtype=ctrpv2_subtypes[[3]])


aa=merge(gcsi_subtypes, ccle_subtypes, by.x="id", by.y="id")
length(which(as.character(aa[,2]) == as.character(aa[,3])))/dim(aa)[1]

aa=merge( ccle_subtypes, ctrpv2_subtypes, by.x="id", by.y="id")
length(which(as.character(aa[,2]) == as.character(aa[,3])))/dim(aa)[1]


aa=merge(gcsi_subtypes,  ctrpv2_subtypes, by.x="id", by.y="id")
length(which(as.character(aa[,2]) == as.character(aa[,3])))/dim(aa)[1]


aa=merge(gdsc_subtypes, ccle_subtypes, by.x="id", by.y="id")
length(which(as.character(aa[,2]) == as.character(aa[,3])))/dim(aa)[1]
# aa[,1][which(as.character(aa[,2]) != as.character(aa[,3]))]
# 
aa=merge(gcsi_subtypes, gdsc_subtypes, by.x="id", by.y="id")
length(which(as.character(aa[,2]) == as.character(aa[,3])))/dim(aa)[1]
# aa[,1][which(as.character(aa[,2]) != as.character(aa[,3]))]
# 
aa=merge(ctrpv2_subtypes, gdsc_subtypes, by.x="id", by.y="id")
length(which(as.character(aa[,2]) == as.character(aa[,3])))/dim(aa)[1]

#########
seq=pam50_subtyping(cohorts$icgc_seq, pam50_centroid)
arr=pam50_subtyping(cohorts$icgc_arr, pam50_centroid)
seq= data.frame(id=seq[[2]], subtype=seq[[3]])
arr= data.frame(id=arr[[2]], subtype=arr[[3]])


aa=merge(seq, arr, by.x="id", by.y="id")
length(which(as.character(aa[,2]) == as.character(aa[,3])))/dim(aa)[1]


