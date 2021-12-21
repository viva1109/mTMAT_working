# data/Lon_data
fol<-"/data/sharedcode/kjkim/Analysis/20190225Combine"

ID_linker<-read.csv("~/analysis/20200928_longi_method/data/linked_IDs_GG_EZ.csv",header=T,stringsAsFactors=F)

rd<-readRDS(file.path(fol,"M678_Feb2019.rds"))
rd1<-rd[[1]]

tax1<-rd1[,1]
otutab<-rd1[,-c(1,dim(rd1)[2])]
library(stringr)
OTU_out1<-which(str_detect(tax1,"New.CleanUp"))
names(rd)

sampleID<-colnames(otutab)

Rcnt<-lapply(rd,function(data){
  otutab<-data[,-c(1,dim(data)[2])]
  Rcnts<-apply(otutab,2,sum)
})

remain_names<-lapply(Rcnt,function(cnts){
  names(cnts)[cnts>3000]
})


ID_linker_left<-c(ID_linker[,1],ID_linker[,1],ID_linker[,1])
ID_linker_right<-c(ID_linker[,2],ID_linker[,3],ID_linker[,4])
getUid<-lapply(remain_names,function(nms){
  ID_linker_left[match(nms,ID_linker_right)]
})


# rd2<-readRDS("rd_content678.rds")

nms_go<-intersect(intersect(getUid[[1]],getUid[[2]]),getUid[[3]])
nms_go<-nms_go[!is.na(nms_go)]
linked_IDs<-ID_linker[match(nms_go,ID_linker[,1]),2:4]       

NewRDs<-lapply(1:length(rd),function(ind_data){
  data<-rd[[ind_data]]
  otutab<-data[,-c(1,dim(data)[2])]
  cbind(data[,1],otutab[,linked_IDs[,ind_data]],data[,dim(data)[2]])
})




otuTab_bf<-lapply(NewRDs,function(data){
  otutab<-data[,-c(1,dim(data)[2])]
})
Info_IDs<-lapply(otuTab_bf,function(otutab){
  Rcnts<-apply(otutab,2,sum)
  timeID<-names(otutab)
  return(list(Rcnts=Rcnts,timeID=timeID))
})

library(parallel)

otu_info<-mclapply(1:length(otuTab_bf),function(ind_data){
  data<-otuTab_bf[[ind_data]]
  zero_prop<-apply(data,1,function(value){
    sum(value!=0)/length(value)
  })
  trc<-Info_IDs[[ind_data]]$Rcnts
  mean_rp<-apply(data,1,function(value){
    mean(value/trc)
  })
  
  list(zero_prop=zero_prop,mean_rp=mean_rp)
},mc.cores=length(otuTab_bf))

otu_info_tax<-mclapply(1:length(otuTab_bf),function(ind_data){
  data<-NewRDs[[ind_data]]
  oid<-as.character(data[,1])
  tax<-as.character(data[,dim(data)[2]])
  list(oid,tax)
},mc.cores=length(otuTab_bf))




otu_info2<-lapply(otu_info,function(info){
  cri1<-which(info[[1]]>0.5)
  cri2<-which(info[[2]]>0.001)
  cri3<-intersect(cri1,cri2)
  list(cri1,cri2,cri3)
})

lapply(otu_info2,sapply,length)

indsx<-lapply(otu_info2,"[[",2)
otu_remain<-intersect(intersect(indsx[[1]],indsx[[2]]),indsx[[3]])

otu_info_tax[[1]][[2]][otu_remain]
otu_info_tax[[2]][[2]][otu_remain]
otu_info_tax[[3]][[2]][otu_remain]

otu_remain2<-which(sapply(str_split(otu_info_tax[[1]][[2]],";"),"[",6)!=" g__")
otu_remiain3<-intersect(otu_remain,otu_remain2)

otu_remiain3

OTU_INFO<-list(otu_info_tax[[1]][[1]][otu_remiain3],otu_info_tax[[1]][[2]][otu_remiain3],sapply(str_split(otu_info_tax[[1]][[2]][otu_remiain3],"; "),"[",6))

OTU_table<-lapply(otuTab_bf,function(data){
  data_out<-t(data[otu_remiain3,])
  rownames(data_out)<-nms_go
  colnames(data_out)<-OTU_INFO[[3]]
  return(data_out)
})


Tax_gg<-fread("/data/sharedcode/scpark/bash/16S_pipeline/file/97_otu_taxonomy.txt",header=F,sep="\t")
  
  