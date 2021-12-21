infol<-"~/analysis/20200928_longi_method/GenusTYPE1_210425P2_ZIBRFIX_n300/"
DVD<-1
HeatSTR<-paste0("GenusTYPE1_210425P2_ZIBRFIX_n300")
library(stringr)

SingleOTU<-3

file_list<-list.files(infol)
file_list2_bf<-file_list[str_detect(file_list,HeatSTR)]

file_list2<-file_list2_bf[order(as.numeric(str_extract(file_list2_bf,"\\d+(?=.rds)")))]
rdrd<-lapply(file.path(infol,file_list2),readRDS)

method_choose<-1:8

Method_data_seq<-c("mTMAT_M","mTMAT_IM","GLMM-MiRKAT","FZINBMM","ZIBR","LMM","ZIBR_logi","ZIBR_beta")[method_choose]
Method_output_seq<-c("mTMAT_IM","mTMAT_M","GLMM-MiRKAT","FZINBMM","ZIBR","LMM","ZIBR_logi","ZIBR_beta")[method_choose]

# infol2<-"/home2/kjkim/analysis/20200928_longi_method/GENRA_diffN_type1"
# HeatSTR2<-paste0("_GENRA_210406_02D1")
# 
# file2_list<-list.files(infol2)
# file2_list2_bf<-file2_list[str_detect(file2_list,HeatSTR2)]
# 
# file2_list2<-file2_list2_bf[order(as.numeric(str_extract(file2_list2_bf,"\\d+(?=.rds)")))]
# 
# rd2rd<-lapply(file.path(infol2,file2_list2),readRDS)


rdrd2<-lapply(rdrd,function(data){
  data[[3]]
})


ttable<-lapply(rdrd,function(data){
  data[[1]]
})


labels<-lapply(rdrd,function(data){
  data[[2]]
})

cutoff<-0.1
cutoff_list<-c(0.1,0.05,0.01,0.005)
TYPE1_table<-lapply(cutoff_list,function(cutoff){
  
  partsTOPLOT<-lapply(1:length(rdrd2),function(ind_data){
    
    dfTMP<-data.frame(Taxa=labels[[ind_data]][[3]][ttable[[ind_data]][,4]],beta=labels[[ind_data]][[4]][ttable[[ind_data]][,5]],CCratio=labels[[ind_data]][[1]][ttable[[ind_data]][,2]],N=labels[[ind_data]][[2]][ttable[[ind_data]][,3]])
    Values<-t(sapply(rdrd2[[ind_data]],function(data){
      apply(data,2,function(value){
        sum(value<cutoff,na.rm=T)/sum(!is.na(value))
      })
    }))
    Values2<-c()
    for( i in 1:dim(Values)[2]){
      Values2<-c(Values2,Values[,i])
    }
    Method<-rep(Method_data_seq,each=dim(Values)[1])
    
    DFDF<-cbind(Values2,dfTMP,Method)
    
  })
  
  partsTOPLOT2_bf2<-do.call("rbind",partsTOPLOT)
  partsTOPLOT2_bf<-partsTOPLOT2_bf2
  # partsTOPLOT2_bf<-partsTOPLOT2_bf2[partsTOPLOT2_bf2$N==50 &partsTOPLOT2_bf2$CCratio==0.25,]
  
  
  
  infofol<-"~/analysis/20200928_longi_method/info/"
  infile_nm<-file.path(infofol,"GENUS_INFO.csv")
  GNR_info<-readRDS(infile_nm)
  
  
  if(SingleOTU==1){
    GENERA_TARGET<-GNR_info[[1]][which(GNR_info[[7]]==1)]
    partsTOPLOT2_bf$Taxa<-as.character(partsTOPLOT2_bf$Taxa)
    partsTOPLOT2<-partsTOPLOT2_bf[partsTOPLOT2_bf$Taxa%in%GENERA_TARGET & partsTOPLOT2_bf$Method!="GLMM-MiRKAT",]
    
    methods_sequence<-setdiff(Method_output_seq,"GLMM-MiRKAT")
    ncol=2
  }else if(SingleOTU==2){
    GENERA_TARGET<-GNR_info[[1]][which(GNR_info[[7]]!=1)]
    partsTOPLOT2_bf$Taxa<-as.character(partsTOPLOT2_bf$Taxa)
    partsTOPLOT2<-partsTOPLOT2_bf[partsTOPLOT2_bf$Taxa%in%GENERA_TARGET,]
    
    methods_sequence<-Method_output_seq
    ncol=3
  }else{
    partsTOPLOT2<-partsTOPLOT2_bf
    methods_sequence<-Method_output_seq
  }
  
  
  Values3<-split(partsTOPLOT2$Values2,paste(partsTOPLOT2$beta,partsTOPLOT2$Method,sep="__"))
  Values4<-sapply(Values3,mean,na.rm=T)
  library(stringr)
  DFPLOT<-data.frame(cbind(Values4,t(sapply(str_split(names(Values4),"__"),"[",1:2))))
  names(DFPLOT)<-c("Value","beta","Method")
  
  library(ggplot2)
  ncol=3
  methods_sequence<-Method_output_seq
  DFPLOT$Method<-factor(DFPLOT$Method, levels = methods_sequence)
  
  
  
  labels_out_bf<-as.character(levels(DFPLOT$Method))
  labels_out<-labels_out_bf
  labels_out[labels_out_bf=="mTMAT_M"]<-expression("mTMAT" ["M"])
  labels_out[labels_out_bf=="mTMAT_IM"]<-expression("mTMAT" ["IM"])
  
  
  ind_methods<-1:length(methods_sequence)
  library(scales)
  ggplotColours <- function(n = 6, h = c(0, 360) + 15){
    if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
    hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
  }
  geom.colors <- ggplotColours(n=length(methods_sequence))
  geom.colors_picked<-geom.colors[ind_methods]
  
  DFPLOT$Value<-as.numeric(as.character(DFPLOT$Value))
  DFPLOT$beta<-as.numeric(as.character(DFPLOT$beta))


  ind_1<-sapply(GNR_info[[2]],mean)<=0.2
  ind_2<-sapply(GNR_info[[2]],mean)<=0.5&sapply(GNR_info[[2]],mean)>0.2
  ind_3<-sapply(GNR_info[[2]],mean)<=0.8&sapply(GNR_info[[2]],mean)>0.5
  ind_4<-sapply(GNR_info[[2]],mean)>0.8

  
  
  DF22<-cbind(GNR_info[[1]],"a")
  DF22[ind_1 ,2]<-"<=20%"
  DF22[ind_2 ,2]<-"20-50%"
  DF22[ind_3 ,2]<-"50-80%"
  DF22[ind_4 ,2]<-">80%"
  
  
  str_SG<-DF22[match(as.character(partsTOPLOT2$Taxa),DF22[,1]),2]
  
  str_SG[str_SG%in% c("50-80%",">80%")]<-">50%"
  
  
  # Values3_2<-split(partsTOPLOT2$Values2,paste(partsTOPLOT2$beta,partsTOPLOT2$Method,str_SG,sep="__"))
  Values3_2<-split(partsTOPLOT2$Values2,paste(partsTOPLOT2$beta,partsTOPLOT2$Method,partsTOPLOT2$CCratio,partsTOPLOT2$N,sep="__"))
  Values4_2<-sapply(Values3_2,mean,na.rm=T)
  
  
  
  library(stringr)
  DFPLOT_2<-data.frame(cbind(Values4_2,t(sapply(str_split(names(Values4_2),"__"),"[",1:4))))
  names(DFPLOT_2)<-c("Value","beta","Method","CCratio","N")
  
  library(ggplot2)
  ncol=4
  methods_sequence<-Method_output_seq
  DFPLOT_2$Method<-factor(DFPLOT_2$Method, levels = methods_sequence)
  
  
  
  labels_out_bf<-as.character(levels(DFPLOT_2$Method))
  labels_out<-labels_out_bf
  labels_out[labels_out_bf=="mTMAT_M"]<-expression("mTMAT" ["M"])
  labels_out[labels_out_bf=="mTMAT_IM"]<-expression("mTMAT" ["IM"])
  
  
  ind_methods<-1:length(methods_sequence)
  library(scales)
  ggplotColours <- function(n = 6, h = c(0, 360) + 15){
    if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
    hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
  }
  geom.colors <- ggplotColours(n=length(methods_sequence))
  geom.colors_picked<-geom.colors[ind_methods]
  
  DFPLOT_2$Value<-as.numeric(as.character(DFPLOT_2$Value))
  DFPLOT_2$beta<-as.numeric(as.character(DFPLOT_2$beta))

  DF_type1<-cbind(DFPLOT_2[DFPLOT_2$beta=="0",],cutoff)
  
})

# new_group<-list(GNR_info[[3]][[1]],GNR_info[[3]][[2]],do.call("c",GNR_info[[3]][3:4]))
# 
# sapply(lapply(new_group,function(ind){
#   GNR_info[[7]][ind]
# }),mean)
# tapply(GNR_info[[3]],,mean)

TYPE1_table2<-do.call("rbind",TYPE1_table)
nm<-paste0("~/TYPE1_",HeatSTR,"_",cutoff,".csv")
write.csv(TYPE1_table2,nm)
system(paste0("scp ",nm," ng:~"))