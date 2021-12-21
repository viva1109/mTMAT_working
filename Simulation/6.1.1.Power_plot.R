Rank<-"Genus"
infol<-"~/analysis/20200928_longi_method/GenusPower_210502_RESULTlogisticTRUE_ZOOMFALSE_logistic_prop_corr_TMATFIX_FILTERoriginal_n100/"
DVD<-1
HeatSTR<-paste0("GenusPower_210502_RESULTlogisticTRUE_ZOOMFALSE_logistic_prop_corr_TMATFIX_FILTERoriginal_n100")
infofol<-"~/analysis/20200928_longi_method/info/"
SingleOTU<-3
library(stringr)

file_list<-list.files(infol)
file_list2_bf<-file_list[str_detect(file_list,HeatSTR)]

file_list2<-file_list2_bf[order(as.numeric(str_extract(file_list2_bf,"\\d+(?=.rds)")))]
rdrd<-lapply(file.path(infol,file_list2),readRDS)

method_choose<-1:5

Method_target<-c("mTMAT_15INDEP","mTMAT_res16INDEP","GLMM_MiRKAT","FZINBMM","ZIBR")
# Method_data_seq<-c("mTMAT_M","mTMAT_IM","GLMM-MiRKAT","FZINBMM","ZIBR","LMM","ZIBR_logi","ZIBR_beta")[method_choose]
# Method_output_seq<-c("mTMAT_IM","mTMAT_M","GLMM-MiRKAT","FZINBMM","ZIBR","LMM","ZIBR_logi","ZIBR_beta")[method_choose]


rdrd2<-lapply(rdrd,function(data){
  data[[3]]
})


ttable<-lapply(rdrd,function(data){
  data[[1]]
})


labels<-lapply(rdrd,function(data){
  data[[2]]
})

standars_critical_val<-readRDS(file.path(infofol,"criticalValues_For_zinb_NBMM.rds"))
###########perform evaluation check
cutoff<-0.1
partsTOPLOT<-lapply(1:length(rdrd2),function(ind_data){
  
  dfTMP<-data.frame(Taxa=labels[[ind_data]][[3]][ttable[[ind_data]][,4]],beta=labels[[ind_data]][[4]][ttable[[ind_data]][,5]],p=c(0,0.5,0.9)[ttable[[ind_data]][,6]])
  Values<-t(sapply(1:length(rdrd2[[ind_data]]),function(ind_data_sub){
    data<-rdrd2[[ind_data]][[ind_data_sub]]
    
    if(is.na(substr(data,1,5)[1])){
      
      ind_separate<-which(colnames(data) %in% c("FZINBMM","NBMM"))
      outval<-c(
        sum(data[,"FZINBMM"]<  standars_critical_val[1,1],na.rm=T)/sum(!is.na(data[,"FZINBMM"])),
        sum(data[,"NBMM"]<  standars_critical_val[1,2],na.rm=T)/sum(!is.na(data[,"NBMM"]))
        ,apply(data[,-ind_separate],2,function(value2){
          sum(value2<cutoff,na.rm=T)/sum(!is.na(value2))
        }))
      
      # outval<-apply(data,2,function(value){
      #   sum(value<cutoff,na.rm=T)/sum(!is.na(value))
      # }) 
    }else{
      if(substr(data,1,5)[1]=="Error"){
        outval<-rep(NA,dim(rdrd2[[ind_data]][[1]])[2])
      }else{
        
        ind_separate<-which(colnames(data) %in% c("FZINBMM","NBMM"))
        outval<-c(
          sum(data[,"FZINBMM"]<  standars_critical_val[1,1],na.rm=T)/sum(!is.na(data[,"FZINBMM"])),
          sum(data[,"NBMM"]<  standars_critical_val[1,2],na.rm=T)/sum(!is.na(data[,"NBMM"]))
          ,apply(data[,-ind_separate],2,function(value2){
            sum(value2<cutoff,na.rm=T)/sum(!is.na(value2))
          }))
        # outval<-apply(data,2,function(value){
        #   sum(value<cutoff,na.rm=T)/sum(!is.na(value))
        # })  
      }
    }
    # print(paste0(ind_data,"_",ind_data_sub))
    names(outval)[1:2]<-c("FZINBMM","NBMM")
    return(outval)
  }))
  Values2<-c()
  golist<-match(Method_target,colnames(Values))
  
  for( i in golist){
    Values2<-c(Values2,Values[,i])
  }
  Method<-rep(Method_target,each=dim(Values)[1])
  print(paste0(ind_data))
  DFDF<-cbind(Values2,dfTMP,Method)
  
})
# lapply(partsTOPLOT,dim)
partsTOPLOT2<-do.call("rbind",partsTOPLOT)




if(Rank=="Family"){
  infile_nm<-file.path(infofol,paste0(Rank,"_INFO.csv"))  
}else{
  infile_nm<-file.path(infofol,paste0(Rank,"_INFO.csv"))
}



GNR_info<-readRDS(infile_nm)

DF22<-cbind(GNR_info[[1]],"a")


DF22[GNR_info[[3]][[1]] ,2]<-"<=20%"
DF22[GNR_info[[3]][[2]] ,2]<-"20-50%"
DF22[GNR_info[[3]][[3]] ,2]<-"50-80%"
DF22[GNR_info[[3]][[4]] ,2]<-">80%"


DF23<-cbind(GNR_info[[1]],"a")
DF23[GNR_info[[4]][[1]] ,2]<-"<=1"
DF23[GNR_info[[4]][[2]] ,2]<-"2"
DF23[GNR_info[[4]][[3]] ,2]<-"3-4"
DF23[GNR_info[[4]][[4]] ,2]<-">4"


str_SG<-DF22[match(as.character(partsTOPLOT2$Taxa),DF22[,1]),2]
str_OG<-DF23[match(as.character(partsTOPLOT2$Taxa),DF23[,1]),2]

str_SG[str_SG%in% c("50-80%",">80%")]<-">50%"



infofol<-"~/analysis/20200928_longi_method/info/"
if(Rank=="Family"){
  infile_nm<-file.path(infofol,paste0(Rank,"_INFO_0.001.csv"))  
}else{
  infile_nm<-file.path(infofol,paste0(Rank,"_INFO_0.001.csv"))
}



GNR_info<-readRDS(infile_nm)

DF22<-cbind(GNR_info[[1]],"a")


DF22[GNR_info[[3]][[1]] ,2]<-"<=20%"
DF22[GNR_info[[3]][[2]] ,2]<-"20-50%"
DF22[GNR_info[[3]][[3]] ,2]<-"50-80%"
DF22[GNR_info[[3]][[4]] ,2]<-">80%"


DF23<-cbind(GNR_info[[1]],"a")
DF23[GNR_info[[4]][[1]] ,2]<-"<=1"
DF23[GNR_info[[4]][[2]] ,2]<-"2"
DF23[GNR_info[[4]][[3]] ,2]<-"3-4"
DF23[GNR_info[[4]][[4]] ,2]<-">4"


str_SG<-DF22[match(as.character(partsTOPLOT2$Taxa),DF22[,1]),2]
str_OG<-DF23[match(as.character(partsTOPLOT2$Taxa),DF23[,1]),2]

str_SG[str_SG%in% c("50-80%",">80%")]<-">50%"


# cbind(GNR_info[[1]],GNR_info[[3]])
methods_sequence<-Method_output_seq
Method_target
# DFPLOT$Method<-as.character(DFPLOT$Method)
# DFPLOT$Method<-factor(DFPLOT$Method, levels = methods_sequence)
if(SingleOTU==1){
  GENERA_TARGET<-GNR_info[[1]][which(GNR_info[[3]]==1)]
  partsTOPLOT2$Taxa<-as.character(partsTOPLOT2$Taxa)
  partsTOPLOT3<-partsTOPLOT2[partsTOPLOT2$Taxa%in%GENERA_TARGET & partsTOPLOT2$Method!="GLMM-MiRKAT",]
  
  ind_methods<-c(1,2,4,5)
  # methods_sequence<-c("mTMAT_IM","mTMAT_M","FZINBMM","ZIBR","LMM")
  ncol=2
}else if(SingleOTU==2){
  GENERA_TARGET<-GNR_info[[1]][which(GNR_info[[3]]!=1)]
  partsTOPLOT2$Taxa<-as.character(partsTOPLOT2$Taxa)
  partsTOPLOT3<-partsTOPLOT2[partsTOPLOT2$Taxa%in%GENERA_TARGET,]
  ind_methods<-1:length(Method_target)
  ncol=3
}else{
  GENERA_TARGET<-GNR_info[[1]]
  partsTOPLOT2$Taxa<-as.character(partsTOPLOT2$Taxa)
  partsTOPLOT3<-partsTOPLOT2[partsTOPLOT2$Taxa%in%GENERA_TARGET & partsTOPLOT2$Method!="GLMM-MiRKAT",]
  
  ind_methods<-c(1,2,4,5)
  ncol=2
}



Values3<-split(partsTOPLOT3$Values2,paste(partsTOPLOT3$beta,partsTOPLOT3$Method,partsTOPLOT3$p,sep="__"))
Values4<-sapply(Values3,mean,na.rm=T)
library(stringr)
DFPLOT<-data.frame(cbind(Values4,t(sapply(str_split(names(Values4),"__"),"[",1:3))))
names(DFPLOT)<-c("Value","beta","Method","p")

library(ggplot2)





Nw1<-sum(GNR_info[[7]]==1)
Nw2<-sum(GNR_info[[7]]!=1)

library(scales)
ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}
geom.colors <- ggplotColours(n=length(Method_target))
geom.colors_picked<-geom.colors[ind_methods]

DFPLOT$Value<-as.numeric(as.character(DFPLOT$Value))
DFPLOT$beta<-as.numeric(as.character(DFPLOT$beta))
DFPLOT2<-DFPLOT[DFPLOT$beta!="0",]

DFPLOT2$Method<-as.character(DFPLOT2$Method)
DFPLOT2$Method<-factor(DFPLOT2$Method, levels = Method_target[ind_methods])

labels_out_bf<-as.character(levels(DFPLOT2$Method))
labels_out<-labels_out_bf
labels_out[labels_out_bf=="mTMAT_15INDEP"]<-expression("mTMAT" ["M"])
labels_out[labels_out_bf=="mTMAT_res16INDEP"]<-expression("mTMAT" ["IM"])
# labels_out[labels_out_bf=="mTMAT_IM"]<-expression("mTMAT" ["IM"])
DFPLOT2_4<-DFPLOT2
p_list<-c(0,0.5,0.9)


DFPLOT2_4$Method

DFPLOT2_5<-DFPLOT2_4
DFPLOT2_5$Method<-factor(DFPLOT2_5$Method, levels = Method_target[ind_methods])
for (z in 1:3){
  p_value_tmp<-p_list[z]
  # DFPLOT3<-DFPLOT2[DFPLOT2$p==p_value_tmp,]
  DFPLOT3<-DFPLOT2_5[DFPLOT2_4$p==p_value_tmp,]
  ggplot(DFPLOT3, aes(x=beta, y=Value, color=Method)) + 
    scale_color_manual(values=geom.colors_picked,labels=labels_out)+
    scale_linetype_manual(values=1:length(Method_target),labels=labels_out)+
    geom_line(aes(linetype=Method), size=1,alpha=0.7)+
    scale_shape_manual(values=c(16,3,5,8,7,15,17,18,19,20)[ind_methods],labels=labels_out) +
    theme_bw()+
    geom_point(aes(shape=Method),alpha=0.7, size=4)+
    guides(linetype = guide_legend("",ncol=ncol,byrow=T),shape = guide_legend("",ncol=ncol,byrow=T),colour = guide_legend("",ncol=ncol,byrow=T),alpha = guide_legend("",ncol=ncol,byrow=T) ,size = guide_legend("none"))+ylab("Power estimates")+xlab("Beta")+theme(legend.position="bottom")+
    ylim(c(0,1))
  
  
  nm_tmp<-paste0(HeatSTR,"_",p_value_tmp,"_Single",SingleOTU,"_",cutoff,".png")
  ggsave(nm_tmp,height=9,width=5,dpi=300)
  system(paste0("scp ",nm_tmp," ng:~"))
  
}



