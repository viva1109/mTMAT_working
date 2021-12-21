Rank<-"Genus"
ddo<-"split"

infol<-"~/analysis/20200928_longi_method/PowerGenus_mTMAT_CORE_GLMM-MiRKAT_ZINB_LMMmissing0.1_210503_FZoomT_n100/"
DVD<-1
HeatSTR<-paste0("PowerGenus_mTMAT_CORE_GLMM-MiRKAT_ZINB_LMMmissing0.1_210503_FZoomT_n100")
infofol<-"~/analysis/20200928_longi_method/info/"
SingleOTU<-3
library(stringr)
forceAlpha<-T
file_list<-list.files(infol)
file_list2_bf<-file_list[str_detect(file_list,HeatSTR)]

file_list2<-file_list2_bf[order(as.numeric(str_extract(file_list2_bf,"\\d+(?=.rds)")))]
rdrd<-lapply(file.path(infol,file_list2),readRDS)

method_choose<-1:7
# c("mTMAT15","mTMAT16","GLMM_MiRKAT","FZINBMM","NBMM","LMM_arcsin","LMM_log","mTMAT_15CS","mTMAT_15AR1","mTMAT_15UN","mTMAT_16CS","mTMAT_16AR1","mTMAT_16UN","TMAT_res15","TMAT_res16","MiSPU_res","MiRKAT_res","Wilcox_res")
# Method_target<-c("mTMAT_15INDEP","mTMAT_res16INDEP","GLMM_MiRKAT","FZINBMM","ZIBR")

Method_target<-c("mTMAT15","mTMAT16","GLMM_MiRKAT","FZINBMM","NBMM","LMM_arcsin","LMM_log")

Method_toplot<-c("mTMAT16","mTMAT15","GLMM_MiRKAT","FZINBMM","LMM_arcsin","LMM_log")

# Method_toplot<-c("mTMAT16","mTMAT15","TMAT_res15","TMAT_res16","MiSPU_res","MiRKAT_res","Wilcox_res")    

# Method_data_seq<-c("mTMAT_M","mTMAT_IM","GLMM-MiRKAT","FZINBMM","ZIBR","LMM","ZIBR_logi","ZIBR_beta")[method_choose]
# Method_output_seq<-c("mTMAT_IM","mTMAT_M","GLMM-MiRKAT","FZINBMM","ZIBR","LMM","ZIBR_logi","ZIBR_beta")[method_choose]


rdrd2<-lapply(rdrd,function(data){
  dd<-data[[3]]
  # lapply(dd,function(data){
  #   colnames(data)<-Method_target
  #   return(data)
  # })
})


ttable<-lapply(rdrd,function(data){
  data[[1]]
})


labels<-lapply(rdrd,function(data){
  data[[2]]
})


# inds_togo_list<-lapply(ttable,function(data){
#   which(data[,2]==2 & data[,3]==2)
# })
# rdrd2_chosen<-lapply(1:length(rdrd2),function(ind_data){
#   rdrd2[[ind_data]][inds_togo_list[[ind_data]]]
# })
# rdrd2_chosen2<-lapply(rdrd2_chosen,function(data1){
#   do.call("rbind",data1)
# })
# gotH0dist<-do.call("rbind",rdrd2_chosen2)
# ALLcriticalValues<-t(sapply(c(0.1,0.05,0.01,0.005),function(cutoff){apply(gotH0dist,2,quantile,cutoff,na.rm=T)} ))
# colnames(ALLcriticalValues)<-Method_target
# saveRDS(ALLcriticalValues,file.path(infofol,"ALLcriticalValues.rds"))


ALLcriticalValues<-readRDS(file.path(infofol,"ALLcriticalValues.rds"))

###########perform evaluation check
cutoff<-0.1


partsTOPLOT<-lapply(1:length(rdrd2),function(ind_data){
  
  dfTMP<-data.frame(Taxa=labels[[ind_data]][[3]][ttable[[ind_data]][,4]],beta=labels[[ind_data]][[4]][ttable[[ind_data]][,5]],p=c(0,0.5,0.9)[ttable[[ind_data]][,6]])
  Values<-t(sapply(1:length(rdrd2[[ind_data]]),function(ind_data_sub){
    data<-apply(rdrd2[[ind_data]][[ind_data_sub]],2,as.numeric)
    
    if(is.na(substr(data,1,5)[1])){
      
      outval<-   sapply(Method_toplot,function(nmsMD){
        if(forceAlpha){
          sum(data[,nmsMD]<  ALLcriticalValues[1,nmsMD],na.rm=T)/sum(!is.na(data[,nmsMD]))
        }else{
          sum(data[,nmsMD]<  cutoff,na.rm=T)/sum(!is.na(data[,nmsMD]))
        }
        
      })
      
      # outval<-apply(data,2,function(value){
      #   sum(value<cutoff,na.rm=T)/sum(!is.na(value))
      # }) 
    }else{
      if(substr(data,1,5)[1]=="Error"){
        outval<-rep(NA,dim(rdrd2[[ind_data]][[1]])[2])
      }else{
        outval<-   sapply(Method_toplot,function(nmsMD){
          if(forceAlpha){
            sum(data[,nmsMD]<  ALLcriticalValues[1,nmsMD],na.rm=T)/sum(!is.na(data[,nmsMD]))
          }else{
            sum(data[,nmsMD]<  cutoff,na.rm=T)/sum(!is.na(data[,nmsMD]))
          }
        })
        # outval<-apply(data,2,function(value){
        #   sum(value<cutoff,na.rm=T)/sum(!is.na(value))
        # })  
      }
    }
    # print(paste0(ind_data,"_",ind_data_sub))
    return(outval)
  }))
  Values2<-c()
  golist<-match(Method_toplot,colnames(Values))
  
  for( i in golist){
    Values2<-c(Values2,Values[,i])
  }
  Method<-rep(Method_toplot,each=dim(Values)[1])
  print(paste0(ind_data))
  DFDF<-cbind(Values2,dfTMP,Method)
  
})
# lapply(partsTOPLOT,dim)
partsTOPLOT2<-do.call("rbind",partsTOPLOT)

partsTOPLOT2_list3<-partsTOPLOT2
partsTOPLOT2_list2<-partsTOPLOT2
partsTOPLOT2_list1<-partsTOPLOT2

partsTOPLOT2_CBD<-rbind(cbind(partsTOPLOT2_list1,Mrate="Balanced"),
cbind(partsTOPLOT2_list2,Mrate="Missing_10%"),
cbind(partsTOPLOT2_list3,Mrate="Missing_20%"))

partsTOPLOT2_CBD2<-partsTOPLOT2_CBD[partsTOPLOT2_CBD$beta==0.04  &partsTOPLOT2_CBD$p==0.5,]

###############
library(ggplot2)
bb<-0.02
pp<-0.5
partsTOPLOT2<-do.call("rbind",partsTOPLOT)
partsTOPLOT2_togo<-partsTOPLOT2[partsTOPLOT2$beta==bb  &partsTOPLOT2$p==pp & partsTOPLOT2$Method%in%c("mTMAT16","FZINBMM","GLMM_MiRKAT","LMM_arcsin"),]

Tax_orderTable<-partsTOPLOT2_togo[partsTOPLOT2_togo$Method=="mTMAT16",]
Tax_order<-Tax_orderTable$Taxa[order(Tax_orderTable$Values2)]

labels_out_bf<-as.character(levels(partsTOPLOT2_togo$Method))
labels_out_bf[labels_out_bf=="mTMAT16"]<-expression("mTMAT" ["IM"])

partsTOPLOT2_togo$Taxa2<-factor(partsTOPLOT2_togo$Taxa,level=Tax_order)
partsTOPLOT2_togo$Method
partsTOPLOT2_togo$Method<-as.character(partsTOPLOT2_togo$Method)
partsTOPLOT2_togo$Method[partsTOPLOT2_togo$Method=="mTMAT16"]<-"mTMAT_IM"

partsTOPLOT2_togo$Method<-factor(partsTOPLOT2_togo$Method,levels=c("mTMAT_IM","GLMM_MiRKAT","FZINBMM","LMM_arcsin"))
# partsTOPLOT2_togo$Method3<-factor(partsTOPLOT2_togo$Method2,levels=labels_out_bf)

ggplot(partsTOPLOT2_togo, aes(x=Taxa2, y=Values2)) +
  geom_point(stat='identity', aes(col=Method), size=3) + geom_hline(yintercept=0.05, linetype='dashed', color='black', size=0.5)+
  # scale_color_manual(name="Mileage (deviation)",
  #                    labels = c("Above Average", "Below Average"),
  #                    values = c("above"="#00ba38", "below"="#0b8fd3")) +
  # geom_segment(aes(y = 0,
  #                  x = CarBrand,
  #                  yend = mpg_z_score,
  #                  xend = CarBrand),
  #              color = "black") +
  # geom_text(color="white", size=2) +
  labs(title="Adjusted")+ylim(0, 1)+xlab("Taxonomy") + ylab("Power estimates")+coord_flip() +theme_bw() +theme(panel.grid.major = element_blank(), panel.grid.minor =element_blank())

nmnm<-"lolpop2_AdjustLine6.png"
ggsave(nmnm,height=13,width=9,dpi=500)
system(paste0("scp ",nmnm," ng:~"))
###############


# partsTOPLOT2[is.nan(partsTOPLOT2$Value)]

# partsTOPLOT2$Taxa
if(Rank=="Family"){
  infile_nm<-file.path(infofol,paste0(Rank,"_INFO.csv"))  
}else{
  infile_nm<-file.path(infofol,paste0(Rank,"_INFO.csv"))
}

# 
# 
# GNR_info<-readRDS(infile_nm)
# 
# DF22<-cbind(GNR_info[[1]],"a")
# 
# 
# DF22[GNR_info[[3]][[1]] ,2]<-"<=20%"
# DF22[GNR_info[[3]][[2]] ,2]<-"20-50%"
# DF22[GNR_info[[3]][[3]] ,2]<-"50-80%"
# DF22[GNR_info[[3]][[4]] ,2]<-">80%"
# 
# 
# DF23<-cbind(GNR_info[[1]],"a")
# DF23[GNR_info[[4]][[1]] ,2]<-"<=1"
# DF23[GNR_info[[4]][[2]] ,2]<-"2"
# DF23[GNR_info[[4]][[3]] ,2]<-"3-4"
# DF23[GNR_info[[4]][[4]] ,2]<-">4"
# 
# 
# str_SG<-DF22[match(as.character(partsTOPLOT2$Taxa),DF22[,1]),2]
# str_OG<-DF23[match(as.character(partsTOPLOT2$Taxa),DF23[,1]),2]
# 
# str_SG[str_SG%in% c("50-80%",">80%")]<-">50%"
# 







infofol<-"~/analysis/20200928_longi_method/info/"
if(Rank=="Family"){
  infile_nm<-file.path(infofol,paste0(Rank,"_INFO_0.001.csv"))  
}else{
  infile_nm<-file.path(infofol,paste0(Rank,"_INFO_0.001.csv"))
}
GNR_info<-readRDS(infile_nm)


# cbind(GNR_info[[1]],GNR_info[[3]])
methods_sequence<-Method_toplot
partsTOPLOT2<-partsTOPLOT2_CBD2
# DFPLOT$Method<-as.character(DFPLOT$Method)
# DFPLOT$Method<-factor(DFPLOT$Method, levels = methods_sequence)
if(SingleOTU==1){
  GENERA_TARGET<-GNR_info[[1]][which(GNR_info[[3]]==1)]
  partsTOPLOT2$Taxa<-as.character(partsTOPLOT2$Taxa)
  partsTOPLOT3<-partsTOPLOT2[partsTOPLOT2$Taxa%in%GENERA_TARGET & partsTOPLOT2$Method!="GLMM_MiRKAT",]
  
  ind_methods<-c(1,2,4,5,6)
  # methods_sequence<-c("mTMAT_IM","mTMAT_M","FZINBMM","ZIBR","LMM")
  ncol=2
}else if(SingleOTU==2){
  GENERA_TARGET<-GNR_info[[1]][which(GNR_info[[3]]!=1)]
  partsTOPLOT2$Taxa<-as.character(partsTOPLOT2$Taxa)
  partsTOPLOT3<-partsTOPLOT2[partsTOPLOT2$Taxa%in%GENERA_TARGET,]
  ind_methods<-c(1,2,3,4,5,6)
  ncol=3
}else{
  GENERA_TARGET<-GNR_info[[1]]
  partsTOPLOT2$Taxa<-as.character(partsTOPLOT2$Taxa)
  partsTOPLOT3<-partsTOPLOT2[partsTOPLOT2$Taxa%in%GENERA_TARGET & partsTOPLOT2$Method!="GLMM_MiRKAT",]
  
  ind_methods<-c(1,2,3,4,5,6)
  ncol=2
}

partsTOPLOT3



DF22<-cbind(GNR_info[[1]],"a")


DF22[sapply(GNR_info[[2]],mean)<=0.2 ,2]<-"<=20%"
DF22[sapply(GNR_info[[2]],mean)<=0.5 & sapply(GNR_info[[2]],mean)>0.2 ,2]<-"20-50%"
DF22[sapply(GNR_info[[2]],mean)<=0.8 & sapply(GNR_info[[2]],mean)>0.5 ,2]<-"50-80%"
DF22[sapply(GNR_info[[2]],mean)>0.8 ,2]<-">80%"


DF23<-cbind(GNR_info[[1]],"a")
DF23[GNR_info[[3]]==1 ,2]<-"=1"
DF23[GNR_info[[3]]>=2 &GNR_info[[3]]<=5 ,2]<-"2-5"
DF23[GNR_info[[3]]<=15 & GNR_info[[3]]>=6,2]<-"6-15"
DF23[GNR_info[[3]]>15 ,2]<-">15"

str_SG<-DF22[match(as.character(partsTOPLOT3$Taxa),DF22[,1]),2]
str_OG<-DF23[match(as.character(partsTOPLOT3$Taxa),DF23[,1]),2]

str_SG[str_SG%in% c("50-80%",">80%")]<-">50%"




# OTU_INFO<-readRDS(file.path(infofol,"OTU_INFO_full_Family.rds"))
# 
# TaxG1<-OTU_INFO[[2]][match(partsTOPLOT3$Taxa,OTU_INFO[[3]])]
# unique(sapply(str_split(TaxG1,"; "),"[",3))

if(ddo=="split"){
  if(Rank=="Family"){
    Values3<-split(partsTOPLOT3$Values2,paste(partsTOPLOT3$beta,partsTOPLOT3$Method,partsTOPLOT3$p,str_OG,sep="__"))
  }else{
    Values3<-split(partsTOPLOT3$Values2,paste(partsTOPLOT3$beta,partsTOPLOT3$Method,partsTOPLOT3$p,str_SG,sep="__"))
  }
  
}else{
  Values3<-split(partsTOPLOT3$Values2,paste(partsTOPLOT3$beta,partsTOPLOT3$Method,partsTOPLOT3$p,sep="__"))
}
Values3<-split(partsTOPLOT3$Values2,paste0(partsTOPLOT3$Method,partsTOPLOT3$Mrate,sep="__"))
Values4<-sapply(Values3,mean,na.rm=T)

if(ddo=="split"){
  DFPLOT<-data.frame(cbind(Values4,t(sapply(str_split(names(Values4),"__"),"[",1:4))))
  if(Rank=="Family"){
    names(DFPLOT)<-c("Value","beta","Method","p","Leafs")
  }else{
    names(DFPLOT)<-c("Value","beta","Method","p","Sparsity")
  }
  
}else{
  DFPLOT<-data.frame(cbind(Values4,t(sapply(str_split(names(Values4),"__"),"[",1:3))))
  names(DFPLOT)<-c("Value","beta","Method","p")
}
DFPLOT<-data.frame(cbind(Values4,t(sapply(str_split(names(Values4),"__"),"[",1:3))))
names(DFPLOT)<-c("Value","Method","Mrate")
DFPLOT$Mrate<-str_sub(DFPLOT$Method,-3,-1)
str(DFPLOT$Method,1,3)
library(stringr)
DFPLOT$Mrate<-as.character(DFPLOT$Mrate)
DFPLOT[DFPLOT$Mrate=="ced"]<-"Balanced"
gsub(DFPLOT$Mrate,"",DFPLOT$Method)
str(DFPLOT)
DFPLOT$Method<-as.character(DFPLOT$Method)
DFPLOT$Method2<-gsub("Missing .+$","",DFPLOT$Method)
DFPLOT$Method3<-gsub("Balanced","",DFPLOT$Method2)
DFPLOT$Mrate2<-gsub("ced","Balanced",DFPLOT$Mrate)
DFPLOT$Method3<-gsub("Missing.+$","",DFPLOT$Method3)
library(ggplot2)





Nw1<-sum(GNR_info[[7]]==1)
Nw2<-sum(GNR_info[[7]]!=1)

library(scales)
ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}
geom.colors <- ggplotColours(n=length(Method_toplot))
geom.colors_picked<-geom.colors[ind_methods]




DFPLOT$Value<-as.numeric(as.character(DFPLOT$Value))
DFPLOT$beta<-as.numeric(as.character(DFPLOT$beta))
DFPLOT2<-DFPLOT[DFPLOT$beta!="0",]

DFPLOT2$Method<-as.character(DFPLOT2$Method)
DFPLOT2$Method<-factor(DFPLOT2$Method, levels = Method_toplot[ind_methods])

labels_out_bf<-as.character(levels(DFPLOT2$Method))
labels_out<-labels_out_bf
labels_out[labels_out_bf=="mTMAT15"]<-expression("mTMAT" ["M"])
labels_out[labels_out_bf=="mTMAT16"]<-expression("mTMAT" ["IM"])
labels_out[labels_out_bf=="TMAT_res15"]<-expression("TMAT" ["M"])
labels_out[labels_out_bf=="TMAT_res16"]<-expression("TMAT" ["IM"])

# labels_out[labels_out_bf=="mTMAT_IM"]<-expression("mTMAT" ["IM"])
DFPLOT2_4<-DFPLOT2
p_list<-c(0,0.5,0.9)


DFPLOT2_4$Method

DFPLOT2_5<-DFPLOT2_4
DFPLOT2_5$Method<-factor(DFPLOT2_5$Method, levels = Method_toplot[ind_methods])
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
  
  
  nm_tmp<-paste0(HeatSTR,"_ForceAl",forceAlpha,"_Single",SingleOTU,"_",p_value_tmp*100,"_",cutoff,".png")
  ggsave(nm_tmp,height=9,width=5,dpi=300)
  system(paste0("scp ",nm_tmp," ng:~"))
  
}

DFPLOT2_6<-DFPLOT2_5[DFPLOT2_5$beta==0.02 & DFPLOT2_5$p==0.5 & complete.cases(DFPLOT2_5),]

byr<-T
if(Rank=="Family"){
  xlabel<-"Number of leaf nodes"
}else{
  xlabel<-"Missing rate"
}
# 
ncol=3

DFPLOT2_6$Sparsity<-factor(DFPLOT2_6$Sparsity,levels=c("<=20%","20-50%",">50%"))

DFPLOT
DFPLOT$Method
DFPLOT$Mrate2<-factor(DFPLOT$Mrate2,levels=c("Balanced","10%","20%"))
DFPLOT$Method3<-factor(DFPLOT$Method3,levels=Method_toplot)

ggplot(data = DFPLOT, aes(x = Mrate2, y = Value, fill = Method3)) +
  scale_fill_manual(values=geom.colors_picked,labels=labels_out)+
  geom_bar(stat="identity", position = position_dodge(), col="black") + # adding bar plo
  # coord_cartesian(ylim=c(1500, 6400)) +
  theme(axis.text = element_text(colour = "black", size = 13),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15))  +
  theme_bw()+
  guides(linetype = guide_legend("",ncol=ncol,byrow=byr),fill = guide_legend("",ncol=ncol,byrow=byr),shape = guide_legend("",ncol=ncol,byrow=byr),colour = guide_legend("",ncol=ncol,byrow=byr),alpha = guide_legend("",ncol=ncol,byrow=byr) ,size = guide_legend("none"))+ylab("Power estimates")+xlab("Beta")+theme(legend.position="bottom")+
  # ylim(c(0,1))+
  xlab(xlabel) + ylab("Power estimates")


nm_tmp<-paste0("POWER_SPLIT_leaf.png")
ggsave(nm_tmp,height=5,width=4,dpi=300)
system(paste0("scp ",nm_tmp," ng:~"))
