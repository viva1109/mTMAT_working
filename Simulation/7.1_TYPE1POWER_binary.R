
What<-"Power"
Filter<-"NewFilter"
Outcome<-"logistic"
type_data<-"prop_corr"
Rank<-"Genus"

n_sim_bf<-100
n_perm<-200
n_core=30

if(Rank=="Family"){
  if(Filter=="original"){
    EEND<-39#17
  }else{
    EEND<-31#13
  }
}else{
  if(Filter=="original"){
    EEND<-73#23
  }else{
    EEND<-60#25
  }
}

HeatSTR<-paste0(Rank,What,"_210501_NEWSIMUL_SCALE_",Outcome,"_",type_data,"_TMATFIX_FILTER",Filter,"_n",n_sim_bf)
part<-1
jump<-4
Done<-50
Start<-1+jump*(part-1)+Done
End<-min(jump*part+Done,EEND)
cat(Start,End,"\n")

library(edgeR)# data/Lon_data
fol<-"/data/sharedcode/kjkim/Analysis/20190225Combine"
outfol<-paste0("~/analysis/20200928_longi_method/",HeatSTR)
infofol<-"~/analysis/20200928_longi_method/info/"
dir.create(outfol)
dir.create(infofol)
source_home<-"http://viva1109.iptime.org/"
source(paste0(source_home,"RFunctions/FunctionsTMAT/functions_tree_down.R"),encoding = "UTF-8")
source(paste0(source_home,"RFunctions/FunctionsTMAT/plot_TMAT.R"),encoding = "UTF-8")
source(paste0(source_home,"RFunctions/FunctionsTMAT/funcs_TMAT.R"),encoding = "UTF-8")
.libPaths("~/tmp")
library(mmmgee)
library(stringr)
library(ape)
library(parallel)
library(data.table)
rd_bf<-data.frame(fread("/data/MICROBIOME/mdprofile/5.ansan/open/ez/otu_table_mc2_new.txt",skip=1),check.names = F)
# ID_linker<-read.csv("~/analysis/20200928_longi_method/data/linked_IDs_GG.csv",header=T,stringsAsFactors=F)
ID_linker<-read.csv("~/analysis/20200928_longi_method/data/linked_IDs_GG_EZ.csv",header=T,stringsAsFactors=F)


rd<-list(rd_bf[,c(names(rd_bf)[1],ID_linker[,2],names(rd_bf)[dim(rd_bf)[2]])],
         rd_bf[,c(names(rd_bf)[1],ID_linker[,3],names(rd_bf)[dim(rd_bf)[2]])],
         rd_bf[,c(names(rd_bf)[1],ID_linker[,4],names(rd_bf)[dim(rd_bf)[2]])]
)


# rd<-readRDS(file.path(fol,"M678_Feb2019.rds"))
# rd1<-rd[[1]]
# 
# tax1<-rd1[,1]
# otutab<-rd1[,-c(1,dim(rd1)[2])]

# names(rd)

# sampleID<-colnames(otutab)

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
linked_IDs<-ID_linker[match(nms_go,ID_linker[,1]),2:4]       

# dir.create("~/analysis/20200928_longi_method/data/")
# write.csv(ID_linker[match(nms_go,ID_linker[,1]),1:4]   ,"~/analysis/20200928_longi_method/data/linked_IDs_GG_EZ.csv",row.names = F)


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



otu_info<-mclapply(1:length(otuTab_bf),function(ind_data){
  data<-otuTab_bf[[ind_data]]
  non_zero_prop<-apply(data,1,function(value){
    sum(value!=0)/length(value)
  })
  trc<-Info_IDs[[ind_data]]$Rcnts
  mean_rp<-apply(data,1,function(value){
    mean(value/trc)
  })
  
  list(non_zero_prop=non_zero_prop,mean_rp=mean_rp)
},mc.cores=length(otuTab_bf))

otuTab_bf_cbd<-do.call("cbind",otuTab_bf)

otu_info_totalinfo<-list(
  non_zero_prop<-apply(otuTab_bf_cbd,1,function(value){
    sum(value!=0)/length(value)
  }),
  total_trc<-c(Info_IDs[[1]]$Rcnts,Info_IDs[[2]]$Rcnts,Info_IDs[[3]]$Rcnts),
  mean_rp<-apply(otuTab_bf_cbd,1,function(value){
    mean(value/total_trc)
  })
)



otu_info_tax<-mclapply(1:length(otuTab_bf),function(ind_data){
  data<-NewRDs[[ind_data]]
  oid<-as.character(data[,1])
  tax<-as.character(data[,dim(data)[2]])
  list(oid,tax)
},mc.cores=length(otuTab_bf))




otu_info2<-lapply(otu_info,function(info){
  cri1<-which(info[[1]]>0.25)
  cri2<-which(info[[2]]>0.001)
  cri3<-intersect(cri1,cri2)
  list(cri1,cri2,cri3)
})

lapply(otu_info2,sapply,length)

# indsx<-lapply(otu_info2,"[[",2)
indsx<-lapply(otu_info2,"[[",1)
# otu_remain<-intersect(intersect(indsx[[1]],indsx[[2]]),indsx[[3]])

if(Filter=="original"){
  filter_cutoff<-0.001
  otu_remain<-otu_info_totalinfo[[3]]>filter_cutoff
}else{
  filter_cutoff<-0.25
  otu_remain_bf0<-otu_info_totalinfo[[1]]>filter_cutoff
  otu_remain_bf1<-substr(otu_info_tax[[1]][[1]],1,3)!="New"
  otu_remain<-otu_remain_bf0&otu_remain_bf1
}




# otu_info_tax[[1]][[2]][otu_remain]
# otu_info_tax[[2]][[2]][otu_remain]
# otu_info_tax[[3]][[2]][otu_remain]

# otu_remain2<-which(sapply(str_split(otu_info_tax[[1]][[2]],";"),"[",6)!=" g__")
# otu_remiain3<-intersect(otu_remain,otu_remain2)
otu_remiain3_bf<-otu_remain
otu_remiain3<-which(otu_remiain3_bf)

if(Rank=="Family"){
  OTU_INFO<-list(
    otu_info_tax[[1]][[1]][otu_remiain3],
    otu_info_tax[[1]][[2]][otu_remiain3],
  
      sapply(str_split(otu_info_tax[[1]][[2]][otu_remiain3],"; "),"[",5),
      sapply(str_split(otu_info_tax[[1]][[2]][otu_remiain3],"; "),"[",7)
  )
}else{
  OTU_INFO<-list(
    otu_info_tax[[1]][[1]][otu_remiain3],
    otu_info_tax[[1]][[2]][otu_remiain3],
    
    sapply(str_split(otu_info_tax[[1]][[2]][otu_remiain3],"; "),"[",6),
    sapply(str_split(otu_info_tax[[1]][[2]][otu_remiain3],"; "),"[",7)  
  )

}


# keep_OTUINFO<-OTU_INFO[[1]]
setdiff(keep_OTUINFO,OTU_INFO[[1]]) #deleted Find single OTU
setdiff(OTU_INFO[[1]],keep_OTUINFO) #newOTU
# intersect(keep_OTUINFO,OTU_INFO[[1]])
list_inve<-c(setdiff(keep_OTUINFO,OTU_INFO[[1]]),setdiff(OTU_INFO[[1]],keep_OTUINFO))

OTU_table<-lapply(otuTab_bf,function(data){
  data_out<-t(data[otu_remiain3,])
  rownames(data_out)<-nms_go
  colnames(data_out)<-OTU_INFO[[1]]
  return(data_out)
})




list_pheno<-file.path("/data/sharedcode/kjkim/data/Lon_data",paste0("HOS_AS",6:8,"_nokor.csv"))
pheD<-lapply(list_pheno,read.csv,stringsAsFactors=F)
sapply(pheD,colnames)

phD2<-lapply(pheD,function(data){
  names(data)<-gsub("^TH._","",toupper(names(data)))
  return(data)
})
sapply(phD2,colnames)





lapply(OTU_table,dim)
# OTU_INFO
# Info_IDs
SAMPLE_UID<-nms_go

tree_EZ<-read.tree("/data/MICROBIOME/mdprofile/5.ansan/open/ez/rep_set.tre")


drop_names<-tree_EZ$tip.label[!tree_EZ$tip.label%in%OTU_INFO[[1]]]
pr_tree<-drop.tip(tree_EZ,drop_names)
genus_otu_cnts<-sort(table(OTU_INFO[[3]]),decreasing = T)
GENERA_INFO<-list(genus_otu_cnts,names(genus_otu_cnts))
GENUS_LIST<-GENERA_INFO[[2]]




# Var_out<-sapply(phD2,function(data){
#   ind_s<-match(SAMPLE_UID,data[,1])
#   as.numeric(data[ind_s,"CREATININE"]>1.15)
# })

sparsities_List<-vector("list",length(GENUS_LIST))

for (genus_ind in 1:length(GENUS_LIST)){
  genus_nm<-GENUS_LIST[genus_ind]
  ind_col<-which(OTU_INFO[[3]]==genus_nm)  
  OTU_table_chosen<-lapply(OTU_table,function(data){
    data[,ind_col]
  })
  if ( length(ind_col)==1){
    sparsities_List[[genus_ind]]<-sapply(OTU_table_chosen,function(data){
      sum(data==0)/length(data)
    })
  }else{
    sparsities_List[[genus_ind]]<-sapply(OTU_table_chosen,function(data){
      apply(data,2,function(data2){
        sum(data2==0)/length(data2)
      })
    })
  }
}

library(edgeR)


mean_List<-vector("list",length(GENUS_LIST))
sd_List<-vector("list",length(GENUS_LIST))
cbe_rcnts<-c(Info_IDs[[1]]$Rcnts,Info_IDs[[2]]$Rcnts,Info_IDs[[3]]$Rcnts)
cbd_OTU_table<-do.call("rbind",OTU_table)
OTU_table_CPM_bf<-t(cpm(t(cbd_OTU_table),log=T,lib.size = cbe_rcnts))
OTU_table_CPM<-OTU_table_CPM_bf-min(OTU_table_CPM_bf)

for (genus_ind in 1:length(GENUS_LIST)){
  genus_nm<-GENUS_LIST[genus_ind]
  ind_col<-which(OTU_INFO[[3]]==genus_nm)  
  OTU_table_CPM_chosen<-OTU_table_CPM[,ind_col]
  
  if ( length(ind_col)==1){
    mean_List[[genus_ind]]<-mean(OTU_table_CPM_chosen[OTU_table_CPM_chosen!=0])
    sd_List[[genus_ind]]<-sd(OTU_table_CPM_chosen[OTU_table_CPM_chosen!=0])
  }else{
    mean_List[[genus_ind]]<- apply(OTU_table_CPM_chosen,2,function(data2){mean(data2[data2!=0])})
    sd_List[[genus_ind]]<- apply(OTU_table_CPM_chosen,2,function(data2){sd(data2[data2!=0])})
  }
}





sapply(sparsities_List,mean)
group1<-which(sapply(sparsities_List,mean)<=0.2)
group2<-which(sapply(sparsities_List,mean)<=0.4 & sapply(sparsities_List,mean)>0.2)
group3<-which(sapply(sparsities_List,mean)<=0.7 & sapply(sparsities_List,mean)>0.5)
group4<-which(sapply(sparsities_List,mean)>0.7)
gROUPs<-list(group1, group2, group3,group4)
gROUPs_otuno<-lapply(gROUPs,function(group_indiv){
  mean(sapply(group_indiv,function(upper_index){
    genus_nm<-GENUS_LIST[upper_index]
    ind_col<-which(OTU_INFO[[3]]==genus_nm)
    length(ind_col)
  }))
})


OTU_no<-sapply(1:length(GENUS_LIST),function(upper_index){
  genus_nm<-GENUS_LIST[upper_index]
  ind_col<-which(OTU_INFO[[3]]==genus_nm)
  length(ind_col)
})

Ogroup1<-which(OTU_no<=1)
Ogroup2<-which(OTU_no<=5 & OTU_no>1)
Ogroup3<-which(OTU_no<=15 & OTU_no>5)
Ogroup4<-which(OTU_no>15)
OgROUPs<-list(Ogroup1, Ogroup2, Ogroup3,Ogroup4)


saveRDS(list(GENUS_LIST,sapply(sparsities_List,mean),gROUPs,OgROUPs,gROUPs_otuno,sparsities_List,OTU_no,mean_List,sd_List),file.path(infofol,paste0(Rank,"_LongINFO",filter_cutoff,".csv")))
valOUT<-list(GENUS_LIST,sparsities_List,OTU_no,mean_List,sd_List,OTU_INFO,pr_tree)

names(valOUT)<-c("GENUS_LIST","sparsities_List","OTU_no","mean_List","sd_List","OTU_INFO","pr_tree")
saveRDS(valOUT,file.path(infofol,paste0(Rank,"_INFO_",filter_cutoff,".csv")))

ind_Cpr<-1
ind_n<-1

if(What=="Power"){
  # if(type_data=="log_prop"){
  #   beta_list<-c(0,10,30,60,120)/3*2.5
  # }else if(type_data=="prop"){
  #   beta_list<-c(0,10,30,60,120)/2.1
  # }else if(type_data=="prop_corr"){
  #   beta_list<-c(0,10,30,60,120)/1.6
  # }else if(type_data=="log_prop_corr"){
  #   beta_list<-c(0,10,30,60,120)/1.8
  # }
  # beta_list<-c(5,15,30,60)
  beta_list<-c(0.05,0.1,0.2,0.4)
  list_n<-c(100)
  list_c_prop<-c(0.25)
  propor_factor_list<-c(0,0.5,0.9)
}else if(What=="Type1"){
  list_n<-c(30,50,100)
  list_c_prop<-c(0.5,0.25)
  beta_list<-c(0)
  propor_factor_list<-c(0)
}

# mean(apply(cbd_OTU_table,2,sd))
#count
# beta_list<-c(0,0.0001,0.0005,0.001,0.005)
##count, corr 0.2
# beta_list<-c(0,0.0001,0.0005,0.001,0.005)
# beta_list<-c(0,0.001,0.005,0.01)
#count log
# beta_list<-c(0,10,30,60,120)/300



GENUS_LIST2<-GENUS_LIST[Start:End]
# GENUS_LIST2<-GENUS_LIST
sn_table<-matrix(nrow=length(list_n)*length(list_c_prop)*length(beta_list)*length(GENUS_LIST2)*length(propor_factor_list),ncol=6)
cnt<-0
sn_table[,1]<-1:dim(sn_table)[1]
colnames(sn_table)<-c("No","Cprop","N","Genera","beta","p")
for (j in 1:length(list_c_prop)){
  for ( i in 1:length(list_n)){
    for ( z in 1:length(GENUS_LIST2)){
      for ( z2 in 1:length(beta_list)){
        for ( z3 in 1:length(propor_factor_list)){
          cnt<-cnt+1
          sn_table[cnt,2]<-j
          sn_table[cnt,3]<-i
          sn_table[cnt,4]<-z
          sn_table[cnt,5]<-z2
          sn_table[cnt,6]<-z3
        }

      }

    }
  }
}
list_labels<-list(list_c_prop,list_n,GENUS_LIST2,beta_list)

# n_sim<-100
# n_perm<-200


####best 2000 / acceptable 1000
n_sim<-n_sim_bf


library(mmmgee)
# str<-"unstructured"
str<-"independence"
# str<-"exchangeable"
# str<-"ar1"
# str<-"m-dependent"

stat<-"score"
# type=16


library(GLMMMiRKAT)
library(CompQuadForm)
library(dirmult)
library(ecodist)
library(GUniFrac)
library(lme4)
library(MASS)
library(Matrix)
library(permute)
library(phyloseq)
library(MiRKAT)

# library(remotes)
# install_github("nyiuab/NBZIMM", force=T, build_vignettes=F)
library(NBZIMM)
library(nlme)
# install_github("chvlyl/ZIBR")
library(ZIBR)
library(mvtnorm)
source("~/analysis/20200928_longi_method/code/Functions_global.R")


output_power<-vector("list",dim(sn_table)[1])
# target_methods<-c("mTMAT15","mTMAT16")
# target_methods<-c("mTMAT15","mTMAT16","GLMM-MiRKAT","ZINB","ZIBR","LMM")
target_methods<-c("mTMAT_CORE","GLMM-MiRKAT","ZINB","ZIBR","LMM")
# target_methods<-c("mTMAT_CORE","mTMAT_ADD","GLMM-MiRKAT","ZINB","ZIBR","LMM")
# target_methods<-c("mTMAT15","mTMAT16","GLMM-MiRKAT","ZINB","ZIBR","LMM")


# target_methods<-c("mTMAT15","mTMAT16","GLMM-MiRKAT","ZINB","ZIBR")
# target_methods<-c("LMM")
# target_methods<-c("mTMAT15","mTMAT16")


b0<-0
sigma<-1
if(type_data=="log_prop"){
  sigma_r<-0 ###correlation 0
}else if(type_data=="prop"){
  sigma_r<-0 ###correlation 0
}else if(type_data=="prop_corr"){
  sigma_r<-0.5 ###correlation 0.2 # 1/3 -> 0.1
}else if(type_data=="log_prop_corr"){
  sigma_r<-0.5 ###correlation 0.2 # 1/3 -> 0.1
}

n_bootstrapping<-10000

if(What=="Type1"){
  lu<-sapply(1:n_bootstrapping,function(o){
    set.seed(o)
    Normi_compo<-rnorm(length(SAMPLE_UID),0,sigma_r)
    Normi<-Normi_compo[match(rep(SAMPLE_UID,3),SAMPLE_UID)]
    Normij<-rnorm(length(SAMPLE_UID)*3,0,sigma)
    pheno<-b0+Normi+Normij
    # cor.test(pheno,apply(Tax_abundance,1,sum))
    return(pheno)
  })
}else{
  
}




for ( ind_table in 1:dim(sn_table)[1]){
  ind_Cpr<-sn_table[ind_table,2]
  ind_n<-sn_table[ind_table,3]
  genus_ind<-sn_table[ind_table,4]
  beta_ind<-sn_table[ind_table,5]
  print(paste(ind_Cpr,ind_n,genus_ind))
  c_prop<-list_c_prop[ind_Cpr]
  n_case<-list_n[ind_n]*c_prop
  n_control<-list_n[ind_n]*(1-c_prop)
  beta_val<-beta_list[beta_ind]
  propor_factor<-propor_factor_list[sn_table[ind_table,6]]

  
  set.seed(1)

  ##########OTU select
  genus_nm<-GENUS_LIST2[genus_ind]
  
  
  ind_col<-which(OTU_INFO[[3]]==genus_nm)
  ind_colTF<-OTU_INFO[[3]]!=genus_nm
  ind_colTF[is.na(ind_colTF)]<-T
  
  if(length(ind_col)==1){
    INDI_SINGLE<-T
  }else{
    INDI_SINGLE<-F
  }
  ###tree process
  otu_id_chosen<-OTU_INFO[[1]][ind_col]
  tree_chosen_bf2<-drop.tip(pr_tree,OTU_INFO[[1]][!OTU_INFO[[1]]%in%otu_id_chosen])
  
  if(!is.rooted(tree_chosen_bf2)){
    tree_chosen_bf<-root(tree_chosen_bf2, 1, r = TRUE)
  }else{
    tree_chosen_bf<-tree_chosen_bf2
  }
  
  if(!is.binary(tree_chosen_bf)){
    tree_chosen<-multi2di(tree_chosen_bf)
  }else{
    tree_chosen<-tree_chosen_bf
  }
  
  whatnode<-sapply((Ntip(tree_chosen)+1):(Ntip(tree_chosen)+Nnode(tree_chosen)),function(data,pt){length(find_childtips(data,pt))},tree_chosen)
  target_tips_info<-get_target_tips(whatnode,"random",tree_chosen)
  ind_which<-target_tips_info[1]
  candi_factor<-find_childtips(ind_which,tree_chosen)
  len_candi_factor<-length(candi_factor)
  no_causal<-len_candi_factor*propor_factor
  if(propor_factor==0){
    hmn_list<-rep(1,n_sim)
  }else{
    no_causal<-max(no_causal,1)
    floor_no_causal<-floor(no_causal)
    hmn_list<-rep(floor_no_causal+1,n_sim)
    hmn_list[1:(n_sim*(1-(no_causal-floor_no_causal)))]<-floor_no_causal
  }
  
  OTU_SELECTED_TMP<-lapply(OTU_table,function(data){
    data[,ind_col]
  })
  
  if(length(ind_col)==1){
    VEC_OTU_SELECTED_TMP<-matrix(c(OTU_SELECTED_TMP[[1]],OTU_SELECTED_TMP[[2]],OTU_SELECTED_TMP[[3]]),ncol=length(ind_col))
  }else{
    VEC_OTU_SELECTED_TMP<-rbind(OTU_SELECTED_TMP[[1]],OTU_SELECTED_TMP[[2]],OTU_SELECTED_TMP[[3]])
  }
  
  
  
  IDsC_TMP<-rep(SAMPLE_UID,3)  
  VEC_total.reads_TMP<-c(Info_IDs[[1]]$Rcnts,Info_IDs[[2]]$Rcnts,Info_IDs[[3]]$Rcnts)
  TIME_TMP<-c(rep(1,length(SAMPLE_UID)),rep(2,length(SAMPLE_UID)),rep(3,length(SAMPLE_UID)))
  
  if(type_data%in%c("log_prop","log_prop_corr")){
    VEC_PROP_OTU_table_TMP<-log((VEC_OTU_SELECTED_TMP/VEC_total.reads_TMP)+1)
  }else if(type_data%in%c("prop","prop_corr")){
    VEC_PROP_OTU_table_TMP<-VEC_OTU_SELECTED_TMP/VEC_total.reads_TMP
  }else{
    print(paste0("type_data FAIL: ",type_data))
  }
  

  

  
  ind_causal_eachSimul<-lapply(1:n_bootstrapping,function(o){
    set.seed(o)
    ind_cs<-sample(candi_factor,sample(hmn_list,1))
    
    otuID_causal<-tree_chosen$tip.label[ind_cs]
    ind_causal<-match(otuID_causal,OTU_INFO[[1]][ind_col])
  })
  
  

  if(What=="Power"){
    beta<-beta_val
    lu<-sapply(1:n_bootstrapping,function(o){
      set.seed(o)
      ind_causal_tmp<-ind_causal_eachSimul[[o]]
      b1_s<-matrix(rep(beta,length(ind_causal_tmp)),ncol=1)
      Tax_abundance<-matrix(scale(VEC_PROP_OTU_table_TMP[,ind_causal_tmp]),ncol=length(ind_causal_tmp))
      Normi_compo<-rnorm(length(SAMPLE_UID),0,sigma_r)
      Normi<-Normi_compo[match(IDsC_TMP,SAMPLE_UID)]
      Normij<-rnorm(dim(VEC_PROP_OTU_table_TMP)[1],0,sigma)
      
      pheno<-b0+Tax_abundance%*%b1_s+Normi+Normij
      
      cor.test(pheno,apply(Tax_abundance,1,sum))
      return(pheno)
    })
  }
  
  Casecutoff<-quantile(lu,0.8)
  PhenoSimul<-lu>Casecutoff
  TMPMat_bf<-apply(PhenoSimul,1,mean)
  TMPMat<-cbind(TMPMat_bf[1:length(SAMPLE_UID)],TMPMat_bf[1:length(SAMPLE_UID)+length(SAMPLE_UID)],TMPMat_bf[1:length(SAMPLE_UID)+length(SAMPLE_UID)*2])
  # cutoff_binary<-quantile(TMPMat,0.8)
  # binary_out<-apply(TMPMat>cutoff_binary,2,as.numeric)
  
  CaseControl<-apply(TMPMat,1,mean)
  cutoff2<-quantile(CaseControl,0.8)
  CaseControl_ind<-CaseControl>cutoff2
  
  base_case<-which(CaseControl_ind)#95
  base_control<-which(!CaseControl_ind)#373
  
  set.seed(1)
  if(length(base_case)<n_case){
    ind_case<-sample(base_case,n_case,replace = T)
  }else{
    ind_case<-sample(base_case,n_case,replace = F)
  }
  
  if(length(base_control)<n_control){
    ind_control<-sample(base_control,n_control,replace = T)
  }else{
    ind_control<-sample(base_control,n_control,replace = F)
  }
  ind_samples<-c(ind_case,ind_control)
  

  #0.005
  


  # mean(apply(lu[ind_case,]>Casecutoff,1,mean))
  # mean(apply(lu[ind_control,]>Casecutoff,1,mean))
  
  OTU_table_SP<-lapply(OTU_table,function(data){
    data[ind_samples,]
  })

 
  # OTU_table_SP_stack<-do.call("rbind",OTU_table_SP)
  # apply(OTU_table_SP_stack,2,sd)
  
  Info_IDs_SP<-lapply(Info_IDs,function(Info_ID){
    Rcnts<-Info_ID[[1]][ind_samples]
    timeID<-Info_ID[[2]][ind_samples]
    return(list(Rcnts=Rcnts,timeID=timeID))
  })


  OTU_table_chosen<-lapply(OTU_table_SP,function(data){
    tmpdf<-matrix(data[,ind_col],ncol=length(ind_col))
    colnames(tmpdf)<-colnames(data)[ind_col]
    rownames(tmpdf)<-rownames(data)
    return(tmpdf)
  })
  
  OTU_table_not_chosen<-lapply(OTU_table_SP,function(data){
    tmpdf<-matrix(data[,ind_colTF],ncol=sum(ind_colTF))
    colnames(tmpdf)<-colnames(data)[ind_colTF]
    rownames(tmpdf)<-rownames(data)
    return(tmpdf)
  })
  
  

  Counts_C<-lapply(1:length(OTU_table_chosen),function(ind_data){
    tt_val<-Info_IDs_SP[[ind_data]]$Rcnts-apply(OTU_table_chosen[[ind_data]],1,sum)
    matrix(tt_val,ncol=1)
  })
  
  
  Counts_C_and_C<-lapply(1:length(OTU_table_chosen),function(ind_data){
    tt_val<-Info_IDs_SP[[ind_data]]$Rcnts-apply(OTU_table_chosen[[ind_data]],1,sum)-apply(OTU_table_not_chosen[[ind_data]],1,sum)
    matrix(tt_val,ncol=1)
  })
  
  OTU_table_not_chosen2<-lapply(1:length(OTU_table_not_chosen),function(ind_data){
    cbind(OTU_table_not_chosen[[ind_data]],Counts_C_and_C[[ind_data]])
  })
  
  
  
  
  
  
  VEC_OTU_table_chosen<-rbind(OTU_table_chosen[[1]],OTU_table_chosen[[2]],OTU_table_chosen[[3]])
  IDsC<-rep(SAMPLE_UID[ind_samples],3)  
  VEC_total.reads<-c(Info_IDs[[1]]$Rcnts[ind_samples],Info_IDs[[2]]$Rcnts[ind_samples],Info_IDs[[3]]$Rcnts[ind_samples])
  
  
  if(type_data=="log_prop"){
    VEC_PROP_OTU_table_chosen<-log((VEC_OTU_table_chosen/VEC_total.reads)+1)
  }else if(type_data=="prop"){
    VEC_PROP_OTU_table_chosen<-VEC_OTU_table_chosen/VEC_total.reads
  }else if(type_data=="prop_corr"){
    VEC_PROP_OTU_table_chosen<-VEC_OTU_table_chosen/VEC_total.reads
  }else if(type_data=="log_prop_corr"){
    VEC_PROP_OTU_table_chosen<-log((VEC_OTU_table_chosen/VEC_total.reads)+1)
  }
  
  
  # VEC_PROP_OTU_table_chosen<-VEC_OTU_table_chosen
  # VEC_PROP_OTU_table_chosen<-log(VEC_OTU_table_chosen+1)
  # Outcome

  if(beta_val==0){

    
    if(INDI_SINGLE){

    }else{
      tmpTable<-cbind(do.call("rbind",OTU_table_chosen),do.call("rbind",Counts_C))
      OTU_table_RF_DF_bf2<-Rarefy(tmpTable,min(rowSums(tmpTable)))$otu.tab.rff[,1:length(ind_col)]
      colnm_tmp<-colnames(tmpTable)[1:length(ind_col)]
      OTU_table_RF_DF_bf<-OTU_table_RF_DF_bf2
      ind_out_Kernel<-!apply(OTU_table_RF_DF_bf,1,sum)==0
      OTU_table_RF_DF<-OTU_table_RF_DF_bf[ind_out_Kernel,]
      Ks <- Kernels(OTU_table_RF_DF, tree_chosen)
      # genera_vec_RF<-apply(OTU_table_RF_DF_bf,1,sum)
    }
  }else{
    
  }
  
  
  p_vals<-mclapply(1:n_sim,function(o){
    set.seed(o)
    vout_bf<-cbind(PhenoSimul[1:length(SAMPLE_UID),o],PhenoSimul[1:length(SAMPLE_UID)+length(SAMPLE_UID),o],PhenoSimul[1:length(SAMPLE_UID)+length(SAMPLE_UID)*2,o])[ind_samples,]
    Var_out2<-apply(vout_bf,2,as.numeric)
    TIME<-c(rep(1,length(SAMPLE_UID[ind_samples])),rep(2,length(SAMPLE_UID[ind_samples])),rep(3,length(SAMPLE_UID[ind_samples])))
    OTU_table_chosen_power<-OTU_table_chosen

    ###KERNEL
    if(beta_val==0){
      
      
    }else{

      if(INDI_SINGLE){
        
      }else{
        
        tmpTable<-cbind(do.call("rbind",OTU_table_chosen_power),do.call("rbind",Counts_C))
        OTU_table_RF_DF_bf2<-Rarefy(tmpTable,min(rowSums(tmpTable)))$otu.tab.rff[,1:length(ind_col)]
        colnm_tmp<-colnames(tmpTable)[1:length(ind_col)]
        
        OTU_table_RF_DF_bf<-OTU_table_RF_DF_bf2
        
        ind_out_Kernel<-!apply(OTU_table_RF_DF_bf,1,sum)==0
        OTU_table_RF_DF<-OTU_table_RF_DF_bf[ind_out_Kernel,]
        Ks <- Kernels(OTU_table_RF_DF, tree_chosen)
        # genera_vec_RF<-apply(OTU_table_RF_DF_bf,1,sum)
      }
    }

    
    
    OTU_table_chosen_power_genera<-lapply(OTU_table_chosen_power,function(data){
      if(INDI_SINGLE){
        data
      }else{
        apply(data,1,sum)
      }
      
    })
    
    y_vec<-c()
    id_vec<-c()
    trc_vec<-c()
    time_vec<-c()
    genera_vec<-c()
    genera_vec_prop<-c()
    
    for ( i in 1:dim(Var_out2)[2]){
      y_vec<-c(y_vec,Var_out2[,i])
      time_vec<-c(time_vec,rep(i,length(Var_out2[,i])))
      id_vec<-c(id_vec,rownames(OTU_table_chosen_power[[1]]))
      trc_vec<-c(trc_vec,Info_IDs_SP[[i]]$Rcnts)
      genera_vec<-c(genera_vec,OTU_table_chosen_power_genera[[i]])
    }
    genera_vec_prop<-genera_vec/trc_vec
    # lapply(list(y_vec,id_vec,trc_vec,genera_vec),length)
    
    
    # vv_OTU_table_chosen_power<-c(OTU_table_chosen_power[[1]],OTU_table_chosen_power[[2]],OTU_table_chosen_power[[3]])
    # cor.test(y_vec,vv_OTU_table_chosen_power)
    # cor.test(tdata$Y.tran,tdata$y_vec)
    
    mTMAT_res15<-NA
    mTMAT_res16<-NA
    mTMAT_res15CS<-NA
    mTMAT_res15AR1<-NA
    mTMAT_res15UN<-NA
    mTMAT_res15_compo<-NA
    mTMAT_res16CS<-NA
    mTMAT_res16AR1<-NA
    mTMAT_res16UN<-NA
    mTMAT_res16_compo<-NA
    mTMAT_res15INDEP<-NA
    mTMAT_res16INDEP<-NA
    if("mTMAT_CORE"%in%target_methods){
      
      str<-"independence"
      mTMAT_res15INDEP<-calculateTMAT(Var_out2,Info_IDs,OTU_table_chosen_power,tree_chosen,INDI_SINGLE,SAMPLE_UID,TMATtype=15,stat_i="score",str_i=str)
      mTMAT_res16INDEP<-calculateTMAT(Var_out2,Info_IDs,OTU_table_chosen_power,tree_chosen,INDI_SINGLE,SAMPLE_UID,TMATtype=16,stat_i="score",str_i=str)
      
 
      str<-"independence"
      mTMAT_res15_compo<-calculateTMAT(Var_out2,Info_IDs,OTU_table_chosen_power,tree_chosen,INDI_SINGLE,SAMPLE_UID,TMATtype=15,stat_i="score",str_i=str,OTU_not_chosen=OTU_table_not_chosen2)
      mTMAT_res16_compo<-calculateTMAT(Var_out2,Info_IDs,OTU_table_chosen_power,tree_chosen,INDI_SINGLE,SAMPLE_UID,TMATtype=16,stat_i="score",str_i=str,OTU_not_chosen=OTU_table_not_chosen2)
      
    }
    
    if("mTMAT15"%in%target_methods){
      mTMAT_res15<-calculateTMAT(Var_out2,Info_IDs,OTU_table_chosen_power,tree_chosen,INDI_SINGLE,SAMPLE_UID,TMATtype=15,stat_i="score",str_i=str)
    }
    if("mTMAT16"%in%target_methods){
      mTMAT_res16<-calculateTMAT(Var_out2,Info_IDs,OTU_table_chosen_power,tree_chosen,INDI_SINGLE,SAMPLE_UID,TMATtype=16,stat_i="score",str_i=str)
    }
    
    if("mTMAT_ADD"%in%target_methods){
      str<-"exchangeable"
      mTMAT_res15CS<-calculateTMAT(Var_out2,Info_IDs,OTU_table_chosen_power,tree_chosen,INDI_SINGLE,SAMPLE_UID,TMATtype=15,stat_i="score",str_i=str)
      mTMAT_res16CS<-calculateTMAT(Var_out2,Info_IDs,OTU_table_chosen_power,tree_chosen,INDI_SINGLE,SAMPLE_UID,TMATtype=16,stat_i="score",str_i=str)
      str<-"ar1"
      mTMAT_res15AR1<-calculateTMAT(Var_out2,Info_IDs,OTU_table_chosen_power,tree_chosen,INDI_SINGLE,SAMPLE_UID,TMATtype=15,stat_i="score",str_i=str)
      mTMAT_res16UN<-calculateTMAT(Var_out2,Info_IDs,OTU_table_chosen_power,tree_chosen,INDI_SINGLE,SAMPLE_UID,TMATtype=16,stat_i="score",str_i=str)
      str<-"unstructured"
      mTMAT_res15UN<-calculateTMAT(Var_out2,Info_IDs,OTU_table_chosen_power,tree_chosen,INDI_SINGLE,SAMPLE_UID,TMATtype=15,stat_i="score",str_i=str)
      mTMAT_res16UN<-calculateTMAT(Var_out2,Info_IDs,OTU_table_chosen_power,tree_chosen,INDI_SINGLE,SAMPLE_UID,TMATtype=16,stat_i="score",str_i=str)
      
    }

    
    output_mirkat_result<-NA
    f4result<-NA
    zibr_result<-NA
    LMM_result<-NA
    zibr_resultLogi<-NA
    zibr_resultBeta<-NA
    LMM_resultlog<-NA
    NBMMresult<-NA
    
    if("GLMM-MiRKAT"%in%target_methods){
      if(INDI_SINGLE){
        output_mirkat_result<-NA
      }else{
        
        if(length(unique(y_vec[ind_out_Kernel]))>10){
          print("level > 10 => continuous response")
          output_mirkat<-try(Fixed_GLMMMiRKAT(y_vec[ind_out_Kernel], cov=NULL, id=id_vec[ind_out_Kernel], Ks=Ks, model="gaussian",n.perm = n_perm))
          # output_mirkat<-try(GLMMMiRKAT::GLMMMiRKAT(y_vec[ind_out_Kernel], cov=NULL, id=id_vec[ind_out_Kernel], Ks=Ks, model="gaussian",n.perm = n_perm))
        }else{
          output_mirkat<-try(Fixed_GLMMMiRKAT(y_vec[ind_out_Kernel], cov=NULL, id=id_vec[ind_out_Kernel], Ks=Ks, model="binomial",n.perm = n_perm))  
        }

        
        
        
        if(class(output_mirkat)!="try-error"){
          output_mirkat_result<-output_mirkat$aGLMMMiRKAT
        }else{
          output_mirkat_result<-NA
        }
      }
    }

    
    # output_cskat<-CSKAT(formula.H0 = y_vec ~ 1 + (1 | id_vec), Ks = Ks,nperm = n_perm)
    # dim(genera_vec_RF)<-NULL
    
    if("ZINB"%in%target_methods){
      zinbDF<-data.frame(Responsezinb=genera_vec,CovIntzinb=y_vec, IDzinb=id_vec,TRCzinb=trc_vec)
      # zinbDF<-data.frame(Response=OTU_table_RF_DF_bf2,CovInt=y_vec, ID=id_vec,TRC=trc_vec)
      f3 = try(glmm.zinb(Responsezinb ~ CovIntzinb + offset(log(TRCzinb)), random = ~ 1|IDzinb,data=zinbDF))
      # f3 = try(glmm.zinb(Response ~ CovInt , random = ~ 1|ID,data=zinbDF))
      if(class(f3)[1]!="try-error"){
        f4<-summary(f3)
        f4result<-f4$tTable[2,'p-value']
      }else{
        f4result<-NA
      }

      NBMM = try(glmm.nb(Responsezinb ~ CovIntzinb + offset(log(TRCzinb)), random = ~1|id_vec, verbose = F,data=zinbDF) )
      
      if(class(NBMM)[1]!="try-error"){
        NBMMresult= summary(NBMM)$tTable["CovIntzinb", ][c(5)]
      }else{
        NBMMresult<-NA
      }
    }
    
    # cor.test(genera_vec/trc_vec,y_vec)
    
    
    if("ZIBR"%in%target_methods){
      zibr.fit <- try(zibr(logistic.cov = matrix(y_vec,ncol=1), beta.cov = matrix(y_vec,ncol=1), Y = genera_vec/trc_vec,subject.ind = id_vec,time.ind =time_vec))
      if(class(zibr.fit)!="try-error"){
        zibr_result<-zibr.fit$joint.p
        zibr_resultLogi<-zibr.fit$logistic.est.table[2,2]
        zibr_resultBeta<-zibr.fit$beta.est.table[2,2]
      }else{
        zibr_result<-NA
        zibr_resultLogi<-NA
        zibr_resultBeta<-NA
      }
    }

    if("LMM"%in%target_methods){
      tdata <- data.frame(Y.tran=asin(sqrt(genera_vec_prop)),y_vec,SID=id_vec)
      lme.fit <- try(lme(Y.tran ~ y_vec,random=~1| SID, data = tdata))
      
      if(class(lme.fit)!="try-error"){
        coef.mat <- summary(lme.fit)$tTable[-1,c(1,5)]
        LMM_result <- coef.mat[2]
      }else{
        LMM_result<-NA
      }
      
      
      # tdata <- data.frame(Y.tran=log(genera_vec_prop+1),y_vec,SID=id_vec)
      # lme.fit <- try(lme(Y.tran ~ y_vec,random=~1| SID, data = tdata))
      # 
      # if(class(lme.fit)!="try-error"){
      #   coef.mat <- summary(lme.fit)$tTable[-1,c(1,5)]
      #   LMM_resultlog <- coef.mat[2]
      # }else{
      #   LMM_resultlog<-NA
      # }
    }
    if(o%%n_core==0) print(o)
    
    resultout<-c(mTMAT_res15,mTMAT_res16,output_mirkat_result,f4result,zibr_result,LMM_result,NBMMresult,zibr_resultLogi,zibr_resultBeta,mTMAT_res15INDEP,mTMAT_res15CS,mTMAT_res15AR1,mTMAT_res15UN,mTMAT_res15_compo,mTMAT_res16INDEP,mTMAT_res16CS,mTMAT_res16AR1,mTMAT_res16UN,mTMAT_res16_compo)
    names(resultout)<-c("mTMAT_15","mTMAT_16","GLMM_MiRKAT","FZINBMM","ZIBR","LMM","NBMM","zibr_resultLogi","zibr_resultBeta","mTMAT_15INDEP","mTMAT_15CS","mTMAT_15AR1","mTMAT_15UN","mTMAT_15_compo","mTMAT_res16INDEP","mTMAT_16CS","mTMAT_16AR1","mTMAT_16UN","mTMAT_16_compo")
    return(resultout)
    
  },mc.cores=n_core)
  
  
  # tmptmp<-do.call("rbind",p_vals)
  # apply(tmptmp<0.1,2,mean,na.rm=T)

  output_power[[ind_table]]<-do.call("rbind",p_vals)
}
output_power[[ind_table]][,c("ZIBR","mTMAT_15INDEP")]
# lapply(output_power[c(7:9,22:24)],function(data){
#   data[,c(1:2,6:7)]
# })

cutoff<-0.1
value<-sapply(output_power,function(value){
  sum(value<0.1,na.rm=T)/length(!is.na(value))
})
value
mean(value)

res<-lapply(c(0.1,0.05,0.01,0.005),function(cutoff){
  value<-sapply(output_power,function(value){
    sum(value<cutoff,na.rm=T)/length(!is.na(value))
  })
  tapply(value,paste0(sn_table[,2],"_",sn_table[,3]),mean)  
})

cutoff<-0.1
Pvalue_bf<-t(sapply(output_power,function(value){
  apply(value,2,function(value2){
    sum(value2<cutoff,na.rm=T)/sum(!is.na(value2))
  })
}))

if(What=="Type1"){
  split_string<-paste0(sn_table[,2],"_",sn_table[,3])  #genera
  # split_string<-paste0(sn_table[,2],"_",sn_table[,3],"_",sn_table[,4])  #genera
}else{
  split_string<-paste0(sn_table[,5],"_",sn_table[,6])
  # split_string<-paste0(sn_table[,4],"_",sn_table[,5],"_",sn_table[,6])  #genera
}



resp<-split(data.frame(Pvalue_bf),split_string)
t(sapply(resp,function(data){
  apply(data,2,mean,na.rm=T)
}))


saveRDS(list(sn_table,list_labels,output_power),paste0(outfol,paste("/",HeatSTR,str,stat,Start,End,sep="_"),".rds"))
system(paste0("scp ",paste0(outfol,paste("/",HeatSTR,str,stat,Start,End,sep="_"),".rds")," ng:~"))
# 
# 
# Pvalue<-sapply(output_power,function(value){
#   apply(value<cutoff,2,sum,na.rm=T)/apply(value,2,function(val){length(!is.na(val))})
# })
# sn_table
# 
# 
# Power_result<-split(data.frame(t(Pvalue)),sn_table[1:17,5])
# sapply(Power_result,apply,2,mean,na.rm=T)
# 
# 
# 
# result2<-do.call("rbind",res)
# write.csv(result2,"type1_longtiTMAT_scoreAR1.csv")
# system("scp type1_longtiTMAT_*.csv ng:~")
# list_n
# 
# 
