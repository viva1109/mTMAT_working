n_sim_bf<-400
DVD<-1
n_core=30

HeatSTR<-paste0("GENRA_210406_2CRSEC_BETA_02D",DVD)
part<-1
jump<-25
Done<-24
Start<-1+jump*(part-1)+Done
End<-min(jump*part+Done,73)
cat(Start,End,"\n")


library(edgeR)# data/Lon_data
fol<-"/data/sharedcode/kjkim/Analysis/20190225Combine"
outfol<-"~/analysis/20200928_longi_method/GENRA_diffN_BETA_2CORSSSEC/"
infofol<-"~/analysis/20200928_longi_method/info/"
dir.create(outfol)
dir.create(infofol)
source_home<-"http://viva1109.iptime.org/"
source(paste0(source_home,"RFunctions/Tree_based_methods.R"),encoding = "UTF-8")
source(paste0(source_home,"RFunctions/FunctionsTMAT/functions_tree_down.R"),encoding = "UTF-8")
source(paste0(source_home,"RFunctions/FunctionsTMAT/plot_TMAT.R"),encoding = "UTF-8")
source(paste0(source_home,"RFunctions/FunctionsTMAT/funcs_TMAT.R"),encoding = "UTF-8")
source(paste0(source_home,"RFunctions/FunctionsTMAT/funcs_MiRKAT.R"),encoding = "UTF-8")
source(paste0(source_home,"RFunctions/FunctionsTMAT/funcs_OMiAT.R"),encoding = "UTF-8")
source(paste0(source_home,"/RFunctions/FunctionsTMAT/ANCOM.R"),encoding = "UTF-8")
source(paste0(source_home,"/RFunctions/FunctionsTMAT/stats_and_funcs.R"),encoding = "UTF-8")
source(paste0(source_home,"/RFunctions/FunctionsTMAT/simulation_funcs.R"),encoding = "UTF-8")
.libPaths("~/tmp")
library(mmmgee)
library(MiSPU)
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
  zero_prop<-apply(data,1,function(value){
    sum(value!=0)/length(value)
  })
  trc<-Info_IDs[[ind_data]]$Rcnts
  mean_rp<-apply(data,1,function(value){
    mean(value/trc)
  })
  
  list(zero_prop=zero_prop,mean_rp=mean_rp)
},mc.cores=length(otuTab_bf))

otuTab_bf_cbd<-do.call("cbind",otuTab_bf)

otu_info_totalinfo<-list(
  zero_prop<-apply(otuTab_bf_cbd,1,function(value){
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
  cri1<-which(info[[1]]>0.5)
  cri2<-which(info[[2]]>0.001)
  cri3<-intersect(cri1,cri2)
  list(cri1,cri2,cri3)
})

lapply(otu_info2,sapply,length)

indsx<-lapply(otu_info2,"[[",2)
# otu_remain<-intersect(intersect(indsx[[1]],indsx[[2]]),indsx[[3]])

otu_remain<-otu_info_totalinfo[[3]]>0.001


otu_info_tax[[1]][[2]][otu_remain]
otu_info_tax[[2]][[2]][otu_remain]
otu_info_tax[[3]][[2]][otu_remain]

# otu_remain2<-which(sapply(str_split(otu_info_tax[[1]][[2]],";"),"[",6)!=" g__")
# otu_remiain3<-intersect(otu_remain,otu_remain2)
otu_remiain3_bf<-otu_remain
otu_remiain3<-which(otu_remiain3_bf)


OTU_INFO<-list(
  otu_info_tax[[1]][[1]][otu_remiain3],
  otu_info_tax[[1]][[2]][otu_remiain3],
  sapply(str_split(otu_info_tax[[1]][[2]][otu_remiain3],"; "),"[",6),
  sapply(str_split(otu_info_tax[[1]][[2]][otu_remiain3],"; "),"[",7)
)


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




Var_out<-sapply(phD2,function(data){
  ind_s<-match(SAMPLE_UID,data[,1])
  as.numeric(data[ind_s,"CREATININE"]>1.15)
})

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
group2<-which(sapply(sparsities_List,mean)<=0.5 & sapply(sparsities_List,mean)>0.2)
group3<-which(sapply(sparsities_List,mean)<=0.8 & sapply(sparsities_List,mean)>0.5)
group4<-which(sapply(sparsities_List,mean)>0.8)
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


saveRDS(list(GENUS_LIST,sapply(sparsities_List,mean),gROUPs,OgROUPs,gROUPs_otuno,sparsities_List,OTU_no,mean_List,sd_List),file.path(infofol,"FAMILY_INFO.csv"))
valOUT<-list(GENUS_LIST,sparsities_List,OTU_no,mean_List,sd_List,OTU_INFO,pr_tree)

names(valOUT)<-c("GENUS_LIST","sparsities_List","OTU_no","mean_List","sd_List","OTU_INFO","pr_tree")
# saveRDS(valOUT,file.path(infofol,"FAMILY_INFO_DASIM.csv"))




ind_Cpr<-1
ind_n<-1


list_n<-c(100)
list_c_prop<-c(0.5)
beta_list<-c(0.01,0.02,0.04)
propor_factor_list<-c(0,0.5,0.9)
# beta_list<-c(0)
GENUS_LIST2<-GENUS_LIST[Start:End]
# GENUS_LIST2<-GENUS_LIST
sn_table<-matrix(nrow=length(list_n)*length(list_c_prop)*length(beta_list)*length(GENUS_LIST2)*length(propor_factor_list),ncol=6)
cnt<-0
sn_table[,1]<-1:dim(sn_table)[1]

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
n_sim<-n_sim_bf/DVD
n_perm<-400/DVD

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

# for ( ind_table in 1:dim(sn_table)[1]){
for ( ind_table in ind_fix){
  
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
  base_case<-which(Var_out[,1]==1)#113
  base_control<-which(Var_out[,1]==0)#363
  set.seed(1)
  ind_case<-sample(base_case,n_case,replace = F)
  ind_control<-sample(base_control,n_control,replace = F)
  ind_samples<-c(ind_case,ind_control)
  
  OTU_table_SP<-lapply(OTU_table,function(data){
    data[ind_samples,]
  })
  
  Info_IDs_SP<-lapply(Info_IDs,function(Info_ID){
    Rcnts<-Info_ID[[1]][ind_samples]
    timeID<-Info_ID[[2]][ind_samples]
    return(list(Rcnts=Rcnts,timeID=timeID))
  })
  genus_nm<-GENUS_LIST2[genus_ind]
  ind_col<-which(OTU_INFO[[3]]==genus_nm)
  if(length(ind_col)==1){
    INDI_SINGLE<-T
  }else{
    INDI_SINGLE<-F
  }
  
  OTU_table_chosen<-lapply(OTU_table_SP,function(data){
    tmpdf<-matrix(data[,ind_col],ncol=length(ind_col))
    colnames(tmpdf)<-colnames(data)[ind_col]
    rownames(tmpdf)<-rownames(data)
    return(tmpdf)
  })
  
  Counts_C<-lapply(1:length(OTU_table_chosen),function(ind_data){
    tt_val<-Info_IDs_SP[[ind_data]]$Rcnts-apply(OTU_table_chosen[[ind_data]],1,sum)
    matrix(tt_val,ncol=1)
  })
  
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
  
  
  
  if(beta_val==0){
    
    
    if(INDI_SINGLE){
      
    }else{
      tmpTable<-cbind(do.call("rbind",OTU_table_chosen),do.call("rbind",Counts_C))
      OTU_table_RF_DF_bf2<-Rarefy(tmpTable,min(rowSums(tmpTable)))$otu.tab.rff[,1:length(ind_col)]
      colnm_tmp<-colnames(tmpTable)[1:length(ind_col)]
      OTU_table_RF_DF_bf<-OTU_table_RF_DF_bf2
      ind_out_Kernel<-!apply(OTU_table_RF_DF_bf,1,sum)==0
      OTU_table_RF_DF<-OTU_table_RF_DF_bf[ind_out_Kernel,]
      Ks <- GLMMMiRKAT::Kernels(OTU_table_RF_DF, tree_chosen)
      
      
      # genera_vec_RF<-apply(OTU_table_RF_DF_bf,1,sum)
    }
  }else{
    
  }
  
  
  
  p_vals<-mclapply(1:n_sim,function(o){
    set.seed(o)
    
    
    # print(ind_table)
    Var_out2_bf<-Var_out[ind_samples,]
    Var_out2<-Var_out2_bf
    Var_out2[,1]<-sample(Var_out2_bf[,1])
    Var_out2[,2]<-sample(Var_out2_bf[,2])
    Var_out2[,3]<-sample(Var_out2_bf[,3])
    
    # OTU_nms<-colnames(OTU_table[[1]])[ind_col]
    
    
    # OTU_table_chosen_remain<-lapply(OTU_table_SP,function(data){
    #   data[,setdiff(1:dim(data)[2],ind_col)]
    # })
    
    
    
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
    ind_cs<-sample(candi_factor,sample(hmn_list,1))
    
    otuID_causal<-tree_chosen$tip.label[ind_cs]
    ind_causal<-match(otuID_causal,colnames(OTU_table_chosen[[1]]))
    
    
    
    
    
    OTU_table_chosen_power_bf<-OTU_table_chosen
    
    for ( col_ind in ind_causal){
      for ( k in 1:3){
        ind_case_tmp<-Var_out2[,k]==1
        core_added<-sd(OTU_table_chosen[[k]][,col_ind])
        OTU_table_chosen_power_bf[[k]][ind_case_tmp,col_ind]<-OTU_table_chosen[[k]][ind_case_tmp,col_ind]+floor(beta_val*core_added)
      } 
      
    }
    # if(beta_val==0){
    OTU_table_chosen_power<-OTU_table_chosen_power_bf
    # }else{
    # OTU_table_chosen_power[[2]]<-floor((OTU_table_chosen_power_bf[[1]]+OTU_table_chosen_power_bf[[2]])/2)
    # OTU_table_chosen_power[[3]]<-floor((OTU_table_chosen_power_bf[[2]]+OTU_table_chosen_power_bf[[3]])/2)
    # }
    
    
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
    
    for ( i in 1:dim(Var_out2)[2]){
      y_vec<-c(y_vec,Var_out2[,i])
      time_vec<-c(time_vec,rep(i,length(Var_out2[,i])))
      id_vec<-c(id_vec,rownames(OTU_table_chosen_power[[1]]))
      trc_vec<-c(trc_vec,Info_IDs_SP[[i]]$Rcnts)
      genera_vec<-c(genera_vec,OTU_table_chosen_power_genera[[i]])
    }
    # lapply(list(y_vec,id_vec,trc_vec,genera_vec),length)

    OTU_table_chosen_power_Vec<-do.call("rbind",OTU_table_chosen_power)
    TMAT_res15<-CrosssecTMAT(y_vec,trc_vec,OTU_table_chosen_power_Vec,tree_chosen,TMATtype=15)$'p-value'
    TMAT_res16<-CrosssecTMAT(y_vec,trc_vec,OTU_table_chosen_power_Vec,tree_chosen,TMATtype=16)$'p-value'
    Wilcox_res<-wilcox.test(genera_vec~y_vec)$p.value
    
    MiSPU_res<-NA
    MiRKAT_res<-NA
    if(INDI_SINGLE){

    }else{
      MiRKAT_res<-try(out_MirKAT_result(OTU_table_RF_DF,tree_chosen,y_vec[ind_out_Kernel],nperm=n_perm,cov=NULL))
      
      if(class(MiRKAT_res)[1]=="try-error"){
        MiRKAT_res<-NA
      }else{
        
      }
      ind_no0<-apply(OTU_table_chosen_power_Vec,1,sum)!=0
      out_mispu <-try(MiSPU(y_vec[ind_no0],OTU_table_chosen_power_Vec[ind_no0,], tree_chosen,model = "binomial", pow = c(2:8, Inf), n.perm = n_perm,cov=NULL))
      if(class(out_mispu)[1]=="try-error"){
        
      }else{
        MiSPU_res<-out_mispu$aMiSPU$pvalue
      }
      
    }

    print(o)
    return(c(TMAT_res15,TMAT_res16,Wilcox_res,MiRKAT_res,MiSPU_res))
    
  },mc.cores=n_core)
  
  output_power[[ind_table]]<-do.call("rbind",p_vals)
}

# ind_fix<-which(substr(sapply(output_power,"[",1),1,3)=="Err")

cutoff<-0.1
value<-sapply(output_power,function(value){
  sum(value<0.1,na.rm=T)/length(!is.na(value))
})
value
mean(value)

res<-lapply(c(0.1,0.05,0.01),function(cutoff){
  value<-sapply(output_power,function(value){
    sum(value<cutoff,na.rm=T)/length(!is.na(value))
  })
  tapply(value,paste0(sn_table[,2],"_",sn_table[,3]),mean)  
})

cutoff<-0.1
Pvalue_bf<-t(sapply(output_power,function(value){
  apply(value,2,function(value2){
    sum(value2<cutoff)/sum(!is.na(value2))
  })
}))

resp<-split(data.frame(Pvalue_bf),paste0(sn_table[,5],"_",sn_table[,6]))
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
