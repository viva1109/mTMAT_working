# n_sim_bf<-100
# DVD<-1
# n_core=30
# 
# 
# HeatSTR<-paste0("DA_210328_17D",DVD)
# part<-1
# jump<-13
# Done<-0
# Start<-1+jump*(part-1)+Done
# End<-min(jump*part+Done,73)
# cat(Start,End,"\n")
# 
# n_sim<-n_sim_bf/DVD
# n_perm<-5000/DVD

str<-"independence"
# str<-"unstructured"
# str<-"exchangeable"
# str<-"ar1"
# str<-"m-dependent"

stat<-"score"
# type=16


.libPaths("~/tmp")

library(mmmgee)
library(stringr)
library(ape)
library(parallel)
library(data.table)
library(edgeR)
source_home<-"http://viva1109.iptime.org/"
source(paste0(source_home,"RFunctions/FunctionsTMAT/functions_tree_down.R"),encoding = "UTF-8")
source(paste0(source_home,"RFunctions/FunctionsTMAT/plot_TMAT.R"),encoding = "UTF-8")
source(paste0(source_home,"RFunctions/FunctionsTMAT/funcs_TMAT.R"),encoding = "UTF-8")
source("~/analysis/20200928_longi_method/code/Functions_global.R")

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("microbiomeDASim")


library(microbiomeDASim)
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
library(microbiomeDASim)



infofol<-"~/analysis/20200928_longi_method/info/"
FAM_info<-readRDS(file.path(infofol,"FAMILY_INFO_DASIM.csv"))


FAMILY_LIST<-FAM_info$GENUS_LIST
OTU_no<-FAM_info$OTU_no
OTU_INFO<-FAM_info$OTU_INFO
pr_tree<-FAM_info$pr_tree
sparsity<-FAM_info$sparsities_List
mean_List<-FAM_info$mean_List
sd_List<-FAM_info$sd_List

Target_otu<-vector("numeric",3)
for ( i in 1:3){
  valtmp<-summary(sapply(sparsity,mean))[c(2,3,5)][i]
  Target_otu[i]<-which.min((sapply(sparsity,mean)-valtmp)^2)
}

ind_togo<-Target_otu

NN<-100
timepoints<-3

ind_togo<-13




str_list<-c("compound","ar1","ind")
rho_list<-c(0.2,0.5,0.8)
# beta_list<-c(0,0.1,0.2)
beta_list<-c(0)
FAMILY_LIST2<-FAMILY_LIST[ind_togo]
sn_table<-matrix(nrow=length(str_list)*length(rho_list)*length(beta_list)*length(FAMILY_LIST2),ncol=5)
cnt<-0
sn_table[,1]<-1:dim(sn_table)[1]

for (j in 1:length(str_list)){
  for ( i in 1:length(rho_list)){
    for ( z in 1:length(FAMILY_LIST2)){
      for ( z2 in 1:length(beta_list)){
        cnt<-cnt+1
        sn_table[cnt,2]<-j
        sn_table[cnt,3]<-i
        sn_table[cnt,4]<-z
        sn_table[cnt,5]<-z2
      }
      
    }
  }
}
list_labels<-list(str_list,rho_list,beta_list,FAMILY_LIST2)
sn_table2<-sn_table[!(sn_table[,2]==3 & sn_table[,3]%in%2:3),]

OUTPUTpheno<-vector("list",length(ind_togo))
OUTPUTdata<-vector("list",length(ind_togo))
OUTPUTdataOther<-vector("list",length(ind_togo))
Tree_list<-vector("list",length(ind_togo))


n_sim<-4000

outputPvalues<-vector("list",dim(sn_table2)[1])

for ( ind_table in 1:dim(sn_table2)[1]){
  i<-ind_togo[sn_table2[ind_table,4]]
  str_i<-str_list[sn_table2[ind_table,2]]
  stat_i<-"score"
  rho_list_i<-rho_list[sn_table2[ind_table,3]]
  beta_i<-beta_list[sn_table2[ind_table,5]]
  
  family_nm<-FAMILY_LIST[i]
  ind_col<-which(OTU_INFO[[3]]==family_nm)
  otu_id_chosen<-OTU_INFO[[1]][ind_col]
  tree_chosen_bf<-drop.tip(pr_tree,OTU_INFO[[1]][!OTU_INFO[[1]]%in%otu_id_chosen])
  
  if(!is.rooted(tree_chosen_bf)){
    tree_chosen<-root(tree_chosen_bf, 1, r = TRUE)
  }else{
    tree_chosen<-tree_chosen_bf
  }
  Tree_list[[cnt]]<-tree_chosen
    zeroprop<-mean(sparsity[[i]])

    meanval<-mean(mean_List[[i]])
    sdval<-mean(sd_List[[i]])
    
    TMATpvalues<-mclapply(1:n_sim,function(o){
      set.seed(o)
      mi_gen <- gen_norm_microbiome(features=length(ind_col)+1, diff_abun_features=length(ind_col), 
                                    n_control=50,n_treat=50, control_mean=meanval, sigma=sdval,
                                    num_timepoints=timepoints, rho=rho_list_i, 
                                    corr_str=str_i,  func_form="linear", 
                                    beta=c(beta_i*sdval, 0),  missing_pct=zeroprop, 
                                    missing_per_subject=2,  miss_val=0)
      Mat_bf<-t(mi_gen[[1]])
      Mattmp<-Mat_bf[,1:length(ind_col)]
      colnames(Mattmp)<-otu_id_chosen
      
      outdata<-list(X_P=log(exp(Mattmp)+1),X_P_comp=matrix(log(exp(Mat_bf[,length(ind_col)+1])+1),ncol=1))
      
      indi_tAnal<-list()
      indi_tAnal$subtree_tmp<-list()
      indi_tAnal$subtree_tmp[[1]]<-tree_chosen
      
      Time<-mi_gen[[2]][,2]
      Ys<-mi_gen[[2]][,3]
      Y2<-vector("numeric",length(Ys))
      Y2[Ys!="Treatment"]<-0
      Y2[Ys=="Treatment"]<-1
      ID<-mi_gen[[2]][,1]
      
      dataSet<-list(y=Y2,simData=outdata,ttr=1000000)
      
      DMsC<-TMAT_func_dm2(dataSet,indi_tAnal,type=16,total.reads=1000000,conti=F,cov=NULL)$dms
      
      
      p_result<-sapply(1:length(DMsC),function(ind_node){
        DATA<-data.frame(response=DMsC[[ind_node]],disease=Y2,subject=ID,TIME=Time)
        DATA_sort <- DATA[order(DATA$subject),]
        formula <- response~disease
        D0<-matrix(c(0,1),nrow=1)
        m2<-geem2(response~disease,id=subject,data=DATA_sort,family=gaussian,corstr=str)
        res<-mmmgee.test(x=m2,L=list(D0),r=list(c(0)),statistic="score",type="quadratic",biascorr=T)
        res$test$global$p
      })
      
      mTMAT_res<-exact_pval(p_result)
    },mc.cores=80)
    
    TMATpvalues2<-do.call("c",TMATpvalues)
    
  # OUTPUTdata[[cnt]]<-matrix(t(mi_gen[[1]])[,1:length(ind_col)],ncol=length(ind_col))
  # colnames(OUTPUTdata[[cnt]])<-otu_id_chosen
  # OUTPUTdataOther[[cnt]]<-matrix(t(mi_gen[[1]])[,length(ind_col)+1],ncol=length(ind_col))
    outputPvalues[[ind_table]]<-TMATpvalues2
    print(ind_table)
}


outputPvalues
list_cutoff<-c(0.1,0.05,0.01,0.005)
outMat<-sapply(list_cutoff,function(cutoff){
  sapply(outputPvalues,function(value){
    sum(value<cutoff)/length(value)
  })
})

sn_table



quantile(OTU_no,0.9)
summary(OTU_no)



