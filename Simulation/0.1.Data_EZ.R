
library(edgeR)# data/Lon_data
fol<-"/data/sharedcode/kjkim/Analysis/20190225Combine"

source_home<-"http://viva1109.iptime.org/"
source(paste0(source_home,"RFunctions/FunctionsTMAT/functions_tree_down.R"),encoding = "UTF-8")
source(paste0(source_home,"RFunctions/FunctionsTMAT/plot_TMAT.R"),encoding = "UTF-8")
source(paste0(source_home,"RFunctions/FunctionsTMAT/funcs_TMAT.R"),encoding = "UTF-8")


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
library(stringr)
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



OTU_INFO<-list(
  otu_info_tax[[1]][[1]][otu_remiain3],
  otu_info_tax[[1]][[2]][otu_remiain3],
  sapply(str_split(otu_info_tax[[1]][[2]][otu_remiain3],"; "),"[",6),
  sapply(str_split(otu_info_tax[[1]][[2]][otu_remiain3],"; "),"[",5)
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
OTU_INFO
Info_IDs
SAMPLE_UID<-nms_go
library(ape)
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



type<-15
conti<-F
p<-1
qq<-1

cTIME<-1



list_n<-c(20,50,100,400)
list_c_prop<-c(0.2,0.5)

set.seed(1)

ind_Cpr<-1
ind_n<-1
c_prop<-list_c_prop[ind_Cpr]
n_case<-list_n[ind_n]*c_prop
n_control<-list_n[ind_n]*(1-c_prop)
base_case<-which(Var_out[,1]==1)#113
base_control<-which(Var_out[,1]==0)#363
ind_case<-sample(base_case,n_case)
ind_control<-sample(base_control,n_control)

ind_samples<-c(ind_case,ind_control)

genus_ind<-1

pvals_bytime<-sapply(1:3,function(cTIME){
  Yvar<-Var_out[ind_samples,cTIME]
  total.reads<-Info_IDs[[cTIME]]$Rcnts[ind_samples]
  out_pval<-sapply(1:length(GENUS_LIST),function(genus_ind){
    print(genus_ind)
    genus_nm<-GENUS_LIST[genus_ind]
    ind_col<-which(OTU_INFO[[3]]==genus_nm)
    
    OTU_table_chosen<-lapply(OTU_table,function(data){
      data[ind_samples,ind_col]
    })
    
    otu_id_chosen<-OTU_INFO[[1]][ind_col]
    tree_chosen<-drop.tip(pr_tree,OTU_INFO[[1]][!OTU_INFO[[1]]%in%otu_id_chosen])
    
    
    
    indi_tAnal<-list()
    indi_tAnal$subtree_tmp<-list()
    indi_tAnal$subtree_tmp[[1]]<-tree_chosen
    X_P<-OTU_table_chosen[[cTIME]]
    dataSet<-simulSet_by_methods_fixed(y=Yvar,dataSet=data.frame(X_P,check.names = F),ttr=total.reads,type=type,method="TMAT") 
    output<-TMAT_func(dataSet,indi_tAnal,type=type,total.reads=total.reads,conti=F,cov=NULL,getBeta=F,DetailResult=F)
    output$pval[1]
  })
  
  
  
})


# sapply(DM,getCores,dataSet$y,cov,conti,DetailResult,contrast=contrast)
.libPaths("~/tmp")
library(mmmgee)

##OTU_ready
genus_ind<-1

print(genus_ind)
genus_nm<-GENUS_LIST[genus_ind]
ind_col<-which(OTU_INFO[[3]]==genus_nm)

OTU_table_chosen<-lapply(OTU_table,function(data){
  data[,ind_col]
})
# library(edgeR)
CPM_table_chosen<-lapply(1:length(OTU_table_chosen),function(ind_data){
  data<-OTU_table_chosen[[ind_data]]
  cpmval<-t(cpm(t(data),log = T,lib.size = Info_IDs[[ind_data]]$Rcnts))
})
RP_table_chosen<-lapply(1:length(OTU_table_chosen),function(ind_data){
  data<-OTU_table_chosen[[ind_data]]
  RPval<-data/Info_IDs[[ind_data]]$Rcnts
})

log_table_chosen<-lapply(1:length(OTU_table_chosen),function(ind_data){
  data<-OTU_table_chosen[[ind_data]]
  RPval<-data/Info_IDs[[ind_data]]$Rcnts
  non_zero<-RPval[RPval!=0]
  log(RPval+(min(non_zero)/2))
})

otu_id_chosen<-OTU_INFO[[1]][ind_col]
tree_chosen<-drop.tip(pr_tree,OTU_INFO[[1]][!OTU_INFO[[1]]%in%otu_id_chosen])

indi_tAnal<-list()
indi_tAnal$subtree_tmp<-list()
indi_tAnal$subtree_tmp[[1]]<-tree_chosen


DMouts<-lapply(1:3,function(cTIME){
    Yvar<-Var_out[,cTIME]
    total.reads<-Info_IDs[[cTIME]]$Rcnts
    X_P<-OTU_table_chosen[[cTIME]]
    dataSet<-simulSet_by_methods_fixed(y=Yvar,dataSet=data.frame(X_P,check.names = F),ttr=total.reads,type=type,method="TMAT") 
    
    # output$pval
    dms<-TMAT_func_dm2(dataSet,indi_tAnal,type=type,total.reads=total.reads,conti=F,cov=NULL)$dms
})

DMsC<-lapply(1:length(DMouts[[1]]),function(indind){
  c(DMouts[[1]][[indind]],DMouts[[2]][[indind]],DMouts[[3]][[indind]])
})

IDsC<-rep(SAMPLE_UID,3)  
TIME<-c(rep(1,length(SAMPLE_UID)),rep(2,length(SAMPLE_UID)),rep(3,length(SAMPLE_UID)))
yC<-c(Var_out[,1],Var_out[,2],Var_out[,3])


# install.packages("geesmv")

p_result<-sapply(1:length(DMsC),function(ind_node){
  DATA<-data.frame(response=DMsC[[ind_node]],disease=yC,subject=IDsC,TIME=TIME)
  DATA_sort <- DATA[order(DATA$subject),]
  formula <- response~disease
  
  D0<-matrix(c(0,1),nrow=1)
  m2<-geem2(response~disease,id=subject,data=DATA_sort,family=gaussian,corstr="exchangeable")
  res<-mmmgee.test(x=m2,L=list(D0),r=list(c(0)),statistic="score",type="quadratic")
  # res$test$global$Chisq
  res$test$global$p
})
exact_pval<-function(tmp_pval){
  1-(1-min(tmp_pval))^length(tmp_pval)
}
mTMAT_res<-exact_pval(p_result)




OTU_table_chosen
tree_chosen
Info_IDs
Var_out

library(lme4)
OTU_table_chosen_
X_stacked<-do.call("rbind",CPM_table_chosen)
X_var<-apply(X_stacked,1,sum)
IDsC<-rep(SAMPLE_UID,3)  
TIME<-c(rep(1,length(SAMPLE_UID)),rep(2,length(SAMPLE_UID)),rep(3,length(SAMPLE_UID)))

DATA_lmer<-data.frame(response=X_var,disease=yC,subject=IDsC,TIME=TIME)
relme_H1<-lmer(X_var~yC+(1|(subject)), data=DATA_lmer,REML=F)
relme_H1_re<-lmer(X_var~yC+(1|(subject)), data=DATA_lmer,REML=T)
sm_re<-summary(relme_H1_re)
coef_res<-coef(sm_re)
re_LR<-drop1(relme_H1,test="Chisq")
unlist(re_LR["yC",4])
       
library(devtools)
install_github("hk1785/GLMM-MiRKAT", force=T)

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

data(bin.phy)
set.seed(100)
bin.phy <- rarefy_even_depth(bin.phy, rngseed=TRUE)

otu.tab <- otu_table(bin.phy)
tree <- phy_tree(bin.phy)
meta <- sample_data(bin.phy)
y <- meta$y
id <- meta$id
x1 <- meta$x1
x2 <- meta$x2
covs <- as.data.frame(cbind(x1,x2))
covs[,2] <- as.factor(covs[,2]) 

Ks <- Kernels(otu.tab, tree)

result_MIRKAT<-GLMMMiRKAT(y, cov=covs, id=id, Ks=Ks, model="binomial",n.perm=1000)
result_MIRKAT$aGLMMMiRKAT
