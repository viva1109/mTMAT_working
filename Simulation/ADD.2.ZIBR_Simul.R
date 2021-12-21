library(parallel)
.libPaths("~/tmp")
library(ZIBR)
source_home<-"http://viva1109.iptime.org/"
source(paste0(source_home,"RFunctions/FunctionsTMAT/functions_tree_down.R"),encoding = "UTF-8")
source(paste0(source_home,"RFunctions/FunctionsTMAT/plot_TMAT.R"),encoding = "UTF-8")
source(paste0(source_home,"RFunctions/FunctionsTMAT/funcs_TMAT.R"),encoding = "UTF-8")
.libPaths("~/tmp")
library(mmmgee)
library(stringr)
library(ape)
library(parallel)
library(edgeR)
library(data.table)
library(NBZIMM)
source("~/analysis/20200928_longi_method/code/Functions_global.R")
TMATtype<-15
INDI_SINGLE<-T
n_core<-60
#beta dispersion 5
# v<-5
#logistic component 1
# s1_val<-1/2
#beta component 0.8
# s2_val<-1/2



# n_sim<-200
# total_result<-lapply(c(2/3,1/2,1/4,0),function(mult){
#   out_sim<-mclapply(1:n_sim,function(o){
#     sim <- simulate_zero_inflated_beta_random_effect_data(
#       subject.n=100,time.n=5,
#       X = as.matrix(c(rep(0,50*5),rep(1,50*5))),
#       Z = as.matrix(c(rep(0,50*5),rep(1,50*5))),
#       alpha = as.matrix(c(-0.5,1)*mult),
#       beta = as.matrix(c(-0.5,0.5)*mult),
#       s1 = s1_val,s2 = s2_val,
#       v = 5,
#       sim.seed=o)
#     
#     readcounts<-10000
#     
# 
#     yC<-sim$X
#     VEC_total.reads<-rep(readcounts,length(yC))
#     VEC_X_P<-round(sim$Y*readcounts)
#     IDsC<-sim$subject.ind
#     TIME<- sim$time.ind
#     
#     
#     LMM_result<-NA
#     f4result<-NA
#     mTMAT_res15<-NA
#     mTMAT_res16<-NA
#     zibr_result<-NA
#     
#     
#     zibr.fit <- zibr(logistic.cov = sim$X, beta.cov = sim$X, Y = sim$Y,
#                      subject.ind = sim$subject.ind,time.ind = sim$time.ind)
#     
#     zibr_result<-zibr.fit$joint.p    
#     
#     
#     
#     zinbDF<-data.frame(Response=VEC_X_P,CovInt=yC, ID=IDsC,TRC=VEC_total.reads)
#     names(zinbDF)[2]<-"CovInt"
#     f3 = try(glmm.zinb(Response ~ CovInt + offset(log(TRC)), random = ~ 1|ID,data=zinbDF))
#     # f3 = try(glmm.zinb(Response ~ CovInt , random = ~ 1|ID,data=zinbDF))
#     if(class(f3)[1]!="try-error"){
#       f4<-summary(f3)
#       f4result<-f4$tTable[2,'p-value']
#     }else{
#       f4result<-NA
#     }
#     
#     
#     tdata <- data.frame(Y.tran=asin(sqrt(sim$Y)),yC,SID=IDsC)
#     lme.fit <- try(lme(Y.tran ~ yC,random=~1| SID, data = tdata))
#     
#     if(class(lme.fit)!="try-error"){
#       coef.mat <- summary(lme.fit)$tTable[-1,c(1,5)]
#       LMM_result <- coef.mat[2]
#     }else{
#       LMM_result<-NA
#     }
#     
#     
#     
#     out_dataSet15<-simulSet_by_methods_fixed(y=yC,dataSet=data.frame(VEC_X_P,check.names = F),ttr=VEC_total.reads,type=15,method="TMAT")
#     out_dataSet16<-simulSet_by_methods_fixed(y=yC,dataSet=data.frame(VEC_X_P,check.names = F),ttr=VEC_total.reads,type=16,method="TMAT") 
#     
#     
#     DMsC_bf15<-TMAT_func_dm2(out_dataSet15,indi_tAnal,type=15,total.reads=total.reads,conti=F,cov=NULL)$dms
#     DMsC_bf16<-TMAT_func_dm2(out_dataSet16,indi_tAnal,type=16,total.reads=total.reads,conti=F,cov=NULL)$dms
#     
#     
#     if(INDI_SINGLE){
#       DMsC15<-list(DMsC_bf15)
#       DMsC16<-list(DMsC_bf16)
#     }else{
#       DMsC15<-DMsC_bf15
#       DMsC16<-DMsC_bf16
#     }
#     
# 
#     
#     
#     p_result<-sapply(1:length(DMsC15),function(ind_node){
#       DATA<-data.frame(response=DMsC15[[ind_node]],disease=yC,subject=IDsC,TIME=TIME)
#       names(DATA)<-c("response","disease","subject","TIME")
#       DATA_sort <- DATA[order(DATA$subject),]
#       formula <- response~disease
#       D0<-matrix(c(0,1),nrow=1)
#       m2<-geem2(response~disease,id=subject,data=DATA_sort,family=gaussian,corstr="independence")
#       res<-mmmgee.test(x=m2,L=list(D0),r=list(c(0)),statistic="score",type="quadratic",biascorr=T)
#       res$test$global$p
#     })
#     mTMAT_res15<-exact_pval(p_result)
#     
#     
#     
#     p_result<-sapply(1:length(DMsC16),function(ind_node){
#       DATA<-data.frame(response=DMsC16[[ind_node]],disease=yC,subject=IDsC,TIME=TIME)
#       names(DATA)<-c("response","disease","subject","TIME")
#       DATA_sort <- DATA[order(DATA$subject),]
#       formula <- response~disease
#       D0<-matrix(c(0,1),nrow=1)
#       m2<-geem2(response~disease,id=subject,data=DATA_sort,family=gaussian,corstr="independence")
#       res<-mmmgee.test(x=m2,L=list(D0),r=list(c(0)),statistic="score",type="quadratic",biascorr=T)
#       res$test$global$p
#     })
#     mTMAT_res16<-exact_pval(p_result)
#     
#     
# 
#     
#     
#     c(zibr.fit$joint.p,mTMAT_res15,mTMAT_res16,f4result,LMM_result)
#   },mc.cores=n_core)
#   print(mult)
#   out_sim2<-do.call("rbind",out_sim)
# })
# 
# cutoff<-0.1
# outResult<-lapply(c(0.1,0.05,0.01,0.005),function(cutoff){
#   total_result2<-t(sapply(total_result,function(values){
#     apply(values<cutoff,2,mean,na.rm=T)
#   }))
# })



set.seed(1)
sparsity<-0.7


n_sim_each<-400
z1 <- rbvunif(n_sim_each,1)
z2 <- rbvunif(n_sim_each,0.5)
z3 <- rbvunif(n_sim_each,0)

Ustart<-0.1
Uend<-1
Ulength<-Uend-Ustart
zNEW1<-z1*Ulength+Ustart
zNEW2<-z2*Ulength+Ustart
zNEW3<-z3*Ulength+Ustart 


sparsity_list<-c(0.1,0.3,0.5,0.7)

out_sim3<-vector("list",length(sparsity_list))
for(k_ind in 1:length(sparsity_list)){
  sparsity<-sparsity_list[k_ind]
  a0<-log(sparsity/(1-sparsity))
  
  compSD1<-runif(n_sim_each,0.1,1)
  compSD2<-runif(n_sim_each,0.1,1)
  Phi<-runif(n_sim_each,2,10)
  
  
  Z1_sp<-cbind(a0,zNEW1[,1],-0.5,zNEW1[,2],compSD1,compSD2,Phi,1)
  Z2_sp<-cbind(a0,zNEW2[,1],-0.5,zNEW2[,2],compSD1,compSD2,Phi,1)
  Z3_sp<-cbind(a0,zNEW3[,1],-0.5,zNEW3[,2],compSD1,compSD2,Phi,1)
  Z4_sp<-cbind(a0,0,-0.5,zNEW1[,2],compSD1,compSD2,Phi,1)
  Z5_sp<-cbind(a0,zNEW1[,1],-0.5,0,compSD1,compSD2,Phi,1)
  Z6_sp<-cbind(a0,0,-0.5,0,compSD1,compSD2,Phi,0)
  
  # exp(-1/5)/(1+exp(-1/5))
  # G1_sp<-cbind(-0.5,runif(100,0.1,1),-0.5,runif(100,0.1,1),runif(100,0.1,1),runif(100,0.1,1),runif(100,2,10),1)
  # G2_sp<-cbind(0.5,runif(100,-1,-0.1),0.5,runif(100,-1,-0.1),runif(100,0.1,1),runif(100,0.1,1),runif(100,2,10),1)
  # G3_sp<-cbind(-0.5,runif(100,0.1,1),0.5,runif(100,-1,-0.1),runif(100,0.1,1),runif(100,0.1,1),runif(100,2,10),1)
  # G4_sp<-cbind(0.5,runif(100,-1,-0.1),-0.5,runif(100,0.1,1),runif(100,0.1,1),runif(100,0.1,1),runif(100,2,10),1)
  # G5_sp<-cbind(rep(0,600),0,-0.5,0,runif(600,0.1,1),runif(600,0.1,1),runif(600,2,10),0)
  
  sim_sena<-rbind(Z1_sp,Z2_sp,Z3_sp,Z4_sp,Z5_sp,Z6_sp)
  
  
  
  n_sim<-n_sim_each*6
  
  return_seq<-function(part,length){
    (length*(part-1)+1):(length*(part))
  }
  # n_sim<-dim(sim_sena)[1]
  
  out_sim<-mclapply(1:n_sim,function(o){
    paras<-sim_sena[o,]
    sim <- simulate_zero_inflated_beta_random_effect_data(
      subject.n=100,time.n=5,
      X = as.matrix(c(rep(0,50*5),rep(1,50*5))),
      Z = as.matrix(c(rep(0,50*5),rep(1,50*5))),
      alpha = as.matrix(c(paras[1],paras[2])),
      beta = as.matrix(c(paras[3],paras[4])),
      s1 = paras[5],s2 = paras[6],
      v = paras[7],
      sim.seed=o)
    
    readcounts<-10000
    
    
    yC<-sim$X
    VEC_total.reads<-rep(readcounts,length(yC))
    VEC_X_P<-round(sim$Y*readcounts)
    IDsC<-sim$subject.ind
    TIME<- sim$time.ind
    
    
    LMM_result<-NA
    f4result<-NA
    mTMAT_res15<-NA
    mTMAT_res16<-NA
    zibr_result<-NA
    
    
    zibr.fit <- zibr(logistic.cov = sim$X, beta.cov = sim$X, Y = sim$Y,
                     subject.ind = sim$subject.ind,time.ind = sim$time.ind)
    
    zibr_result<-zibr.fit$joint.p    
    
    
    
    zinbDF<-data.frame(Response=VEC_X_P,CovInt=yC, ID=IDsC,TRC=VEC_total.reads)
    names(zinbDF)[2]<-"CovInt"
    f3 = try(glmm.zinb(Response ~ CovInt + offset(log(TRC)), random = ~ 1|ID,data=zinbDF))
    # f3 = try(glmm.zinb(Response ~ CovInt , random = ~ 1|ID,data=zinbDF))
    if(class(f3)[1]!="try-error"){
      f4<-summary(f3)
      f4result<-f4$tTable[2,'p-value']
    }else{
      f4result<-NA
    }
    
    
    tdata <- data.frame(Y.tran=asin(sqrt(sim$Y)),yC,SID=IDsC)
    names(tdata)<-c("Y.tran","yC","SID")
    lme.fit <- try(lme(Y.tran ~ yC,random=~1| SID, data = tdata))
    
    if(class(lme.fit)!="try-error"){
      coef.mat <- summary(lme.fit)$tTable[-1,c(1,5)]
      LMM_result <- coef.mat[2]
    }else{
      LMM_result<-NA
    }
    
    
    
    out_dataSet15<-simulSet_by_methods_fixed(y=yC,dataSet=data.frame(VEC_X_P,check.names = F),ttr=VEC_total.reads,type=15,method="TMAT")
    out_dataSet16<-simulSet_by_methods_fixed(y=yC,dataSet=data.frame(VEC_X_P,check.names = F),ttr=VEC_total.reads,type=16,method="TMAT") 
    
    
    DMsC_bf15<-TMAT_func_dm2(out_dataSet15,indi_tAnal,type=15,total.reads=total.reads,conti=F,cov=NULL)$dms
    DMsC_bf16<-TMAT_func_dm2(out_dataSet16,indi_tAnal,type=16,total.reads=total.reads,conti=F,cov=NULL)$dms
    
    
    if(INDI_SINGLE){
      DMsC15<-list(DMsC_bf15)
      DMsC16<-list(DMsC_bf16)
    }else{
      DMsC15<-DMsC_bf15
      DMsC16<-DMsC_bf16
    }
    
    
    
    
    p_result<-sapply(1:length(DMsC15),function(ind_node){
      DATA<-data.frame(response=DMsC15[[ind_node]],disease=yC,subject=IDsC,TIME=TIME)
      names(DATA)<-c("response","disease","subject","TIME")
      DATA_sort <- DATA[order(DATA$subject),]
      formula <- response~disease
      D0<-matrix(c(0,1),nrow=1)
      m2<-geem2(response~disease,id=subject,data=DATA_sort,family=gaussian,corstr="independence")
      res<-mmmgee.test(x=m2,L=list(D0),r=list(c(0)),statistic="score",type="quadratic",biascorr=T)
      res$test$global$p
    })
    mTMAT_res15<-exact_pval(p_result)
    
    
    
    p_result<-sapply(1:length(DMsC16),function(ind_node){
      DATA<-data.frame(response=DMsC16[[ind_node]],disease=yC,subject=IDsC,TIME=TIME)
      names(DATA)<-c("response","disease","subject","TIME")
      DATA_sort <- DATA[order(DATA$subject),]
      formula <- response~disease
      D0<-matrix(c(0,1),nrow=1)
      m2<-geem2(response~disease,id=subject,data=DATA_sort,family=gaussian,corstr="independence")
      res<-mmmgee.test(x=m2,L=list(D0),r=list(c(0)),statistic="score",type="quadratic",biascorr=T)
      res$test$global$p
    })
    mTMAT_res16<-exact_pval(p_result)
    
    
    p_result<-sapply(1:length(DMsC16),function(ind_node){
      DATA<-data.frame(response=DMsC16[[ind_node]],disease=yC,subject=IDsC,TIME=TIME)
      names(DATA)<-c("response","disease","subject","TIME")
      DATA_sort <- DATA[order(DATA$subject),]
      formula <- response~disease
      D0<-matrix(c(0,1),nrow=1)
      m2<-geem2(response~disease,id=subject,data=DATA_sort,family=gaussian,corstr="unstructured")
      res<-mmmgee.test(x=m2,L=list(D0),r=list(c(0)),statistic="score",type="quadratic",biascorr=T)
      res$test$global$p
    })
    mTMAT_res16_uc<-exact_pval(p_result)
    
    
    p_result<-sapply(1:length(DMsC16),function(ind_node){
      DATA<-data.frame(response=DMsC16[[ind_node]],disease=yC,subject=IDsC,TIME=TIME)
      names(DATA)<-c("response","disease","subject","TIME")
      DATA_sort <- DATA[order(DATA$subject),]
      formula <- response~disease
      D0<-matrix(c(0,1),nrow=1)
      m2<-geem2(response~disease,id=subject,data=DATA_sort,family=gaussian,corstr="independence")
      res<-mmmgee.test(x=m2,L=list(D0),r=list(c(0)),statistic="wald",type="quadratic",biascorr=T)
      res$test$global$p
    })
    mTMAT_res16_wald<-exact_pval(p_result)
    
    
    p_result<-sapply(1:length(DMsC16),function(ind_node){
      DATA<-data.frame(response=DMsC16[[ind_node]],disease=yC,subject=IDsC,TIME=TIME)
      names(DATA)<-c("response","disease","subject","TIME")
      DATA_sort <- DATA[order(DATA$subject),]
      formula <- response~disease
      D0<-matrix(c(0,1),nrow=1)
      m2<-geem2(response~disease,id=subject,data=DATA_sort,family=gaussian,corstr="independence")
      res<-mmmgee.test(x=m2,L=list(D0),r=list(c(0)),statistic="wald",type="quadratic",biascorr=F)
      res$test$global$p
    })
    mTMAT_res16_nocor<-exact_pval(p_result)
    
    
    if(o%%100==0){
      print(o)
    }
    
    
    c(zibr.fit$joint.p,mTMAT_res15,mTMAT_res16,mTMAT_res16_uc,mTMAT_res16_wald,mTMAT_res16_nocor,f4result,LMM_result)
  },mc.cores=n_core)
  out_sim3[[k_ind]]<-do.call("rbind",out_sim)
}




# saveRDS(out_sim3,"~/out_sim3.rds")

cutoffs<-c(0.1,0.05,0.01,0.005)
sparsity_list
true_y<-sim_sena[,8]
cutoffs<-c(0.1)
out_sim4<-lapply(1:length(sparsity_list),function(ind_sim){
  out_sim2<-out_sim3[[ind_sim]]
  resfifi<-lapply(cutoffs,function(cutoff){
    return_output<-cbind(data.frame(t(sapply(1:6,function(part){
      apply(out_sim2[return_seq(part,n_sim_each),]<cutoff,2,mean)
    }))),TYPE=c("Cor1","Cor0.5","Indep","onlyBeta","onlyAlpha","NULL"),cutoff)
    names(return_output)<-c("ZIBR","mTMAT_15","mTMAT_16","mTMAT_16_uc","mTMAT_16_wald","mTMAT_16_nocor","FZINBMM","LMM","TYPE")
    return(return_output)
  })
  resfifi2<-do.call("rbind",resfifi)
  cbind(resfifi2,1-sparsity_list[ind_sim])
})






resfifi2

# apply(out_sim2[301:400,]<cutoff,2,mean)


#https://wernerantweiler.ca/blog.php?item=2020-07-05 rbvunif


