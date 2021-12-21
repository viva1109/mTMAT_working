.libPaths("~/tmp")
library(phyloseq)
library(NBZIMM)
library(nlme)
source_home<-"http://viva1109.iptime.org/"
source(paste0(source_home,"RFunctions/FunctionsTMAT/functions_tree_down.R"),encoding = "UTF-8")
source(paste0(source_home,"RFunctions/FunctionsTMAT/plot_TMAT.R"),encoding = "UTF-8")
source(paste0(source_home,"RFunctions/FunctionsTMAT/funcs_TMAT.R"),encoding = "UTF-8")
library(stringr)
library(raster, quietly = TRUE)
suppressWarnings(suppressMessages(library(ape, quietly = TRUE)))
library(edgeR, quietly = TRUE)
library(ggplot2)
library(grid)
library(edgeR)
library(parallel)
setwd("/home2/kjkim/analysis/20200928_longi_method/RealData/data")
clinical_bf <- read.csv("temporal_spatial_clinical_merged.csv",stringsAsFactors=F)
otu_bf <- read.csv("temporal_spatial_otu.csv", check.names = F,stringsAsFactors=F)
taxonomy <- read.csv("temporal_spatial_taxonomy.csv",stringsAsFactors=F)
# clusters <- read.csv("temporal_clust.csv")

# clinical <- merge(clinical_bf, clusters, by.x = "X", by.y = "X")
########### otu data #################
otu_bf[,1]
dim(otu_bf)
dim(clinical_bf)
colnames(otu_bf)
rownames((otu_bf))
clinical_bf$SampleID==otu_bf[,1]
# spl_remains<-clinical_bf$NumReads>=3000
spl_remains<-rep(T,length(clinical_bf$NumReads))
otu_bf2 <- otu_bf[spl_remains, ]
clinical<-clinical_bf[spl_remains,]
  
library(GUniFrac)
rf_table<-Rarefy(otu_bf2[,-1])

library(ape)
library(GUniFrac)
otu2<-otu_bf2[,-1]
library(vegan)
colnames(rf_table$otu.tab.rff)<-colnames(otu_bf2[,-1])
rftable_filtered<-rf_table$otu.tab.rff[,apply(rf_table$otu.tab.rff,2,sum)!=0]

braydist<-vegdist(t(rftable_filtered), method="bray")


# dotu<-as.dist(otu)

arbol <- nj(as.dist(braydist))
# braydist2<-as.dist(braydist)

# arbol<-readRDS("nj.rds")
if(!is.rooted(arbol)){
  tree_ori<-root(arbol, 1, r = TRUE)
}else{
  tree_ori<-arbol
}


# saveRDS(braydist,"braydist.rds")
# saveRDS(arbol,"nj.rds")


############### clinical data ############################
dim(clinical)
colnames(clinical)

sites = clinical$BodySite  # four different body sites
table(sites)
clinical <- clinical[clinical$BodySite == "Vaginal_Swab", ]

# braydist


# Samples collected before delivery
Preg <- clinical$Preg
table(sites, Preg)
PrePreg <- clinical$PrePreg
table(Preg, PrePreg)
ind_out<-Preg==F&PrePreg==T
# Samples collected after delivery
Post <- clinical$PostPreg
table(Preg, Post)
clinical2<-clinical[!ind_out,]
#### keep only samples collected during pregnancy
##### Keep only samples in vaginal swab ####
# clinical <- clinical[clinical$Outcome != "Marginal", ] ### Same as in the paper
# clinical <- clinical[clinical$V1 == "4", ] ### Same as in the paper
# 
# clinical <- clinical[clinical$Preg == "TRUE", ]
# clinical <- clinical[clinical$BodySite == "Vaginal_Swab", ] ### To compare with CST analysis in the paper



library(pldist)

sp_clinical<-split(clinical2,clinical2$SubjectID)
sp_clinical_1<-sp_clinical[sapply(sp_clinical,dim)[1,]!=1]
sp_clinical2<-lapply(sp_clinical_1,function(data){
  data2<-cbind(data[order(data$GDColl),],1:dim(data)[1])
  names(data2)<-c(names(data),"Time")
  return(data2)
})
cbd_clinical<-do.call("rbind",sp_clinical2)


ind_forOTU<-match(cbd_clinical$SampleID,otu_bf2$Sample_ID)

SampleIDs_mat<-otu_bf2$Sample_ID[ind_forOTU]
otu_togo<-otu2[ind_forOTU,]
rownames(otu_togo)<-SampleIDs_mat
rf_table2<-rf_table$otu.tab.rff[ind_forOTU,]
rownames(rf_table2)<-SampleIDs_mat
paired.meta2<-cbd_clinical[,c("SubjectID","SampleID","Time")]
names(paired.meta2)<-c("subjID","sampID","time")




paired.bray= pldist(rf_table2, paired.meta2, paired = F, binary = FALSE, method = "bray")$D
paired.jaccard = pldist(rf_table2, paired.meta2, paired = F, binary = FALSE, method = "jaccard")$D
paired.kulczynski = pldist(rf_table2, paired.meta2, paired = F, binary = FALSE, method = "kulczynski")$D
paired.gower = pldist(rf_table2, paired.meta2, paired = F, binary = FALSE, method = "gower")$D

cbd_clinical$Time==1

combined<-data.frame(cbind(cbd_clinical[,c("GDColl","Birthweight_kg","Length_at_birth","Ethnicity","Gender_Baby","Race","Labor_Initiation","Indication","Operator","SES_maternal_education","SES_household_income","Hypertensive.Disorder","Preg","CSection","Outcome","Term")]))
Cov<-combined[cbd_clinical$Time==1,]
outDF<-data.frame(matrix(nrow=dim(Cov)[2],ncol=3))
for( i in 1:dim(Cov)[2]){
  Target_Var<-Cov[,i]
  outDF[i,1]<-colnames(Cov)[i]
  ind_remain<-!is.na(Target_Var)
  paired.bray2<-as.dist(paired.bray[ind_remain,ind_remain])
  if(length(unique(Target_Var[ind_remain]))==1){
    print(paste0(colnames(Cov)[i]," has only one group level"))
    outDF[i,2]<-NA
    outDF[i,3]<-NA
  }else{
    permanova <- adonis(paired.bray2 ~ Target_Var[ind_remain], permutations=10000)
    outDF[i,2]<-permanova$aov.tab$R2[1]
    outDF[i,3]<-permanova$aov.tab$Pr[1]
  }
  
  cat(i)
  cat(" ")
}

bf<-outDF
bf$X1<-factor(bf$X1,bf$X1[order(bf$X2)])
ggplot(bf, aes(x=X1, y=X2)) +
  geom_bar(stat="identity", col="black") +
  # geom_point(stat="identity",stroke=0.001,shape=8,size=2.5,aes(alpha=c(0,1)[bf$FDR_factor]),show.legend = F) +
  ylab(expression("r"^2))+
  xlab("Phenotype")+
  theme(
    axis.text.y = element_text(),
    axis.text.x = element_text())+
  coord_flip()+
  theme_classic(base_size=13)+
  #theme options
  theme(
    strip.background = element_rect(fill="white"),
    axis.text.y = element_text(face="bold",size=11),
    axis.text.x = element_text(face="bold", size=11),
    #bold font for both axis text
    axis.text=element_text(face="bold"),
    #set thickness of axis ticks
    axis.ticks=element_line(size=0.4),
    #remove plot background
    plot.background=element_blank(),
    #remove plot border
    panel.border=element_blank())+
  ylim(c(0,max(bf$X2)))
ggsave(paste0("output.png"),width=5,height=8,dpi=300)

system("scp output.png ng:~")



rf_table3<-rf_table2[cbd_clinical$Time==1,]
braydist<-vegdist(rf_table3, method="bray")



Cov<-combined[cbd_clinical$Time==1,]
outDF<-data.frame(matrix(nrow=dim(Cov)[2],ncol=3))
for( i in 1:dim(Cov)[2]){
  Target_Var<-Cov[,i]
  outDF[i,1]<-colnames(Cov)[i]
  ind_remain<-!is.na(Target_Var)
  paired.bray2<-as.dist(braydist)
  if(length(unique(Target_Var[ind_remain]))==1){
    print(paste0(colnames(Cov)[i]," has only one group level"))
    outDF[i,2]<-NA
    outDF[i,3]<-NA
  }else{
    permanova <- adonis(paired.bray2 ~ Target_Var[ind_remain], permutations=10000)
    outDF[i,2]<-permanova$aov.tab$R2[1]
    outDF[i,3]<-permanova$aov.tab$Pr[1]
  }
  
  
  cat(i)
  cat(" ")
}

bf<-outDF
bf$X1<-factor(bf$X1,bf$X1[order(bf$X2)])
ggplot(bf, aes(x=X1, y=X2)) +
  geom_bar(stat="identity", col="black") +
  # geom_point(stat="identity",stroke=0.001,shape=8,size=2.5,aes(alpha=c(0,1)[bf$FDR_factor]),show.legend = F) +
  ylab(expression("r"^2))+
  xlab("Phenotype")+
  theme(
    axis.text.y = element_text(),
    axis.text.x = element_text())+
  coord_flip()+
  theme_classic(base_size=13)+
  #theme options
  theme(
    strip.background = element_rect(fill="white"),
    axis.text.y = element_text(face="bold",size=11),
    axis.text.x = element_text(face="bold", size=11),
    #bold font for both axis text
    axis.text=element_text(face="bold"),
    #set thickness of axis ticks
    axis.ticks=element_line(size=0.4),
    #remove plot background
    plot.background=element_blank(),
    #remove plot border
    panel.border=element_blank())+
  ylim(c(0,max(bf$X2)))
ggsave(paste0("output_Time1.png"),width=5,height=8,dpi=300)

system("scp output_Time1.png ng:~")




rf_table3<-rf_table2[cbd_clinical$Time==2,]
braydist<-vegdist(rf_table3, method="bray")



Cov<-combined[cbd_clinical$Time==2,]
outDF<-data.frame(matrix(nrow=dim(Cov)[2],ncol=3))
for( i in 1:dim(Cov)[2]){
  Target_Var<-Cov[,i]
  outDF[i,1]<-colnames(Cov)[i]
  ind_remain<-!is.na(Target_Var)
  paired.bray2<-as.dist(braydist)
  if(length(unique(Target_Var[ind_remain]))==1){
    print(paste0(colnames(Cov)[i]," has only one group level"))
    outDF[i,2]<-NA
    outDF[i,3]<-NA
  }else{
    permanova <- adonis(paired.bray2 ~ Target_Var[ind_remain], permutations=10000)
    outDF[i,2]<-permanova$aov.tab$R2[1]
    outDF[i,3]<-permanova$aov.tab$Pr[1]
  }
  
  
  cat(i)
  cat(" ")
}

bf<-outDF
bf$X1<-factor(bf$X1,bf$X1[order(bf$X2)])
ggplot(bf, aes(x=X1, y=X2)) +
  geom_bar(stat="identity", col="black") +
  # geom_point(stat="identity",stroke=0.001,shape=8,size=2.5,aes(alpha=c(0,1)[bf$FDR_factor]),show.legend = F) +
  ylab(expression("r"^2))+
  xlab("Phenotype")+
  theme(
    axis.text.y = element_text(),
    axis.text.x = element_text())+
  coord_flip()+
  theme_classic(base_size=13)+
  #theme options
  theme(
    strip.background = element_rect(fill="white"),
    axis.text.y = element_text(face="bold",size=11),
    axis.text.x = element_text(face="bold", size=11),
    #bold font for both axis text
    axis.text=element_text(face="bold"),
    #set thickness of axis ticks
    axis.ticks=element_line(size=0.4),
    #remove plot background
    plot.background=element_blank(),
    #remove plot border
    panel.border=element_blank())+
  ylim(c(0,max(bf$X2)))
ggsave(paste0("output_Time2.png"),width=5,height=8,dpi=300)

system("scp output_Time2.png ng:~")



yy = as.matrix(otu_togo[,-1])
# yy = ifelse(is.na(yy), 0, yy)

zero.p = apply(yy, 2, function(x) {length(x[x != 0])/length(x)} )
zero.p = sort(zero.p, decreasing = T)
zero.p = data.frame(zero.p)
zero.p$id = rownames(zero.p)
zero.p = data.frame(zero.p[zero.p$zero.p>0.25, ]) ## same as the proportion of non zero used in the paper for CST = 4

tree_chosen<-drop.tip(tree_ori,setdiff(tree_ori$tip.label,rownames(zero.p)))


yy = yy[, rownames(zero.p)]
rownames(yy) = otu$Sample_ID

sampleID = cbd_clinical$SampleID # different samples for one subject
SubjectID = as.factor(cbd_clinical$SubjectID)
outcome = cbd_clinical$Outcome # preterm vs term groups
table(outcome)
term_label = ifelse(outcome == "Term", 0, 1)
preg_label<-ifelse(cbd_clinical$Preg,0,1)

weeks = as.numeric(log(cbd_clinical$GDColl)) # Days for samples collected
summary(weeks)
days<-cbd_clinical$GDColl
N = cbd_clinical$NumReads  # total reads
mean(N); sd(N)
mean(log(N)); sd(log(N))

##########################
## LMM model
### compare random intercept model with random slope (dropping the correlation between random intercept and slope) model
f1 = f2 = f3 = f4 = f5 = f6 = f7 = f8 = list()
out1 = out2 = out3 = out4 = out5 = out6 = out7 = out8 = out9 = out10 = out11 = out12 = matrix(NA, ncol(yy), 1)

for (j in 1:ncol(yy)){
  tryCatch({
    y = yy[, j]
    y0 = asin(sqrt(y/N))
    f1 = lme(y0 ~ preg_label, random = ~ 1|SubjectID)
    out1[j, ] = summary(f1)$tTable[2, 5]
    
    f2 = lme(y0 ~ preg_label, random = list(SubjectID = pdDiag(~weeks)))
    out3[j, ] = summary(f2)$tTable[2, 5]
    
    f3 = lme(y0 ~ preg_label*weeks, random = ~ 1|SubjectID)
    out5[j, ] = summary(f3)$tTable[2, 5]
    out7[j, ] = summary(f3)$tTable[4, 5]
    
    f4 = lme(y0 ~ preg_label*weeks, random = list(SubjectID = pdDiag(~weeks)))
    out9[j, ] = summary(f4)$tTable[2, 5]
    out11[j, ] = summary(f4)$tTable[4, 5]
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
out_lmm <- cbind(out1, out3, out5, out7, out9, out11)
apply(out_lmm,2,p.adjust,method="fdr")<0.05

##########################
## NBMM model
### compare random intercept model with random slope (dropping the correlation between random intercept and slope) model
for (j in 1:ncol(yy)){
  tryCatch({
    y = as.numeric(yy[, j])
    #y0 = log2(y + 1)#asin(sqrt(y/N))
    f5 = glmm.nb(y ~ preg_label + offset(log(N)), random = ~ 1|SubjectID)
    out2[j, ] = summary(f5)$tTable[2, 5]
    #
    f6 = glmm.nb(y ~ preg_label + offset(log(N)), random = list(SubjectID = pdDiag(~weeks)))
    out4[j, ] = summary(f6)$tTable[2, 5]
    
    f7 = glmm.nb(y ~ preg_label*weeks + offset(log(N)), random = ~ 1|SubjectID)
    out6[j, ] = summary(f7)$tTable[2, 5]
    out8[j, ] = summary(f7)$tTable[4, 5]
    #
    f8 = glmm.nb(y ~ preg_label*weeks + offset(log(N)), random = list(SubjectID = pdDiag(~weeks)))
    out10[j, ] = summary(f8)$tTable[2, 5]
    out12[j, ] = summary(f8)$tTable[4, 5]
    
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
out_nbmm <- cbind(out2, out4, out6, out8, out10, out12)

out <- cbind(out1, out2, out3,  out4, out5, out6, out7, out8, out9, out10, out11, out12)
rownames(out) = colnames(yy)
est = out2[complete.cases(out2), 1]
est = est
res = sim.out(coefs.p = t(out), coefs.est = t(est), alpha = c(0.05, 0.01, 0.005, 0.001))
res

p.adjust(out2[,1],method="fdr")

res2 = sim.out(coefs.p = t(p.adjust), coefs.est = t(est), alpha = c(0.05, 0.01, 0.005, 0.001))
res2

cbind(est,row.names(zero.p))


# 
# SubjectID
# 
# SubjectID
library(mmmgee)
library(stringr)
library(ape)
library(parallel)
library(data.table)
library(edgeR)# data/Lon_data
source_home<-"http://viva1109.iptime.org/"
source(paste0(source_home,"RFunctions/FunctionsTMAT/functions_tree_down.R"),encoding = "UTF-8")
source(paste0(source_home,"RFunctions/FunctionsTMAT/plot_TMAT.R"),encoding = "UTF-8")
source(paste0(source_home,"RFunctions/FunctionsTMAT/funcs_TMAT.R"),encoding = "UTF-8")
OTU_id_togo<-row.names(zero.p)
Taxo_togo<-taxonomy[match(row.names(zero.p),taxonomy[,1]),]
cnt_genus<-table(Taxo_togo$Genus)
genus_togo1<-names(sort(cnt_genus,decreasing = T))
which(Taxo_togo$Genus==u_Taxo_togo[1])




# data_TMAT<-lapply(1:length(genus_togo1),function(j){
#   target_genus<-genus_togo1[j]
#   col_ind<-which(Taxo_togo$Genus==target_genus)
#   Yvar<-cbd_clinical$Preg
#   total.reads<-N
#   X_P<-matrix(yy[, col_ind],ncol=length(col_ind))
#   colnames(X_P)<-colnames(yy)[col_ind]
#   # tree_chosen
#   
#   
#   if(is.binary.tree(tree_chosen)){
#     tree_chosen_1<-tree_chosen
#   }else{
#     tree_chosen_1<-multi2di(tree_chosen)
#   }
#   
#   tree_chosen2<-drop.tip(tree_chosen_1,setdiff(tree_chosen_1$tip.label,colnames(X_P)))
#   
#   
#   
#   indi_tAnal<-list()
#   indi_tAnal$subtree_tmp<-list()
#   
#   if(!is.rooted(tree_chosen2)){
#     indi_tAnal$subtree_tmp[[1]]<-root(tree_chosen2, 1, r = TRUE)
#   }else{
#     indi_tAnal$subtree_tmp[[1]]<-tree_chosen2
#   }
#   
#   dataSet<-simulSet_by_methods_fixed(y=Yvar,dataSet=data.frame(X_P,check.names = F),ttr=total.reads,type=16,method="TMAT") 
#   
# })


# TMAT_pval_bygenus_bf<-mclapply(1:length(data_TMAT),function(ind_target){
#   target<-data_TMAT[[ind_target]]
#   pruned.tree_perphy_omiat<-drop.tip(tree_chosen,tree_chosen$tip.label[!tree_chosen$tip.label%in%colnames(target$simData$X_P)])
#   indi_tAnal<-list()
#   indi_tAnal$subtree_tmp<-list()
#   indi_tAnal$subtree_tmp[[1]]<-pruned.tree_perphy_omiat
#   TMAT_func(target,indi_tAnal,type=16,total.reads=Data_Final$totalcounts,conti=F,cov=NULL,getBeta=T,DetailResult=T)
# },mc.cores = 80)


outTMATp<-mclapply(1:length(genus_togo1),function(j){
  target_genus<-genus_togo1[j]
  col_ind<-which(Taxo_togo$Genus==target_genus)
  Yvar<-cbd_clinical$Preg
  total.reads<-N
  X_P<-matrix(yy[, col_ind],ncol=length(col_ind))
  colnames(X_P)<-colnames(yy)[col_ind]
  # tree_chosen
  
  
  if(is.binary.tree(tree_chosen)){
    tree_chosen_1<-tree_chosen
  }else{
    tree_chosen_1<-multi2di(tree_chosen)
  }
 
  tree_chosen2<-drop.tip(tree_chosen_1,setdiff(tree_chosen_1$tip.label,colnames(X_P)))
  
  
  
  indi_tAnal<-list()
  indi_tAnal$subtree_tmp<-list()
  
  if(!is.rooted(tree_chosen2)){
    indi_tAnal$subtree_tmp[[1]]<-root(tree_chosen2, 1, r = TRUE)
  }else{
    indi_tAnal$subtree_tmp[[1]]<-tree_chosen2
  }
  
  dataSet<-simulSet_by_methods_fixed(y=Yvar,dataSet=data.frame(X_P,check.names = F),ttr=total.reads,type=16,method="TMAT") 
  
  # TMAT_pval_bygenus_bf<-TMAT_func(dataSet,indi_tAnal,type=16,total.reads=total.reads,conti=F,cov=NULL,getBeta=T,DetailResult=T)
  
  # output$pval
  if(dim(dataSet[[2]][[1]])[2]==1){
    DMsC<-TMAT_func_dm2(dataSet,indi_tAnal,type=16,total.reads=total.reads,conti=F,cov=NULL,getBeta=T,DetailResult=F)
    DMsC2<-DMsC
    DMsC2$dms<-list(DMsC2$dms)
  }else{
    DMsC<-TMAT_func_dm2(dataSet,indi_tAnal,type=16,total.reads=total.reads,conti=F,cov=NULL,getBeta=T,DetailResult=F)
    DMsC2<-DMsC
  }
  
  # str<-"unstructured"
  str<-"independence"
  # str<-"exchangeable"
  # str<-"ar1"
  # str<-"m-dependent"

  
  p_result<-t(sapply(1:length(DMsC2$dms),function(ind_node){
    DATA<-data.frame(response=DMsC2$dms[[ind_node]],disease=Yvar,subject=SubjectID,weeks=weeks)
    DATA_sort <- DATA[order(DATA$subject),]
    formula <- response~disease+weeks
    D0<-matrix(c(0,1,0),nrow=1)
    # D1<-matrix(c(0,0,0,1),nrow=1)
    m2<-geem2(formula,id=subject,data=DATA_sort,family=gaussian,corstr=str)
    res1<-mmmgee.test(x=m2,L=list(D0),r=list(c(0)),statistic="score",type="quadratic",biascorr=T)
    # res2<-mmmgee.test(x=m2,L=list(D1),r=list(c(0)),statistic="score",type="quadratic",biascorr=T)
    c(res1$test$global$p,m2$beta[2])
    # c(res1$test$global$p,res2$test$global$p)
  }))
  print(j)
  mTMAT_res<-exact_pval(p_result[,1])
  list(mTMAT_res,pvals_betas=cbind(p_result[,1],p_result[,2]),output=DMsC,dataSet=dataSet)
},mc.cores=80)
outTMATp2<-lapply(outTMATp,function(data){
  data[[1]]
})

outTMATbeta<-lapply(outTMATp,function(data){
  data[[2]]
})
outTMATpoutput<-lapply(outTMATp,function(data){
  data[[3]]
})

outTMATdataset<-lapply(outTMATp,function(data){
  data[[4]]
})



outTMATp2<-do.call("rbind",outTMATp2)
adj_Ps1<-p.adjust(outTMATp2[,1],method="fdr")
# adj_Ps2<-p.adjust(outTMATp2[,2],method="fdr")
outTB<-cbind(data.frame(genus_togo1,stringsAsFactors=F),outTMATp2,adj_Ps1,stringsAsFactors=F)
ind_sig<-which(outTB$adj_Ps1<0.05)
outtogo<-outTB[ind_sig,]

outTMATp[adj_Ps1<0.05]

ind_tmp<-1
ind_direct_genera<-ind_sig[ind_tmp]
outTB

# plot_TMAT(outTB[ind_sig[ind_tmp],1],target_ori=Data_Final,data_TMAT=data_TMAT,TMAT_pval_bygenus_bf=TMAT_pval_bygenus_bf,conti=conti)



dev.off()
png("TMATPlot0406.png" ,width=1200*1.5,height=700*1.5)

x_axis_label<-"log of days from first visit"
genus_index<-NULL
target_genus<-outTB[ind_direct_genera,1]
output<-outTMATpoutput[[ind_direct_genera]]
target<-outTMATdataset[[ind_direct_genera]]

cexSize<-7
nodeCex<-1.09
# browser()
type<-16
# notready<-target_ori$tree





# margins<-c(0.2,2,2,2)
margins<-c(0.2,1.5,1.5,1.5)
# par(mfrow=c(1,1), mar=margins, oma=c(0, 4, 0, 4))
par(mfrow=c(1,1), mar=margins, oma=c(0, 0, 0, 0))
# layout(matrix(c(1,1,2,2), 2, 2, byrow = TRUE),widths = lcm(c(8, 4)),heights = lcm(c(4, 8)))
# layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE),widths = c(3, 2),heights = c(3, 2))
layout(matrix(c(1,2,3,4,5,6), 2, 3, byrow = TRUE),widths = c(2,10, 8),heights = c(3,2))


col_ind<-which(Taxo_togo$Genus==target_genus)
X_P<-matrix(yy[, col_ind],ncol=length(col_ind))
colnames(X_P)<-colnames(yy)[col_ind]

contrast<-c(1,2)
DetailResult<-F




if(is.binary.tree(tree_chosen)){
  tree_chosen_1<-tree_chosen
}else{
  tree_chosen_1<-multi2di(tree_chosen)
}

tree_chosen2_bf<-drop.tip(tree_chosen_1,setdiff(tree_chosen_1$tip.label,colnames(X_P)))
tree_chosen2<-tree_chosen2_bf
tree_chosen2$edge.length[tree_chosen2$edge.length<0.1]<-runif(1,0.07,0.12)
tiplabels_ready<-Taxo_togo[match(tree_chosen2$tip.label,OTU_id_togo),"Species"]
ind_uncultered<-str_detect(tiplabels_ready,"uncultured")
tiplabels_ready[ind_uncultered]<- paste0(tiplabels_ready," (",tree_chosen2$tip.label,")")[ind_uncultered]
tiplabels_ready<-c("Other Genera",tiplabels_ready)
tiplabels_ready<-c(tiplabels_ready,paste0("[Genus] ",target_genus))


indi_tAnal<-list()
indi_tAnal$subtree_tmp<-list()

if(!is.rooted(tree_chosen2)){
  indi_tAnal$subtree_tmp[[1]]<-root(tree_chosen2, 1, r = TRUE)
}else{
  indi_tAnal$subtree_tmp[[1]]<-tree_chosen2
}

betah<-outTMATbeta[[ind_direct_genera]][,2]
vvaall<-exp(betah)
colfunc <- colorRampPalette(c("blue","white","red"))
ind_col<-ceiling(vvaall*199.5-100)
ind_col[ind_col>200]<-200
ind_col[ind_col<0]<-1
col_picked<-colfunc(200)
col_nodes<-col_picked[ind_col]
col_nodes2<-col_nodes[c(length(col_nodes),1:(length(col_nodes)-1))]

x <- rtree(2, tip.label = LETTERS[1:2])
x$edge.length[1]<-0.02
x$edge.length[2]<-0
tree.Groups_toplot<-bind.tree(x,tree_chosen2,where = 1)


n_tips<-Ntip(tree.Groups_toplot)

xl <- 1
yb <- 1
xr <- 1.5
yt <- 2
ncol_picked<-20
ncol_texts<-10
col_picked<-colfunc(ncol_picked)
# ?plot
plot(NA,type="n",xlim=c(1,2),ylim=c(1,2),xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
mtext(expression(exp(hat(beta))), side=3, adj=0, line=-1, cex=0.8, font=2);
# mtext("exp(beta)", side=4, adj=0, line=0, cex=0.8, font=2);
rect(
  xl,
  head(seq(yb,yt,(yt-yb)/ncol_picked),-1),
  xr,
  tail(seq(yb,yt,(yt-yb)/ncol_picked),-1),
  col=col_picked
)



# ttt<-round(seq(0.5,2,length=ncol_texts)*10)/10
ttt_bf<-seq(1,200,length=ncol_texts)
# ttt<-1
# ttt_bf<-200
ttt<-round(((ttt_bf-1)+100)/199.5*10)/10
mtext(ttt,side=2,at=tail(seq(yb,yt,(yt-yb)/ncol_texts),-1)-0.05,las=2,cex=0.7)
p_phy<-plot_phylo(tree.Groups_toplot,main=target_genus,direction="downwards",show.tip.label = FALSE,edge.width=2)



p_vals<-outTMATbeta[[ind_direct_genera]][,1]

if(!min(p_vals) %in% p_vals[-length(p_vals)]){
  ind_minP_node<- 0
  col_level<-rep(1,n_tips)
  if(as.logical(output$dir[which.min(p_vals)])){
    
    col_level[-n_tips]<-2
    col_level[n_tips]<-3
    
    
  }else{
    
    col_level[-n_tips]<-3
    col_level[n_tips]<-2
    
  }
  # tree.Groups_toplot$node.label
  makedm<-function(ind_sig_node){
    C1<-apply(target$simData$X_P,1,sum)
    C2<-target$simData$X_P_comp
    dm<- type_to_dm(C1=C1,C2=C2,type=type,total.reads=target_ori$totalcounts)[,1]
  }
  
  nodes_left<-1:length(p_vals)
  nodes_right<-0
  
}else{
  ind_minP_node<- which(p_vals[-length(p_vals)]==min(p_vals))
  
  # title("Title text", adj = 0, line = 0)
  
  target_plot<-output$testnodes[[which.min(p_vals)]] 
  col_level<-rep(1,n_tips)
  if(as.logical(output$dir[which.min(p_vals)])){
    col_level[which(target_plot$group1)]<-2
    col_level[which(target_plot$group2)]<-3
  }else{
    col_level[which(target_plot$group1)]<-3
    col_level[which(target_plot$group2)]<-2
  }
  makedm<-function(ind_sig_node){
    C1<-target$simData$X_P%*%ttttt$testnodes[[ind_sig_node]]$group1
    C2<-target$simData$X_P%*%ttttt$testnodes[[ind_sig_node]]$group2
    dm<- type_to_dm(C1=C1,C2=C2,type=type,total.reads=target_ori$totalcounts)[,1]
  }
  
  nodes_left<-which(target_plot$group1)
  nodes_right<-which(target_plot$group2)
  
}


col_level2<-col_level[c(length(col_level),1:(length(col_level)-1))]
co_alpha<-0.7
co_black<-function(co_alpha){rgb(0, 0, 0, alpha=co_alpha)}
co_red<-function(co_alpha){rgb(1, 0, 0, alpha=co_alpha)}
co_blue<-function(co_alpha){rgb(0, 0, 1, alpha=co_alpha)}
coco<-c("white","#ff756b","#00bec6")
# mycol<-c("black","red","blue")[col_level]
co_alpha<-1
col_nodes_picked<-sapply(data.frame(col2rgb(col_nodes2)),function(data){
  rgb(data[1],data[2],data[3],alpha=co_alpha,maxColorValue=255)
})
mycol<-coco[col_level2]
pchs<-rep(21,Nnode(tree.Groups_toplot))
# pchs[which.min(output$min_index)]<-25
# ?nodelabels
nodelabels(pch=pchs, col="black", bg=col_nodes2, cex=cexSize,srt=0, frame = "none")
nodelabels(paste0("k=",1:Nnode(tree.Groups_toplot)-1),col="black" ,bg=col_nodes2,cex=nodeCex,srt=0, frame = "none")

p_vals2_bf<-p.adjust(p_vals[c(length(p_vals),1:(length(p_vals)-1))],method="fdr")
p_vals2<-round(p_vals2_bf*1000000)/1000000


tt_targetX<-p_phy[[1]]
targetX<-tt_targetX[(n_tips+1):length(tt_targetX)]
tt_targetY<-p_phy[[2]]
targetY<-tt_targetY[(n_tips+1):length(tt_targetY)]

text(targetX,targetY-max(tt_targetY)/15,paste0("P-value: ",p_vals2),col="black" ,bg=col_nodes2,cex=1,srt=0)


if(n_tips==2){
  tmp1<-mycol[1]
  tmp2<-mycol[2]
  mycol<-c(tmp2,tmp1)
  tiplabels(pch=22, col="black", bg=mycol, cex=cexSize,srt=0,frame="none")
  tiplabels(paste0("m=",c(1,0)),pch="", col="black", bg=co_black(co_alpha), cex=nodeCex,srt=0,frame="none")  
}else{
  tiplabels(pch=22, col="black", bg=mycol, cex=cexSize,srt=0,frame="none")
  tiplabels(paste0("m=",c(1:n_tips-1)),pch="", col="black", bg=co_black(co_alpha), cex=nodeCex,srt=0,frame="none")  
}


# tiplabels(tiplabels_ready,pch="", col="white", adj=0, bg=mycol, cex=1)
nlabel<-vector("character",Nnode(tree.Groups_toplot))


# dftmp_bf<-data.frame(target$simData$X_P[,pruned.tree_perphy_omiat$tip.label])
# dftmp<-cbind(dftmp_bf,target$simData$X_P_comp)
# df_orinm<-pruned.tree_perphy_omiat$tip.label
# names(dftmp)<-c(paste0("T",1:length(df_orinm)),"T0")



dftmp_bf<-data.frame(target$simData$X_P[,tree_chosen2$tip.label])
sumVal<-apply(dftmp_bf,1,sum)
dftmp<-cbind(dftmp_bf,sumVal)
df_orinm<-tree_chosen2$tip.label

if(length(df_orinm)>5){
  str_m<-paste0(1,",...,",length(df_orinm))
}else{
  str_m<-paste0(1:length(df_orinm),collapse=",")
}
names(dftmp)<-c(paste0("m=",1:length(df_orinm)),paste0("m=",str_m))

stacked<-stack(dftmp) 



groups_for_plot_bf<-c("Control","Case")[Yvar+1]
groups_for_plot_bf<-factor(groups_for_plot_bf)

groups_for_plot<-factor(groups_for_plot_bf,levels = c("Case","Control"))
cbd_stacked<-cbind(stacked,groups_for_plot)
names(cbd_stacked)<-c("LogCPM","Leaf_nodes","groups")



# 
# png("lala0406.png")
# p1.0 <- ggplot(outToplot, aes(x=x, y=y,color=groups,group=groups,fill=groups)) +geom_line()
# p1.2 <- p1.0 + geom_ribbon(aes(ymax=upper, ymin=lower,color=groups,fill=groups), alpha=1/5)
# p1.3 <- p1.2 +  geom_point(aes(y=y_points))
# p1.4 <- p1.3 + theme_classic()
# p1.4
# 
# dev.off()
# system("scp lala0406.png ng:~")
# df_dm<-data.frame(do.call("cbind",lala))

# cbd_stacked













# p <- ggplot(data = cbd_stacked, aes(x = Leaf_nodes, y = LogCPM)) + 
#   geom_boxplot(aes(fill = groups), width = 0.8) + theme_bw()+ theme(legend.position="bottom",plot.title = element_text(hjust = 0.5))+labs(title="Log counts per million by groups",y="Log counts per million",x="Leaf nodes")



df_dm_tmp_tips<-cbind(cbd_stacked,week=weeks)
loess_bf_tip<-data.frame(LogCPM=df_dm_tmp_tips$LogCPM , y_points=df_dm_tmp_tips$LogCPM, Leaf_nodes=df_dm_tmp_tips$Leaf_nodes, groups=df_dm_tmp_tips$groups,week=df_dm_tmp_tips$week)

spdf_tip<-split(loess_bf_tip,paste0(as.character(loess_bf_tip$Leaf_nodes),"_",loess_bf_tip$groups))
sp_les_tip<-lapply(spdf_tip,function(data){
  loess_output<-loess.sd(x=data$week, y =data$LogCPM, nsigma = 1,span=0.75)
  cbind(data.frame(x=loess_output$x,y=loess_output$y,upper=loess_output$upper,lower=loess_output$lower),data)
})

outToplot_tips<-do.call("rbind",sp_les_tip)


if(length(unique(outToplot$x))>3){
  
  
# png("lala0406_tips.png",height=1024,width=1024)
  
p1.0 <- ggplot(outToplot_tips, aes(x=x, y=y,color=groups,group=groups,fill=groups)) +geom_line()
p1.2 <- p1.0 + geom_ribbon(aes(ymax=upper, ymin=lower,color=groups,fill=groups), alpha=1/5)
p1.3 <- p1.2 +  geom_point(aes(y=y_points))
p1.4 <- p1.3 + theme_classic()+ theme(legend.position="bottom",plot.title = element_text(hjust = 0.5))+labs(title="logCPM for each leaf node",y="logCPM",x=x_axis_label)
p<-p1.4+facet_wrap(~Leaf_nodes, scales="free_y", nrow=1)

# p
# dev.off()
# system("scp lala0406_tips.png ng:~")




}else{
  
}



vp <- viewport(height = unit(0.4,"npc"), width=unit(0.6, "npc"), 
               just = c("left","top"),   y = 0.4, x = 0)
print(p, vp = vp)

AllNodesPlot<-F

fontsize=min((21-n_tips)/2,5)
text<-paste(paste0("m=",c((1:n_tips-1),str_m),": ",tiplabels_ready),collapse="\n")






sp <- ggplot(NULL, aes(0, 1, label = text))+geom_point(pch="")
p2<-sp + geom_text(hjust=0,size=fontsize) + theme_bw()+ xlim(c(0, 1))+theme(axis.line=element_blank(),
                                                                            axis.text.x=element_blank(),
                                                                            axis.text.y=element_blank(),
                                                                            axis.ticks=element_blank(),
                                                                            axis.title.x=element_blank(),
                                                                            axis.title.y=element_blank(),
                                                                            legend.position="none",
                                                                            panel.background=element_blank(),
                                                                            panel.border=element_blank(),
                                                                            panel.grid.major=element_blank(),
                                                                            panel.grid.minor=element_blank(),
                                                                            plot.background=element_blank())



vp <- viewport(height = unit(0.4,"npc"), width=unit(0.6, "npc"), 
               just = c("left","top"),   y = 0.4, x = 0.6)
print(p2, vp = vp)


# ind_sig_node<-which.min(p_vals2)

ind_sig_node<-ind_minP_node



if(AllNodesPlot){
  lala<-lapply(2:(length(p_vals)),makedm)
}else{
  lala<-lapply(ind_minP_node,makedm)
}

# install.packages("statVisual")
library(statVisual)



# statVisual(type = 'Box', 
#            data = df_dm_tmp, 
#            x = 'week', 
#            y = 'LogCPM', 
#            group = 'groups',
#            title = "Boxplots across time") 
# install.packges("msir")
library(msir)

# 



df_dm<-data.frame(do.call("cbind",lala))



if(dim(df_dm)[2]!=1){
  stacked_dm<-stack(df_dm) 
  names(stacked_dm)<-c("LogCPM","Nodes")
  names(stacked_dm)[1]<-"LogCPM"
  df_dm<-cbind(stacked_dm,data.frame(groups=groups_for_plot))
  p3 <- ggplot(data = df_dm, aes(x = Nodes, y = LogCPM)) + 
    geom_boxplot(aes(fill = groups), width = 0.8) + theme_bw()+ theme(legend.position="bottom",plot.title = element_text(hjust = 0.5))+labs(title=paste("Log ratio of logCPMs of the each nodes"),y="Log ratio of log CPM")
}else{
  stacked_dm<-df_dm
  names(stacked_dm)[1]<-"LogCPM"
  df_dm<-cbind(stacked_dm,data.frame(groups=groups_for_plot))

  df_dm_tmp<-cbind(stacked_dm,data.frame(groups=groups_for_plot,week=weeks))
  loess_bf<-data.frame(LogCPM=df_dm_tmp$LogCPM , y_points=df_dm_tmp$LogCPM,groups=df_dm_tmp$groups,week=df_dm_tmp$week)
  
  spdf<-split(loess_bf,loess_bf$groups)
  sp_les<-lapply(spdf,function(data){
    loess_output<-loess.sd(x=data$week, y =data$LogCPM, nsigma = 1,span=0.75)
    cbind(data.frame(x=loess_output$x,y=loess_output$y,upper=loess_output$upper,lower=loess_output$lower),data)
  })
  
  outToplot<-do.call("rbind",sp_les)
  # loess_output<-loess.sd(x=weeks, y =df_dm_tmp$LogCPM, nsigma = 1)
  


  if(length(unique(outToplot$x))>3){
    p1.0 <- ggplot(outToplot, aes(x=x, y=y,color=groups,group=groups,fill=groups)) +geom_line()
    p1.2 <- p1.0 + geom_ribbon(aes(ymax=upper, ymin=lower,color=groups,fill=groups), alpha=1/5)
    p1.3 <- p1.2 +  geom_point(aes(y=y_points))
    p3 <- p1.3 + theme_classic()+ theme(legend.position="bottom",plot.title = element_text(hjust = 0.5))+labs(title=substitute(log~(frac(d,l))~d2, list(d2=paste0(" for the most significant node k=",ind_sig_node),d=paste0("logCPM of m=",paste0(nodes_left,collapse = ",")),l=paste0("logCPM of m=",paste0(nodes_right,collapse = ","))) ),y="Log ratio of logCPM")
    
  }else{
    
  }
  

  # p3<-ggplot(data=df_dm, aes(x = groups, y = LogCPM, fill = groups)) +
  #   geom_boxplot() + theme_bw()+ theme(legend.position="bottom",plot.title = element_text(hjust = 0.5))+labs(title=paste0("Log ratio of log CPMs of test node k=",ind_sig_node),y="Log ratio of log CPM")
  
}

vp <- viewport(height = unit(0.6,"npc"), width=unit(0.4, "npc"), 
               just = c("left","top"),   y = 1, x = 0.6)
print(p3, vp = vp)


dev.off()
system("scp TMATPlot0406.png ng:~")

