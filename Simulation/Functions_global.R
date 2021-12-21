calculateTMAT<-function(Var_out2,Info_IDs,OTU_table_chosen,tree_chosen,INDI_SINGLE,SAMPLE_UID,TMATtype,str_i="independence",stat_i="score",OTU_not_chosen=NULL){
  # browser()
  indi_tAnal<-list()
  indi_tAnal$subtree_tmp<-list()
  
  if(!is.rooted(tree_chosen)){
    indi_tAnal$subtree_tmp[[1]]<-root(tree_chosen, 1, r = TRUE)
  }else{
    indi_tAnal$subtree_tmp[[1]]<-tree_chosen
  }
  # 
  # TMATtype<-15
  # OTU_not_chosen<-OTU_table_not_chosen2#NO!!!!
  
  
  # DMouts<-lapply(1:3,function(cTIME){
  #   Yvar<-Var_out2[,cTIME]
  #   total.reads<-Info_IDs[[cTIME]]$Rcnts[ind_samples]
  #   X_P<-OTU_table_chosen[[cTIME]]
  #   dataSet<-simulSet_by_methods_fixed(y=Yvar,dataSet=data.frame(X_P,check.names = F),ttr=total.reads,type=TMATtype,method="TMAT") 
  #   # output$pval
  #   TMAT_func_dm2(dataSet,indi_tAnal,type=TMATtype,total.reads=total.reads,conti=F,cov=NULL)$dms
  # })
  
  # VEC_Yvar<-c(Var_out2[,1],Var_out2[,2],Var_out2[,3])
  yC<-c(Var_out2[,1],Var_out2[,2],Var_out2[,3])
  VEC_total.reads<-c(Info_IDs[[1]]$Rcnts[ind_samples],Info_IDs[[2]]$Rcnts[ind_samples],Info_IDs[[3]]$Rcnts[ind_samples])
  VEC_X_P<-rbind(OTU_table_chosen[[1]],OTU_table_chosen[[2]],OTU_table_chosen[[3]])

  
  out_dataSet<-simulSet_by_methods_fixed(y=yC,dataSet=data.frame(VEC_X_P,check.names = F),ttr=VEC_total.reads,type=TMATtype,method="TMAT") 

  
  if(is.null(OTU_not_chosen)){
    DMsC_bf<-TMAT_func_dm2(out_dataSet,indi_tAnal,type=TMATtype,total.reads=total.reads,conti=F,cov=NULL)$dms
    if(INDI_SINGLE){
      DMsC<-list(DMsC_bf)
    }else{
      DMsC<-DMsC_bf
    }
  }else{
    
    # OTU_not_chosenFIX<-lapply(OTU_not_chosen,function(data){
    #   data[,-dim(data)[2]]
    # })
    # VEC_OTU_not_chosen<-rbind(OTU_not_chosenFIX[[1]],OTU_not_chosenFIX[[2]],OTU_not_chosenFIX[[3]])
    # VEC_GENUS<-apply(out_dataSet$simData$X_P,1,sum)
    # cpm_tmp<-t(cpm(t(VEC_OTU_not_chosen),log=T,lib.size=VEC_total.reads))
    # cpm_corrected<-log(exp(cpm_tmp)+1)
    # cpm_geo_togo<-cbind(VEC_GENUS,cpm_corrected)
    # geomean<-apply(cpm_geo_togo,1,function(value){
    #   (prod(value))^(1/length(value))
    # })
    # 
    # cpm_corrected_sum<-apply(cpm_corrected,1,sum)

    
    VEC_OTU_not_chosen<-rbind(OTU_not_chosen[[1]],OTU_not_chosen[[2]],OTU_not_chosen[[3]])
    VEC_GENUS<-apply(out_dataSet$simData$X_P,1,sum)
    cpm_tmp<-t(cpm(t(VEC_OTU_not_chosen),log=T,lib.size=VEC_total.reads))
    cpm_corrected<-log(exp(cpm_tmp)+1)
    cpm_geo_togo<-cbind(VEC_GENUS,cpm_corrected)
    geomean<-apply(cpm_geo_togo,1,function(value){
      (prod(value))^(1/length(value))
    })
    
    DMsC_bf<-TMAT_func_dm2(out_dataSet,indi_tAnal,type=TMATtype,total.reads=total.reads,conti=F,cov=NULL)$dms
    if(INDI_SINGLE){
      DMsC<-list(DMsC_bf)
    }else{
      DMsC<-DMsC_bf
    }
    
    if(TMATtype==16){
      DMsC[[length(DMsC)]]<-inv_norm_trans(log(VEC_GENUS/geomean))
    }else{
      DMsC[[length(DMsC)]]<-log(VEC_GENUS/geomean)
    }
    
    
    # DMsC_bf[[length(DMsC_bf)]]<-log(VEC_GENUS/cpm_corrected_sum)
    
    # cor(VEC_GENUS/cpm_corrected_sum,VEC_GENUS/geomean)
    # cor(VEC_GENUS,log(VEC_GENUS/geomean))
    # cor(VEC_GENUS,log(VEC_GENUS/cpm_corrected_sum))
    # sd(geomean)
    # sd(geomean2)
    # cbind(geomean,geomean2)
  }
  # 
  # cor.test(VEC_GENUS,yC)
  # cor.test(VEC_GENUS,yC,method="spearman")
  # 
  # 
  # cbind(VEC_X_P,out_dataSet[[2]][[1]])
  # 
  # cor.test(yC,VEC_X_P)
  # cor.test(yC,out_dataSet[[2]][[1]])
  # cor.test(yC,DMsC_bf)
  # cor.test(yC,log(out_dataSet[[2]][[1]]/out_dataSet[[2]][[2]]))
  # out_dataSet[[1]]
  # out_dataSet[[2]]
  
  # DMouts<-lapply(1:3,function(cTIME){
  #   Yvar<-Var_out2[,cTIME]
  #   total.reads<-Info_IDs[[cTIME]]$Rcnts[ind_samples]
  #   X_P<-OTU_table_chosen[[cTIME]]
  #   dataSet<-simulSet_by_methods_fixed(y=Yvar,dataSet=data.frame(X_P,check.names = F),ttr=total.reads,type=TMATtype,method="TMAT") 
  #   # output$pval
  #   TMAT_func_dm2(dataSet,indi_tAnal,type=TMATtype,total.reads=total.reads,conti=F,cov=NULL)$dms
  # })
  
  
  

  
  IDsC<-rep(SAMPLE_UID[ind_samples],3)  
  TIME<-c(rep(1,length(SAMPLE_UID[ind_samples])),rep(2,length(SAMPLE_UID[ind_samples])),rep(3,length(SAMPLE_UID[ind_samples])))
  # cbind(do.call("cbind",DMsC),yC)
  
  
  # cor.test(VEC_GENUS,yC)
  # cor.test(VEC_GENUS,yC,method="spearman")
  # cor.test(DMsC[[1]],yC)
  # cor.test(DMsC[[1]],yC,method="spearman")
  # split(DMsC[[1]],yC)
  # 
  # cbind(DMsC[[1]],yC)
  # cor.test(DMsC[[ind_node]],yC)
  p_result<-sapply(1:length(DMsC),function(ind_node){
    DATA<-data.frame(response=DMsC[[ind_node]],disease=yC,subject=IDsC,TIME=TIME)
    DATA_sort <- DATA[order(DATA$subject),]
    formula <- response~disease
    D0<-matrix(c(0,1),nrow=1)
    m2<-geem2(response~disease,id=subject,data=DATA_sort,family=gaussian,corstr=str_i)
    res<-mmmgee.test(x=m2,L=list(D0),r=list(c(0)),statistic=stat_i,type="quadratic",biascorr=F)
    res$test$global$p
  })
  mTMAT_res<-exact_pval(p_result)
  
  return(mTMAT_res)
}


estimate_rho_TMAT<-function(Var_out2,Info_IDs,OTU_table_chosen,tree_chosen,INDI_SINGLE,SAMPLE_UID,TMATtype,str_i="independence",stat_i="score"){
  indi_tAnal<-list()
  indi_tAnal$subtree_tmp<-list()
  
  if(!is.rooted(tree_chosen)){
    indi_tAnal$subtree_tmp[[1]]<-root(tree_chosen, 1, r = TRUE)
  }else{
    indi_tAnal$subtree_tmp[[1]]<-tree_chosen
  }
  
  
  
  
  DMouts<-lapply(1:3,function(cTIME){
    Yvar<-Var_out2[,cTIME]
    total.reads<-Info_IDs[[cTIME]]$Rcnts[ind_samples]
    X_P<-OTU_table_chosen[[cTIME]]
    dataSet<-simulSet_by_methods_fixed(y=Yvar,dataSet=data.frame(X_P,check.names = F),ttr=total.reads,type=TMATtype,method="TMAT") 
    # output$pval
    TMAT_func_dm2(dataSet,indi_tAnal,type=TMATtype,total.reads=total.reads,conti=F,cov=NULL)$dms
  })
  
  # VEC_Yvar<-c(Var_out2[,1],Var_out2[,2],Var_out2[,3])
  yC<-c(Var_out2[,1],Var_out2[,2],Var_out2[,3])
  VEC_total.reads<-c(Info_IDs[[1]]$Rcnts[ind_samples],Info_IDs[[2]]$Rcnts[ind_samples],Info_IDs[[3]]$Rcnts[ind_samples])
  VEC_X_P<-rbind(OTU_table_chosen[[1]],OTU_table_chosen[[2]],OTU_table_chosen[[3]])
  
  out_dataSet<-simulSet_by_methods_fixed(y=yC,dataSet=data.frame(VEC_X_P,check.names = F),ttr=VEC_total.reads,type=TMATtype,method="TMAT") 
  DMouts<-TMAT_func_dm2(out_dataSet,indi_tAnal,type=TMATtype,total.reads=total.reads,conti=F,cov=NULL)$dms
  
  # DMouts<-lapply(1:3,function(cTIME){
  #   Yvar<-Var_out2[,cTIME]
  #   total.reads<-Info_IDs[[cTIME]]$Rcnts[ind_samples]
  #   X_P<-OTU_table_chosen[[cTIME]]
  #   dataSet<-simulSet_by_methods_fixed(y=Yvar,dataSet=data.frame(X_P,check.names = F),ttr=total.reads,type=TMATtype,method="TMAT") 
  #   # output$pval
  #   TMAT_func_dm2(dataSet,indi_tAnal,type=TMATtype,total.reads=total.reads,conti=F,cov=NULL)$dms
  # })
  
  
  
  # if(INDI_SINGLE){
  #   DMsC<-list(c(DMouts[[1]],DMouts[[2]],DMouts[[3]]))
  # }else{
  #   DMsC<-lapply(1:length(DMouts[[1]]),function(indind){
  #     c(DMouts[[1]][[indind]],DMouts[[2]][[indind]],DMouts[[3]][[indind]])
  #   })
  # }
  
  IDsC<-rep(SAMPLE_UID[ind_samples],3)  
  TIME<-c(rep(1,length(SAMPLE_UID[ind_samples])),rep(2,length(SAMPLE_UID[ind_samples])),rep(3,length(SAMPLE_UID[ind_samples])))
  
  
  rho_result<-sapply(1:length(DMsC),function(ind_node){
    DATA<-data.frame(response=DMsC[[ind_node]],disease=yC,subject=IDsC,TIME=TIME)
    DATA_sort <- DATA[order(DATA$subject),]
    formula <- response~disease
    D0<-matrix(c(0,1),nrow=1)
    m2<-geem2(response~disease,id=subject,data=DATA_sort,family=gaussian,corstr=str_i)
    m2$alpha
  })
  
  return(rho_result)
}


exact_pval<-function(tmp_pval){
  1-(1-min(tmp_pval))^length(tmp_pval)
}


MakeMat1<-function(inputVec,colnm){
  valMat<-matrix(inputVec,ncol=1)
  colnames(valMat)<-colnm
  return(valMat)
}




CrosssecTMAT<-function(pheno_vec,total_reads,OTU_table_chosen_vec,tree_chosen,TMATtype,str_i="independence",stat_i="score",TestType="Chi-squared"){
  indi_tAnal<-list()
  indi_tAnal$subtree_tmp<-list()
  
  if(!is.rooted(tree_chosen)){
    indi_tAnal$subtree_tmp[[1]]<-root(tree_chosen, 1, r = TRUE)
  }else{
    indi_tAnal$subtree_tmp[[1]]<-tree_chosen
  }
  
  
  
  

  Yvar<-pheno_vec
  total.reads<-total_reads
  X_P<-OTU_table_chosen_vec
  dataSet<-simulSet_by_methods_fixed(y=Yvar,dataSet=data.frame(X_P,check.names = F),ttr=total.reads,type=TMATtype,method="TMAT") 
    # output$pval
  TMAT_pval_bygenus_bf<-list(TMAT_func(dataSet,indi_tAnal,type=TMATtype,total.reads=total.reads,conti=F,cov=NULL,getBeta=T,DetailResult=F))

  TMAT_pval_bygenus_bf2<-lapply(TMAT_pval_bygenus_bf,function(data){
    data$pval
    # print("data")
    # print(data)
    
  })
  
  
  TMAT_pval_bygenus<-do.call("rbind",TMAT_pval_bygenus_bf2)
  genus<-"blank"
  TMAT_results<-cbind(data.frame(genus,stringsAsFactors=F),TMAT_pval_bygenus)
  if(TestType=="Chi-squared"){
    result_output<-cbind(data.frame(TMAT_results[,1],stringsAsFactors=F),TMAT_results[,3],p.adjust(TMAT_results[,3],method="fdr"))
  }else if(TestType=="F-test"){
    result_output<-cbind(data.frame(TMAT_results[,1],stringsAsFactors=F),TMAT_results[,2],p.adjust(TMAT_results[,2],method="fdr"))
  }else{
    print('TestType has to be "Chi-squared" or "F-test"')
  }
  names(result_output)<-c("Genus","p-value","FDR")

  return(result_output)
}



Fixed_GLMMMiRKAT<-function (y, cov = NULL, id, time.pt = NULL, Ks, model, slope = FALSE, 
                            n.perm = 5000) 
{
  #browser()
  if (is.null(time.pt) & slope) {
    stop("time.pt is required for the random slope model")
  }
  if (is.null(time.pt)) {
    id <- as.character(id)
    ind <- order(id)
    id <- id[ind]
    y <- y[ind]
    if (!is.null(cov)) {
      if (is.matrix(cov) | is.data.frame(cov)) {
        cov <- as.data.frame(cov)[ind, ]
      }
      else {
        cov <- as.data.frame(cov[ind])
      }
    }
    for (i in 1:length(Ks)) {
      Ks[[i]] <- Ks[[i]][ind, ind]
    }
  }
  else {
    ind <- order(as.character(time.pt))
    id <- as.character(id)
    id <- id[ind]
    y <- y[ind]
    if (!is.null(cov)) {
      if (is.matrix(cov) | is.data.frame(cov)) {
        cov <- as.data.frame(cov)[ind, ]
      }
      else {
        cov <- as.data.frame(cov[ind])
      }
    }
    time.pt <- time.pt[ind]
    for (i in 1:length(Ks)) {
      Ks[[i]] <- Ks[[i]][ind, ind]
    }
    ind <- order(id)
    id <- id[ind]
    y <- y[ind]
    if (!is.null(cov)) {
      if (is.matrix(cov) | is.data.frame(cov)) {
        cov <- as.data.frame(cov)[ind, ]
      }
      else {
        cov <- as.data.frame(cov[ind])
      }
    }
    time.pt <- time.pt[ind]
    for (i in 1:length(Ks)) {
      Ks[[i]] <- Ks[[i]][ind, ind]
    }
  }
  if (model == "gaussian") {
    y <- as.numeric(scale(y))
    if (is.null(cov)) {
      if (!is.null(time.pt) & slope) {
        fit <- lmer(y ~ (time.pt | id))
      }
      else {
        fit <- lmer(y ~ (1 | id))
      }
    }
    else {
      cov <- as.data.frame(cov)
      if (ncol(cov) == 1) {
        colnames(cov) <- "cov1"
        if (!is.null(time.pt) & slope) {
          formula.0 <- as.formula("y ~ (time.pt | id) + cov1")
          fit <- lmer(formula.0, data = cov)
        }
        else {
          formula.0 <- as.formula("y ~ (1 | id) + cov1")
          fit <- lmer(formula.0, data = cov)
        }
      }
      else {
        colnames(cov) <- paste("cov", 1:ncol(cov), sep = "")
        if (!is.null(time.pt) & slope) {
          formula.0 <- as.formula(paste(paste(c("y ~ (time.pt | id)", 
                                                paste("cov", 1:ncol(cov), sep = "")), collapse = " + "), 
                                        collapse = " + "))
          fit <- lmer(formula.0, data = cov)
        }
        else {
          formula.0 <- as.formula(paste(paste(c("y ~ (1 | id)", 
                                                paste("cov", 1:ncol(cov), sep = "")), collapse = " + "), 
                                        collapse = " + "))
          fit <- lmer(formula.0, data = cov)
        }
      }
    }
  }
  else {
    if (is.null(cov)) {
      if (!is.null(time.pt) & slope) {
        fit <- glmer(y ~ (time.pt | id), family = model)
      }
      else {
        fit <- glmer(y ~ (1 | id), family = model)
      }
    }
    else {
      cov <- as.data.frame(cov)
      if (ncol(cov) == 1) {
        colnames(cov) <- "cov1"
        if (!is.null(time.pt) & slope) {
          formula.0 <- as.formula("y ~ (time.pt | id) + cov1")
          fit <- glmer(formula.0, family = model, data = cov)
        }
        else {
          formula.0 <- as.formula("y ~ (1 | id) + cov1")
          fit <- glmer(formula.0, family = model, data = cov)
        }
      }
      else {
        colnames(cov) <- paste("cov", 1:ncol(cov), sep = "")
        if (!is.null(time.pt) & slope) {
          formula.0 <- as.formula(paste(paste(c("y ~ (time.pt | id)", 
                                                paste("cov", 1:ncol(cov), sep = "")), collapse = " + "), 
                                        collapse = " + "))
          fit <- glmer(formula.0, family = model, data = cov)
        }
        else {
          formula.0 <- as.formula(paste(paste(c("y ~ (1 | id)", 
                                                paste("cov", 1:ncol(cov), sep = "")), collapse = " + "), 
                                        collapse = " + "))
          fit <- glmer(formula.0, family = model, data = cov)
        }
      }
    }
  }
  if(!is.null(fit@optinfo$conv$lme4$message)){
    if(fit@optinfo$conv$lme4$message=="boundary (singular) fit: see ?isSingular"){
      stop("Boundary fit for lme4. Singular problem")
    }
  }
  
  
  fe <- as.numeric(getME(fit, "X") %*% getME(fit, "fixef"))
  re <- as.numeric(getME(fit, "Z") %*% getME(fit, "b"))
  if (model == "gaussian") {
    mu <- fe + re
    inv.de.mu <- diag(length(y))
  }
  if (model == "binomial") {
    mu <- exp(fe + re)/(1 + exp(fe + re))
    inv.de.mu <- diag(as.numeric(mu * (1 - mu)))
  }
  if (model == "poisson") {
    mu <- exp(fe + re)
    inv.de.mu <- diag(as.numeric(mu))
  }
  y.star <- as.numeric(fe + re + ginv(inv.de.mu) %*% (y - mu))
  r <- y.star - fe
  r.cov <- as.matrix((getME(fit, "Z") %*% Matrix::crossprod(getME(fit, 
                                                                  "Lambdat")) %*% getME(fit, "Zt")))
  if (model == "gaussian") {
    e.var <- diag(rep(var(mu), length(y)))
  }
  if (model == "binomial") {
    e.var <- diag(rep(1, length(y)))
  }
  if (model == "poisson") {
    e.var <- diag(rep(1, length(y)))
  }
  v.cov <- r.cov + e.var
  inv.vcov <- ginv(v.cov)
  r.s <- list()
  if (!slope) {
    id.time.pt <- id
    p.ind <- as.matrix(getPermuteMatrix(n.perm, length(r), 
                                        strata = id.time.pt))
    clust.sizes <- as.numeric(names(table(table(id.time.pt))))
    block.r.ids <- list()
    for (l in 1:n.perm) {
      block.r.id <- list()
      for (j in 1:length(clust.sizes)) {
        block.r.id.clust <- list()
        u.id.clust <- names(which(table(id.time.pt) == 
                                    clust.sizes[j]))
        r.u.id.clust <- u.id.clust[shuffle(length(u.id.clust))]
        for (k in 1:length(u.id.clust)) {
          block.r.id.clust[[k]] <- which(id.time.pt == 
                                           r.u.id.clust[k])
        }
        block.r.id[[j]] <- unlist(block.r.id.clust)
      }
      unlist.block.r.id <- rep(NA, length(id.time.pt))
      for (j in 1:length(clust.sizes)) {
        unlist.block.r.id[id.time.pt %in% names(which(table(id.time.pt) == 
                                                        clust.sizes[j]))] <- block.r.id[[j]]
      }
      block.r.ids[[l]] <- unlist.block.r.id
    }
    for (l in 1:n.perm) {
      p.ind[l, ] <- p.ind[l, block.r.ids[[l]]]
    }
    for (j in 1:n.perm) {
      r.s[[j]] <- r[p.ind[j, ]]
    }
  }
  if (!is.null(time.pt) & slope) {
    id.names <- names(table(id))
    id.time.pt <- rep(NA, length(id))
    for (j in 1:length(id.names)) {
      ind <- which(id == id.names[j])
      id.time.pt[ind] <- paste(c(id.names[j], as.character(time.pt[ind])), 
                               collapse = ".")
    }
    noid.time.pt <- rep(NA, length(id))
    for (j in 1:length(id.names)) {
      ind <- which(id == id.names[j])
      noid.time.pt[ind] <- paste(as.character(time.pt[ind]), 
                                 collapse = ".")
    }
    ex.clust <- names(table(noid.time.pt))
    ids <- rep(NA, length(id))
    for (j in 1:length(id.names)) {
      ind <- which(id == id.names[j])
      ids[ind] <- id.names[j]
    }
    if (length(names(which(table(noid.time.pt) == 1))) != 
        0) {
      warn.ids <- unique(ids[noid.time.pt %in% names(which(table(noid.time.pt) == 
                                                             1))])
      if (length(warn.ids) == 1) {
        warning("The cluster id (", warn.ids, ") is a silent block which is not exchangeable with any other blocks. It did not contribute to the analysis.")
      }
      if (length(warn.ids) > 1) {
        warning("The cluster ids (", paste(warn.ids, 
                                           collapse = ", "), ") are silent blocks which are not exchangeable with any other blocks. They did not contribute to the analysis.")
      }
    }
    block.r.ids <- list()
    for (l in 1:n.perm) {
      block.r.id <- list()
      for (j in 1:length(ex.clust)) {
        block.r.id.clust <- list()
        u.id.clust <- unique(id.time.pt[which(noid.time.pt == 
                                                ex.clust[j])])
        r.u.id.clust <- u.id.clust[shuffle(length(u.id.clust))]
        for (k in 1:length(u.id.clust)) {
          block.r.id.clust[[k]] <- which(id.time.pt == 
                                           r.u.id.clust[k])
        }
        block.r.id[[j]] <- unlist(block.r.id.clust)
      }
      unlist.block.r.id <- rep(NA, length(id.time.pt))
      for (j in 1:length(ex.clust)) {
        unlist.block.r.id[id.time.pt %in% unique(id.time.pt[which(noid.time.pt == 
                                                                    ex.clust[j])])] <- block.r.id[[j]]
      }
      block.r.ids[[l]] <- unlist.block.r.id
    }
    for (j in 1:n.perm) {
      r.s[[j]] <- r[block.r.ids[[j]]]
    }
  }
  Qs <- rep(NA, length(Ks))
  for (j in 1:length(Ks)) {
    Qs[j] <- (t(r) %*% inv.vcov %*% Ks[[j]] %*% inv.vcov %*% 
                r)/(t(r) %*% inv.vcov %*% r)
  }
  Q0s <- list()
  for (j in 1:length(Ks)) {
    Q0s.inv <- rep(NA, n.perm)
    for (k in 1:n.perm) {
      Q0s.inv[k] <- (t(r.s[[k]]) %*% inv.vcov %*% Ks[[j]] %*% 
                       inv.vcov %*% r.s[[k]])/(t(r.s[[k]]) %*% inv.vcov %*% 
                                                 r.s[[k]])
    }
    Q0s[[j]] <- Q0s.inv
  }
  pvs <- rep(NA, length(Ks))
  for (j in 1:length(Ks)) {
    pvs[j] <- (length(which(Q0s[[j]] > Qs[[j]])) + 1)/(n.perm + 
                                                         1)
  }
  T <- min(pvs)
  T0 <- rep(NA, n.perm)
  for (l in 1:n.perm) {
    T0.s.n <- list()
    for (m in 1:length(Ks)) {
      T0.s.n[[m]] <- Q0s[[m]][-l]
    }
    a.Ts <- unlist(lapply(Ks, function(x) return((t(r.s[[l]]) %*% 
                                                    inv.vcov %*% x %*% inv.vcov %*% r.s[[l]])/(t(r.s[[l]]) %*% 
                                                                                                 inv.vcov %*% r.s[[l]]))))
    a.pvs <- unlist(mapply(function(x, y) (length(which(x > 
                                                          y)) + 1)/n.perm, T0.s.n, a.Ts))
    T0[l] <- min(a.pvs)
  }
  pv.opt <- (length(which(T0 < T)) + 1)/(n.perm + 1)
  names(pvs) <- names(Ks)
  return(list(ItembyItem = pvs, aGLMMMiRKAT = pv.opt))
}



rbvunif <- function(n,rho) {
  x <- runif(n)
  if ((rho > 1.0) || (rho < -1.0)) {
    stop('rbvunif::rho not in [-1,+1]')
  }
  else if (rho==1.0) xy <- cbind(x,x)
  else if (rho==-1.0) xy <- cbind(x,1-x)
  else if (rho==0.0) xy <- cbind(x,runif(n))
  else {
    a <- (sqrt((49+rho)/(1+rho))-5)/2
    u <- rbeta(n,a,1.0)
    y <- runif(n)
    y <- ifelse(y<0.5 ,abs(u-x), 1-abs(1-u-x) )
    xy <- cbind(x,y)
  }
  return(xy)
}



