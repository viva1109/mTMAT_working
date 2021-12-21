
sapply(dms,getCores,dataSet$y,NULL,conti,DetailResult,contrast=contrast)

#t_X<-dm
#t_y_bf<-y
# y<-Yvar

dm<-dms[[1]]

betahat<-get_betahat(y,dm,cov=cov,conti=conti)
n_indv<-length(y)
if(conti){
  q<-1
}else{
  q<-length(unique(y))-1
}
if(!is.null(cov)){
  t_re_1<-getT2(y,dm,cov=cov,conti=conti,DetailResult =DetailResult)
  p<-dim(cov)[2]
}else{
  t_re_1<-getT2(y,dm,conti=conti,DetailResult =DetailResult)
  p<-1
}  


t_y_bf<-y
t_X<-dm

N<- length(t_y_bf)

v_one_n<-rep(1,N)
if(!is.null(cov)){
  t_A <- diag(N) - (cov %*% solve(t(cov)%*%cov) %*% t(cov))
  p<-dim(cov)[2]
}else{
  t_A <- diag(N)-v_one_n%*%solve(t(v_one_n)%*%v_one_n)%*%t(v_one_n)
  p<-1
}
if(!is.null(cov)){
  ready_tt<-as.matrix(t(t_X)%*%t_A%*%t_X,ncol=1)
  t_SIGMA <- (ready_tt/(N-dim(cov)[2]))[1,1]
  # t_SIGMA <- cov(as.matrix(t_X,ncol=1))[1,1]
}else{
  t_SIGMA <- cov(as.matrix(t_X,ncol=1))[1,1]
}

if(conti){
  t_y<-t_y_bf
  qq<-1
}else{
  t_y<-c()
  if(length(unique(t_y_bf))>2){
    for( i in 1:length(unique(t_y_bf))){
      if(i!=1){
        t_y<-cbind(t_y,as.numeric(t_y_bf==unique(t_y_bf)[i]))
      }
    }
    qq<-dim(t_y)[2]
  }else{
    t_y<-t_y_bf
    qq<-1
  }
}


y_sq<-(t(t_y)%*%t_y)
t_S <- t(t_y)%*%t_A%*%t_X
t_yAy <- t(t_y)%*%t_A%*%t_y
t_varS <- t_SIGMA*t_yAy
# print("t_varS")
# print(t_varS)
t_T<-t(t_S)%*%solve(t_varS)%*%t_S
betahat<-solve(y_sq)%*%t_S

t_T

####t_T is the objective

#
1-pchisq(t_re_1,df=1)
#
exact_pval<-function(tmp_pval){
  1-(1-min(tmp_pval))^length(tmp_pval)
}


