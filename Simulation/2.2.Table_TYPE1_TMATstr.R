infol<-"~/analysis/20200928_longi_method/GENRA_diffN_type1_2CORSSSEC/"
DVD<-1
HeatSTR<-paste("GENRA_210406_2CRSEC",sep="_")

SingleOTU<-TRUE
library(stringr)

file_list<-list.files(infol)
file_list2_bf<-file_list[str_detect(file_list,HeatSTR)]

file_list2<-file_list2_bf[order(as.numeric(str_extract(file_list2_bf,"\\d+(?=.rds)")))]
rdrd<-lapply(file.path(infol,file_list2),readRDS)
rdrd2<-lapply(rdrd,function(data){
  data[[3]]
})


ttable<-lapply(rdrd,function(data){
  data[[1]]
})


labels<-lapply(rdrd,function(data){
  data[[2]]
})

cutoff<-0.1
table_TYPE1_list<-lapply(c(0.1,0.05,0.01,0.005),function(cutoff){
  partsTOPLOT<-lapply(1:length(rdrd2),function(ind_data){
    
    dfTMP<-data.frame(Taxa=labels[[ind_data]][[3]][ttable[[ind_data]][,4]],beta=labels[[ind_data]][[4]][ttable[[ind_data]][,5]])
    Values<-t(sapply(rdrd2[[ind_data]],function(data){
      apply(data,2,function(value){
        sum(value<cutoff,na.rm=T)/sum(!is.na(value))
      })
    }))
    # Values2<-c()
    # for( i in 1:dim(Values)[2]){
    #   Values2<-c(Values2,Values[,i])
    # }
    # Method<-rep(c("mTMAT_M","mTMAT_IM","GLMM-MiRKAT","FZINBMM","ZIBR"),each=dim(Values)[1])
    # 
    # DFDF<-cbind(Values2,dfTMP,Method)
    
  })
  lapply(partsTOPLOT,dim)
  partsTOPLOT2<-do.call("rbind",partsTOPLOT)
  
  
  
  sn_table<-do.call("rbind",ttable)
  
  
  list_n<-c(30,50,100)
  list_c_prop<-c(0.25,0.5)

  
  
  resp<-split(data.frame(partsTOPLOT2),paste(sn_table[,2],sn_table[,3]))
  table_TYPE1<-sapply(resp,function(data){
    apply(data,2,mean,na.rm=T)
  })
  cbind(table_TYPE1,cutoff)
})
table_TYPE1_out<-do.call("rbind",table_TYPE1_list)


nm<-"table_TYPE1_othermethodsCross.csv"
write.csv(table_TYPE1_out,nm)

system(paste0("scp ",nm," ng:~"))
