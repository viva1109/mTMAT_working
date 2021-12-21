infol<-"~/analysis/20200928_longi_method/output/"
file_list<-list.files(infol)
file_list2<-file_list[c(1,9,2:8,10)]
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
partsTOPLOT<-lapply(1:length(rdrd2),function(ind_data){
  
  dfTMP<-data.frame(Taxa=labels[[ind_data]][[3]][ttable[[ind_data]][,4]],beta=labels[[ind_data]][[4]][ttable[[ind_data]][,5]])
  Values<-t(sapply(rdrd2[[ind_data]],function(data){
    apply(data,2,function(value){
      sum(value<cutoff)/length(value)
    })
  }))
  Values2<-c()
  for( i in 1:dim(Values)[2]){
    Values2<-c(Values2,Values[,i])
  }
  Method<-rep(c("mTMAT_M","mTMAT_IM","GLMM-MiRKAT","FZINBMM","FZINBMM"),each=dim(Values)[1])
  
  DFDF<-cbind(Values2,dfTMP,Method)
  
})

partsTOPLOT2<-do.call("rbind",partsTOPLOT)

Values3<-split(partsTOPLOT2$Values2,paste(partsTOPLOT2$beta,partsTOPLOT2$Method,sep="__"))
Values4<-sapply(Values3,mean)
library(stringr)
DFPLOT<-data.frame(cbind(Values4,t(sapply(str_split(names(Values4),"__"),"[",1:2))))
names(DFPLOT)<-c("Value","beta","Method")

library(ggplot2)
ncol=4
methods_sequence<-c("mTMAT_M","mTMAT_IM","GLMM-MiRKAT","FZINBMM","FZINBMM")
DFPLOT$Method<-factor(DFPLOT$Method, levels = methods_sequence)



labels_out_bf<-as.character(levels(DFPLOT$Method))
labels_out<-labels_out_bf
labels_out[labels_out_bf=="mTMAT_M"]<-expression("mTMAT" ["M"])
labels_out[labels_out_bf=="mTMAT_IM"]<-expression("mTMAT" ["IM"])


ind_methods<-1:length(methods_sequence)
library(scales)
ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}
geom.colors <- ggplotColours(n=length(methods_sequence))
geom.colors_picked<-geom.colors[ind_methods]

DFPLOT$Value<-as.numeric(as.character(DFPLOT$Value))
DFPLOT$beta<-as.numeric(as.character(DFPLOT$beta))

ggplot(DFPLOT, aes(x=beta, y=Value, color=Method)) + 
  scale_color_manual(values=geom.colors_picked,labels=labels_out)+
  scale_linetype_manual(values=1:length(methods_sequence),labels=labels_out)+
  geom_line(aes(linetype=Method), size=1,alpha=0.7)+
  scale_shape_manual(values=c(16,3,5,8,7,15,17,18,19,20)[ind_methods],labels=labels_out) +
  geom_point(aes(shape=Method),alpha=0.7, size=4)+
  guides(linetype = guide_legend("",ncol=ncol,byrow=T),shape = guide_legend("",ncol=ncol,byrow=T),colour = guide_legend("",ncol=ncol,byrow=T),alpha = guide_legend("",ncol=ncol,byrow=T) ,size = guide_legend("none"))+ylab("Power estimates")+xlab("Beta")+theme(legend.position="bottom")+
  ylim(c(0,1))








