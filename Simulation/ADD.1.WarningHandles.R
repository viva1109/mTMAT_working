print("warning handles")


tryCatch( { result <- log("not a number"); print(res) }
          , warning = function(e) {an.error.occured <<- TRUE})
print(an.error.occured)

tryCatch( { warning("wow") }
          , warning = function(e) { an.warning <<- TRUE})


Fixed_GLMMMiRKAT(y_vec[ind_out_Kernel], cov=NULL, id=id_vec[ind_out_Kernel], Ks=Ks, model="gaussian",n.perm = n_perm)

duplist<-id[duplicated(id)]

newy<-y[which(id%in%duplist)]
newid<-id[which(id%in%duplist)]

fit <- lmer(newy ~ (1 | newid))