library(MRInstruments)
data(gwas_catalog)
head(gwas_catalog)
data(metab_qtls)
head(metab_qtls)

MET2<-metab_qtls[metab_qtls$phenotype=="AcAce",]
MET2$