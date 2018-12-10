
source("analysisFunctions.Rscript")

 

res<- lapply(c(1:length(protnamTb)),function(x){ 
tmp1<-globalPlotOnlyDta(protnamTb[x],c(1:9),c(5:10),0)
tmp<-summary(fitMe(tmp1,"V1","V2"))$coef
obssl<-as.numeric(tmp[2])
obserr<-as.numeric(tmp[4])
tmp1<-globalPlotOnlyDta(protnamTa[x],c(1:9),c(5:10),0)
tmp<-summary(fitMe(tmp1,"V1","V2"))$coef
ransl<-as.numeric(tmp[2])
ranerr<-as.numeric(tmp[4])
c(obssl,obserr,ransl,ranerr)
})

tmp<-do.call(rbind,res)
tmp<-as.data.frame(tmp)


nams<-lapply(protnamTb, function(x) { tmp <- sub("/home/dfc/seqsPTb/","",x); tmp <- sub("ForPtTb.txt","",tmp);tmp     })
nams<-unlist(nams)
tmp<-cbind(nams,tmp)
colnames(tmp) <-c("species","slope","error","slopeRan","errorRan")
write.table(file="slopes5To10P.txt",tmp,sep=",",row.names=F,quote=F)
