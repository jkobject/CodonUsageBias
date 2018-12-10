
source("analysisFunctions.Rscript")


fnamT<-dir("/home/dfc/seqs",full.names=T)
fnamTb<-dir("/home/dfc/seqsTb",full.names=T)
fnamTa<-dir("/home/dfc/seqsTa",full.names=T)
fnamTab<-dir("/home/dfc/seqsTab",full.names=T)


res<- lapply(c(1:length(fnamTb)),function(x){ 
tmp1<-globalPlotOnlyDta(fnamTb[x],c(1:9),c(5:10),3)
tmp<-summary(fitMe(tmp1,"V1","V2"))$coef
obssl<-as.numeric(tmp[2])
obserr<-as.numeric(tmp[4])
tmp1<-globalPlotOnlyDta(fnamTab[x],c(1:9),c(5:10),3)
tmp<-summary(fitMe(tmp1,"V1","V2"))$coef
ransl<-as.numeric(tmp[2])
ranerr<-as.numeric(tmp[4])
c(obssl,obserr,ransl,ranerr)
})

tmp<-do.call(rbind,res)
tmp<-as.data.frame(tmp)


nams<-lapply(fnamTb, function(x) { tmp <- sub("/home/dfc/seqsTb/","",x); tmp <- sub("ForTb.txt","",tmp);tmp     })
nams<-unlist(nams)
tmp<-cbind(nams,tmp)
colnames(tmp) <-c("species","slope","error","slopeRan","errorRan")
write.table(file="slopes5To10Smaller3.txt",tmp,sep=",",row.names=F,quote=F)
