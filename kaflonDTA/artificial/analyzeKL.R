
source("analysisFunctions.Rscript")


fnamT<-dir("/home/dfc/seqs",full.names=T)
fnamTb<-dir("/home/dfc/seqsTb",full.names=T)
fnamTa<-dir("/home/dfc/seqsTa",full.names=T)
fnamTab<-dir("/home/dfc/seqsTab",full.names=T)



tmpnam<- lapply(fnamTb,function(x) sub("ForTb.txt","",strsplit(x,"/")[[1]][5]))
tmpnam<-unlist(tmpnam)

#tmp<-getKLValues(c(1:length(fnamTb)),c(1:18),c(5:15))
kl<-read.csv("KLValues.txt.gz",header=T)
klmeans<-aggregate(kl$KLNorm,list(kl$species),mean)
klsd<-aggregate(kl$KLNorm,list(kl$species),sd)
klstats<- cbind(tmpnam,klmeans,klsd[[2]])
colnames(klstats)<- c("name","species","KLmean","Klsd")
tmp<-order(klstats$KLmean)
klSorted<- klstats[tmp,]
