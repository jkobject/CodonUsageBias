
source("analysisFunctions.Rscript")


fnamT<-dir("/home/dfc/seqs",full.names=T)
fnamTb<-dir("/home/dfc/seqsTb",full.names=T)
fnamTa<-dir("/home/dfc/seqsTa",full.names=T)
fnamTab<-dir("/home/dfc/seqsTab",full.names=T)

res<- lapply(fnamT[],function(x) calcSanov(x,c(1:18),c(5:20),F))
tmp<-do.call(rbind.data.frame,res)
write.table(file="sanovDetail.txt",tmp,sep=",",row.names=F,quote=F)
