
source("analysisFunctions.Rscript")


fnamT<-dir("/home/dfc/seqs",full.names=T)
fnamTb<-dir("/home/dfc/seqsTb",full.names=T)
fnamTa<-dir("/home/dfc/seqsTa",full.names=T)
fnamTab<-dir("/home/dfc/seqsTab",full.names=T)
res<-compareSanov(c(1:length(fnamTb)),c(1:9),c(5:10))
write.table(file="sanovRatios.txt",res,sep=",",row.names=F,quote=F)
