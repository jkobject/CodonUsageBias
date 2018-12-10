
source("analysisFunctions.Rscript")


fnamT<-dir("/home/dfc/seqs",full.names=T)
fnamTb<-dir("/home/dfc/seqsTb",full.names=T)
fnamTa<-dir("/home/dfc/seqsTa",full.names=T)
fnamTab<-dir("/home/dfc/seqsTab",full.names=T)


tmp<-getKLValues(c(1:length(fnamTb)),c(1:18),c(5:15))
