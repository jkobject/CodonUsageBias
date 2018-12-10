
#This is a variant of mkDistrib and makes it possible.
#To look at the distribution of individual amino acids inbteractively

##USE interactiveFitting.mw with maple to visualize the fit!
##The maple files will read the correct dta.csv from /tmp

#"E" "H" "Q" "F" "Y" "C" "N" "K" "D" 

system("mv nonlinfit.txt /tmp")
source("analysisFunctions.Rscript")
#Finally write the first line of the file
#com4<-'echo "species,aa,a,b,c,residual" > nonlinfit.txt '
#system(com4)

len<-15
spec<-33
bst<-"F"
#now produce the correct maple file by setting in the correct length
com3<-paste('sed "s/xxx/N:=',len,';/g" fitModel.xxx > fitModel.mpl',sep="")
system(com3) 

print(spec)
a<-mkPlotC(bacnamTb[spec],bst,len,T)
dta<-cbind(a[["k"]],a[["Freq"]])
write.table(dta,file="/tmp/dta.csv",quote=F,row.names=F,sep=",",col.names=F)
sink("nonlinfit.txt", append=TRUE)
cat(spec,",",bst,",")
sink()
system("maple fitModel.mpl")
filnam<- paste("nonlinfit",len,".txt",sep="")
com2 <- paste("mv nonlinfit.txt",filnam,sep=" ")
system(com2)



#cleaning up errors in the files
system("bash cleannonlinfiles.sh");
a<-read.csv("/tmp/dta.csv")
a;
plot(a)