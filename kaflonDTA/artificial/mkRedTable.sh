#This file takes an Yun-style synotable
#produces a synotable with only the first 9 codons (where codon number is 2)
#outputs the cub where the number is smaller than 0.5
#the output file is tmp, so you need to rename
echo "E,H,Q,F,Y,C,N,K,D" > tmpcub
tail -n +2 $1 > tmp
cut -f1,2,4,6,8,10,12,14,16,18 -d"," tmp  >> tmpcub
echo "species, E,H,Q,F,Y,C,N,K,D" > tmp
r -e 'a <- read.table("tmpcub",sep=",",as.is=T) ; b<-apply(a,2,function(x) ifelse(x<0.5,x,1-x)) ;write.table(file="tmp",b,col.names=F,row.names=T,quote=F,sep=",",append=T)'
