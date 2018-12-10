%%%%%compare difference between hist pi or pref
syno=6;
NN=74;
itr=100000;
pmax=Efor(syno,NN);

P=zeros(1,syno);
P(:,:)=1/syno;
R=mnrnd(NN,P,itr);

for i=1:itr
    mnp(i)=mnpdf(R(i,:),P);
end

[NHIST2,edges2]=histcounts(mnp,'BinWidth',0.00000001,'Normalization','cdf');
x2=log(pmax./edges2(2:end))/NN;
y2=NHIST2;
figure
plot(x2,y2,'.');

pref2=log(pmax./mnp)/NN;
[NHIST2,edges2]=histcounts(pref2,'BinWidth',0.001,'Normalization','cdf');
idbegin2=find(edges2>=0,1);
edges2New=edges2(idbegin2:end);
NHIST2new=NHIST2(idbegin2:end);

figure
plot(edges2New(2:end),NHIST2new,'.');
