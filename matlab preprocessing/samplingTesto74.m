%%distribution exploration: histogram(pi) first VS histogram(Si);
%%Si=log(pmax./pi)/NN;
syno=6;
NN=74;
itr=100000;

P=zeros(1,syno);
P(:,:)=1/syno;
R=mnrnd(NN,P,itr);

pref1=getPref(NN,syno);
[NHIST,edges]=histcounts(pref1,'BinWidth',0.0001);  %%this bin width should agree with function 'getTheoryPlot'
lhist=length(edges);
if edges(1)<0
    NHISTnew(1)=sum(NHIST(1:2)); %%previously delete all the negative edges;but not good
    NHISTnew(2:(lhist-2))=NHIST(3:(lhist-1));
    edgesNew(1:(lhist-2))=edges(2:(lhist-1)); %% since -0 take up big portion and should consider as 0;-0 due to accuration loss like 4/3*3-1~=0;
    %%% equal to find the midpoint corresponding to NHIST
else
    NHISTnew=NHIST;
    edgesNew(1:(lhist-1))=edges(1:(lhist-1));
end

pmax=Efor(syno,NN);
PiOriginal=pmax./exp(edgesNew.*NN); %% find corresponding orignial Pi value
NHISTfinal=NHISTnew.*PiOriginal;

for i=1:itr
    mnp(i)=mnpdf(R(i,:),P);
end

pref2=log(pmax./mnp)/NN;
[NHIST2,edges2]=histcounts(pref2,'BinWidth',0.0001,'Normalization','probability');
lhist2=length(edges2);
if edges2(1)<0          %%%%%let edgesNew and NHIST match each other
                NHIST2new(1)=sum(NHIST2(1:2));
                NHIST2new(2:(lhist2-2))=NHIST2(3:(lhist2-1));
                edges2New(1:(lhist2-2))=edges2(2:(lhist2-1));
            else
                NHIST2new=NHIST2;
                edges2New(1:(lhist2-1))=edges2(1:(lhist2-1));
end

x=min(length(NHIST2new),length(NHISTfinal));
sumId=find((NHIST2new(1:x).*NHISTfinal(1:x))~=0);
selectP=sum(NHIST2new(sumId).*log(NHIST2new(sumId)./NHISTfinal(sumId))); %%%%Kullback Leibler divergence

figure
plot(edgesNew,NHISTfinal,'.');
title('distribution of partitions of length 74');
xlabel('S value');
ylabel('probability');
figure
plot(edges2New,NHIST2new,'.');
