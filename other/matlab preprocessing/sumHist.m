function SumHist = sumHist(SumHist,pref,syno,NN)

if syno==6 && NN>100
    [NHIST,edges]=histcounts(pref,'BinWidth',0.001,'Normalization','probability');
    idbegin=find(edges>=0,1);
    edgesNew=edges(idbegin:end);
    NHISTnew=NHIST(idbegin:end);
    
else
    [NHIST,edges]=histcounts(pref,'BinWidth',0.001);  %%this bin width should agree with function 'getTheoryPlot'
    indbegin=find(edges>=0,1);  %%delete all the negative edges
    edgesNew=edges(indbegin:end);
    pmax=Efor(syno,NN);
    PiOriginal=pmax./exp(edgesNew(2:end).*NN); %% find corresponding orignial Pi value
    NHISTnew=NHIST(indbegin:end).*PiOriginal;
    
end
lnew=length(NHISTnew);  %%compare length for sum
lold=length(SumHist);
if lnew>lold
    SumHist(1:lold)=NHISTnew(1:lold)+SumHist(1:lold);
    SumHist((lold+1):lnew)=NHISTnew((lold+1):lnew);
else
    SumHist(1:lnew)=NHISTnew(1:lnew)+SumHist(1:lnew);
    SumHist((lnew+1):lold)=SumHist((lnew+1):lold);
end

end

