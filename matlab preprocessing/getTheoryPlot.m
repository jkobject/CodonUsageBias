function [SumHistF,edgesF,ThoeryPercent]=getTheoryPlot(SumHist)
SumHistF=SumHist./sum(SumHist); %%SumHistF final data for plot
edgesF=0.0005:0.001:0.001*length(SumHistF); %% bin width should agree with function 'sumHist'
CumSum=cumsum(SumHistF);%%%%%find 99% percentil
ThoeryPercent=edgesF(find(CumSum>0.99,1));
end