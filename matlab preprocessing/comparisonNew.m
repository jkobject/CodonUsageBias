% setSynonymousCodonTable;
% load 'pasteCodon11.mat';
% pasteCodon=pasteCodon11;
% global ctG
% for i=1:length(pasteCodon) %% note: pasteCodon is column vector
% [NNg,~,~]=GlyAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA
% if ~(isnan(NNg)||NNg>300)
% InGd=replaceE(NNg,ctG);
% [~,~,Gp]=GlyAminoAcidH(InGd); %%get Ygg for each mRNA
% G(i)=-log(Gp)/NNg;
% else
% G(i)=NaN;
% end
% end
% [NHISTrealRandom,EdgesRealRandom]=getRealPlot(G);
% CumSum=cumsum(SumHistF);%%%%%find 95% percentil
% ThoeryPercent=edgesF(find(CumSum>0.99,1));
% 
figure
plot(edgesF,SumHistF,'.');
plot(edgesF,SumHistF,'g.');
hold on
plot(EdgesReal,NHISTreal,'r.');
plot(ThoeryPercent,0,'b*');
hold off
title('distribution of miexed lengths for Gly in sah6');
xlabel('entropy value');
ylabel('probability');
legend('theoretical','real','99%percentil','random generated 10 times');