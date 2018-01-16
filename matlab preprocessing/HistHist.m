% setSynonymousCodonTable;
% load('pasteCodon12.mat')
% pasteCodon=pasteCodon12;
% 
% 
% % 6 syno
% %%%%%%%%%%%%%%theoretical reference calculation
% SumHist=zeros(1,1);
% for i=1:length(NNr12)
%     NN=NNr12(i);%
%     if NN<=400
%         pref=getPref(NN,6); %% remember to change syno value;  read all the Entropy value for possible configurations--pref
%         SumHist=sumHist(SumHist,pref);%% SumHist resore sum value
%     end
% end
% [SumHistF,edgesF,ThoeryPercent]=getTheoryPlot(SumHist);
% 
% %%%%%%%%%%%%%real entropy distribution calculation
% for i=1:length(pasteCodon) %% note: pasteCodon is column vector
%     [NNr,~,Rp]=ArgAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA
%     if NNr<=400
%         R(i)=-log(Rp)/NNr;
%     else
%         R(i)=NaN;
%     end
% end
% [NHISTreal,EdgesReal]=getRealPlot(R);
% 
% figure
% plot(edgesF,SumHistF,'g.');
% hold on
% plot(EdgesReal,NHISTreal,'r.');
% plot(ThoeryPercent,0,'b*');
% hold off
% title('distribution of miexed lengths for Arg in saeu');
% xlabel('entropy value');
% ylabel('probability');
% legend('theoretical','real','99%percentil');
% 
% save 'ArgAccurate12.mat' SumHistF edgesF NHISTreal EdgesReal



%4syno
% %%%%%%%%%%%%%%theoretical reference calculation
% SumHist=zeros(1,1);
% for i=1:length(NNg)
%     NN=NNg(i);
%     if NN<=300
%     pref=getPref(NN,4); %%read all the Entropy value for possible configurations--pref 
%     SumHist=sumHist(SumHist,pref);%% SumHist resore sum value
%     end
% end
% [SumHistF,edgesF,ThoeryPercent]=getTheoryPlot(SumHist);

%%%%%%%%%%%%%real entropy distribution calculation
% for i=1:length(pasteCodon) %% note: pasteCodon is column vector
% [NNgs,~,Gp]=GlyAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA
% if NNgs<=300
% G(i)=-log(Gp)/NNgs;
% else
% G(i)=NaN;
% end
% end
% [NHISTreal,EdgesReal]=getRealPlot(G);
% 
% figure
% plot(edgesF,SumHistF,'g.');
% hold on
% plot(EdgesReal,NHISTreal,'r.');
% plot(ThoeryPercent,0,'b*');
% hold off
% title('distribution of miexed lengths for Gly in saeu');
% xlabel('entropy value');
% ylabel('probability');
% legend('theoretical','real','99%percentil','random generated 10 times');
% 
% save 'GlyAccurate12.mat' SumHistF edgesF NHISTreal EdgesReal
%%%%%%generate random sequences for comparison
% % % global ctG
% % % for i=1:length(pasteCodon) %% note: pasteCodon is column vector
% % % [NNg,~,~]=GlyAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA
% % % if ~(isnan(NNg)||NNg>300)
% % % InGd=replaceE(NNg,ctG);
% % % [~,~,Gp]=GlyAminoAcidH(InGd); %%get Ygg for each mRNA
% % % G(i)=-log(Gp)/NNg;
% % % else
% % % G(i)=NaN;
% % % end
% % % end
% % % [NHISTrealRandom,EdgesRealRandom]=getRealPlot(G);
% % % 
% % % plot(EdgesRealRandom,NHISTrealRandom,'y.');
%



% 3 syno
%%%%%%%%%%%%%%theoretical reference calculation
SumHist=zeros(1,1);
for i=1:length(NNi12)
    NN=NNi12(i);
    if NN<=500
    pref=getPref(NN,3); %% remember to change syno value;  read all the Entropy value for possible configurations--pref 
    SumHist=sumHist(SumHist,pref);%% SumHist resore sum value
    end
end
[SumHistF,edgesF,ThoeryPercent]=getTheoryPlot(SumHist);

%%%%%%%%real entropy distribution calculation
for i=1:length(pasteCodon) %% note: pasteCodon is column vector  
[NNi,~,Ip]=IleAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA
if NNi<=500
I(i)=-log(Ip)/NNi;
else
I(i)=NaN;
end
end
[NHISTreal,EdgesReal]=getRealPlot(I);

figure
plot(edgesF,SumHistF,'g.');
hold on
plot(EdgesReal,NHISTreal,'r.');
plot(ThoeryPercent,0,'b*');
hold off
title('distribution of miexed lengths for Ile in saeu');
xlabel('entropy value');
ylabel('probability');
legend('theoretical','real','99%percentil');

save 'IleAccurate12.mat' SumHistF edgesF NHISTreal EdgesReal



%2syno
%%%%%%%%%%%%%%theoretical reference calculation
SumHist=zeros(1,1);
for i=1:length(NNd)
     NN=NNd(i);
     if NN<=800
     pref=getPref(NN,2); %% function 'getPfre', read all the Entropy value for possible configurations--pref 
     SumHist=sumHist(SumHist,pref);%% function 'sumHist' resore sum value
     end
end
[SumHistF,edgesF,ThoeryPercent]=getTheoryPlot(SumHist); %%function 'getTheoryPlot' calculate data for plot

%%%%%%%%%%%%real entropy distribution calculation
for i=1:length(pasteCodon) %% note: pasteCodon is column vector
[NNds,~,Dp]=AspAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA
if NNds<=800
D(i)=-log(Dp)/NNds;
else
D(i)=NaN;
end
end
[NHISTreal,EdgesReal]=getRealPlot(D);  %% function 'getRealPlot' calculate data for plot

figure
plot(edgesF,SumHistF,'g.');
hold on
plot(EdgesReal,NHISTreal,'r.');
plot(ThoeryPercent,0,'b*');
hold off
title('distribution of miexed lengths for Asp in saeu');
xlabel('entropy value');
ylabel('probability');
legend('theoretical','real','99%percentil');

save 'AspAccurate12.mat' SumHistF edgesF NHISTreal EdgesReal