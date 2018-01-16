setSynonymousCodonTable

load pasteCodon11.mat
%%real Gly11
for i=1:length(pasteCodon11) %% note: pasteCodon is column vector
    
[NNg,~,Gp]=GlyAminoAcidH(pasteCodon11{1,i}); %%get Ygg for each mRNA

if NNg<=300
G(i)=-log(Gp)/NNg;
else
G(i)=NaN;
end
end
Gdone1=G(~isnan(G));
Gdone2=Gdone1(Gdone1>0);
% figure
% histogram(Gdone2,'BinWidth',0.01,'Normalization','Probability');
% title('Observed Entropy Distribution for Gly in Sah6');
% xlabel('Entropy');
% ylabel('probability');

% %%real Asp11
%  for i=1:length(pasteCodon11) %% note: pasteCodon is column vector
%     
% [NNd,~,Dp]=GlyAminoAcidH(pasteCodon11{1,i}); %%get Ygg for each mRNA
% 
% if NNd<=800
% D(i)=-log(Dp)/NNd;
% else
% D(i)=NaN;
% end
% end
% Ddone1=D(~isnan(D));
% Ddone2=Ddone1(Ddone1>0);
% figure
% histogram(Ddone2,'BinWidth',0.01,'Normalization','Probability');
% title('Observed Entropy Distribution for Asp in Sah6');
% xlabel('Entropy');
% ylabel('probability');

%%real Ile11
 for i=1:length(pasteCodon11) %% note: pasteCodon is column vector
    
[NNi,~,Ip]=IleAminoAcidH(pasteCodon11{1,i}); %%get Ygg for each mRNA
NN=NNi;
if NNi<=500
I(i)=-log(Ip)/NN;
else
I(i)=NaN;
end
end
Ddone1=I(~isnan(I));
Ddone2=Ddone1(Ddone1>0);
figure
histogram(Ddone2,'BinWidth',0.001,'Normalization','Probability');
title('Observed Entropy Distribution for Ile in Sah6');
xlabel('Entropy');
ylabel('probability');

[NHISTreal,EdgesRealp]=histcounts(Ddone2,'BinWidth',0.001,'Normalization','Probability'); 
EdgesReal=EdgesRealp(2:end);
%clear I Ddone1 Ddone2