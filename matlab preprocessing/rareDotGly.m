% %%Analysis for Gly 4 synonymous codons
% setSynonymousCodonTable
% 
% global ctG cfG
% % % % 
% clear I Io Ir IDf InRi NN %% important for accuracy
% % 
% % % % % generate and store entropy value first
% load('pasteCodon11.mat');
% load('AveEntropy4f.mat');
% load('sequenceName11.mat');
% 
% 
% 
% for i=1:length(pasteCodon11)
%     
% [NNg,~,Gop]=GlyAminoAcidH(pasteCodon11{1,i}); %%get Ygg for each mRNA
% 
% if ~(isnan(NNg)||NNg>700) %% store original entropy Io, lengthNN
% Grp=AveEntropy4(NNg);
% Gr(i,:)=-Grp/NNg;
% Go(i,:)=-log(Gop)/NNg;
% NN(i,:)=NNg;
% else 
% Go(i,:)=NaN;
% Gr(i,:)=NaN;
% NN(i,:)=NaN;  %% only 3 sequences length over 400
% end
% end
% figure
% plot(Go,Gr,'.y');
% hold on
% % % 
% save 'GlyComparison1.mat' Go Gr NN %% store the entropy value
% % % 
% % % 
% filename='GlyComparison1.mat';
% m=matfile(filename);
% Go=m.Go;
% NN=m.NN;
% % % %%%generate random sequences for x times
% for x=1:10
%     clear G Gr NN InRg%% important for accuracy
%     
% for i=1:length(pasteCodon11)
%     
% [NNg,~,~]=GlyAminoAcidH(pasteCodon11{1,i}); %%get Ygg for each mRNA
% 
% if ~(isnan(NNg)||NNg>700)
% 
% InRg=replaceE(NNg,ctG);%%input random sequence for G
% [~,~,Gp]=GlyAminoAcidH(InRg);
% 
% Grp=AveEntropy4(NNg);
% Gr(i,:)=-Grp/NNg;
% G(i,:)=-log(Gp)/NNg;
% NN(i,:)=NNg;
% else 
% G(i,:)=NaN; 
% Gr(i,:)=NaN;
% NN(i,:)=NaN;
% end
% end
% plot(G,Gr,'.r');
% end
% hold off
% xlabel('per codon entropy of each codon sequence');
% ylabel('per codon expected entropy of sequences after equal codon replacement');
% title(' Gly in Sah6');
% legend('real observed sequence','10 times of randomly generated sequence');
% xlim([0 1.4]);

%%find particular genes

% aimGene='>EJS44836';
% for w=1:length(sequenceName11)
%     tempGene=sequenceName11{w};
% if strcmp(aimGene,tempGene);
%     id=w;
%     break;
% end
% end


l=NN(w);
pmax=Efor(4,l);
phis=pf;
GlyHist=(-log(phis./pmax))/l;
figure
histogram(GlyHist,'Normalization','probability');
hold on
x=Go(w);
y=0;
plot(x,y,'*r');

title('distribution of entropy for 4 synonymous codons in Sah6 at length 31');
xlabel('entropy value for length 31');
ylabel('frequency');
legend('','GeneName:EJS44836');