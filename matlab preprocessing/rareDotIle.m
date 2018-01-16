% setSynonymousCodonTable
% 
% global ctI cfI
% % % % 
% clear I Io Ir IDf InRi NN %% important for accuracy
% % 
% % % % % generate and store entropy value first
% load('pasteCodon11.mat');
% load('AveEntropy3f.mat');
% load('sequenceName11.mat');
% load('AveEntropy3f1t200.mat');
% load('AveEntropy3f201t400.mat');
% 
% 
% for i=1:length(pasteCodon11)
%     
% [NNi,~,Iop]=IleAminoAcidH(pasteCodon11{1,i}); %%get Ygg for each mRNA
% 
% if ~(isnan(NNi)||NNi>1000) %% store original entropy Io, lengthNN
% Irp=AveEntropy3(NNi);
% Ir(i,:)=-Irp/NNi;
% Io(i,:)=-log(Iop)/NNi;
% NN(i,:)=NNi;
% else 
% Io(i,:)=NaN;
% Ir(i,:)=NaN;
% NN(i,:)=NaN;
% end
% end
% figure
% plot(Io,Ir,'.y');
% hold on
% % % 
% save 'IleComparison1.mat' Io Ir NN %% store the entropy value
% % % 
% % % 
% filename='IleComparison1.mat';
% m=matfile(filename);
% Io=m.Io;
% NN=m.NN;
% % % %%%generate random sequences for x times
% for x=1:10
%     clear I Ir IDf InRi NN %% important for accuracy
%     
% for i=1:length(pasteCodon11)
%     
% [NNi,~,~]=IleAminoAcidH(pasteCodon11{1,i}); %%get Ygg for each mRNA
% 
% if ~(isnan(NNi)||NNi>1000)
% 
% InRi=replaceE(NNi,ctI);%%input random sequence for D
% [~,~,Ip]=IleAminoAcidH(InRi);
% 
% Irp=AveEntropy3(NNi);
% Ir(i,:)=-Irp/NNi;
% I(i,:)=-log(Ip)/NNi;
% NN(i,:)=NNi;
% else 
% I(i,:)=NaN; 
% Ir(i,:)=NaN;
% NN(i,:)=NaN;
% end
% end
% plot(I,Ir,'.r');
% end
% hold off
% xlabel('per codon entropy of each codon sequence');
% ylabel('per codon expected entropy of sequences after equal codon replacement');
% title(' Ile in Sah6');
% legend('real observed sequence','10 times of randomly generated sequence');

clear NNhist
% NNhist=NN(NN<20);

% for lhist=1:length(NNhist)
 
clear pmax phis ileHis idSpecial VSpecial idInterest GeneNameInterest x y
Nhist=32;
% Nhist=NNhist(lhist);  %% Nhist: find sequence of particular length whether it is small probability
pmax=Efor(3,Nhist);
if Nhist>200
phis=pf3f201t400{1,Nhist};
else
phis=pf3f1t200{1,Nhist};
end
ileHist=(-log(phis./pmax))/Nhist;

idSpecial=find(NN==Nhist);
VSpecial=Io(idSpecial); 
x=max(VSpecial);
if x>0.7           %%bug: same length will plot the highly selected each time, code here need modificaiton
figure
histogram(ileHist,'Normalization','probability');
% xlim([0 1.2]);
title(['distribution of entropy for 3 synonymous codons in Sah6 at length',num2str(Nhist)]);
xlabel(['entropy value for length',num2str(Nhist)]);
ylabel('frequency');
hold on
y=zeros(length(x));
plot(x,y,'*r');  %% mark the highly selected species
hold off
idInterest=find(Io==x);
GeneNameInterest=sequenceName11(idInterest);
legend('',sprintf('%s',GeneNameInterest{1}));
end

% end