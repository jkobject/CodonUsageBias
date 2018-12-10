function [Sp,SelWeakId] = SelectionPressure(pasteCodon)  %% Sp: portion of real sequence entropy larger, Sid corresponding index of codon sequence 

global ctG cfG ctE cfE ctA cfA ctP cfP ctD cfD

load('AveEntropy2f.mat')

for i=1:length(pasteCodon) %% note: pasteCodon is column vector
    
disp([num2str(i), 'codon sequence calculation begins']);
    
% [NNe,~,E(i,:)]=GluAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA

[NNa,~,A]=AlaAminoAcidH(pasteCodon{1,i});


if (isnan(NNa) || NNa>700)
% if  isnan(NNa)
    continue;
else   
%     RGT = replaceE(NNg,ctG);  %% replace codon according to fraction of codon table 
%    [~,~,Grp]=GlyAminoAcidH(RGT);
   
   Ao(i,:)=-log(A)/NNa;
%    Gor(i,:)=-log(Grp)/NNg;
   Ar(i,:)=-AveEntropy4(NNa)/NNa;
end
       
end 


% 2 synonymous codons
% global ctD cfD 
% 
% load AveEntropy2f.mat
%  
% for i=1:length(pasteCodon) %% note: pasteCodon is column vector
%     
% disp([num2str(i), 'codon sequence calculation begins']);
%     
% [NNd,~,D(i,:)]=AspAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA
% 
% if isnan(NNd) 
%     continue;
% else   
% Do(i,:)=log(D(i,:))/NNd;    %% o: original; r:replacement
% Dr(i,:)=AveEntropy2(NNd)/NNd;
% end
% 
% end



%% figure and table
figure    
%     
plot(Ao,Ar,'.');
% 
% xlim([-0.9,0]);
% % 
% hold on
%  
% x=-0.1:0.05:0;
% 
% y=x;
% 
% plot(x,y);
% 
% hold off
% 
% set(gca,'xaxislocation','bottom','yaxislocation','left','xdir','reverse','ydir','reverse')
% 
% xlabel('per codon entropy of original sequence');
%  
% ylabel('per codon average entropy of sequence after equal codon replacement');
% 
% title(['Species--sah6: ','Amino Acid--Asp(4 synonymous)']);

SelDif=Ao-Ar;                       %% SelDif: selection difference; SelWeakId: weak selection dots ID 
SelWeakId=find(SelDif>0);
Sp=length(SelWeakId)/length(SelDif);  %% Sp: weak selection portion
    
    
end