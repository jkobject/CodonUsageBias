% % % % Equally Random sequence entropy comparison vs Expected entropy
% % % % Not complete code
% 
% setSynonymousCodonTable
% 
% global ctD cfD ctR cfR
% 
% 
% load('AveEntropy2f.mat');
% load('AveEntropy6.mat');
% load('pasteCodon21.mat'); %% species Sp
% load('pasteCodon32.mat'); %% species afl
 %load('pasteCodon12.mat'); %% species saeu

clear D Dr InRd  %% important for accuracy

% % % generate and store entropy value first
for i=1:length(pasteCodon32)
    
[NNd,~,~]=AspAminoAcidH(pasteCodon32{1,i}); %%get Ygg for each mRNA

if ~(isnan(NNd)||NNd>1500)

InRd=replaceE(NNd,ctD);%%input random sequence for D
[~,~,Dp]=AspAminoAcidH(InRd);

Drp=AveEntropy2(NNd);
Dr(i,:)=-Drp/NNd;
D(i,:)=-log(Dp)/NNd;

else 
D(i,:)=NaN; 
Dr(i,:)=NaN;
end
end 

save 'AspErComparison10.mat' D Dr %% store the entropy value


clear R Rr InRr  
for i=1:length(pasteCodon32)
    
[NNr,~,~]=ArgAminoAcidH(pasteCodon32{1,i}); %%get Ygg for each mRNA

if ~(isnan(NNr)||NNr>500)

InRr=replaceE(NNr,ctR);%%input random sequence for D
[~,~,Rp]=ArgAminoAcidH(InRr);

Rrp=AveEntropy6(NNr);
Rr(i,:)=-Rrp/NNr;
R(i,:)=-log(Rp)/NNr;

else 
R(i,:)=NaN; 
Rr(i,:)=NaN;
end
end 

save 'ArgErComparison32.mat' R Rr

% plots for Equally Random artificial sequence: Asp Arg in saeu, afl:
fileName1='AspErComparison6.mat';
m1=matfile(fileName1);
D1=m1.D;
Dr1=m1.Dr;


fileName2='AspErComparison10.mat';
m2=matfile(fileName2);
D2=m2.D;
Dr2=m2.Dr;

fileName3='ArgErComparison2.mat';
m3=matfile(fileName3);
R1=m3.R;
Rr1=m3.Rr;


fileName4='ArgErComparison10.mat';
m4=matfile(fileName4);
R2=m4.R;
Rr2=m4.Rr;

% % Asp for Sp
Ddis1=(D1(~isnan(D1)))./(Dr1(~isnan(Dr1)));
n1bins=floor(length(D1)/15);

figure
histogram(Ddis1,n1bins);
xlabel('ratio between artificial sequence entropy and expected entropy');
ylabel('frequency');
title('Distribution: Equal Random Artificial Sequence: Schizzosaccharomyces pombe--Asp(2sysnonymous)');

%Asp in afl
Ddis2=(D2(~isnan(D2)))./(Dr2(~isnan(Dr2)));
n2bins=floor(length(D2)/15);

figure 
histogram(Ddis2,n2bins);
xlabel('ratio between artificial sequence entropy and expected entropy');
ylabel('frequency');
title('Distribution: Equal Random Artificial Sequence: Aspergillus flavus--Asp(2sysnonymous)');

%histogram Arg in saeu
Rdis1=(R1(~isnan(R1)))./(Rr1(~isnan(Rr1)));
n3bins=floor(length(R1)/15);

figure
histogram(Rdis1,n3bins);
xlim([0 7]);
xlabel('ratio between artificial sequence entropy and expected entropy');
ylabel('frequency');
title('Distribution: Equal Random Artificial Sequence: Saccharomyces eubayanus--Arg(6sysnonymous)');

% Arg in afl
Rdis2=(R2(~isnan(R2)))./(Rr2(~isnan(Rr2)));
n4bins=floor(length(R2)/15);

% % entropy histogram plot
figure 
histogram(Rdis2,n4bins);
xlim([0 7]);
xlabel('ratio between artificial sequence entropy and expected entropy');
ylabel('frequency');
title('Distribution: Equal Random Artificial Sequence: Aspergillus flavus--Arg(6sysnonymous)');


%%entropy comparison plot
figure    
    
plot(R2,Rr2,'.');

hold on
 
x=(0:0.1:0.3); 
y=x;

plot(x,y);

hold off

xlim([0 7]);


xlabel('per codon entropy of artificial sequence');
 
ylabel('per codon expected entropy of sequence after equal codon replacement');

title('Equal Random Artificial Sequence: Aspergillus flavus--Arg(6sysnonymous)');
