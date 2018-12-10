%%to explain the 'v' patten in entropy comparison plot, short randomly
%%sequence 10000 times 

setSynonymousCodonTable

global ctD cfD
fileName1='AspComparison6.mat';
m1=matfile(fileName1);
D1=m1.D;
Dr1=m1.Dr;
figure
plot(D1,Dr1,'.y')
hold on


load('AveEntropy2f.mat');
NNd=20;
clear D Dr InRd
for u=1:10000
InRd=replaceE(NNd,ctD);%%input random sequence for D
[~,~,Dp]=AspAminoAcidH(InRd);

Drp=AveEntropy2(NNd);
D(u,:)=-log(Dp)/NNd;
Dr(u,:)=-Drp/NNd;
end

plot(D,Dr,'.r')

NNd=10;
clear D Dr InRd
for u=1:10000
InRd=replaceE(NNd,ctD);%%input random sequence for D
[~,~,Dp]=AspAminoAcidH(InRd);

Drp=AveEntropy2(NNd);
D(u,:)=-log(Dp)/NNd;
Dr(u,:)=-Drp/NNd;
end
plot(D,Dr,'.g')

NNd=5;
clear D Dr InRd
for u=1:10000
InRd=replaceE(NNd,ctD);%%input random sequence for D
[~,~,Dp]=AspAminoAcidH(InRd);

Drp=AveEntropy2(NNd);
D(u,:)=-log(Dp)/NNd;
Dr(u,:)=-Drp/NNd;
end
plot(D,Dr,'.b')
