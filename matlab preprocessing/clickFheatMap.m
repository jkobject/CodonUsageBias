% % click button get 18 amino acids 'selection pressure' for one species
% % cd '/Users/yd46/Documents/Matlab'
setSynonymousCodonTable

global cfG cfGT ctG cfD cfDT ctD cfE cfET ctE cfV cfVT ctV cfA cfAT ctA cfR cfRT ctR cfS cfST ctS cfK cfKT ctK cfN cfNT ctN cfM cfMT ctM cfI cfIT ctI cfT cfTT ctT cfW cfWT ctW cfZ cfZT ctZ cfC cfCT ctC cfY cfYT  ctY cfL cfLT ctL cfF cfFT ctF cfQ cfQT ctQ cfH cfHT ctH cfP cfPT ctP;

load('AveEntropy2f.mat');
load('AveEntropy3f.mat');
load('AveEntropy4f.mat');
load('AveEntropy6.mat');

m=matfile('speciesName.mat');
speciesName=m.speciesName; %% get all the species name, benefit further clear expression

for tP=1:length(speciesName)
    
fileName=sprintf('%s',[speciesName{tP},'.fa']);

pasteCodon=getCodonSequence(fileName);
%%%%%%Ile begins
% for i=1:length(pasteCodon) %% note: pasteCodon is column vector        
% [NNi,~,Ip]=IleAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA
% if ~(isnan(NNi)||NNi>1000) 
% Irp=AveEntropy3(NNi);
% Ir(i,:)=-Irp/NNi;
% I(i,:)=-log(Ip)/NNi;
% else 
% I(i,:)=NaN;   
% Ir(i,:)=NaN;
% end 
% end
% 
% SelDif=I(~isnan(I))-Ir(~isnan(Ir));                       %% SelDif: selection difference; SelWeakId: weak selection dots ID 
% SelLid=find(SelDif<0);
% Ltest=abs(sum(SelDif(SelLid)));
% SelRid=find(SelDif>0);
% Rtest=sum(SelDif(SelRid));
% Sptest=(Ltest-Rtest)/length(SelDif);   
% 
% fileID1=fopen('Ile.txt','a');
% fprintf(fileID1,'%s,%d\n',speciesName{tP},Sptest);
% fclose(fileID1);
% 
% clear I Ir SelDif SelLid SelRid

%%%%% Asp begins 
for i=1:length(pasteCodon)
[NNd,~,Dp]=AspAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA
if ~(isnan(NNd)||NNd>1500) 
Drp=AveEntropy2(NNd);
Dr(i,:)=-Drp/NNd;
D(i,:)=-log(Dp)/NNd;
else 
D(i,:)=NaN;   
Dr(i,:)=NaN;
end 
end
SelDif=D(~isnan(D))-Dr(~isnan(Dr));                       %% SelDif: selection difference; SelWeakId: weak selection dots ID 
SelLid=find(SelDif<0);
Ltest=abs(sum(SelDif(SelLid)));
SelRid=find(SelDif>0);
Rtest=sum(SelDif(SelRid));
Sptest=(Ltest-Rtest)/length(SelDif); 

fileID2=fopen('Asp.txt','a');
fprintf(fileID2,'%s,%d\n',speciesName{tP},Sptest);
fclose(fileID2);
clear D Dr SelDif SelLid SelRid

%%%%%%%Glu begins
for i=1:length(pasteCodon) %% note: pasteCodon is column vector        
[NNe,~,Ep]=GluAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA
if ~(isnan(NNe)||NNe>1500) 
Erp=AveEntropy2(NNe);
Er(i,:)=-Erp/NNe;
E(i,:)=-log(Ep)/NNe;
else 
E(i,:)=NaN;   
Er(i,:)=NaN;
end 
end

SelDif=E(~isnan(E))-Er(~isnan(Er));                       %% SelDif: selection difference; SelWeakId: weak selection dots ID 
SelLid=find(SelDif<0);
Ltest=abs(sum(SelDif(SelLid)));
SelRid=find(SelDif>0);
Rtest=sum(SelDif(SelRid));
Sptest=(Ltest-Rtest)/length(SelDif); 

fileID3=fopen('Glu.txt','a');
fprintf(fileID3,'%s,%d\n',speciesName{tP},Sptest);
fclose(fileID3);
clear E Er SelDif SelLid SelRid

%%% His begins
for i=1:length(pasteCodon)
[NNh,~,Hp]=HisAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA
if ~(isnan(NNh)||NNh>1500)
Hrp=AveEntropy2(NNh);
Hr(i,:)=-Hrp/NNh;
H(i,:)=-log(Hp)/NNh;
else 
H(i,:)=NaN;    
Hr(i,:)=NaN;
end 
end

SelDif=H(~isnan(H))-Hr(~isnan(Hr));                       %% SelDif: selection difference; SelWeakId: weak selection dots ID 
SelLid=find(SelDif<0);
Ltest=abs(sum(SelDif(SelLid)));
SelRid=find(SelDif>0);
Rtest=sum(SelDif(SelRid));
Sptest=(Ltest-Rtest)/length(SelDif);

fileID4=fopen('His.txt','a');
fprintf(fileID4,'%s,%d\n',speciesName{tP},Sptest);
fclose(fileID4);
clear H Hr SelDif SelLid SelRid


%%% Pro begins
for i=1:length(pasteCodon)
[NNp,~,Pp]=ProAminoAcidH(pasteCodon{1,i});
if ~(isnan(NNp)||NNp>700) 
Prp=AveEntropy4(NNp);
Pr(i,:)=-Prp/NNp;
P(i,:)=-log(Pp)/NNp;
else 
P(i,:)=NaN;   
Pr(i,:)=NaN;
end 
end

SelDif=P(~isnan(P))-Pr(~isnan(Pr));                       %% SelDif: selection difference; SelWeakId: weak selection dots ID 
SelLid=find(SelDif<0);
Ltest=abs(sum(SelDif(SelLid)));
SelRid=find(SelDif>0);
Rtest=sum(SelDif(SelRid));
Sptest=(Ltest-Rtest)/length(SelDif);  

fileID5=fopen('Pro.txt','a');
fprintf(fileID5,'%s,%d\n',speciesName{tP},Sptest);
fclose(fileID5);
clear P Pr SelDif SelLid SelRid


%%% Gly begins
for i=1:length(pasteCodon)  
[NNg,~,Gp]=GlyAminoAcidH(pasteCodon{1,i}); 
if ~(isnan(NNg)||NNg>700)
Grp=AveEntropy4(NNg);
Gr(i,:)=-Grp/NNg;
G(i,:)=-log(Gp)/NNg;
else 
G(i,:)=NaN;    
Gr(i,:)=NaN;
end 
end

SelDif=G(~isnan(G))-Gr(~isnan(Gr));                       %% SelDif: selection difference; SelWeakId: weak selection dots ID 
SelLid=find(SelDif<0);
Ltest=abs(sum(SelDif(SelLid)));
SelRid=find(SelDif>0);
Rtest=sum(SelDif(SelRid));
Sptest=(Ltest-Rtest)/length(SelDif);  

fileID6=fopen('Gly.txt','a');
fprintf(fileID6,'%s,%d\n',speciesName{tP},Sptest);
fclose(fileID6);
clear G Gr SelDif SelLid SelRid

%%% Val begins
for i=1:length(pasteCodon) %% note: pasteCodon is column vector        
[NNv,~,Vp]=ValAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA
if ~(isnan(NNv)||NNv>700) 
Vrp=AveEntropy4(NNv);
Vr(i,:)=-Vrp/NNv;
V(i,:)=-log(Vp)/NNv;
else 
V(i,:)=NaN;   
Vr(i,:)=NaN;
end 
end

SelDif=V(~isnan(V))-Vr(~isnan(Vr));                       %% SelDif: selection difference; SelWeakId: weak selection dots ID 
SelLid=find(SelDif<0);
Ltest=abs(sum(SelDif(SelLid)));
SelRid=find(SelDif>0);
Rtest=sum(SelDif(SelRid));
Sptest=(Ltest-Rtest)/length(SelDif);   

fileID7=fopen('Val.txt','a');
fprintf(fileID7,'%s,%d\n',speciesName{tP},Sptest);
fclose(fileID7);
clear V Vr SelDif SelLid SelRid


%%% Lys begins
for i=1:length(pasteCodon) %% note: pasteCodon is column vector        
[NNk,~,Kp]=LysAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA
if ~(isnan(NNk)||NNk>1500) 
Krp=AveEntropy2(NNk);
Kr(i,:)=-Krp/NNk;
K(i,:)=-log(Kp)/NNk;
else 
K(i,:)=NaN;   
Kr(i,:)=NaN;
end 
end

SelDif=K(~isnan(K))-Kr(~isnan(Kr));                       %% SelDif: selection difference; SelWeakId: weak selection dots ID 
SelLid=find(SelDif<0);
Ltest=abs(sum(SelDif(SelLid)));
SelRid=find(SelDif>0);
Rtest=sum(SelDif(SelRid));
Sptest=(Ltest-Rtest)/length(SelDif);   

fileID8=fopen('Lys.txt','a');
fprintf(fileID8,'%s,%d\n',speciesName{tP},Sptest);
fclose(fileID8);
clear K Kr SelDif SelLid SelRid


%%% Asn begins
for i=1:length(pasteCodon) %% note: pasteCodon is column vector   
[NNn,~,Np]=AsnAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA
if ~(isnan(NNn)||NNn>1500) 
Nrp=AveEntropy2(NNn);
Nr(i,:)=-Nrp/NNn;
N(i,:)=-log(Np)/NNn;
else 
N(i,:)=NaN;   
Nr(i,:)=NaN;
end 
end

SelDif=N(~isnan(N))-Nr(~isnan(Nr));                       %% SelDif: selection difference; SelWeakId: weak selection dots ID 
SelLid=find(SelDif<0);
Ltest=abs(sum(SelDif(SelLid)));
SelRid=find(SelDif>0);
Rtest=sum(SelDif(SelRid));
Sptest=(Ltest-Rtest)/length(SelDif);   

fileID9=fopen('Asn.txt','a');
fprintf(fileID9,'%s,%d\n',speciesName{tP},Sptest);
fclose(fileID9);
clear N Nr SelDif SelLid SelRid


%%%%% Thr begins
for i=1:length(pasteCodon) %% note: pasteCodon is column vector     
[NNt,~,Tp]=ThrAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA
if ~(isnan(NNt)||NNt>700) 
Trp=AveEntropy4(NNt);
Tr(i,:)=-Trp/NNt;
T(i,:)=-log(Tp)/NNt;
else 
T(i,:)=NaN;   
Tr(i,:)=NaN;
end 
end

SelDif=T(~isnan(T))-Tr(~isnan(Tr));                       %% SelDif: selection difference; SelWeakId: weak selection dots ID 
SelLid=find(SelDif<0);
Ltest=abs(sum(SelDif(SelLid)));
SelRid=find(SelDif>0);
Rtest=sum(SelDif(SelRid));
Sptest=(Ltest-Rtest)/length(SelDif); 

fileID10=fopen('Thr.txt','a');
fprintf(fileID10,'%s,%d\n',speciesName{tP},Sptest);
fclose(fileID10);
clear T Tr SelDif SelLid SelRid


%%%% Cys begins
for i=1:length(pasteCodon) %% note: pasteCodon is column vector      
[NNc,~,Cp]=CysAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA
if ~(isnan(NNc)||NNc>1500) 
Crp=AveEntropy2(NNc);
Cr(i,:)=-Crp/NNc;
C(i,:)=-log(Cp)/NNc;
else 
C(i,:)=NaN;   
Cr(i,:)=NaN;
end 
end

SelDif=C(~isnan(C))-Cr(~isnan(Cr));                       %% SelDif: selection difference; SelWeakId: weak selection dots ID 
SelLid=find(SelDif<0);
Ltest=abs(sum(SelDif(SelLid)));
SelRid=find(SelDif>0);
Rtest=sum(SelDif(SelRid));
Sptest=(Ltest-Rtest)/length(SelDif);   

fileID11=fopen('Cys.txt','a');
fprintf(fileID11,'%s,%d\n',speciesName{tP},Sptest);
fclose(fileID11);
clear C Cr SelDif SelLid SelRid


%%%%% Tyr beins
for i=1:length(pasteCodon) %% note: pasteCodon is column vector   
[NNy,~,Yp]=TyrAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA
if ~(isnan(NNy)||NNy>1500) 
Yrp=AveEntropy2(NNy);
Yr(i,:)=-Yrp/NNy;
Y(i,:)=-log(Yp)/NNy;
else 
Y(i,:)=NaN;   
Yr(i,:)=NaN;
end 
end

SelDif=Y(~isnan(Y))-Yr(~isnan(Yr));                       %% SelDif: selection difference; SelWeakId: weak selection dots ID 
SelLid=find(SelDif<0);
Ltest=abs(sum(SelDif(SelLid)));
SelRid=find(SelDif>0);
Rtest=sum(SelDif(SelRid));
Sptest=(Ltest-Rtest)/length(SelDif);   

fileID12=fopen('Tyr.txt','a');
fprintf(fileID12,'%s,%d\n',speciesName{tP},Sptest);
fclose(fileID12);
clear Y Yr SelDif SelLid SelRid


%%%%%Phe begins
for i=1:length(pasteCodon) %% note: pasteCodon is column vector    
[NNf,~,Fp]=PheAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA
if ~(isnan(NNf)||NNf>1500) 
Frp=AveEntropy2(NNf);
Fr(i,:)=-Frp/NNf;
F(i,:)=-log(Fp)/NNf;
else 
F(i,:)=NaN;   
Fr(i,:)=NaN;
end 
end

SelDif=F(~isnan(F))-Fr(~isnan(Fr));                       %% SelDif: selection difference; SelWeakId: weak selection dots ID 
SelLid=find(SelDif<0);
Ltest=abs(sum(SelDif(SelLid)));
SelRid=find(SelDif>0);
Rtest=sum(SelDif(SelRid));
Sptest=(Ltest-Rtest)/length(SelDif);   

fileID13=fopen('Phe.txt','a');
fprintf(fileID13,'%s,%d\n',speciesName{tP},Sptest);
fclose(fileID13);
clear F Fr SelDif SelLid SelRid



%%%%%Gln begins
for i=1:length(pasteCodon) %% note: pasteCodon is column vector  
[NNq,~,Qp]=GlnAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA
if ~(isnan(NNq)||NNq>1500) 
Qrp=AveEntropy2(NNq);
Qr(i,:)=-Qrp/NNq;
Q(i,:)=-log(Qp)/NNq;
else 
Q(i,:)=NaN;   
Qr(i,:)=NaN;
end 
end
SelDif=Q(~isnan(Q))-Qr(~isnan(Qr));                       %% SelDif: selection difference; SelWeakId: weak selection dots ID 
SelLid=find(SelDif<0);
Ltest=abs(sum(SelDif(SelLid)));
SelRid=find(SelDif>0);
Rtest=sum(SelDif(SelRid));
Sptest=(Ltest-Rtest)/length(SelDif); 

fileID14=fopen('Gln.txt','a');
fprintf(fileID14,'%s,%d\n',speciesName{tP},Sptest);
fclose(fileID14);
clear Q Qr SelDif SelLid SelRid


%%%% Ala begins
for i=1:length(pasteCodon) %% note: pasteCodon is column vector   
[NNa,~,Ap]=AlaAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA
if ~(isnan(NNa)||NNa>700) 
Arp=AveEntropy4(NNa);
Ar(i,:)=-Arp/NNa;
A(i,:)=-log(Ap)/NNa;
else 
A(i,:)=NaN;   
Ar(i,:)=NaN;
end 
end

SelDif=A(~isnan(A))-Ar(~isnan(Ar));                       %% SelDif: selection difference; SelWeakId: weak selection dots ID 
SelLid=find(SelDif<0);
Ltest=abs(sum(SelDif(SelLid)));
SelRid=find(SelDif>0);
Rtest=sum(SelDif(SelRid));
Sptest=(Ltest-Rtest)/length(SelDif); 

fileID15=fopen('Ala.txt','a');
fprintf(fileID15,'%s,%d\n',speciesName{tP},Sptest);
fclose(fileID15);
clear A Ar SelDif SelLid SelRid


%%%%%Arg begins
for i=1:length(pasteCodon) %% note: pasteCodon is column vector  
[NNr,~,Rp]=ArgAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA
if ~(isnan(NNr)||NNr>500) 
Rrp=AveEntropy6(NNr);
Rr(i,:)=-Rrp/NNr;
R(i,:)=-log(Rp)/NNr;
else 
R(i,:)=NaN;   
Rr(i,:)=NaN;
end 
end

SelDif=R(~isnan(R))-Rr(~isnan(Rr));                       %% SelDif: selection difference; SelWeakId: weak selection dots ID 
SelLid=find(SelDif<0);
Ltest=abs(sum(SelDif(SelLid)));
SelRid=find(SelDif>0);
Rtest=sum(SelDif(SelRid));
Sptest=(Ltest-Rtest)/length(SelDif);  

fileID16=fopen('Arg.txt','a');
fprintf(fileID16,'%s,%d\n',speciesName{tP},Sptest);
fclose(fileID16);
clear R Rr SelDif SelLid SelRid


%%%%Leu begins
for i=1:length(pasteCodon) %% note: pasteCodon is column vector      
[NNl,~,Lp]=LeuAminoAcidH(pasteCodon{1,i}); 
if ~(isnan(NNl)||NNl>500) 
Lrp=AveEntropy6(NNl);
Lr(i,:)=-Lrp/NNl;
L(i,:)=-log(Lp)/NNl;
else 
L(i,:)=NaN;   
Lr(i,:)=NaN;
end 
end

SelDif=L(~isnan(L))-Lr(~isnan(Lr));                       %% SelDif: selection difference; SelWeakId: weak selection dots ID 
SelLid=find(SelDif<0);
Ltest=abs(sum(SelDif(SelLid)));
SelRid=find(SelDif>0);
Rtest=sum(SelDif(SelRid));
Sptest=(Ltest-Rtest)/length(SelDif);  

fileID17=fopen('Leu.txt','a');
fprintf(fileID17,'%s,%d\n',speciesName{tP},Sptest);
fclose(fileID17);
clear L Lr SelDif SelLid SelRid


%%%%%%Ser begins
for i=1:length(pasteCodon) %% note: pasteCodon is column vector   
[NNs,~,Sp]=SerAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA
if ~(isnan(NNs)||NNs>500) 
Srp=AveEntropy6(NNs);
Sr(i,:)=-Srp/NNs;
S(i,:)=-log(Sp)/NNs;
else 
S(i,:)=NaN;   
Sr(i,:)=NaN;
end 
end

SelDif=S(~isnan(S))-Sr(~isnan(Sr));                       %% SelDif: selection difference; SelWeakId: weak selection dots ID 
SelLid=find(SelDif<0);
Ltest=abs(sum(SelDif(SelLid)));
SelRid=find(SelDif>0);
Rtest=sum(SelDif(SelRid));
Sptest=(Ltest-Rtest)/length(SelDif);  

fileID18=fopen('Ser.txt','a');
fprintf(fileID18,'%s,%d\n',speciesName{tP},Sptest);
fclose(fileID18);
clear S Sr SelDif SelLid SelRid

clear fileName pasteCodon
end
