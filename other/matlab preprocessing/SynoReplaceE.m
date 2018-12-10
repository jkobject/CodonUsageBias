setSynonymousCodonTable

global cfG cfGT ctG cfD cfDT ctD cfE cfET ctE cfV cfVT ctV cfA cfAT ctA cfR cfRT ctR cfS cfST ctS cfK cfKT ctK cfN cfNT ctN cfM cfMT ctM cfI cfIT ctI cfT cfTT ctT cfW cfWT ctW cfZ cfZT ctZ cfC cfCT ctC cfY cfYT  ctY cfL cfLT ctL cfF cfFT ctF cfQ cfQT ctQ cfH cfHT ctH cfP cfPT ctP;

load('AveEntropy2f.mat');
load('AveEntropy3f.mat');
load('AveEntropy4f.mat');
load('AveEntropy6.mat');

FileName={'pasteCodon11.mat','pasteCodon12.mat','pasteCodon13.mat','pasteCodon14.mat','pasteCodon15.mat',...
    'pasteCodon21.mat','pasteCodon22.mat','pasteCodon23.mat',...
    'pasteCodon31.mat','pasteCodon32.mat','pasteCodon33.mat','pasteCodon34.mat','pasteCodon35.mat','pasteCodon36.mat','pasteCodon37.mat',...
    'pasteCodon41.mat','pasteCodon42.mat','pasteCodon43.mat','pasteCodon44.mat','pasteCodon45.mat'};


PasteCodon=cell(1,20);

for r=1:20

m=matfile(FileName{r});

v=who(m);

vn=v{1};
    
PasteCodon{r}=m.(vn);

end

SpeciesName={'sah6','saeu','saku','Ag','Yl','sp','so','sj','ac','afl','ans','anr','ao','at','afu','ff','fg','fo','fv','tr'};

load('AveEntropy2f.mat');
load('AveEntropy3f.mat');
load('AveEntropy4f.mat');
load('AveEntropy6.mat');

% %%calculation for Ile begins

for j=1:20

pasteCodon=PasteCodon{j};

for i=1:length(pasteCodon) %% note: pasteCodon is column vector
        
[NNi,~,Ip]=IleAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA

if ~(isnan(NNi)||NNi>1000) 

Irp=AveEntropy3(NNi);
Ir(i,:)=-Irp/NNi;
I(i,:)=-log(Ip)/NNi;

else 
I(i,:)=NaN;   
Ir(i,:)=NaN;

end 
end

SelDif=I(~isnan(I))-Ir(~isnan(Ir));                       %% SelDif: selection difference; SelWeakId: weak selection dots ID 
SelLid=find(SelDif<0);
Ltest=abs(sum(SelDif(SelLid)));
SelRid=find(SelDif>0);
Rtest=sum(SelDif(SelRid));
Sptest(j)=(Ltest-Rtest)/length(SelDif);  

clear I Ir SelDif SelLid SelRid

end
tblIle= table((SpeciesName)',(Sptest)','VariableNames',{'SpeciesName' 'SelectionPressure'});
writetable(tblIle);


% %%calculation for Asp begins


for j=1:20

pasteCodon=PasteCodon{j};

for i=1:length(pasteCodon) %% note: pasteCodon is column vector
        
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
% fileName=['AspComparison',num2str(j),'.mat'];
% save(fileName,'D','Dr');

SelDif=D(~isnan(D))-Dr(~isnan(Dr));                       %% SelDif: selection difference; SelWeakId: weak selection dots ID 
SelLid=find(SelDif<0);
Ltest=abs(sum(SelDif(SelLid)));
SelRid=find(SelDif>0);
Rtest=sum(SelDif(SelRid));
Sptest(j)=(Ltest-Rtest)/length(SelDif);  
clear D Dr SelDif SelLid SelRid
end
tblAsp= table((SpeciesName)',(Sptest)','VariableNames',{'SpeciesName' 'SelectionPressure'});
writetable(tblAsp);


% %%calculation for Glu begins

for j=1:20

pasteCodon=PasteCodon{j};

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
Sptest(j)=(Ltest-Rtest)/length(SelDif); 
clear E Er SelDif SelLid SelRid
end
tblGlu= table((SpeciesName)',(Sptest)','VariableNames',{'SpeciesName' 'SelectionPressure'});
writetable(tblGlu);

%% His calculation begins

for j=1:20

pasteCodon=PasteCodon{j};

for i=1:length(pasteCodon) %% note: pasteCodon is column vector
    
% disp([num2str(i), 'codon sequence calculation begins']);
    
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
% 
% % % figure    
% % %     
% % % plot(H,Hr,'.');
% % % 
% % %  hold on
% % %  
% % % x=(0:0.05:0.1); 
% % % y=x;
% % % 
% % % plot(x,y);
% % % 
% % % hold off
% % % 
% % % xlim([0 0.7]);
% % % 
% % % % set(gca,'xaxislocation','bottom','yaxislocation','left','xdir','reverse','ydir','reverse')
% % % 
% % % xlabel('per codon entropy of original observed sequence');
% % %  
% % % ylabel('per codon expected entropy of sequence after equal codon replacement');
% % % 
% % % title([SpeciesName{j},' Amino Acid--His(2sysnonymous)']);
% % % 
% % % fileName=['HisComparison',num2str(j),'.mat'];
% % % 
% % % save(fileName,'H','Hr');
% % 
SelDif=H(~isnan(H))-Hr(~isnan(Hr));                       %% SelDif: selection difference; SelWeakId: weak selection dots ID 
SelLid=find(SelDif<0);
Ltest=abs(sum(SelDif(SelLid)));
SelRid=find(SelDif>0);
Rtest=sum(SelDif(SelRid));
Sptest(j)=(Ltest-Rtest)/length(SelDif);  
clear H Hr SelDif SelLid SelRid
end
tblHis = table((SpeciesName)',(Sptest)','VariableNames',{'SpeciesName' 'SelectionPressure'});
writetable(tblHis);

% % Pro calculation begins

for j=1:20

pasteCodon=PasteCodon{j};

for i=1:length(pasteCodon) %% note: pasteCodon is column vector
    
% disp([num2str(i), 'codon sequence calculation begins']);
    
[NNp,~,Pp]=ProAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA

if ~(isnan(NNp)||NNp>700) 

Prp=AveEntropy4(NNp);

Pr(i,:)=-Prp/NNp;

P(i,:)=-log(Pp)/NNp;

else 
P(i,:)=NaN;   
Pr(i,:)=NaN;
end 
end
% % % 
% % % figure    
% % %     
% % % plot(P,Pr,'.');
% % % 
% % % hold on
% % %  
% % % x=(0:0.05:0.25); 
% % % y=x;
% % % 
% % % plot(x,y);
% % % 
% % % hold off
% % % 
% % % % set(gca,'xaxislocation','bottom','yaxislocation','left','xdir','reverse','ydir','reverse')
% % % xlim([0 1.2]);
% % % 
% % % xlabel('per codon entropy of original observed sequence');
% % %  
% % % ylabel('per codon expected entropy of sequence after equal codon replacement');
% % % 
% % % title([SpeciesName{j},' Amino Acid--Pro(4sysnonymous)']);
% % % 
% % % fileName=['ProComparison',num2str(j),'.mat'];
% % % 
% % % save(fileName,'P','Pr');
% % % 
SelDif=P(~isnan(P))-Pr(~isnan(Pr));                       %% SelDif: selection difference; SelWeakId: weak selection dots ID 
SelLid=find(SelDif<0);
Ltest=abs(sum(SelDif(SelLid)));
SelRid=find(SelDif>0);
Rtest=sum(SelDif(SelRid));
Sptest(j)=(Ltest-Rtest)/length(SelDif);  
clear P Pr SelDif SelLid SelRid
end
tblPro = table((SpeciesName)',(Sptest)','VariableNames',{'SpeciesName' 'SelectionPressure'});
writetable(tblPro);

%%calculation for Gly begins

for j=1:20

pasteCodon=PasteCodon{j};

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
Sptest(j)=(Ltest-Rtest)/length(SelDif);  
clear G Gr SelDif SelLid SelRid
end
tblGly = table((SpeciesName)',(Sptest)','VariableNames',{'SpeciesName' 'SelectionPressure'});
writetable(tblGly);

%%calculation for Val begins

for j=1:20

pasteCodon=PasteCodon{j};

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
Sptest(j)=(Ltest-Rtest)/length(SelDif);   
clear V Vr SelDif SelLid SelRid
end
tblVal = table((SpeciesName)',(Sptest)','VariableNames',{'SpeciesName' 'SelectionPressure'});
writetable(tblVal);

%%calculation for Lys begins
for j=1:20

pasteCodon=PasteCodon{j};

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
Sptest(j)=(Ltest-Rtest)/length(SelDif); 
clear K Kr SelDif SelLid SelRid
end
tblLys= table((SpeciesName)',(Sptest)','VariableNames',{'SpeciesName' 'SelectionPressure'});
writetable(tblLys);

%%calculation for Asn begins

for j=1:20

pasteCodon=PasteCodon{j};

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
Sptest(j)=(Ltest-Rtest)/length(SelDif);  
clear N Nr SelDif SelLid SelRid
end
tblAsn= table((SpeciesName)',(Sptest)','VariableNames',{'SpeciesName' 'SelectionPressure'});
writetable(tblAsn);

%%calculation for Thr begins

for j=1:20

pasteCodon=PasteCodon{j};

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
Sptest(j)=(Ltest-Rtest)/length(SelDif); 
clear T Tr SelDif SelLid SelRid
end
tblThr= table((SpeciesName)',(Sptest)','VariableNames',{'SpeciesName' 'SelectionPressure'});
writetable(tblThr);

%%calculation for Cys begins

for j=1:20

pasteCodon=PasteCodon{j};

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
Sptest(j)=(Ltest-Rtest)/length(SelDif);  
clear C Cr SelDif SelLid SelRid
end
tblCys= table((SpeciesName)',(Sptest)','VariableNames',{'SpeciesName' 'SelectionPressure'});
writetable(tblCys);


%%calculation for Tyr begins
for j=1:20

pasteCodon=PasteCodon{j};

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
Sptest(j)=(Ltest-Rtest)/length(SelDif);  
clear Y Yr SelDif SelLid SelRid
end
tblTyr= table((SpeciesName)',(Sptest)','VariableNames',{'SpeciesName' 'SelectionPressure'});
writetable(tblTyr);

%%calculation for Phe begins
for j=1:20

pasteCodon=PasteCodon{j};

for i=1:length(pasteCodon) %% note: pasteCodon is column vector
        
[NNf,~,Fp]=PheAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA

if ~(isnan(NNf)||NNf>1500) 

Frp=AveEntropy2(NNf);
Fr(i,:)=-Frp/NNf;
F(i,:)=-log(Fp)/NNf;

else 
F(i,:)=NaN;   %%calculate again  
Fr(i,:)=NaN;

end 
end

SelDif=F(~isnan(F))-Fr(~isnan(Fr));                       %% SelDif: selection difference; SelWeakId: weak selection dots ID 
SelLid=find(SelDif<0);
Ltest=abs(sum(SelDif(SelLid)));
SelRid=find(SelDif>0);
Rtest=sum(SelDif(SelRid));
Sptest(j)=(Ltest-Rtest)/length(SelDif);  
clear F Fr SelDif SelLid SelRid
end
tblPhe= table((SpeciesName)',(Sptest)','VariableNames',{'SpeciesName' 'SelectionPressure'});
writetable(tblPhe);


%%calculation for Gln begins
for j=1:20

pasteCodon=PasteCodon{j};

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
Sptest(j)=(Ltest-Rtest)/length(SelDif);  
clear Q Qr SelDif SelLid SelRid
end
tblGln= table((SpeciesName)',(Sptest)','VariableNames',{'SpeciesName' 'SelectionPressure'});
writetable(tblGln);

%%calculation for Ala begins

for j=1:20

pasteCodon=PasteCodon{j};

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
Sptest(j)=(Ltest-Rtest)/length(SelDif);  
clear A Ar SelDif SelLid SelRid
end
tblAla= table((SpeciesName)',(Sptest)','VariableNames',{'SpeciesName' 'SelectionPressure'});
writetable(tblAla);

%% calculation for Arg begins

for j=1:20

pasteCodon=PasteCodon{j};

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
Sptest(j)=(Ltest-Rtest)/length(SelDif);  
clear R Rr SelDif SelLid SelRid
end

tblArg = table((SpeciesName)',(Sptest)','VariableNames',{'SpeciesName' 'SelectionPressure'});
writetable(tblArg);

%%calculation for Leu begins
for j=1:20

pasteCodon=PasteCodon{j};

for i=1:length(pasteCodon) %% note: pasteCodon is column vector
        
[NNl,~,Lp]=LeuAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA

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
Sptest(j)=(Ltest-Rtest)/length(SelDif);  
clear L Lr SelDif SelLid SelRid
end

tblLeu = table((SpeciesName)',(Sptest)','VariableNames',{'SpeciesName' 'SelectionPressure'});
writetable(tblLeu);


%%calculation for Ser begins
for j=1:20

pasteCodon=PasteCodon{j};

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
Sptest(j)=(Ltest-Rtest)/length(SelDif);  
clear S Sr SelDif SelLid SelRid
end

tblSer = table((SpeciesName)',(Sptest)','VariableNames',{'SpeciesName' 'SelectionPressure'});
writetable(tblSer);