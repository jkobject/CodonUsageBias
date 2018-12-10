setSynonymousCodonTable
%%randomly sequence entropy cost comparison

load('AveEntropy2f.mat');
load('AveEntropy6.mat');

pasteCodon=PasteCodon{1};

%%calculation for Asp begins  
clear SelDif SelLid Ltest SelRid Rtest

for i=1:length(pasteCodon)
    
[NNd,~,~]=AspAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA

if ~(isnan(NNd)||NNd>1500)

InRd=replaceE(NNd,ctD);%%input random sequence for D
[~,~,Dp]=AspAminoAcidH(InRd);

InRTd=replaceT(NNd,ctD,cfDT);%%input random sequence for D
[~,~,DpT]=AspAminoAcidH(InRTd);
    
Drp=AveEntropy2(NNd);
Dr(i,:)=-Drp/NNd;
Dt(i,:)=-log(DpT)/NNd;
D(i,:)=-log(Dp)/NNd;

else 
D(i,:)=NaN; 
Dt(i,:)=NaN;
Dr(i,:)=NaN;

end
end 
 
SelDif=D(~isnan(D))-Dr(~isnan(Dr));                       %% SelDif: selection difference; SelWeakId: weak selection dots ID 
SelLid=find(SelDif<0);
Ltest=abs(sum(SelDif(SelLid)));
SelRid=find(SelDif>0);
Rtest=sum(SelDif(SelRid));
Sp1=(Ltest-Rtest)/length(SelDif);   

SelDif=Dt(~isnan(Dt))-Dr(~isnan(Dr));                       %% SelDif: selection difference; SelWeakId: weak selection dots ID 
SelLid=find(SelDif<0);
Ltest=abs(sum(SelDif(SelLid)));
SelRid=find(SelDif>0);
Rtest=sum(SelDif(SelRid));
Sp1t=(Ltest-Rtest)/length(SelDif);  

clear SelDif SelLid Ltest SelRid Rtest

%%calculation for Glu begines
for i=1:length(pasteCodon) %% note: pasteCodon is column vector
    
[NNe,~,~]=GluAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA

if ~(isnan(NNe)||NNe>1500) 

InRe=replaceE(NNe,ctE);
[~,~,Ep]=GluAminoAcidH(InRe);

InRTe=replaceT(NNe,ctE,cfET);
[~,~,EpT]=GluAminoAcidH(InRTe);

Erp=AveEntropy2(NNe);
Er(i,:)=-Erp/NNe;
E(i,:)=-log(Ep)/NNe;
Et(i,:)=-log(EpT)/NNe;

else 
E(i,:)=NaN;   
Er(i,:)=NaN;
Et(i,:)=NaN;

end 
end

SelDif=E(~isnan(E))-Er(~isnan(Er));                       %% SelDif: selection difference; SelWeakId: weak selection dots ID 
SelLid=find(SelDif<0);
Ltest=abs(sum(SelDif(SelLid)));
SelRid=find(SelDif>0);
Rtest=sum(SelDif(SelRid));
Sp2=(Ltest-Rtest)/length(SelDif);   

SelDif=Et(~isnan(Et))-Er(~isnan(Er));                       %% SelDif: selection difference; SelWeakId: weak selection dots ID 
SelLid=find(SelDif<0);
Ltest=abs(sum(SelDif(SelLid)));
SelRid=find(SelDif>0);
Rtest=sum(SelDif(SelRid));
Sp2t=(Ltest-Rtest)/length(SelDif);   

clear SelDif SelLid Ltest SelRid Rtest

%calculation for Ile begin

for i=1:length(pasteCodon) %% note: pasteCodon is column vector
        
[NNi,~,~]=IleAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA

if ~(isnan(NNi)||NNi>1000) 
    
InRi=replaceE(NNi,ctI);%%input random sequence for D
[~,~,Ip]=IleAminoAcidH(InRi);

InRTi=replaceT(NNi,ctI,cfIT);%%input random sequence for D
[~,~,IpT]=IleAminoAcidH(InRTi);

Irp=AveEntropy3(NNi);
Ir(i,:)=-Irp/NNi;
I(i,:)=-log(Ip)/NNi;
It(i,:)=-log(IpT)/NNi;

else
I(i,:)=NaN;   
Ir(i,:)=NaN;
It(i,:)=NaN;
end

end 


SelDif=I(~isnan(I))-Ir(~isnan(Ir));                   
SelLid=find(SelDif<0);
Ltest=abs(sum(SelDif(SelLid)));
SelRid=find(SelDif>0);
Rtest=sum(SelDif(SelRid));
Sp3=(Ltest-Rtest)/length(SelDif); 

SelDif=It(~isnan(It))-Ir(~isnan(Ir));                       
SelLid=find(SelDif<0);
Ltest=abs(sum(SelDif(SelLid)));
SelRid=find(SelDif>0);
Rtest=sum(SelDif(SelRid));
Sp3t=(Ltest-Rtest)/length(SelDif); 

clear SelDif SelLid Ltest SelRid Rtest
%calculation for Pro begins

for i=1:length(pasteCodon) %% note: pasteCodon is column vector
        
[NNp,~,~]=ProAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA

if ~(isnan(NNp)||NNp>700) 
    
InRp=replaceE(NNp,ctP);%%input random sequence for D
[~,~,Pp]=ProAminoAcidH(InRp);

InRTp=replaceT(NNp,ctP,cfPT);%%input random sequence for D
[~,~,PpT]=ProAminoAcidH(InRTp);


Prp=AveEntropy4(NNp);
Pr(i,:)=-Prp/NNi;
P(i,:)=-log(Pp)/NNp;
Pt(i,:)=-log(PpT)/NNp;

else 
P(i,:)=NaN;   
Pr(i,:)=NaN;
Pt(i,:)=NaN;

end 
end

SelDif=P(~isnan(P))-Pr(~isnan(Pr));                       %% SelDif: selection difference; SelWeakId: weak selection dots ID 
SelLid=find(SelDif<0);
Ltest=abs(sum(SelDif(SelLid)));
SelRid=find(SelDif>0);
Rtest=sum(SelDif(SelRid));
Sp4=(Ltest-Rtest)/length(SelDif); 

SelDif=Pt(~isnan(Pt))-Pr(~isnan(Pr));                       %% SelDif: selection difference; SelWeakId: weak selection dots ID 
SelLid=find(SelDif<0);
Ltest=abs(sum(SelDif(SelLid)));
SelRid=find(SelDif>0);
Rtest=sum(SelDif(SelRid));
Sp4t=(Ltest-Rtest)/length(SelDif); 
clear SelDif SelLid Ltest SelRid Rtest

%calculation for Ser begins
for i=1:length(pasteCodon) %% note: pasteCodon is column vector
        
[NNs,~,~]=SerAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA

if ~(isnan(NNs)||NNs>500) 
    
InRs=replaceE(NNs,ctS);%%input random sequence for D
[~,~,Sp]=SerAminoAcidH(InRs);

InRTs=replaceT(NNs,ctS,cfST);%%input random sequence for D
[~,~,SpT]=SerAminoAcidH(InRTs);


Srp=AveEntropy6(NNs);
Sr(i,:)=-Srp/NNs;
S(i,:)=-log(Sp)/NNs;
St(i,:)=-log(SpT)/NNs;

else 
S(i,:)=NaN;   
Sr(i,:)=NaN;
St(i,:)=NaN;

end 
end

SelDif=S(~isnan(S))-Sr(~isnan(Sr));                       %% SelDif: selection difference; SelWeakId: weak selection dots ID 
SelLid=find(SelDif<0);
Ltest=abs(sum(SelDif(SelLid)));
SelRid=find(SelDif>0);
Rtest=sum(SelDif(SelRid));
Sp5=(Ltest-Rtest)/length(SelDif); 

SelDif=St(~isnan(St))-Sr(~isnan(Sr));                       %% SelDif: selection difference; SelWeakId: weak selection dots ID 
SelLid=find(SelDif<0);
Ltest=abs(sum(SelDif(SelLid)));
SelRid=find(SelDif>0);
Rtest=sum(SelDif(SelRid));
Sp5t=(Ltest-Rtest)/length(SelDif); 

clear SelDif SelLid Ltest SelRid Rtest

%calculation for His begins

for i=1:length(pasteCodon) %% note: pasteCodon is column vector
        
[NNh,~,~]=HisAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA

if ~(isnan(NNh)||NNh>1500) 
    
InRh=replaceE(NNh,ctH);%%input random sequence for D
[~,~,Hp]=HisAminoAcidH(InRh);

InRTh=replaceT(NNh,ctH,cfHT);%%input random sequence for D
[~,~,HpT]=HisAminoAcidH(InRTh);

Hrp=AveEntropy2(NNh);
Hr(i,:)=-Hrp/NNh;
H(i,:)=-log(Hp)/NNh;
Ht(i,:)=-log(HpT)/NNh;
else 
    
H(i,:)=NaN;   
Hr(i,:)=NaN;
Ht(i,:)=NaN;

end 
end

SelDif=H(~isnan(H))-Hr(~isnan(Hr));                       %% SelDif: selection difference; SelWeakId: weak selection dots ID 
SelLid=find(SelDif<0);
Ltest=abs(sum(SelDif(SelLid)));
SelRid=find(SelDif>0);
Rtest=sum(SelDif(SelRid));
Sp6=(Ltest-Rtest)/length(SelDif); 

SelDif=Ht(~isnan(Ht))-Hr(~isnan(Hr));                       %% SelDif: selection difference; SelWeakId: weak selection dots ID 
SelLid=find(SelDif<0);
Ltest=abs(sum(SelDif(SelLid)));
SelRid=find(SelDif>0);
Rtest=sum(SelDif(SelRid));
Sp6t=(Ltest-Rtest)/length(SelDif); 
clear SelDif SelLid Ltest SelRid Rtest


%calculation for Leu begins

for i=1:length(pasteCodon) %% note: pasteCodon is column vector
        
[NNl,~,~]=LeuAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA

if ~(isnan(NNl)||NNl>500) 
    
InRl=replaceE(NNl,ctL);%%input random sequence for D
[~,~,Lp]=LeuAminoAcidH(InRl);

InRTl=replaceT(NNl,ctL,cfLT);%%input random sequence for D
[~,~,LpT]=LeuAminoAcidH(InRTl);

Lrp=AveEntropy6(NNl);
Lr(i,:)=-Lrp/NNl;
L(i,:)=-log(Lp)/NNl;
Lt(i,:)=-log(LpT)/NNl;
else 
    
L(i,:)=NaN;   
Lr(i,:)=NaN;
Lt(i,:)=NaN;

end 
end

SelDif=L(~isnan(L))-Lr(~isnan(Lr));                       %% SelDif: selection difference; SelWeakId: weak selection dots ID 
SelLid=find(SelDif<0);
Ltest=abs(sum(SelDif(SelLid)));
SelRid=find(SelDif>0);
Rtest=sum(SelDif(SelRid));
Sp7=(Ltest-Rtest)/length(SelDif);

SelDif=Lt(~isnan(Lt))-Lr(~isnan(Lr));                       %% SelDif: selection difference; SelWeakId: weak selection dots ID 
SelLid=find(SelDif<0);
Ltest=abs(sum(SelDif(SelLid)));
SelRid=find(SelDif>0);
Rtest=sum(SelDif(SelRid));
Sp7t=(Ltest-Rtest)/length(SelDif);

calculation for Gly begins
clear SelDif SelLid Ltest SelRid Rtest

for i=1:length(pasteCodon) %% note: pasteCodon is column vector
        
[NNg,~,~]=GlyAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA

if ~(isnan(NNg)||NNg>700) 
    
InRg=replaceE(NNg,ctG);%%input random sequence for D
[~,~,Gp]=GlyAminoAcidH(InRg);

InRTg=replaceT(NNg,ctG,cfGT);%%input random sequence for D
[~,~,GpT]=GlyAminoAcidH(InRTg);

Grp=AveEntropy4(NNg);
Gr(i,:)=-Grp/NNg;
G(i,:)=-log(Gp)/NNg;
Gt(i,:)=-log(GpT)/NNg;
else 
    
G(i,:)=NaN;   
Gr(i,:)=NaN;
Gt(i,:)=NaN;

end 
end

SelDif=G(~isnan(G))-Gr(~isnan(Gr));                       %% SelDif: selection difference; SelWeakId: weak selection dots ID 
SelLid=find(SelDif<0);
Ltest=abs(sum(SelDif(SelLid)));
SelRid=find(SelDif>0);
Rtest=sum(SelDif(SelRid));
Sp8=(Ltest-Rtest)/length(SelDif);

SelDif=Gt(~isnan(Gt))-Gr(~isnan(Gr));                       %% SelDif: selection difference; SelWeakId: weak selection dots ID 
SelLid=find(SelDif<0);
Ltest=abs(sum(SelDif(SelLid)));
SelRid=find(SelDif>0);
Rtest=sum(SelDif(SelRid));
Sp8t=(Ltest-Rtest)/length(SelDif);

calculation for Val begins

for i=1:length(pasteCodon) %% note: pasteCodon is column vector
        
[NNv,~,~]=ValAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA

if ~(isnan(NNv)||NNv>700) 
    
InRv=replaceE(NNv,ctV);%%input random sequence for D
[~,~,Vp]=ValAminoAcidH(InRv);

InRTv=replaceT(NNv,ctV,cfVT);%%input random sequence for D
[~,~,VpT]=ValAminoAcidH(InRTv);

Vrp=AveEntropy4(NNv);
Vr(i,:)=-Vrp/NNv;
V(i,:)=-log(Vp)/NNv;
Vt(i,:)=-log(VpT)/NNv;
else 
    
V(i,:)=NaN;   
Vr(i,:)=NaN;
Vt(i,:)=NaN;

end 
end

SelDif=V(~isnan(V))-Vr(~isnan(Vr));                       %% SelDif: selection difference; SelWeakId: weak selection dots ID 
SelLid=find(SelDif<0);
Ltest=abs(sum(SelDif(SelLid)));
SelRid=find(SelDif>0);
Rtest=sum(SelDif(SelRid));
Sp9=(Ltest-Rtest)/length(SelDif);

SelDif=Vt(~isnan(Vt))-Vr(~isnan(Vr));                       %% SelDif: selection difference; SelWeakId: weak selection dots ID 
SelLid=find(SelDif<0);
Ltest=abs(sum(SelDif(SelLid)));
SelRid=find(SelDif>0);
Rtest=sum(SelDif(SelRid));
Sp9t=(Ltest-Rtest)/length(SelDif);

%%calculation for Ala begins
clear SelDif SelLid Ltest SelRid Rtest

for i=1:length(pasteCodon) %% note: pasteCodon is column vector
        
[NNa,~,~]=AlaAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA

if ~(isnan(NNa)||NNa>700) 
    
InRa=replaceE(NNa,ctA);%%input random sequence for D
[~,~,Ap]=AlaAminoAcidH(InRa);

InRTa=replaceT(NNa,ctA,cfAT);%%input random sequence for D
[~,~,ApT]=AlaAminoAcidH(InRTa);

Arp=AveEntropy6(NNa);
Ar(i,:)=-Arp/NNa;
A(i,:)=-log(Ap)/NNa;
At(i,:)=-log(ApT)/NNa;
else 
    
A(i,:)=NaN;   
Ar(i,:)=NaN;
At(i,:)=NaN;

end 
end

SelDif=A(~isnan(A))-Ar(~isnan(Ar));                       %% SelDif: selection difference; SelWeakId: weak selection dots ID 
SelLid=find(SelDif<0);
Ltest=abs(sum(SelDif(SelLid)));
SelRid=find(SelDif>0);
Rtest=sum(SelDif(SelRid));
Sp10=(Ltest-Rtest)/length(SelDif);

SelDif=At(~isnan(At))-Ar(~isnan(Ar));                       %% SelDif: selection difference; SelWeakId: weak selection dots ID 
SelLid=find(SelDif<0);
Ltest=abs(sum(SelDif(SelLid)));
SelRid=find(SelDif>0);
Rtest=sum(SelDif(SelRid));
Sp10t=(Ltest-Rtest)/length(SelDif);

%%calculation for Arg begins
clear SelDif SelLid Ltest SelRid Rtest

for i=1:length(pasteCodon) %% note: pasteCodon is column vector
        
[NNr,~,~]=ArgAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA

if ~(isnan(NNr)||NNr>500) 
    
InRr=replaceE(NNr,ctR);%%input random sequence for D
[~,~,Rp]=ArgAminoAcidH(InRr);

InRTr=replaceT(NNr,ctR,cfRT);%%input random sequence for D
[~,~,RpT]=ArgAminoAcidH(InRTr);

Rrp=AveEntropy6(NNr);
Rr(i,:)=-Rrp/NNr;
R(i,:)=-log(Rp)/NNr;
Rt(i,:)=-log(RpT)/NNr;
else 
    
R(i,:)=NaN;   
Rr(i,:)=NaN;
Rt(i,:)=NaN;

end 
end

SelDif=R(~isnan(R))-Rr(~isnan(Rr));                       %% SelDif: selection difference; SelWeakId: weak selection dots ID 
SelLid=find(SelDif<0);
Ltest=abs(sum(SelDif(SelLid)));
SelRid=find(SelDif>0);
Rtest=sum(SelDif(SelRid));
Sp11=(Ltest-Rtest)/length(SelDif);

SelDif=Rt(~isnan(Rt))-Rr(~isnan(Rr));                       %% SelDif: selection difference; SelWeakId: weak selection dots ID 
SelLid=find(SelDif<0);
Ltest=abs(sum(SelDif(SelLid)));
SelRid=find(SelDif>0);
Rtest=sum(SelDif(SelRid));
Sp11t=(Ltest-Rtest)/length(SelDif);


calculation for Lys begins
clear SelDif SelLid Ltest SelRid Rtest

for i=1:length(pasteCodon) %% note: pasteCodon is column vector
        
[NNk,~,~]=LysAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA

if ~(isnan(NNk)||NNk>1500) 
    
InRk=replaceE(NNk,ctK);%%input random sequence for D
[~,~,Kp]=LysAminoAcidH(InRk);

InRTk=replaceT(NNk,ctK,cfKT);%%input random sequence for D
[~,~,KpT]=LysAminoAcidH(InRTk);

Krp=AveEntropy2(NNk);
Kr(i,:)=-Krp/NNk;
K(i,:)=-log(Kp)/NNk;
Kt(i,:)=-log(KpT)/NNk;
else 
    
K(i,:)=NaN;   
Kr(i,:)=NaN;
Kt(i,:)=NaN;

end 
end

SelDif=K(~isnan(K))-Kr(~isnan(Kr));                       %% SelDif: selection difference; SelWeakId: weak selection dots ID 
SelLid=find(SelDif<0);
Ltest=abs(sum(SelDif(SelLid)));
SelRid=find(SelDif>0);
Rtest=sum(SelDif(SelRid));
Sp12=(Ltest-Rtest)/length(SelDif);

SelDif=Kt(~isnan(Kt))-Kr(~isnan(Kr));                       %% SelDif: selection difference; SelWeakId: weak selection dots ID 
SelLid=find(SelDif<0);
Ltest=abs(sum(SelDif(SelLid)));
SelRid=find(SelDif>0);
Rtest=sum(SelDif(SelRid));
Sp12t=(Ltest-Rtest)/length(SelDif);

calculation for Asn begins
clear SelDif SelLid Ltest SelRid Rtest

for i=1:length(pasteCodon) %% note: pasteCodon is column vector
        
[NNn,~,~]=AsnAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA

if ~(isnan(NNn)||NNn>1500) 
    
InRn=replaceE(NNn,ctN);%%input random sequence for D
[~,~,Np]=AsnAminoAcidH(InRn);

InRTn=replaceT(NNn,ctN,cfNT);%%input random sequence for D
[~,~,NpT]=AsnAminoAcidH(InRTn);

Nrp=AveEntropy2(NNn);
Nr(i,:)=-Nrp/NNn;
N(i,:)=-log(Np)/NNn;
Nt(i,:)=-log(NpT)/NNn;
else 
    
N(i,:)=NaN;   
Nr(i,:)=NaN;
Nt(i,:)=NaN;

end 
end

SelDif=N(~isnan(N))-Nr(~isnan(Nr));                       %% SelDif: selection difference; SelWeakId: weak selection dots ID 
SelLid=find(SelDif<0);
Ltest=abs(sum(SelDif(SelLid)));
SelRid=find(SelDif>0);
Rtest=sum(SelDif(SelRid));
Sp13=(Ltest-Rtest)/length(SelDif);

SelDif=Nt(~isnan(Nt))-Nr(~isnan(Nr));                       %% SelDif: selection difference; SelWeakId: weak selection dots ID 
SelLid=find(SelDif<0);
Ltest=abs(sum(SelDif(SelLid)));
SelRid=find(SelDif>0);
Rtest=sum(SelDif(SelRid));
Sp13t=(Ltest-Rtest)/length(SelDif);

calculation for Cys begins
clear SelDif SelLid Ltest SelRid Rtest

for i=1:length(pasteCodon) %% note: pasteCodon is column vector
        
[NNc,~,~]=CysAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA

if ~(isnan(NNc)||NNc>1500) 
    
InRc=replaceE(NNc,ctC);%%input random sequence for D
[~,~,Cp]=CysAminoAcidH(InRc);

InRTc=replaceT(NNc,ctC,cfCT);%%input random sequence for D
[~,~,CpT]=CysAminoAcidH(InRTc);

Crp=AveEntropy2(NNc);
Cr(i,:)=-Crp/NNc;
C(i,:)=-log(Cp)/NNc;
Ct(i,:)=-log(CpT)/NNc;
else 
    
C(i,:)=NaN;   
Cr(i,:)=NaN;
Ct(i,:)=NaN;

end 
end

SelDif=C(~isnan(C))-Cr(~isnan(Cr));                       %% SelDif: selection difference; SelWeakId: weak selection dots ID 
SelLid=find(SelDif<0);
Ltest=abs(sum(SelDif(SelLid)));
SelRid=find(SelDif>0);
Rtest=sum(SelDif(SelRid));
Sp14=(Ltest-Rtest)/length(SelDif);

SelDif=Ct(~isnan(Ct))-Cr(~isnan(Cr));                       %% SelDif: selection difference; SelWeakId: weak selection dots ID 
SelLid=find(SelDif<0);
Ltest=abs(sum(SelDif(SelLid)));
SelRid=find(SelDif>0);
Rtest=sum(SelDif(SelRid));
Sp14t=(Ltest-Rtest)/length(SelDif);

calculation for Thr begins
clear SelDif SelLid Ltest SelRid Rtest

for i=1:length(pasteCodon) %% note: pasteCodon is column vector
        
[NNt,~,~]=ThrAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA

if ~(isnan(NNt)||NNt>700) 
    
InRt=replaceE(NNt,ctT);%%input random sequence for D
[~,~,Tp]=ThrAminoAcidH(InRt);

InRTt=replaceT(NNt,ctT,cfTT);%%input random sequence for D
[~,~,TpT]=ThrAminoAcidH(InRTt);

Trp=AveEntropy4(NNt);
Tr(i,:)=-Trp/NNt;
T(i,:)=-log(Tp)/NNt;
Tt(i,:)=-log(TpT)/NNt;
else 
    
T(i,:)=NaN;   
Tr(i,:)=NaN;
Tt(i,:)=NaN;

end 
end

SelDif=T(~isnan(T))-Tr(~isnan(Tr));                       %% SelDif: selection difference; SelWeakId: weak selection dots ID 
SelLid=find(SelDif<0);
Ltest=abs(sum(SelDif(SelLid)));
SelRid=find(SelDif>0);
Rtest=sum(SelDif(SelRid));
Sp15=(Ltest-Rtest)/length(SelDif);

SelDif=Tt(~isnan(Tt))-Tr(~isnan(Tr));                       %% SelDif: selection difference; SelWeakId: weak selection dots ID 
SelLid=find(SelDif<0);
Ltest=abs(sum(SelDif(SelLid)));
SelRid=find(SelDif>0);
Rtest=sum(SelDif(SelRid));
Sp15t=(Ltest-Rtest)/length(SelDif);

calculation for Tyr begins
clear SelDif SelLid Ltest SelRid Rtest

for i=1:length(pasteCodon) %% note: pasteCodon is column vector
        
[NNy,~,~]=TyrAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA

if ~(isnan(NNy)||NNy>1500) 
    
InRy=replaceE(NNy,ctY);%%input random sequence for D
[~,~,Yp]=TyrAminoAcidH(InRy);

InRTy=replaceT(NNy,ctY,cfYT);%%input random sequence for D
[~,~,YpT]=TyrAminoAcidH(InRTy);

Yrp=AveEntropy2(NNy);
Yr(i,:)=-Yrp/NNy;
Y(i,:)=-log(Yp)/NNy;
Yt(i,:)=-log(YpT)/NNy;
else 
    
Y(i,:)=NaN;   
Yr(i,:)=NaN;
Yt(i,:)=NaN;

end 
end

SelDif=Y(~isnan(Y))-Yr(~isnan(Yr));                       %% SelDif: selection difference; SelWeakId: weak selection dots ID 
SelLid=find(SelDif<0);
Ltest=abs(sum(SelDif(SelLid)));
SelRid=find(SelDif>0);
Rtest=sum(SelDif(SelRid));
Sp16=(Ltest-Rtest)/length(SelDif);

SelDif=Yt(~isnan(Yt))-Yr(~isnan(Yr));                       %% SelDif: selection difference; SelWeakId: weak selection dots ID 
SelLid=find(SelDif<0);
Ltest=abs(sum(SelDif(SelLid)));
SelRid=find(SelDif>0);
Rtest=sum(SelDif(SelRid));
Sp16t=(Ltest-Rtest)/length(SelDif);

calculation for Phe begins
clear SelDif SelLid Ltest SelRid Rtest

for i=1:length(pasteCodon) %% note: pasteCodon is column vector
        
[NNf,~,~]=PheAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA

if ~(isnan(NNf)||NNf>1500) 
    
InRf=replaceE(NNf,ctF);%%input random sequence for D
[~,~,Fp]=PheAminoAcidH(InRf);

InRTf=replaceT(NNf,ctF,cfFT);%%input random sequence for D
[~,~,FpT]=PheAminoAcidH(InRTf);

Frp=AveEntropy2(NNf);
Fr(i,:)=-Frp/NNf;
F(i,:)=-log(Fp)/NNf;
Ft(i,:)=-log(FpT)/NNf;
else 
    
F(i,:)=NaN;   
Fr(i,:)=NaN;
Ft(i,:)=NaN;

end 
end

SelDif=F(~isnan(F))-Fr(~isnan(Fr));                       %% SelDif: selection difference; SelWeakId: weak selection dots ID 
SelLid=find(SelDif<0);
Ltest=abs(sum(SelDif(SelLid)));
SelRid=find(SelDif>0);
Rtest=sum(SelDif(SelRid));
Sp17=(Ltest-Rtest)/length(SelDif);

SelDif=Ft(~isnan(Ft))-Fr(~isnan(Fr));                       %% SelDif: selection difference; SelWeakId: weak selection dots ID 
SelLid=find(SelDif<0);
Ltest=abs(sum(SelDif(SelLid)));
SelRid=find(SelDif>0);
Rtest=sum(SelDif(SelRid));
Sp17t=(Ltest-Rtest)/length(SelDif);

calculation for Gln begins
clear SelDif SelLid Ltest SelRid Rtest

for i=1:length(pasteCodon) %% note: pasteCodon is column vector
        
[NNq,~,~]=GlnAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA

if ~(isnan(NNq)||NNq>1500) 
    
InRq=replaceE(NNq,ctQ);%%input random sequence for D
[~,~,Qp]=GlnAminoAcidH(InRq);

InRTq=replaceT(NNq,ctQ,cfQT);%%input random sequence for D
[~,~,QpT]=GlnAminoAcidH(InRTq);

Qrp=AveEntropy2(NNq);
Qr(i,:)=-Qrp/NNq;
Q(i,:)=-log(Qp)/NNq;
Qt(i,:)=-log(QpT)/NNq;
else 
    
Q(i,:)=NaN;   
Qr(i,:)=NaN;
Qt(i,:)=NaN;

end 
end

SelDif=Q(~isnan(Q))-Qr(~isnan(Qr));                       %% SelDif: selection difference; SelWeakId: weak selection dots ID 
SelLid=find(SelDif<0);
Ltest=abs(sum(SelDif(SelLid)));
SelRid=find(SelDif>0);
Rtest=sum(SelDif(SelRid));
Sp18=(Ltest-Rtest)/length(SelDif);

SelDif=Qt(~isnan(Qt))-Qr(~isnan(Qr));                       %% SelDif: selection difference; SelWeakId: weak selection dots ID 
SelLid=find(SelDif<0);
Ltest=abs(sum(SelDif(SelLid)));
SelRid=find(SelDif>0);
Rtest=sum(SelDif(SelRid));
Sp18t=(Ltest-Rtest)/length(SelDif);


X=[Sp1,Sp2,Sp3,Sp4,Sp5,Sp6,Sp7,Sp8,Sp9,Sp10,Sp11,Sp12,Sp13,Sp14,Sp15,Sp16,Sp17,Sp18,Sp19;...
Sp1t,Sp2t,Sp3t,Sp4t,Sp5t,Sp6t,Sp7t,Sp8t,Sp9t,Sp10t,Sp11t,Sp12t,Sp13t,Sp14t,Sp15t,Sp16t,Sp17t,Sp18t,Sp19t];

dlmwrite('RandomFsah6.txt',X);

