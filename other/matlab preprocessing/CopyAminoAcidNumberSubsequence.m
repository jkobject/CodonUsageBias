function Y = AminoAcidNumberSubsequence(codonSequence)

disp('A new calculation begins');

codonFrequency=zeros(1,length(codonSequence));

ctG={'GGG','GGA','GGT','GGC'}; %%Gly
%cfG=[0.12,0.23,0.45,0.20];
cfG=[1/4,1/4,1/4,1/4];

ctD={'GAT','GAC'};  %%Asp
%cfD=[0.65,0.35];
cfD=[1/2,1/2];

ctE={'GAG','GAA'};  %%Glu
%cfE=[0.30,0.70];
cfE=[1/2,1/2];

ctV={'GTG','GTA','GTT','GTC'}; %%Val
%cfV=[0.20,0.21,0.39,0.20];%%change 0.22
cfV=[1/4,1/4,1/4,1/4];

ctA={'GCG','GCA','GCT','GCC'}; %%Ala
%cfA=[0.11,0.30,0.37,0.22];
cfA=[1/4,1/4,1/4,1/4];

ctR={'AGG','AGA','CGG','CGA','CGT','CGC'};   %%Arg
%cfR=[0.22,0.47,0.04,0.07,0.14,0.06];%%change 0.21
cfR=[1/6,1/6,1/6,1/6,1/6,1/6];

ctS={'AGT','AGC','TCG','TCA','TCT','TCC'};   %%Ser
%cfS=[0.16,0.11,0.10,0.21,0.26,0.16];
cfS=[1/6,1/6,1/6,1/6,1/6,1/6];

ctK={'AAG','AAA'};   %%Lys
%cfK=[0.42,0.58];
cfK=[1/2,1/2];

ctN={'AAT','AAC'};   %%Asn
%cfN=[0.60,0.40];
cfN=[1/2,1/2];

ctM={'ATG'};  %%Met
%cfM=[1.00];
cfM=[1.00];

ctI={'ATA','ATT','ATC'}; %%Ile
%cfI=[0.28,0.46,0.26];
cfI=[1/3,1/3,1/3];

ctT={'ACG','ACA','ACT','ACC'};  %%Thr
%cfT=[0.14,0.31,0.34,0.21];
cfT=[1/4,1/4,1/4,1/4];

ctW={'TGG'};   %%Trp
%cfW=[1.00];   
cfW=[1.00];

ctZ={'TGA','TAG','TAA'};  %%End
%cfZ=[0.30,0.23,0.47];
cfZ=[1/3,1/3,1/3];

ctC={'TGT','TGC'};  %%Cys
%cfC=[0.62,0.38];
cfC=[1/2,1/2];

ctY={'TAT','TAC'};   %%Tyr
%cfY=[0.57,0.43];
cfY=[1/2,1/2];

ctL={'TTG','TTA','CTG','CTA','CTT','CTC'};  %%Leu
%cfL=[0.28,0.28,0.11,0.14,0.13,0.06];
cfL=[1/6,1/6,1/6,1/6,1/6,1/6];

ctF={'TTT','TTC'}; %%Phe
%cfF=[0.59,0.41];
cfF=[1/2,1/2];

ctQ={'CAG','CAA'};  %%Gln
%cfQ=[0.32,0.68];
cfQ=[1/2,1/2];

ctH={'CAT','CAC'};  %%His
%cfH=[0.64,0.36];
cfH=[1/2,1/2];

ctP={'CCG','CCA','CCT','CCC'};  %%Pro
%cfP=[0.12,0.41,0.31,0.16];
cfP=[1/4,1/4,1/4,1/4];
                                        %% naming principle : 'p'--primilinary; 'c'-- codon or cell;'i'--index;'f'--frequency;
disp('codon subsequence coding Gly number calculation begins');
Gia=find(ismember(codonSequence, 'GGG'));      %%subsequence for Gly begins
if (Gia~=0)
    codonFrequency(Gia)=cfG(1,1);%% Gia, Gib, Gic, Gid -- indices for codons coding Gly 
    NumberCodonGia=length(Gia)
else NumberCodonGia=0
end


Gib=find(ismember(codonSequence,'GGA'));  %% find relevant codon and assign corresponding probability
if (Gib~=0)
    codonFrequency(Gib)=cfG(1,2);
    NumberCodonGib=length(Gib)
else NumberCodonGib=0
end
 

Gic=find(ismember(codonSequence,'GGT'));
if (Gic~=0)
    codonFrequency(Gic)=cfG(1,3);
    NumberCodonGic=length(Gic)
else NumberCodonGic=0
end


Gid=find(ismember(codonSequence,'GGC'));
if (Gid~=0)
    codonFrequency(Gid)=cfG(1,4);
    NumberCodonGid=length(Gid)
else NumberCodonGid=0
end

Gip=[(Gia)',(Gib)',(Gic)',(Gid)'];
Gi=find(Gip);

if (Gi~=0)
Gig=Gip(Gi);   %% Gig -- remove all the zeros among Gia, Gib, Gic, Gid
codonSequenceR=(codonSequence)';  
subsequenceG=codonSequenceR((Gig)'); %%subsequenceG -- codon subsequence coding Gly
Gf=codonFrequency(Gig);       %% remove all the zeros within codonFrequency
disp(subsequenceG);
NNg=length(subsequenceG)

Xg=[NumberCodonGia,NumberCodonGib,NumberCodonGic,NumberCodonGid] %% if subsequence coding for Gly exists, calculate Yg
Pg=cfG          %% cfG is corresponding synonymous fraction values from 'yeast codon table'
Yg=mnpdf(Xg,Pg)  %% by function 'mnpdf' calculate Yg, which is the probability of 'comination of occurrence'(Xg)
Eg=Efor4(NNg)    %% Eg, maximum probability, call function Efor4
Ygg=log(Yg/Eg);  %% Ygg, used for final plotting

else NNg=0   %% if no subsequece coding for Gly, assign relevant values of such subsequence as '0', which for further convenient operation 
    Yg=0;
    Ygg=0;
    disp('codons coding Gly are not existed');     %% subsequence for Gly ends
end



disp('codon subsequence coding Asp number calculation begins'); %%subsequence for Asp begins
Dia=find(ismember(codonSequence, 'GAT'));  
if (Dia~=0)
    codonFrequency(Dia)=cfD(1,1); 
    NumberCodonDia=length(Dia)
    else NumberCodonDia=0
end

Dib=find(ismember(codonSequence,'GAC'));
if (Dib~=0)
    codonFrequency(Dib)=cfD(1,2);
    NumberCodonDib=length(Dib)
    else NumberCodonDib=0
end

Dip=[(Dia)',(Dib)'];
Di=find(Dip);

if (Di~=0)
   Dig=Dip(Di);
   codonSequenceR=(codonSequence)';  
   subsequenceD=codonSequenceR((Dig)'); 
   Df=codonFrequency(Dig);       
   disp(subsequenceD);
   NNd=length(subsequenceD)
   
Xd=[NumberCodonDia,NumberCodonDib]
Pd=cfD
Yd=mnpdf(Xd,Pd)
Ed=Efor2(NNd)
Ydd=log(Yd/Ed);

else NNd=0
    Yd=0;
    Ydd=0;
    disp('codons coding Asp are not existed');     %% subsequence for Asp ends
end


disp('codon subsequence coding Glu number calculation begins');
Eia=find(ismember(codonSequence,'GAG'));
if (Eia~=0)
    codonFrequency(Eia)=cfE(1,1);
    NumberCodonEia=length(Eia)
    else NumberCodonEia=0
end

Eib=find(ismember(codonSequence,'GAA'));
if (Eib~=0)
    codonFrequency(Eib)=cfE(1,2);
    NumberCodonEib=length(Eib)
    else NumberCodonEib=0
end

Eip=[(Eia)',(Eib)'];
Eipi=find(Eip);

if (Eipi~=0)
Eig=Eip(Eipi);
codonSequenceR=(codonSequence)';  
subsequenceE=codonSequenceR((Eig)'); 
Ef=codonFrequency(Eig);       
disp(subsequenceE);
NNe=length(subsequenceE)

Xe=[NumberCodonEia,NumberCodonEib]
Pe=cfE
Ye=mnpdf(Xe,Pe)
Ee=Efor2(NNe)
Yee=log(Ye/Ee);

else NNe=0
    Ye=0;
    Yee=0;
    disp('codons coding Glu are not existed');     %% subsequence for Glu ends
end


disp('codon subsequence coding Val number calculation begins'); %%subsequence for Val begins
Via=find(ismember(codonSequence, 'GTG'));  
if (Via~=0)
    codonFrequency(Via)=cfV(1,1); 
    NumberCodonVia=length(Via)
    else NumberCodonVia=0
end

Vib=find(ismember(codonSequence, 'GTA'));  
if (Vib~=0)
    codonFrequency(Vib)=cfV(1,2); 
    NumberCodonVib=length(Vib)
    else NumberCodonVib=0
end

Vic=find(ismember(codonSequence, 'GTT'));  
if (Vic~=0)
    codonFrequency(Vic)=cfV(1,3);
    NumberCodonVic=length(Vic)
    else NumberCodonVic=0
end

Vid=find(ismember(codonSequence, 'GTC'));  
if (Vid~=0)
    codonFrequency(Vid)=cfV(1,4); 
    NumberCodonVid=length(Vid)
    else NumberCodonVid=0
end

Vip=[(Via)',(Vib)',(Vic)',(Vid)'];
Vi=find(Vip);

if (Vi~=0)
    Vig=Vip(Vi);
    codonSequenceR=(codonSequence)';  
    subsequenceV=codonSequenceR((Vig)'); 
    Vf=codonFrequency(Vig);       
    disp(subsequenceV);
    NNv=length(subsequenceV)
    
Xv=[NumberCodonVia,NumberCodonVib,NumberCodonVic,NumberCodonVid]
Pv=cfV
Yv=mnpdf(Xv,Pv)
Ev=Efor4(NNv)
Yvv=log(Yv/Ev);
else
    NNv=0
    Yv=0;
    Yvv=0;
    disp('codons coding Val are not existed');     %% subsequence for Val ends
end


disp('codon subsequence coding Ala number calculation begins');
Aia=find(ismember(codonSequence, 'GCG'));  
if (Aia~=0)
    codonFrequency(Aia)=cfA(1,1);
    NumberCodonAia=length(Aia)
    else NumberCodonAia=0
end

Aib=find(ismember(codonSequence, 'GCA'));  
if (Aib~=0)
    codonFrequency(Aib)=cfA(1,2);
    NumberCodonAib=length(Aib)
    else NumberCodonAib=0
end

Aic=find(ismember(codonSequence, 'GCT'));  
if (Aic~=0)
    codonFrequency(Aic)=cfA(1,3);
    NumberCodonAic=length(Aic)
    else NumberCodonAic=0
end

Aid=find(ismember(codonSequence, 'GCC'));  
if (Aid~=0)
    codonFrequency(Aid)=cfA(1,4); 
    NumberCodonAid=length(Aid)
    else NumberCodonAid=0
end

Aip=[(Aia)',(Aib)',(Aic)',(Aid)'];
Ai=find(Aip);

if (Ai~=0)
    Aig=Aip(Ai);
    codonSequenceR=(codonSequence)';  
    subsequenceA=codonSequenceR((Aig)'); 
    Af=codonFrequency(Aig);       
    disp(subsequenceA);
    NNa=length(subsequenceA)
    
Xa=[NumberCodonAia,NumberCodonAib,NumberCodonAic,NumberCodonAid]
Pa=cfA
Ya=mnpdf(Xa,Pa)
Ea=Efor4(NNa)
Yaa=log(Ya/Ea);
else
    NNa=0
    Ya=0;
    Yaa=0;
    disp('codons coding Val are not existed');     %% subsequence for Ala ends
end

disp('codon subsequence coding Arg number calculation begins');
Ria=find(ismember(codonSequence, 'AGG'));      %%subsequence for Arg begins
if (Ria~=0)
    codonFrequency(Ria)=cfR(1,1);%% Ria, Rib, Ric, Rid -- indices for codons coding Arg 
    NumberCodonRia=length(Ria)
else NumberCodonRia=0
end


Rib=find(ismember(codonSequence,'AGA'));  %% find relevant codon and assign corresponding probability
if (Rib~=0)
    codonFrequency(Rib)=cfR(1,2);
    NumberCodonRib=length(Rib)
else NumberCodonRib=0
end
 

Ric=find(ismember(codonSequence,'CGG'));
if (Ric~=0)
    codonFrequency(Ric)=cfR(1,3);
    NumberCodonRic=length(Ric)
else NumberCodonRic=0
end


Rid=find(ismember(codonSequence,'CGA'));
if (Rid~=0)
    codonFrequency(Rid)=cfR(1,4);
    NumberCodonRid=length(Rid)
else NumberCodonRid=0
end

Rie=find(ismember(codonSequence,'CGT'));
if (Rie~=0)
    codonFrequency(Rie)=cfR(1,5);
    NumberCodonRie=length(Rie)
else NumberCodonRie=0
end

Rif=find(ismember(codonSequence,'CGC'));
if (Rif~=0)
    codonFrequency(Rif)=cfR(1,6);
    NumberCodonRif=length(Rif)
else NumberCodonRif=0
end
Rip=[(Ria)',(Rib)',(Ric)',(Rid)',(Rie)',(Rif)'];
Ri=find(Rip);

if (Ri~=0)
Rig=Rip(Ri);   %% Rig -- remove all the zeros
codonSequenceR=(codonSequence)';  
subsequenceR=codonSequenceR((Rig)'); %%subsequenceG -- codon subsequence coding Gly
Rf=codonFrequency(Rig);       %% remove all the zeros within codonFrequency
disp(subsequenceR);
NNr=length(subsequenceR)

Xr=[NumberCodonRia,NumberCodonRib,NumberCodonRic,NumberCodonRid,NumberCodonRie,NumberCodonRif]
Pr=cfR
Yr=mnpdf(Xr,Pr)
Er=Efor6(NNr)
Yrr=log(Yr/Er);

else NNr=0
    Yr=0;
    Yrr=0;
    disp('codons coding Arg are not existed');     %% subsequence for Arg ends
end


disp('codon subsequence coding Ser number calculation begins'); %%subsequence for Ser begins
Sia=find(ismember(codonSequence, 'AGT'));  
if (Sia~=0)
    codonFrequency(Sia)=cfS(1,1); 
    NumberCodonSia=length(Sia)
    else NumberCodonSia=0
end

Sib=find(ismember(codonSequence,'AGC'));
if (Sib~=0)
    codonFrequency(Sib)=cfS(1,2);
    NumberCodonSib=length(Sib)
    else NumberCodonSib=0
end

Sic=find(ismember(codonSequence,'TCG'));
if (Sic~=0)
    codonFrequency(Sic)=cfS(1,3);
    NumberCodonSic=length(Sic)
    else NumberCodonSic=0
end

Sid=find(ismember(codonSequence,'TCA'));
if (Sid~=0)
    codonFrequency(Sid)=cfS(1,4);
    NumberCodonSid=length(Sid)
    else NumberCodonSid=0
end

Sie=find(ismember(codonSequence,'TCT'));
if (Sie~=0)
    codonFrequency(Sie)=cfS(1,5);
    NumberCodonSie=length(Sie)
    else NumberCodonSie=0
end

Sif=find(ismember(codonSequence,'TCC'));
if (Sif~=0)
    codonFrequency(Sif)=cfS(1,6);
    NumberCodonSif=length(Sif)
    else NumberCodonSif=0
end

Sip=[(Sia)',(Sib)',(Sic)',(Sid)',(Sie)',(Sif)'];
Si=find(Sip);

if (Si~=0)
   Sig=Sip(Si);
   codonSequenceR=(codonSequence)';  
   subsequenceS=codonSequenceR((Sig)'); 
   Df=codonFrequency(Sig);       
   disp(subsequenceS);
   NNs=length(subsequenceS)
   
Xs=[NumberCodonSia,NumberCodonSib,NumberCodonSic,NumberCodonSid,NumberCodonSie,NumberCodonSif]
Ps=cfS
Ys=mnpdf(Xs,Ps)
Es=Efor6(NNs)
Yss=log(Ys/Es);

else NNs=0
    Ys=0;
    Yss=0;
    disp('codons coding Ser are not existed');     %% subsequence for Ser ends
end


disp('codon subsequence coding Lys number calculation begins');
Kia=find(ismember(codonSequence,'AAG'));
if (Kia~=0)
    codonFrequency(Kia)=cfK(1,1);
    NumberCodonKia=length(Kia)
    else NumberCodonKia=0
end

Kib=find(ismember(codonSequence,'AAA'));
if (Kib~=0)
    codonFrequency(Kib)=cfK(1,2);
    NumberCodonKib=length(Kib)
    else NumberCodonKib=0
end

Kip=[(Kia)',(Kib)'];
Ki=find(Kip);

if (Ki~=0)
Kig=Kip(Ki);
codonSequenceR=(codonSequence)';  
subsequenceK=codonSequenceR((Kig)'); 
Kf=codonFrequency(Kig);       
disp(subsequenceK);
NNk=length(subsequenceK)

Xk=[NumberCodonKia,NumberCodonKib]
Pk=cfK
Yk=mnpdf(Xk,Pk)
Ek=Efor2(NNk)
Ykk=log(Yk/Ek);

else NNk=0
    Yk=0;
    Ykk=0;
    disp('codons coding Lys are not existed');     %% subsequence for Lys ends
end


disp('codon subsequence coding Asn number calculation begins');
Nia=find(ismember(codonSequence,'AAT'));
if (Nia~=0)
    codonFrequency(Nia)=cfN(1,1);
    NumberCodonNia=length(Nia)
    else NumberCodonNia=0
end

Nib=find(ismember(codonSequence,'AAC'));
if (Nib~=0)
    codonFrequency(Nib)=cfN(1,2);
    NumberCodonNib=length(Nib)
    else NumberCodonNib=0
end

Nip=[(Nia)',(Nib)'];
Ni=find(Nip);

if (Ni~=0)
Nig=Nip(Ni);
codonSequenceR=(codonSequence)';  
subsequenceN=codonSequenceR((Nig)'); 
Nf=codonFrequency(Nig);       
disp(subsequenceN);
NNn=length(subsequenceN)

Xn=[NumberCodonNia,NumberCodonNib]
Pn=cfN
Yn=mnpdf(Xn,Pn)
En=Efor2(NNn)
Ynn=log(Yn/En);

else NNn=0
    Yn=0;
    Ynn=0;
    disp('codons coding Asn are not existed');     %% subsequence for Asn ends
end


disp('codon subsequence coding Met number calculation begins');
Mia=find(ismember(codonSequence,'ATG'));
if (Mia~=0)
    codonFrequency(Mia)=cfM(1,1);
    NumberCodonMia=length(Mia)
    else NumberCodonMia=0
end

Mip=(Mia)';
Mi=find(Mip);

if (Mi~=0)
Mig=Mip(Mi);
codonSequenceR=(codonSequence)';  
subsequenceM=codonSequenceR((Mig)'); 
Mf=codonFrequency(Mig);       
disp(subsequenceM);
NNm=length(subsequenceM)

Xm=[NumberCodonMia]
Pm=cfM
Ym=mnpdf(Xm,Pm)
Em=Efor1(NNm)
Ymm=log(Ym/Em);

else NNm=0
    Ym=0;
    Ymm=0;
    disp('codons coding Met are not existed');     
end

disp('codon subsequence coding Ile number calculation begins');
Iia=find(ismember(codonSequence,'ATA'));
if (Iia~=0)
    codonFrequency(Iia)=cfI(1,1);
    NumberCodonIia=length(Iia)
    else NumberCodonIia=0
end

Iib=find(ismember(codonSequence,'ATT'));
if (Iib~=0)
    codonFrequency(Iib)=cfI(1,2);
    NumberCodonIib=length(Iib)
    else NumberCodonIib=0
end

Iic=find(ismember(codonSequence,'ATC'));
if (Iic~=0)
    codonFrequency(Iic)=cfI(1,3);
    NumberCodonIic=length(Iic)
    else NumberCodonIic=0
end

Iip=[(Iia)',(Iib)',(Iic)'];
Ii=find(Iip);

if (Ii~=0)
Iig=Iip(Ii);
codonSequenceR=(codonSequence)';  
subsequenceI=codonSequenceR((Iig)'); 
If=codonFrequency(Iig);       
disp(subsequenceI);
NNi=length(subsequenceI)

Xi=[NumberCodonIia,NumberCodonIib,NumberCodonIic]
Pi=cfI
Yii=mnpdf(Xi,Pi)
Ei=Efor3(NNi)
Yiii=log(Yii/Ei);

else NNi=0
    Yi=0;
    Yiii=0;
    disp('codons coding Ile are not existed');     
end

disp('codon subsequence coding Thr number calculation begins');
Tia=find(ismember(codonSequence,'ACG'));
if (Tia~=0)
    codonFrequency(Tia)=cfT(1,1);
    NumberCodonTia=length(Tia)
    else NumberCodonTia=0
end

Tib=find(ismember(codonSequence,'ACA'));
if (Tib~=0)
    codonFrequency(Tib)=cfT(1,2);
    NumberCodonTib=length(Tib)
    else NumberCodonTib=0
end

Tic=find(ismember(codonSequence,'ACT'));
if (Tic~=0)
    codonFrequency(Tic)=cfT(1,3);
    NumberCodonTic=length(Tic)
    else NumberCodonTic=0
end

Tid=find(ismember(codonSequence,'ACC'));
if (Tid~=0)
    codonFrequency(Tid)=cfT(1,4);
    NumberCodonTid=length(Tid)
    else NumberCodonTid=0
end

Tip=[(Tia)',(Tib)',(Tic)',(Tid)'];
Ti=find(Tip);

if (Ti~=0)
Tig=Tip(Ti);
codonSequenceR=(codonSequence)';  
subsequenceT=codonSequenceR((Tig)'); 
Tf=codonFrequency(Tig);       
disp(subsequenceT);
NNt=length(subsequenceT)

Xt=[NumberCodonTia,NumberCodonTib,NumberCodonTic,NumberCodonTid]
Pt=cfT
Yt=mnpdf(Xt,Pt)
Et=Efor4(NNt)
Ytt=log(Yt/Et);

else NNt=0
    Yt=0;
    Ytt=0;
    disp('codons coding Thr are not existed');     
end

disp('codon subsequence coding Trp number calculation begins');
Wia=find(ismember(codonSequence,'TGG'));
if (Wia~=0)
    codonFrequency(Wia)=cfW(1,1);
    NumberCodonWia=length(Wia)
    else NumberCodonWia=0
end

Wip=(Wia)';
Wi=find(Wip);

if (Wi~=0)
Wig=Wip(Wi);
codonSequenceR=(codonSequence)';  
subsequenceW=codonSequenceR((Wig)'); 
Wf=codonFrequency(Wig);       
disp(subsequenceW);
NNw=length(subsequenceW)

Xw=[NumberCodonWia]
Pw=cfW
Yw=mnpdf(Xw,Pw)
Ew=Efor1(NNw)
Yww=log(Yw/Ew);

else NNw=0
    Yw=0;
    Yww=0;
    disp('codons coding Trp are not existed');     
end

disp('codon subsequence coding End number calculation begins');
Zia=find(ismember(codonSequence,'TGA'));
if (Zia~=0)
    codonFrequency(Zia)=cfZ(1,1);
    NumberCodonZia=length(Zia)
    else NumberCodonZia=0
end

Zib=find(ismember(codonSequence,'TAG'));
if (Zib~=0)
    codonFrequency(Zib)=cfZ(1,2);
    NumberCodonZib=length(Zib)
    else NumberCodonZib=0
end

Zic=find(ismember(codonSequence,'TAA'));
if (Zic~=0)
    codonFrequency(Zic)=cfZ(1,3);
    NumberCodonZic=length(Zic)
    else NumberCodonZic=0
end

Zip=[(Zia)',(Zib)',(Zic)'];
Zi=find(Zip);

if (Zi~=0)
Zig=Zip(Zi);
codonSequenceR=(codonSequence)';  
subsequenceZ=codonSequenceR((Zig)'); 
Zf=codonFrequency(Zig);       
disp(subsequenceZ);
NNz=length(subsequenceZ)

Xz=[NumberCodonZia,NumberCodonZib,NumberCodonZic]
Pz=cfZ
Yz=mnpdf(Xz,Pz)
Ez=Efor3(NNz)
Yzz=log(Yz/Ez);

else NNz=0
    Yz=0;
    Yzz=0;
    disp('codons coding End are not existed');     
end

disp('codon subsequence coding Cys number calculation begins');
Cia=find(ismember(codonSequence,'TGT'));
if (Cia~=0)
    codonFrequency(Cia)=cfC(1,1);
    NumberCodonCia=length(Cia)
    else NumberCodonCia=0
end

Cib=find(ismember(codonSequence,'TGC'));
if (Cib~=0)
    codonFrequency(Cib)=cfC(1,2);
    NumberCodonCib=length(Cib)
    else NumberCodonCib=0
end

Cip=[(Cia)',(Cib)'];
Ci=find(Cip);

if (Ci~=0)
Cig=Cip(Ci);
codonSequenceR=(codonSequence)';  
subsequenceC=codonSequenceR((Cig)'); 
Cf=codonFrequency(Cig);       
disp(subsequenceC);
NNc=length(subsequenceC)

Xc=[NumberCodonCia,NumberCodonCib]
Pc=cfC
Yc=mnpdf(Xc,Pc)
Ec=Efor2(NNc)
Ycc=log(Yc/Ec);

else NNc=0
    Yc=0;
    Ycc=0;
    disp('codons coding Cys are not existed');     
end


disp('codon subsequence coding Tyr number calculation begins');
Yia=find(ismember(codonSequence,'TAT'));
if (Yia~=0)
    codonFrequency(Yia)=cfY(1,1);
    NumberCodonYia=length(Yia)
    else NumberCodonYia=0
end

Yib=find(ismember(codonSequence,'TAC'));
if (Yib~=0)
    codonFrequency(Yib)=cfY(1,2);
    NumberCodonYib=length(Yib)
    else NumberCodonYib=0
end

Yip=[(Yia)',(Yib)'];
Yi=find(Yip);

if (Yi~=0)
Yig=Yip(Yi);
codonSequenceR=(codonSequence)';  
subsequenceY=codonSequenceR((Yig)'); 
Yf=codonFrequency(Yig);       
disp(subsequenceY);
NNy=length(subsequenceY)

Xy=[NumberCodonYia,NumberCodonYib]
Py=cfY
Yy=mnpdf(Xy,Py)
Ey=Efor2(NNy)
Yyy=log(Yy/Ey);

else NNy=0
    Yy=0;
    Yyy=0;
    disp('codons coding Tyr are not existed');     
end


disp('codon subsequence coding Leu number calculation begins');
Lia=find(ismember(codonSequence,'TTG'));
if (Lia~=0)
    codonFrequency(Lia)=cfL(1,1);
    NumberCodonLia=length(Lia)
    else NumberCodonLia=0
end

Lib=find(ismember(codonSequence,'TTA'));
if (Lib~=0)
    codonFrequency(Lib)=cfL(1,2);
    NumberCodonLib=length(Lib)
    else NumberCodonLib=0
end

Lic=find(ismember(codonSequence,'CTG'));
if (Lic~=0)
    codonFrequency(Lic)=cfL(1,3);
    NumberCodonLic=length(Lic)
    else NumberCodonLic=0
end

Lid=find(ismember(codonSequence,'CTA'));
if (Lid~=0)
    codonFrequency(Lid)=cfL(1,4);
    NumberCodonLid=length(Lid)
    else NumberCodonLid=0
end

Lie=find(ismember(codonSequence,'CTT'));
if (Lie~=0)
    codonFrequency(Lie)=cfL(1,5);
    NumberCodonLie=length(Lie)
    else NumberCodonLie=0
end

Lif=find(ismember(codonSequence,'CTC'));
if (Lif~=0)
    codonFrequency(Lif)=cfL(1,6);
    NumberCodonLif=length(Lif)
    else NumberCodonLif=0
end


Lip=[(Lia)',(Lib)',(Lic)',(Lid)',(Lie)',(Lif)'];
Li=find(Lip);

if (Li~=0)
Lig=Lip(Li);
codonSequenceR=(codonSequence)';  
subsequenceL=codonSequenceR((Lig)'); 
Lf=codonFrequency(Lig);       
disp(subsequenceL);
NNl=length(subsequenceL)

Xl=[NumberCodonLia,NumberCodonLib,NumberCodonLic,NumberCodonLid,NumberCodonLie,NumberCodonLif]
Pl=cfL
Yl=mnpdf(Xl,Pl)
El=Efor6(NNl)
Yll=log(Yl/El);

else NNl=0
    Yl=0;
    Yll=0;
    disp('codons coding Leu are not existed');     
end


disp('codon subsequence coding Phe number calculation begins');
Fia=find(ismember(codonSequence,'TTT'));
if (Fia~=0)
    codonFrequency(Fia)=cfF(1,1);
    NumberCodonFia=length(Fia)
    else NumberCodonFia=0
end

Fib=find(ismember(codonSequence,'TTC'));
if (Fib~=0)
    codonFrequency(Fib)=cfF(1,2);
    NumberCodonFib=length(Fib)
    else NumberCodonFib=0
end

Fip=[(Fia)',(Fib)'];
Fi=find(Fip);

if (Fi~=0)
Fig=Fip(Fi);
codonSequenceR=(codonSequence)';  
subsequenceF=codonSequenceR((Fig)'); 
Ff=codonFrequency(Fig);       
disp(subsequenceF);
NNf=length(subsequenceF)

Xf=[NumberCodonFia,NumberCodonFib]
Pf=cfF
Yf=mnpdf(Xf,Pf)
Ef=Efor2(NNf)
Yff=log(Yf/Ef);

else NNf=0
    Yf=0;
    Yff=0;
    disp('codons coding Phe are not existed');     
end


disp('codon subsequence coding Gln number calculation begins');
Qia=find(ismember(codonSequence,'CAG'));
if (Qia~=0)
    codonFrequency(Qia)=cfQ(1,1);
    NumberCodonQia=length(Qia)
    else NumberCodonQia=0
end

Qib=find(ismember(codonSequence,'CAA'));
if (Qib~=0)
    codonFrequency(Qib)=cfQ(1,2);
    NumberCodonQib=length(Qib)
    else NumberCodonQib=0
end

Qip=[(Qia)',(Qib)'];
Qi=find(Qip);

if (Qi~=0)
Qig=Qip(Qi);
codonSequenceR=(codonSequence)';  
subsequenceQ=codonSequenceR((Qig)'); 
Qf=codonFrequency(Qig);       
disp(subsequenceQ);
NNq=length(subsequenceQ)

Xq=[NumberCodonQia,NumberCodonQib]
Pq=cfQ
Yq=mnpdf(Xq,Pq)
Eq=Efor2(NNq)
Yqq=log(Yq/Eq);

else NNq=0
    Yq=0;
    Yqq=0;
    disp('codons coding Gln are not existed');     
end


disp('codon subsequence coding His number calculation begins');
Hia=find(ismember(codonSequence,'CAT'));
if (Hia~=0)
    codonFrequency(Hia)=cfH(1,1);
    NumberCodonHia=length(Hia)
    else NumberCodonHia=0
end

Hib=find(ismember(codonSequence,'CAC'));
if (Hib~=0)
    codonFrequency(Hib)=cfH(1,2);
    NumberCodonHib=length(Hib)
    else NumberCodonHib=0
end

Hip=[(Hia)',(Hib)'];
Hi=find(Hip);

if (Hi~=0)
Hig=Hip(Hi);
codonSequenceR=(codonSequence)';  
subsequenceH=codonSequenceR((Hig)'); 
Hf=codonFrequency(Hig);       
disp(subsequenceH);
NNh=length(subsequenceH)

Xh=[NumberCodonHia,NumberCodonHib]
Ph=cfH
Yh=mnpdf(Xh,Ph)
Eh=Efor2(NNh)
Yhh=log(Yh/Eh);

else NNh=0
    Yh=0;
    Yhh=0;
    disp('codons coding His are not existed');     
end


disp('codon subsequence coding Pro number calculation begins');
Pia=find(ismember(codonSequence,'CCG'));
if (Pia~=0)
    codonFrequency(Pia)=cfP(1,1);
    NumberCodonPia=length(Pia)
    else NumberCodonPia=0
end

Pib=find(ismember(codonSequence,'CCA'));
if (Pib~=0)
    codonFrequency(Pib)=cfP(1,2);
    NumberCodonPib=length(Pib)
    else NumberCodonPib=0
end

Pic=find(ismember(codonSequence,'CCT'));
if (Pic~=0)
    codonFrequency(Pic)=cfP(1,3);
    NumberCodonPic=length(Pic)
    else NumberCodonPic=0
end

Pid=find(ismember(codonSequence,'CCC'));
if (Pid~=0)
    codonFrequency(Pid)=cfP(1,4);
    NumberCodonPid=length(Pid)
    else NumberCodonPid=0
end

Pip=[(Pia)',(Pib)',(Pic)',(Pid)'];
Pi=find(Pip);

if (Pi~=0)
Pig=Pip(Pi);
codonSequenceR=(codonSequence)';  
subsequenceP=codonSequenceR((Pig)'); 
Pf=codonFrequency(Pig);       
disp(subsequenceP);
NNp=length(subsequenceP)

Xp=[NumberCodonPia,NumberCodonPib,NumberCodonPic,NumberCodonPid]
Pp=cfP
Yp=mnpdf(Xp,Pp)
Ep=Efor4(NNp)
Ypp=log(Yp/Ep);

else NNp=0
    Yp=0;
    Ypp=0;
    disp('codons coding Pro are not existed');     
end

Y=[Ygg,Ydd,Yee,Yvv,Yaa,Yrr,Yss,Ykk,Ynn,Ymm,Yiii,Ytt,Yww,Yzz,Ycc,Yyy,Yll,Yff,Yqq,Yhh,Ypp];

disp('End of number calculation of subsequences coding for : Gly,Asp,Glu,Val,Ala,Arg,Ser,Lys,Asn,Met,Ile,Thr,Trp,End,Cys,Tyr,Leu,Phe,Gln,His,Pro');

end

