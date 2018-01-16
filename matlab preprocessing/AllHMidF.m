function MidFq = AllHMidF(pasteCodon)   %%plot ratio of probability distributions (histograms) for21 amino acids &
%%(bin midpoint,log frequency) 
%%among all the mRNAS

global cfG cfGT ctG cfD cfDT ctD cfE cfET ctE cfV cfVT ctV cfA cfAT ctA cfR cfRT ctR cfS cfST ctS cfK cfKT ctK cfN cfNT ctN cfM cfMT ctM cfI cfIT ctI cfT cfTT ctT cfW cfWT ctW cfZ cfZT ctZ cfC cfCT ctC cfY cfYT  ctY cfL cfLT ctL cfF cfFT ctF cfQ cfQT ctQ cfH cfHT ctH cfP cfPT ctP;


% pasteCodon=getCodonSequence; %% call function getCodonSequence to creat input 

NumSynonymous={'4','2','4','4','4','6','6','2','2','3','4','2','2','6','2','2','2','4'};

NameSynonymous={'Gly(4synonymous)','Asp(2synonymous)','Glu(4synonymous)','Val(4synonymous)','Ala(4synonymous)','Arg(6synonymous)','Ser(6synonymous)','Lys(2synonymous)','Asn(2synonymous)','Ile(3synonymous)','Thr(4synonymous)','Cys(2synonymous)','Tyr(2synonymous)','Leu(6synonymous)','Phe(2synonymous)','Gln(2synonymous)','His(2synonymous)','Pro(4synonymous)'};

NameAminoAcid={'Gly','Asp','Glu','Val','Ala','Arg','Ser','Lys','Asn','Ile','Thr','Cys','Tyr','Leu','Phe','Gln','His','Pro'};

% NumSynonymous={'4','2'};  %% used for test
% NameSynonymous={'Gly(4synonymous)','Asp(2synonymous)'};
% NameAminoAcid={'Gly','Asp'};


for i=1:length(pasteCodon) %% note: pasteCodon is column vector
    
disp([num2str(i), 'codon sequence calculation begins']);
G=zeros(1,length(pasteCodon{1,i}));
[NNg,~,G(i)]=GlyAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA
if isnan(NNg)
    continue;
else
    Gr=zeros(1,NNg);

%     RG = replaceE(NNg,ctG); %% Replaced sequence G
    
%     [~,~,Gr(i,:)]=GlyAminoAcidH(RG); %% obtain new ratio matix Gr 


RGT = replaceT(NNg,ctG,cfGT);  %% replace codon according to fraction of codon table 


[~,~,Gr(i)]=GlyAminoAcidH(RGT);

clear G Gr RGT
end


D=zeros(1,length(pasteCodon{1,i}));
[NNd,~,D(i)]=AspAminoAcidH(pasteCodon{1,i}); %%get Ydd for each mRNA
if isnan(NNd)
    continue;
else
%     RD = replaceE(NNd,ctD);
%     
%     [~,~,Dr(i)]=AspAminoAcidH(RD);    


RDT = replaceT(NNd,ctD,cfDT);  

[~,~,Dr(i)]=AspAminoAcidH(RDT);

clear D Dr RDT
end


E=zeros(1,length(pasteCodon{1,i}));
[NNe,~,E(i)]=GluAminoAcidH(pasteCodon{1,i});
if isnan(NNe)
    continue;
else
%     RE = replaceE(NNe,ctE);
%     
%     [~,~,Er(i)]=GluAminoAcidH(RE);    

RET = replaceT(NNe,ctE,cfET);  

[~,~,Er(i)]=GluAminoAcidH(RET);
clear E Er RET
end

V=zeros(1,length(pasteCodon{1,i}));
[NNv,~,V(i)]=ValAminoAcidH(pasteCodon{1,i});
if isnan(NNv)
    continue;
else
%     RV = replaceE(NNv,ctV);
%     
%     [~,~,Vr(i)]=ValAminoAcidH(RV); 

RVT = replaceT(NNv,ctV,cfVT);  

[~,~,Vr(i)]=ValAminoAcidH(RVT);
clear V Vr RVT
end


A=zeros(1,length(pasteCodon{1,i}));
[NNa,~,A(i)]=AlaAminoAcidH(pasteCodon{1,i});
if isnan(NNa)
    continue;
else
%     RA = replaceE(NNa,ctA);
%     
%     [~,~,Ar(i)]=AlaAminoAcidH(RA);   

RAT = replaceT(NNa,ctA,cfAT);  

[~,~,Ar(i)]=AlaAminoAcidH(RAT);
clear A Ar RAT
end


R=zeros(1,length(pasteCodon{1,i}));
[NNr,~,R(i)]=ArgAminoAcidH(pasteCodon{1,i});
if isnan(NNr)
    continue;
else
%     RR = replaceE(NNr,ctR);
%     
%     [~,~,Rr(i)]=ArgAminoAcidH(RR);

RRT = replaceT(NNr,ctR,cfRT);  

[~,~,Rr(i)]=ArgAminoAcidH(RRT);
clear R Rr RRT
end 


S=zeros(1,length(pasteCodon{1,i}));
[NNs,~,S(i)]=SerAminoAcidH(pasteCodon{1,i});

if isnan(NNs)
    continue;
else
%     RS = replaceE(NNs,ctS);
%     
%     [~,~,Sr(i)]=SerAminoAcidH(RS); 

RST = replaceT(NNs,ctS,cfST);  

[~,~,Sr(i)]=SerAminoAcidH(RST);    
clear S Sr RST
end


K=zeros(1,length(pasteCodon{1,i}));
[NNk,~,K(i)]=LysAminoAcidH(pasteCodon{1,i});
if isnan(NNk)
    continue;
else
%     RK = replaceE(NNk,ctK);
%     
%     [~,~,Kr(i)]=LysAminoAcidH(RK); 
RKT = replaceT(NNk,ctK,cfKT);  
[~,~,Kr(i)]=LysAminoAcidH(RKT);  
clear K Kr RKT
end
 

N=zeros(1,length(pasteCodon{1,i}));
[NNn,~,N(i)]=AsnAminoAcidH(pasteCodon{1,i});
if isnan(NNn)
    continue;
else
%     RN = replaceE(NNn,ctN);
%     
%     [~,~,Nr(i)]=AsnAminoAcidH(RN);  
RNT = replaceT(NNn,ctN,cfNT);  
[~,~,Nr(i)]=AsnAminoAcidH(RNT);
clear N Nr RNT
end


I=zeros(1,length(pasteCodon{1,i}));
[NNi,~,I(i)]=IleAminoAcidH(pasteCodon{1,i});
if isnan(NNi)
    continue;
else
%     RI = replaceE(NNi,ctI);
%     
%     [~,~,Ir(i)]=IleAminoAcidH(RI); 
RIT = replaceT(NNi,ctI,cfIT);  
[~,~,Ir(i)]=IleAminoAcidH(RIT);   
clear I Ir RIT
end 


T=zeros(1,length(pasteCodon{1,i}));
[NNt,~,T(i)]=ThrAminoAcidH(pasteCodon{1,i});
if isnan(NNt)
    continue;
else
%     RT = replaceE(NNt,ctT);
%     
%     [~,~,Tr(i)]=ThrAminoAcidH(RT); 
RTT = replaceT(NNt,ctT,cfTT);  
[~,~,Tr(i)]=ThrAminoAcidH(RTT); 
clear T Tr RTT
end


C=zeros(1,length(pasteCodon{1,i}));
[NNc,~,C(i)]=CysAminoAcidH(pasteCodon{1,i});
if isnan(NNc)
    continue;
else
%     RC = replaceE(NNc,ctC);
%     
%     [~,~,Cr(i)]=CysAminoAcidH(RC); 
RCT = replaceT(NNc,ctC,cfCT);  
[~,~,Cr(i)]=CysAminoAcidH(RCT); 
clear C Cr RCT
end


Y=zeros(1,length(pasteCodon{1,i}));
[NNy,~,Y(i,:)]=TyrAminoAcidH(pasteCodon{1,i});
if isnan(NNy)
    continue;
else
%     RY = replaceE(NNy,ctY);
%     
%     [~,~,Yr(i)]=TyrAminoAcidH(RY); 
RYT = replaceT(NNy,ctY,cfYT);  
[~,~,Yr(i)]=TyrAminoAcidH(RYT); 
clear Y Yr RYT
end


L=zeros(1,length(pasteCodon{1,i}));
[NNl,~,L(i)]=LeuAminoAcidH(pasteCodon{1,i});
if isnan(NNl)
    continue;
else
%     RL = replaceE(NNl,ctL);
%     
%     [~,~,Lr(i)]=LeuAminoAcidH(RL);  
RLT = replaceT(NNl,ctL,cfLT);  
[~,~,Lr(i)]=LeuAminoAcidH(RLT);
clear L Lr RLT
end


F=zeros(1,length(pasteCodon{1,i}));
[NNf,~,F(i)]=PheAminoAcidH(pasteCodon{1,i});
if isnan(NNf)
    continue;
else
%     RF = replaceE(NNf,ctF);
%     
%     [~,~,Fr(i)]=PheAminoAcidH(RF); 
RFT = replaceT(NNf,ctF,cfFT);  
[~,~,Fr(i)]=PheAminoAcidH(RFT);
clear F Fr RFT
end


Q=zeros(1,length(pasteCodon{1,i}));
[NNq,~,Q(i)]=GlnAminoAcidH(pasteCodon{1,i});
if isnan(NNq)
    continue;
else
%     RQ = replaceE(NNq,ctQ);
%     
%     [~,~,Qr(i)]=GlnAminoAcidH(RQ);  
RQT = replaceT(NNq,ctQ,cfQT);  
[~,~,Qr(i)]=GlnAminoAcidH(RQT);
clear Q Qr RQT
end


H=zeros(1,length(pasteCodon{1,i}));
[NNh,~,H(i)]=HisAminoAcidH(pasteCodon{1,i});
if isnan(NNh)
    continue;
else
%     RH = replaceE(NNh,ctH);
%     
%     [~,~,Hr(i)]=HisAminoAcidH(RH); 
RHT = replaceT(NNh,ctH,cfHT);  
[~,~,Hr(i)]=HisAminoAcidH(RHT);
clear H Hr RHT
end


P=zeros(1,length(pasteCodon{1,i}));
[NNp,~,P(i)]=ProAminoAcidH(pasteCodon{1,i});
if isnan(NNp)
    continue;
else
%     RP = replaceE(NNp,ctP);
%     
%     [~,~,Pr(i)]=ProAminoAcidH(RP);  
RPT = replaceT(NNp,ctP,cfPT);  
[~,~,Pr(i)]=ProAminoAcidH(RPT);
clear P Pr RPT
end
 
end

% X={G,D,E,V,A,R,S,K,N,I,T,C,Y,L,F,Q,H,P};
% X={Gr,Dr,Er,Vr,Ar,Rr,Sr,Kr,Nr,Ir,Tr,Cr,Yr,Lr,Fr,Qr,Hr,Pr};
% X={log(G),log(D),log(E),log(V),log(A),log(R),log(S),log(K),log(N),log(I),log(T),log(C),log(Y),log(L),log(F),log(Q),log(H),log(P)};
X={log(Gr),log(Dr),log(Er),log(Vr),log(Ar),log(Rr),log(Sr),log(Kr),log(Nr),log(Ir),log(Tr),log(Cr),log(Yr),log(Lr),log(Fr),log(Qr),log(Hr),log(Pr)};

for j=1:length(X)
% Xpr=zeros(length(X{j}),1);
% Xpr=X{j};    
% Xok=Xpr(~(isnan(Xpr))); %% for mean calculation matrix
% Xst=sort(Xok);
% xlen=abs(Xst(length(Xst)-Xst(1)));
% mean=sum(Xok)/xlen;  

% figure 

% h=histogram(X{j},'BinWidth',0.01); %% normalize bin size
% 
% xlim([-180,0]);

% h=histogram(X{j});
% 
% xlabel([NameAminoAcid(j),'ratio of probability']);
% 
% ylabel('frequency');
% 
% title([NameSynonymous(j), 'Histogram in all mRNAs']);


[Nbin,edges]=histcounts(X{j}); %% prepare to plot midpoint,slope.

k=1;
mid=zeros(1,(length(edges)-1)); %%get (mid) midpoints vector
for k=1:(length(edges)-1)
%      if Nbin(k)==0    %%discard outliers
%         mid(k)=NaN;
%         Nbin=NaN;
%     else
        mid(k)=(edges(k)+edges(k+1))/2;
%     end
end


% pointPosition=[mid',Nbin'];
% 
% [hn,cn]=size(pointPosition);

% p=1;
% for p=1:hn
%     pointLabel=['(',num2str(mid(p)),',',num2str(Nbin(p)),')']; %%label points
%     text(mid(p),Nbin(p)+0.1,pointLabel);
% end
Good = find(~(isinf(log(Nbin)))); %%discard NaN values for further plot
midG=mid(Good);
NbinG=Nbin(Good); %% NbinG is the histogram value excluding NaN
NgbinG=log(NbinG);

MidFq{1,j}={(midG)',(NgbinG)'}; 

% plot(mid,log(Nbin),'o');
% hold on;
% % 
% Good = ~(isnan(mid) | isnan(log(Nbin))) ; %%discard NaN values for further plot
% midG=mid(Good);
% NbinG=Nbin(Good); %% NbinG is the histogram value excluding NaN
% NgbinG=log(NbinG); %%NgbinG is the log value 
% 
% % mean=sum(midG.*NbinG)./sum(NbinG);
% 
% slopeOneC = polyfit(midG,NgbinG,1); %%coefficiencies of fitted line; must delet NaN one otherwise polyfit doesn't work
% % 
% % filin=polyval(slopeOneC,midG);
% 
% % yresid=filin-NgbinG;     %% calculate fitting error
% % SSresid=sum(yresid.^2);
% % SStotal=(length(NgbinG)-1)*var(NgbinG);
% % rsq=1-SSresid/SStotal;
% 
% slopeOne=slopeOneC(1);
% % 
% % %plot(midG,filin,'g-',midG,filin+2*delta,'r:',midG,filin-2*delta,'r:');
% % 
% % plot(midG,filin);
% % 
% % legend(['slope:',num2str(slopeOne)]);
% % % 
xlabel('midpoint of bin');
 
ylabel('log of frequency');

title([NameAminoAcid(j),'bin plot']);

mdl = fitlm((midG)',(log(NgbinG))');

figure
plotResiduals(mdl,'probability');
title(NameAminoAcid(j));


meanAll(1,j)=mean;
RsqAll(1,j)=mdl.Rsquared.Ordinary;
SlopeAll(1,j)=slopeOne;

end

tbl = table((NumSynonymous)',(SlopeAll)',(RsqAll)','VariableNames',{'NumSynonymous' 'Slope' 'Rsq'},'RowNames',(NameAminoAcid)');
    
Tbl1 = sortrows(tbl,{'NumSynonymous'},'ascend');

end

