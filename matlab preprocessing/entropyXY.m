global cfG cfGT ctG cfD cfDT ctD cfE cfET ctE cfV cfVT ctV cfA cfAT ctA cfR cfRT ctR cfS cfST ctS cfK cfKT ctK cfN cfNT ctN cfM cfMT ctM cfI cfIT ctI cfT cfTT ctT cfW cfWT ctW cfZ cfZT ctZ cfC cfCT ctC cfY cfYT  ctY cfL cfLT ctL cfF cfFT ctF cfQ cfQT ctQ cfH cfHT ctH cfP cfPT ctP;


% pasteCodon=getCodonSequence; %% call function getCodonSequence to creat input 

NumSynonymous={'4','2','4','4','4','6','6','2','2','3','4','2','2','6','2','2','2','4'};

NameSynonymous={'Gly(4synonymous)','Asp(2synonymous)','Glu(4synonymous)','Val(4synonymous)','Ala(4synonymous)','Arg(6synonymous)','Ser(6synonymous)','Lys(2synonymous)','Asn(2synonymous)','Ile(3synonymous)','Thr(4synonymous)','Cys(2synonymous)','Tyr(2synonymous)','Leu(6synonymous)','Phe(2synonymous)','Gln(2synonymous)','His(2synonymous)','Pro(4synonymous)'};

NameAminoAcid={'Gly','Asp','Glu','Val','Ala','Arg','Ser','Lys','Asn','Ile','Thr','Cys','Tyr','Leu','Phe','Gln','His','Pro'};

% NumSynonymous={'4','2'};  %% used for test
% NameSynonymous={'Gly(4synonymous)','Asp(2synonymous)'};
% NameAminoAcid={'Gly','Asp'};


for i=1:length(pasteCodon6) %% note: pasteCodon is column vector
    
disp([num2str(i), 'codon sequence calculation begins']);
 
[NNg,~,G(i,:)]=GlyAminoAcidH(pasteCodon6{1,i}); %%get Ygg for each mRNA
if isnan(NNg)
%     continue;
Gr(i,:)=NaN;
else
    
    rn=1;
    for rn=1:20
    RG = replaceE(NNg,ctG); %% Replaced sequence G
    end
     [~,~,Gr(i,:)]=GlyAminoAcidH(RG); %% obtain new ratio matix Gr 


%RGT = replaceT(NNg,ctG,cfGT);  %% replace codon according to fraction of codon table 

%[~,~,Gr(i,:)]=GlyAminoAcidH(RGT);

end



[NNd,~,D(i,:)]=AspAminoAcidH(pasteCodon6{1,i}); %%get Ydd for each mRNA
if isnan(NNd)
%     continue;
Dr(i,:)=NaN;
else
    rn=1;
    for rn=1:20
    RD = replaceE(NNd,ctD);
    end
    [~,~,Dr(i,:)]=AspAminoAcidH(RD);    


% RDT = replaceT(NNd,ctD,cfDT);  
% [~,~,Dr(i,:)]=AspAminoAcidH(RDT);

end



[NNe,~,E(i,:)]=GluAminoAcidH(pasteCodon6{1,i});
if isnan(NNe)
%     continue;
Er(i,:)=NaN;
else
    RE = replaceE(NNe,ctE);
    
    [~,~,Er(i,:)]=GluAminoAcidH(RE);    

% RET = replaceT(NNe,ctE,cfET);  
% 
% [~,~,Er(i,:)]=GluAminoAcidH(RET);

end


[NNv,~,V(i,:)]=ValAminoAcidH(pasteCodon6{1,i});
if isnan(NNv)
%     continue;
Vr(i,:)=NaN;
else
    RV = replaceE(NNv,ctV);    
    [~,~,Vr(i,:)]=ValAminoAcidH(RV); 

% RVT = replaceT(NNv,ctV,cfVT);  
% [~,~,Vr(i,:)]=ValAminoAcidH(RVT);

end


[NNa,~,A(i,:)]=AlaAminoAcidH(pasteCodon6{1,i});
if isnan(NNa)
%     continue;
Ar(i,:)=NaN;
else
    RA = replaceE(NNa,ctA);    
    [~,~,Ar(i,:)]=AlaAminoAcidH(RA);   

% RAT = replaceT(NNa,ctA,cfAT);  
% [~,~,Ar(i,:)]=AlaAminoAcidH(RAT);
end


[NNr,~,R(i,:)]=ArgAminoAcidH(pasteCodon6{1,i});
if isnan(NNr)
%     continue;
Rr(i,:)=NaN;
else
    RR = replaceE(NNr,ctR);    
    [~,~,Rr(i,:)]=ArgAminoAcidH(RR);

% RRT = replaceT(NNr,ctR,cfRT);  
% [~,~,Rr(i,:)]=ArgAminoAcidH(RRT);
end 



[NNs,~,S(i,:)]=SerAminoAcidH(pasteCodon6{1,i});

if isnan(NNs)
%     continue;
Sr(i,:)=NaN;
else
    RS = replaceE(NNs,ctS);     
    [~,~,Sr(i,:)]=SerAminoAcidH(RS); 

% RST = replaceT(NNs,ctS,cfST);  
% [~,~,Sr(i,:)]=SerAminoAcidH(RST);    
end 


[NNk,~,K(i,:)]=LysAminoAcidH(pasteCodon6{1,i});
if isnan(NNk)
%     continue;
Kr(i,:)=NaN;
else
    RK = replaceE(NNk,ctK);  
    [~,~,Kr(i,:)]=LysAminoAcidH(RK); 
    
% RKT = replaceT(NNk,ctK,cfKT);  
% [~,~,Kr(i,:)]=LysAminoAcidH(RKT);     
end
 

[NNn,~,N(i,:)]=AsnAminoAcidH(pasteCodon6{1,i});
if isnan(NNn)
%     continue;
Nr(i,:)=NaN;
else
    RN = replaceE(NNn,ctN);    
    [~,~,Nr(i,:)]=AsnAminoAcidH(RN); 
    
% RNT = replaceT(NNn,ctN,cfNT);  
% [~,~,Nr(i,:)]=AsnAminoAcidH(RNT);   
end


[NNi,~,I(i,:)]=IleAminoAcidH(pasteCodon6{1,i});
if isnan(NNi)
%     continue;
Ir(i,:)=NaN;
else
    RI = replaceE(NNi,ctI);     
    [~,~,Ir(i,:)]=IleAminoAcidH(RI); 
    
% RIT = replaceT(NNi,ctI,cfIT);  
% [~,~,Ir(i,:)]=IleAminoAcidH(RIT);   

end 


[NNt,~,T(i,:)]=ThrAminoAcidH(pasteCodon6{1,i});
if isnan(NNt)
%     continue;
Tr(i,:)=NaN;
else
    RT = replaceE(NNt,ctT);     
    [~,~,Tr(i,:)]=ThrAminoAcidH(RT); 
    
% RTT = replaceT(NNt,ctT,cfTT);  
% [~,~,Tr(i,:)]=ThrAminoAcidH(RTT); 
end



[NNc,~,C(i,:)]=CysAminoAcidH(pasteCodon6{1,i});
if isnan(NNc)
%     continue;
Cr(i,:)=NaN;
else
    RC = replaceE(NNc,ctC);     
    [~,~,Cr(i,:)]=CysAminoAcidH(RC); 
    
   
% RCT = replaceT(NNc,ctC,cfCT);  
% [~,~,Cr(i,:)]=CysAminoAcidH(RCT); 

end


[NNy,~,Y(i,:)]=TyrAminoAcidH(pasteCodon6{1,i});
if isnan(NNy)
%     continue;
Yr(i,:)=NaN;
else
    RY = replaceE(NNy,ctY);     
    [~,~,Yr(i,:)]=TyrAminoAcidH(RY); 

% RYT = replaceT(NNy,ctY,cfYT);  
% [~,~,Yr(i,:)]=TyrAminoAcidH(RYT); 
end



[NNl,~,L(i,:)]=LeuAminoAcidH(pasteCodon6{1,i});
if isnan(NNl)
%     continue;
Lr(i,:)=NaN;
else
    RL = replaceE(NNl,ctL);     
    [~,~,Lr(i,:)]=LeuAminoAcidH(RL); 

% RLT = replaceT(NNl,ctL,cfLT);  
% [~,~,Lr(i,:)]=LeuAminoAcidH(RLT);

end



[NNf,~,F(i,:)]=PheAminoAcidH(pasteCodon6{1,i});
if isnan(NNf)
%     continue;
Fr(i,:)=NaN;
else
    RF = replaceE(NNf,ctF);  
    [~,~,Fr(i,:)]=PheAminoAcidH(RF); 
    
% RFT = replaceT(NNf,ctF,cfFT);  
% [~,~,Fr(i,:)]=PheAminoAcidH(RFT);
end



[NNq,~,Q(i,:)]=GlnAminoAcidH(pasteCodon6{1,i});
if isnan(NNq)
%     continue;
Qr(i,:)=NaN;
else
    RQ = replaceE(NNq,ctQ);  
    [~,~,Qr(i,:)]=GlnAminoAcidH(RQ);
    
% RQT = replaceT(NNq,ctQ,cfQT);  
% [~,~,Qr(i,:)]=GlnAminoAcidH(RQT);
end



[NNh,~,H(i,:)]=HisAminoAcidH(pasteCodon6{1,i});
if isnan(NNh)
%     continue;
Hr(i,:)=NaN;
else
    RH = replaceE(NNh,ctH);     
    [~,~,Hr(i,:)]=HisAminoAcidH(RH); 
    
% RHT = replaceT(NNh,ctH,cfHT);  
% [~,~,Hr(i,:)]=HisAminoAcidH(RHT);
end



[NNp,~,P(i,:)]=ProAminoAcidH(pasteCodon6{1,i});
if isnan(NNp)
%     continue;
Pr(i,:)=NaN;
else
    RP = replaceE(NNp,ctP);     
    [~,~,Pr(i,:)]=ProAminoAcidH(RP);  
    
% RPT = replaceT(NNp,ctP,cfPT);  
% [~,~,Pr(i,:)]=ProAminoAcidH(RPT);
end
 
end

% X={G,D,E,V,A,R,S,K,N,I,T,C,Y,L,F,Q,H,P};
% X={Gr,Dr,Er,Vr,Ar,Rr,Sr,Kr,Nr,Ir,Tr,Cr,Yr,Lr,Fr,Qr,Hr,Pr};
X1={log(G),log(D),log(E),log(V),log(A),log(R),log(S),log(K),log(N),log(I),log(T),log(C),log(Y),log(L),log(F),log(Q),log(H),log(P)};
X2={log(Gr),log(Dr),log(Er),log(Vr),log(Ar),log(Rr),log(Sr),log(Kr),log(Nr),log(Ir),log(Tr),log(Cr),log(Yr),log(Lr),log(Fr),log(Qr),log(Hr),log(Pr)};

for j=1:18
    
figure    
    
plot(X1{j},X2{j},'o');

xlabel('entropy of original sequence');
 
ylabel('entropy of sequence after codon replacement');

title(['sj: ',NameSynonymous(j)]);
    
end