%%plot ratio of probability distributions (histograms) for21 amino acids &
%%(bin midpoint,log frequency) 
%%among all the mRNAS
clc 
clear


pasteCodon=getCodonSequence; %% call function getCodonSequence to creat input 

NameSynonymous={'Gly(4synonymous)','Asp(2synonymous)','Glu(4synonymous)','Val(4synonymous)','Ala(4synonymous)','Arg(6synonymous)','Ser(6synonymous)','Lys(2synonymous)','Asn(2synonymous)','Met(1synonymous)','Ile(3synonymous)','Thr(4synonymous)','Trp(1synonymous)','End(3synonymous)','Cys(2synonymous)','Tyr(2synonymous)','Leu(6synonymous)','Phe(2synonymous)','Gln(2synonymous)','His(2synonymous)','Pro(4synonymous)'};

NameAminoAcid={'Gly','Asp','Glu','Val','Ala','Arg','Ser','Lys','Asn','Met','Ile','Thr','Trp','End','Cys','Tyr','Leu','Phe','Gln','His','Pro'};


for i=1:length(pasteCodon) %% note: pasteCodon is column vector
    
G(i,:)=GlyAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA

D(i,:)=AspAminoAcidH(pasteCodon{1,i}); %%get Ydd for each mRNA

E(i,:)=GluAminoAcidH(pasteCodon{1,i});

V(i,:)=ValAminoAcidH(pasteCodon{1,i});

A(i,:)=AlaAminoAcidH(pasteCodon{1,i});

R(i,:)=ArgAminoAcidH(pasteCodon{1,i});

S(i,:)=SerAminoAcidH(pasteCodon{1,i});

K(i,:)=LysAminoAcidH(pasteCodon{1,i});

N(i,:)=AsnAminoAcidH(pasteCodon{1,i});

M(i,:)=MetAminoAcidH(pasteCodon{1,i});

I(i,:)=IleAminoAcidH(pasteCodon{1,i});

T(i,:)=ThrAminoAcidH(pasteCodon{1,i});

W(i,:)=TrpAminoAcidH(pasteCodon{1,i});

Z(i,:)=EndAminoAcidH(pasteCodon{1,i});

C(i,:)=CysAminoAcidH(pasteCodon{1,i});

Y(i,:)=TyrAminoAcidH(pasteCodon{1,i});

L(i,:)=LeuAminoAcidH(pasteCodon{1,i});

F(i,:)=PheAminoAcidH(pasteCodon{1,i});

Q(i,:)=GlnAminoAcidH(pasteCodon{1,i});

H(i,:)=HisAminoAcidH(pasteCodon{1,i});

P(i,:)=ProAminoAcidH(pasteCodon{1,i});

end

X={G,D,E,V,A,R,S,K,N,M,I,T,W,Z,C,Y,L,F,Q,H,P};

for j=1:21

figure 

histogram(X{j},'BinWidth',3);

xlim([-170,0]);

[Nbin,edges]=histcounts(X{j});

k=1;
mid=zeros(1,(length(edges)-1)); %%get (mid) midpoints vector
for k=1:(length(edges)-1)
mid(k)=(edges(k)+edges(k+1))/2;
end

pointPosition=[mid',Nbin'];

[hn,cn]=size(pointPosition);

p=1;
for p=1:hn
    pointLabel=['(',num2str(mid(p)),',',num2str(Nbin(p)),')']; %%label points
    text(mid(p),Nbin(p)+0.1,pointLabel);
end

xlabel([NameAminoAcid(j),'ratio of probability']);

ylabel('frequency');

title([NameSynonymous(j), 'Histogram in all mRNAs']);

figure 

plot(mid,log(Nbin),'o');  %% plot (midpoint,log frequency)

xlabel('midpoint of bin');

ylabel('log of frequency');

title([NameAminoAcid(j),'bin plot']);

end

