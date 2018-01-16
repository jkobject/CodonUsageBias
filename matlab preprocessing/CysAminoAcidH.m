function [NNc,Xc,Y] = CysAminoAcidH(codonSequence)

global cfC ctC;

codonFrequency=zeros(1,length(codonSequence));

% disp('codon subsequence coding Cys number calculation begins');
Cia=find(ismember(codonSequence,'TGT'));
if (Cia~=0)
    codonFrequency(Cia)=cfC(1,1);
    NumberCodonCia=length(Cia);
else NumberCodonCia=0;
end

Cib=find(ismember(codonSequence,'TGC'));
if (Cib~=0)
    codonFrequency(Cib)=cfC(1,2);
    NumberCodonCib=length(Cib);
else NumberCodonCib=0;
end

Cip=[(Cia)',(Cib)'];
Ci=find(Cip);

if (Ci~=0)
    Cig=Cip(Ci);
    codonSequenceR=(codonSequence)';
    subsequenceC=codonSequenceR((Cig)');
    Cf=codonFrequency(Cig);
    NNc=length(subsequenceC);
    
    Xc=[NumberCodonCia,NumberCodonCib];
    Pc=cfC;
    Yc=mnpdf(Xc,Pc);
    Ec=Efor(length(ctC),NNc);
    Ycc=Yc/Ec;
    
else NNc=NaN;
    Yc=NaN;
    Ycc=NaN;
    Xc=NaN;
    subsequenceC=NaN;
%     disp('codons coding Cys are not existed');
end

Y=Ycc;
end