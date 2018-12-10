function [NNa,Xa,Y] = AlaAminoAcidH(codonSequence)

global cfA ctA;

% disp('codon subsequence coding Ala number calculation begins');

codonFrequency=zeros(1,length(codonSequence));

Aia=find(ismember(codonSequence, 'GCG'));
if (Aia~=0)
    codonFrequency(Aia)=cfA(1,1);
    NumberCodonAia=length(Aia);
else NumberCodonAia=0;
end

Aib=find(ismember(codonSequence, 'GCA'));
if (Aib~=0)
    codonFrequency(Aib)=cfA(1,2);
    NumberCodonAib=length(Aib);
else NumberCodonAib=0;
end

Aic=find(ismember(codonSequence, 'GCT'));
if (Aic~=0)
    codonFrequency(Aic)=cfA(1,3);
    NumberCodonAic=length(Aic);
else NumberCodonAic=0;
end

Aid=find(ismember(codonSequence, 'GCC'));
if (Aid~=0)
    codonFrequency(Aid)=cfA(1,4);
    NumberCodonAid=length(Aid);
else NumberCodonAid=0;
end

Aip=[(Aia)',(Aib)',(Aic)',(Aid)'];
Ai=find(Aip);

if (Ai~=0)
    Aig=Aip(Ai);
    codonSequenceR=(codonSequence)';
    subsequenceA=codonSequenceR((Aig)');
    Af=codonFrequency(Aig);
    NNa=length(subsequenceA);
    
    Xa=[NumberCodonAia,NumberCodonAib,NumberCodonAic,NumberCodonAid];
    Pa=cfA;
    Ya=mnpdf(Xa,Pa);
    Ea=Efor(length(ctA),NNa);
    Yaa=Ya/Ea;
else
    NNa=NaN;
    Xa=NaN;
    Yaa=NaN;
%     disp('codons coding Val are not existed');     %% subsequence for Ala ends
end

Y=Yaa;

end