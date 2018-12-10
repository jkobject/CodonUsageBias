function [NNh,Xh,Y] = HisAminoAcidH(codonSequence)

global cfH ctH;

codonFrequency=zeros(1,length(codonSequence));

% disp('codon subsequence coding His number calculation begins');
Hia=find(ismember(codonSequence,'CAT'));
if (Hia~=0)
    codonFrequency(Hia)=cfH(1,1);
    NumberCodonHia=length(Hia);
else NumberCodonHia=0;
end

Hib=find(ismember(codonSequence,'CAC'));
if (Hib~=0)
    codonFrequency(Hib)=cfH(1,2);
    NumberCodonHib=length(Hib);
else NumberCodonHib=0;
end

Hip=[(Hia)',(Hib)'];
Hi=find(Hip);

if (Hi~=0)
    Hig=Hip(Hi);
    codonSequenceR=(codonSequence)';
    subsequenceH=codonSequenceR((Hig)');
    Hf=codonFrequency(Hig);
    NNh=length(subsequenceH);
    
    Xh=[NumberCodonHia,NumberCodonHib];
    Ph=cfH;
    Yh=mnpdf(Xh,Ph);
    Eh=Efor(length(ctH),NNh);
    Yhh=Yh/Eh;
    
else NNh=NaN;
    Yh=NaN;
    Yhh=NaN;
    Xh=NaN;
    subsequenceH=NaN;
%     disp('codons coding His are not existed');
end

Y=Yhh;
end