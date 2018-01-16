function Y = MetAminoAcidH(codonSequence)

global cfM ctM;

codonFrequency=zeros(1,length(codonSequence));

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
    Em=Efor(length(ctM),NNm)
    Ymm=Ym/Em;
    
else NNm=0
    Ym=0;
    Ymm=NaN;
    disp('codons coding Met are not existed');
end

Y=Ymm;
end