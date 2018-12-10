function Y = TrpAminoAcidH(codonSequence)

global cfW ctW;

codonFrequency=zeros(1,length(codonSequence));

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
    Ew=Efor(length(ctW),NNw)
    Yww=Yw/Ew;
    
else NNw=0
    Yw=0;
    Yww=NaN;
    disp('codons coding Trp are not existed');
end

Y=Yww;