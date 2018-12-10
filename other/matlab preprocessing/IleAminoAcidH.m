function [NNi,Xi,Y] = IleAminoAcidH(codonSequence)

global cfI ctI;

codonFrequency=zeros(1,length(codonSequence));

% disp('codon subsequence coding Ile number calculation begins');
Iia=find(ismember(codonSequence,'ATA'));
if (Iia~=0)
    codonFrequency(Iia)=cfI(1,1);
    NumberCodonIia=length(Iia);
else NumberCodonIia=0;
end

Iib=find(ismember(codonSequence,'ATT'));
if (Iib~=0)
    codonFrequency(Iib)=cfI(1,2);
    NumberCodonIib=length(Iib);
else NumberCodonIib=0;
end

Iic=find(ismember(codonSequence,'ATC'));
if (Iic~=0)
    codonFrequency(Iic)=cfI(1,3);
    NumberCodonIic=length(Iic);
else NumberCodonIic=0;
end

Iip=[(Iia)',(Iib)',(Iic)'];
Ii=find(Iip);

if (Ii~=0)
    Iig=Iip(Ii);
    codonSequenceR=(codonSequence)';
    subsequenceI=codonSequenceR((Iig)');
    If=codonFrequency(Iig);
    NNi=length(subsequenceI);
    
    Xi=[NumberCodonIia,NumberCodonIib,NumberCodonIic];
    Pi=cfI;
    Yii=mnpdf(Xi,Pi);
    Ei=Efor(length(ctI),NNi);
    Yiii=Yii/Ei;
    
else NNi=NaN;
%     Yi=NaN;
    Yiii=NaN;
    Xi=NaN;
%     subsequenceI=NaN;
%     disp('codons coding Ile are not existed');
end

Y=Yiii;
end