function [NNy,Xy,Y] = TyrAminoAcidH(codonSequence)

global cfY ctY;

codonFrequency=zeros(1,length(codonSequence));

% disp('codon subsequence coding Tyr number calculation begins');
Yia=find(ismember(codonSequence,'TAT'));
if (Yia~=0)
    codonFrequency(Yia)=cfY(1,1);
    NumberCodonYia=length(Yia);
else NumberCodonYia=0;
end

Yib=find(ismember(codonSequence,'TAC'));
if (Yib~=0)
    codonFrequency(Yib)=cfY(1,2);
    NumberCodonYib=length(Yib);
else NumberCodonYib=0;
end

Yip=[(Yia)',(Yib)'];
Yi=find(Yip);

if (Yi~=0)
    Yig=Yip(Yi);
    codonSequenceR=(codonSequence)';
    subsequenceY=codonSequenceR((Yig)');
    Yf=codonFrequency(Yig);
    NNy=length(subsequenceY);
    
    Xy=[NumberCodonYia,NumberCodonYib];
    Py=cfY;
    Yy=mnpdf(Xy,Py);
    Ey=Efor(length(ctY),NNy);
    Yyy=Yy/Ey;
    
else NNy=NaN;
    Yy=NaN;
    Yyy=NaN;
    Xy=NaN;
    subsequenceY=NaN;
%     disp('codons coding Tyr are not existed');
end

Y=Yyy;
end