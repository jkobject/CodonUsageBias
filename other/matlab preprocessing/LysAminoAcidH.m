function [NNk,Xk,Y] = LysAminoAcidH(codonSequence)

global cfK ctK;

codonFrequency=zeros(1,length(codonSequence));

% disp('codon subsequence coding Lys number calculation begins');
Kia=find(ismember(codonSequence,'AAG'));
if (Kia~=0)
    codonFrequency(Kia)=cfK(1,1);
    NumberCodonKia=length(Kia);
else NumberCodonKia=0;
end

Kib=find(ismember(codonSequence,'AAA'));
if (Kib~=0)
    codonFrequency(Kib)=cfK(1,2);
    NumberCodonKib=length(Kib);
else NumberCodonKib=0;
end

Kip=[(Kia)',(Kib)'];
Ki=find(Kip);

if (Ki~=0)
    Kig=Kip(Ki);
    codonSequenceR=(codonSequence)';
    subsequenceK=codonSequenceR((Kig)');
    Kf=codonFrequency(Kig);
    NNk=length(subsequenceK);
    
    Xk=[NumberCodonKia,NumberCodonKib];
    Pk=cfK;
    Yk=mnpdf(Xk,Pk);
    Ek=Efor(length(ctK),NNk);
    Ykk=Yk/Ek;
    
else NNk=NaN;
    Yk=NaN;
    Ykk=NaN;
    Xk=NaN;
    subsequenceK=NaN;
%     disp('codons coding Lys are not existed');     %% subsequence for Lys ends
end

Y=Ykk;
end
