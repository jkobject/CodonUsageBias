function [NNf,Xf,Y] = PheAminoAcidH(codonSequence)

global cfF ctF;

codonFrequency=zeros(1,length(codonSequence));

% disp('codon subsequence coding Phe number calculation begins');
Fia=find(ismember(codonSequence,'TTT'));
if (Fia~=0)
    codonFrequency(Fia)=cfF(1,1);
    NumberCodonFia=length(Fia);
else NumberCodonFia=0;
end

Fib=find(ismember(codonSequence,'TTC'));
if (Fib~=0)
    codonFrequency(Fib)=cfF(1,2);
    NumberCodonFib=length(Fib);
else NumberCodonFib=0;
end

Fip=[(Fia)',(Fib)'];
Fi=find(Fip);

if (Fi~=0)
    Fig=Fip(Fi);
    codonSequenceR=(codonSequence)';
    subsequenceF=codonSequenceR((Fig)');
    Ff=codonFrequency(Fig);
    NNf=length(subsequenceF);
    
    Xf=[NumberCodonFia,NumberCodonFib];
    Pf=cfF;
    Yf=mnpdf(Xf,Pf);
    Ef=Efor(length(ctF),NNf);
    Yff=Yf/Ef;
    
else NNf=NaN;
    Yf=NaN;
    Yff=NaN;
    Xf=NaN;
    subsequenceF=NaN;
%     disp('codons coding Phe are not existed');
end

Y=Yff;
end