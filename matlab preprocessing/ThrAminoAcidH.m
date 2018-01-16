function [NNt,Xt,Y] = ThrAminoAcidH(codonSequence)

global cfT ctT;

codonFrequency=zeros(1,length(codonSequence));

% disp('codon subsequence coding Thr number calculation begins');
Tia=find(ismember(codonSequence,'ACG'));
if (Tia~=0)
    codonFrequency(Tia)=cfT(1,1);
    NumberCodonTia=length(Tia);
else NumberCodonTia=0;
end

Tib=find(ismember(codonSequence,'ACA'));
if (Tib~=0)
    codonFrequency(Tib)=cfT(1,2);
    NumberCodonTib=length(Tib);
else NumberCodonTib=0;
end

Tic=find(ismember(codonSequence,'ACT'));
if (Tic~=0)
    codonFrequency(Tic)=cfT(1,3);
    NumberCodonTic=length(Tic);
else NumberCodonTic=0;
end

Tid=find(ismember(codonSequence,'ACC'));
if (Tid~=0)
    codonFrequency(Tid)=cfT(1,4);
    NumberCodonTid=length(Tid);
else NumberCodonTid=0;
end

Tip=[(Tia)',(Tib)',(Tic)',(Tid)'];
Ti=find(Tip);

if (Ti~=0)
    Tig=Tip(Ti);
    codonSequenceR=(codonSequence)';
    subsequenceT=codonSequenceR((Tig)');
%     Tf=codonFrequency(Tig);
    NNt=length(subsequenceT);
    
    Xt=[NumberCodonTia,NumberCodonTib,NumberCodonTic,NumberCodonTid];
    Pt=cfT;
    Yt=mnpdf(Xt,Pt);
    Et=Efor(length(ctT),NNt);
    Ytt=Yt/Et;
    
else NNt=NaN;
    Yt=NaN;
    Ytt=NaN;
    Xt=NaN;
    subsequenceT=NaN;
%     disp('codons coding Thr are not existed');
end

Y=Ytt;
end