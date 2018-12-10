function [NNs,Xs,Y] = SerAminoAcidH(codonSequence)

global cfS ctS;

codonFrequency=zeros(1,length(codonSequence));


% disp('codon subsequence coding Ser number calculation begins'); %%subsequence for Ser begins
Sia=find(ismember(codonSequence, 'AGT'));
if (Sia~=0)
    codonFrequency(Sia)=cfS(1,1);
    NumberCodonSia=length(Sia);
else NumberCodonSia=0;
end

Sib=find(ismember(codonSequence,'AGC'));
if (Sib~=0)
    codonFrequency(Sib)=cfS(1,2);
    NumberCodonSib=length(Sib);
else NumberCodonSib=0;
end

Sic=find(ismember(codonSequence,'TCG'));
if (Sic~=0)
    codonFrequency(Sic)=cfS(1,3);
    NumberCodonSic=length(Sic);
else NumberCodonSic=0;
end

Sid=find(ismember(codonSequence,'TCA'));
if (Sid~=0)
    codonFrequency(Sid)=cfS(1,4);
    NumberCodonSid=length(Sid);
else NumberCodonSid=0;
end

Sie=find(ismember(codonSequence,'TCT'));
if (Sie~=0)
    codonFrequency(Sie)=cfS(1,5);
    NumberCodonSie=length(Sie);
else NumberCodonSie=0;
end

Sif=find(ismember(codonSequence,'TCC'));
if (Sif~=0)
    codonFrequency(Sif)=cfS(1,6);
    NumberCodonSif=length(Sif);
else NumberCodonSif=0;
end

Sip=[(Sia)',(Sib)',(Sic)',(Sid)',(Sie)',(Sif)'];
Si=find(Sip);

if (Si~=0)
    Sig=Sip(Si);
    codonSequenceR=(codonSequence)';
    subsequenceS=codonSequenceR((Sig)');
%     Df=codonFrequency(Sig);
    NNs=length(subsequenceS);
     
Xs=[NumberCodonSia,NumberCodonSib,NumberCodonSic,NumberCodonSid,NumberCodonSie,NumberCodonSif];
    Ps=cfS;
    Ys=mnpdf(Xs,Ps);
    Es=Efor(length(ctS),NNs);
    Yss=Ys/Es;
    
else NNs=NaN;
    Ys=NaN;
    Yss=NaN;
    Xs=NaN;
    subsequenceS=NaN;
%     disp('codons coding Ser are not existed');     %% subsequence for Ser ends
end

 Y=Yss;
end