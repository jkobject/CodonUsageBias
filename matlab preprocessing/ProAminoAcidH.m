function [NNp,Xp,Y] = ProAminoAcidH(codonSequence)

global cfP ctP;

codonFrequency=zeros(1,length(codonSequence));

% disp('codon subsequence coding Pro number calculation begins');
Pia=find(ismember(codonSequence,'CCG'));
if (Pia~=0)
    codonFrequency(Pia)=cfP(1,1);
    NumberCodonPia=length(Pia);
else NumberCodonPia=0;
end

Pib=find(ismember(codonSequence,'CCA'));
if (Pib~=0)
    codonFrequency(Pib)=cfP(1,2);
    NumberCodonPib=length(Pib);
else NumberCodonPib=0;
end

Pic=find(ismember(codonSequence,'CCT'));
if (Pic~=0)
    codonFrequency(Pic)=cfP(1,3);
    NumberCodonPic=length(Pic);
else NumberCodonPic=0;
end

Pid=find(ismember(codonSequence,'CCC'));
if (Pid~=0)
    codonFrequency(Pid)=cfP(1,4);
    NumberCodonPid=length(Pid);
else NumberCodonPid=0;
end

Pip=[(Pia)',(Pib)',(Pic)',(Pid)'];
Pi=find(Pip);

if (Pi~=0)
    Pig=Pip(Pi);
    codonSequenceR=(codonSequence)';
    subsequenceP=codonSequenceR((Pig)');
    Pf=codonFrequency(Pig);
    NNp=length(subsequenceP);
    
    Xp=[NumberCodonPia,NumberCodonPib,NumberCodonPic,NumberCodonPid];
    Pp=cfP;
    Yp=mnpdf(Xp,Pp);
    Ep=Efor(length(ctP),NNp);
    Ypp=Yp/Ep;
    
else NNp=NaN;
    Yp=NaN;
    Ypp=NaN;
    Xp=NaN;
    subsequenceP=NaN;
%     disp('codons coding Pro are not existed');
end

Y=Ypp;
end