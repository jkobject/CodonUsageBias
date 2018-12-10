function [NNl,Xl,Y] = LeuAminoAcidH(codonSequence)

global cfL ctL;

% codonFrequency=zeros(1,length(codonSequence));

% disp('codon subsequence coding Leu number calculation begins');
Lia=find(ismember(codonSequence,'TTG'));
if (Lia~=0)
    codonFrequency(Lia)=cfL(1,1);
    NumberCodonLia=length(Lia);
else NumberCodonLia=0;
end

Lib=find(ismember(codonSequence,'TTA'));
if (Lib~=0)
    codonFrequency(Lib)=cfL(1,2);
    NumberCodonLib=length(Lib);
else NumberCodonLib=0;
end

Lic=find(ismember(codonSequence,'CTG'));
if (Lic~=0)
    codonFrequency(Lic)=cfL(1,3);
    NumberCodonLic=length(Lic);
else NumberCodonLic=0;
end

Lid=find(ismember(codonSequence,'CTA'));
if (Lid~=0)
    codonFrequency(Lid)=cfL(1,4);
    NumberCodonLid=length(Lid);
else NumberCodonLid=0;
end

Lie=find(ismember(codonSequence,'CTT'));
if (Lie~=0)
    codonFrequency(Lie)=cfL(1,5);
    NumberCodonLie=length(Lie);
else NumberCodonLie=0;
end

Lif=find(ismember(codonSequence,'CTC'));
if (Lif~=0)
    codonFrequency(Lif)=cfL(1,6);
    NumberCodonLif=length(Lif);
else NumberCodonLif=0;
end


Lip=[(Lia)',(Lib)',(Lic)',(Lid)',(Lie)',(Lif)'];
Li=find(Lip);

if (Li~=0)
    Lig=Lip(Li);
    codonSequenceR=(codonSequence)';
    subsequenceL=codonSequenceR((Lig)');
%     Lf=codonFrequency(Lig);
    NNl=length(subsequenceL);
    
    Xl=[NumberCodonLia,NumberCodonLib,NumberCodonLic,NumberCodonLid,NumberCodonLie,NumberCodonLif];
    Pl=cfL;
    Yl=mnpdf(Xl,Pl);
    El=Efor(length(ctL),NNl);
    Yll=Yl/El;
    
else NNl=NaN;
    Yl=NaN;
    Yll=NaN;
    Xl=NaN;
    subsequenceL=NaN;
%     disp('codons coding Leu are not existed');
end

 Y=Yll;
end