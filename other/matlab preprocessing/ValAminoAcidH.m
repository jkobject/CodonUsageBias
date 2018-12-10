function [NNv,Xv,Y] = ValAminoAcidH(codonSequence)

global cfV ctV;

% disp('A new calculation begins');

codonFrequency=zeros(1,length(codonSequence));


%% naming principle : 'p'--primilinary; 'c'-- codon or cell;'i'--index;'f'--frequency;
% disp('codon subsequence coding Val number calculation begins'); %%subsequence for Val begins
Via=find(ismember(codonSequence, 'GTG'));
if (Via~=0)
    codonFrequency(Via)=cfV(1,1);
    NumberCodonVia=length(Via);
else NumberCodonVia=0;
end

Vib=find(ismember(codonSequence, 'GTA'));
if (Vib~=0)
    codonFrequency(Vib)=cfV(1,2);
    NumberCodonVib=length(Vib);
else NumberCodonVib=0;
end

Vic=find(ismember(codonSequence, 'GTT'));
if (Vic~=0)
    codonFrequency(Vic)=cfV(1,3);
    NumberCodonVic=length(Vic);
else NumberCodonVic=0;
end

Vid=find(ismember(codonSequence, 'GTC'));
if (Vid~=0)
    codonFrequency(Vid)=cfV(1,4);
    NumberCodonVid=length(Vid);
else NumberCodonVid=0;
end

Vip=[(Via)',(Vib)',(Vic)',(Vid)'];
Vi=find(Vip);

if (Vi~=0)
    Vig=Vip(Vi);
    codonSequenceR=(codonSequence)';
    subsequenceV=codonSequenceR((Vig)');
    Vf=codonFrequency(Vig);
    NNv=length(subsequenceV);
    
    Xv=[NumberCodonVia,NumberCodonVib,NumberCodonVic,NumberCodonVid];
    Pv=cfV;
    Yv=mnpdf(Xv,Pv);
    Ev=Efor(length(ctV),NNv);
    Yvv=Yv/Ev;
else
    NNv=NaN;
    Yv=NaN;
    Xv=NaN;
    subsequenceV=NaN;
    Yvv=NaN;
%     disp('codons coding Val are not existed');     %% subsequence for Val ends
end

Y=Yvv;
end