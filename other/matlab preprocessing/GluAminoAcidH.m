function [NNe,Xe,Y] = GluAminoAcidH(codonSequence)

global cfE ctE;

%% naming principle : 'p'--primilinary; 'c'-- codon or cell;'i'--index;'f'--frequency;
% disp('codon subsequence coding Glu number calculation begins');
Eia=find(ismember(codonSequence,'GAG'));
if (Eia~=0)
    NumberCodonEia=length(Eia);
else NumberCodonEia=0; 
end

Eib=find(ismember(codonSequence,'GAA'));
if (Eib~=0)
    NumberCodonEib=length(Eib);
else NumberCodonEib=0;
end

Eip=[(Eia)',(Eib)'];
Eipi=find(Eip);

if (Eipi~=0)
    Eig=Eip(Eipi);
    codonSequenceR=(codonSequence)';
    subsequenceE=codonSequenceR((Eig)');
    NNe=length(subsequenceE);
    
    Xe=[NumberCodonEia,NumberCodonEib];
    Pe=cfE;
    Ye=mnpdf(Xe,Pe);
    Ee=Efor(length(ctE),NNe);
    Yee=Ye/Ee;
    
else NNe=NaN;
    Ye=NaN;
    Xe=NaN;
    subsequenceE=NaN;
    Yee=NaN;
%     disp('codons coding Glu are not existed');     %% subsequence for Glu ends
end

Y=Yee;
end