function[NNd,Xd,Y]= AspAminoAcidH(codonSequence)

global cfD ctD;

% disp('A new calculation begins');

codonFrequency=zeros(1,length(codonSequence));


%% naming principle : 'p'--primilinary; 'c'-- codon or cell;'i'--index;'f'--frequency;
% disp('codon subsequence coding Asp number calculation begins'); %%subsequence for Asp begins
Dia=find(ismember(codonSequence, 'GAT'));
if (Dia~=0)
    codonFrequency(Dia)=cfD(1,1);
    NumberCodonDia=length(Dia);
else NumberCodonDia=0;
end

Dib=find(ismember(codonSequence,'GAC'));
if (Dib~=0)
    codonFrequency(Dib)=cfD(1,2);
    NumberCodonDib=length(Dib);
else NumberCodonDib=0;
end

Dip=[(Dia)',(Dib)'];
Di=find(Dip);

if (Di~=0)
    Dig=Dip(Di);
    codonSequenceR=(codonSequence)';
    subsequenceD=codonSequenceR((Dig)');
    Df=codonFrequency(Dig);
    NNd=length(subsequenceD);
    
    Xd=[NumberCodonDia,NumberCodonDib];
    Pd=cfD;
    Yd=mnpdf(Xd,Pd);
    Ed=Efor(length(ctD),NNd);
    Ydd=Yd/Ed;
    
else
    NNd=NaN;
    Xd=NaN;
%     subsequenceD=NaN;
%     Yd=NaN;
    Ydd=NaN;
%     disp('codons coding Asp are not existed');     %% subsequence for Asp ends
end

Y=Ydd;
end