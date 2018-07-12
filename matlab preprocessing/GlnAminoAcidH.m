function [NNq,Xq,Y] = GlnAminoAcidH(codonSequence)

global cfQ ctQ;

codonFrequency=zeros(1,length(codonSequence));
% a matrix of same size as the codon sequence

% disp('codon subsequence coding Gln number calculation begins');
Qia=find(ismember(codonSequence,'CAG')); % a matrix of position of CAG in the codon sequence
if (Qia~=0)
    codonFrequency(Qia)=cfQ(1,1); % at the position where there is 
                                % CAG we set the values CfQ which are average frequencies
    NumberCodonQia=length(Qia);
else NumberCodonQia=0;
end

Qib=find(ismember(codonSequence,'CAA'));
if (Qib~=0)
    codonFrequency(Qib)=cfQ(1,2);
    NumberCodonQib=length(Qib);
else NumberCodonQib=0;
end

Qip=[(Qia)',(Qib)'];
Qi=find(Qip);

if (Qi~=0)
    Qig=Qip(Qi); 
    codonSequenceR=(codonSequence)';
    subsequenceQ=codonSequenceR((Qig)');
    Qf=codonFrequency(Qig);
    NNq=length(subsequenceQ);
    
    Xq=[NumberCodonQia,NumberCodonQib];
    Pq=cfQ;
    Yq=mnpdf(Xq,Pq);
    Eq=Efor(length(ctQ),NNq);
    Yqq=Yq/Eq;
    
else NNq=NaN;
    Yq=NaN;
    Yqq=NaN;
    Xq=NaN;
    subsequenceQ=NaN;
%    disp('codons coding Gln are not existed');
end

Y=Yqq;
end