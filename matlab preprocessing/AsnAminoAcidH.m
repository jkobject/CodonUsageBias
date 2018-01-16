function [NNn,Xn,Y] = AsnAminoAcidH(codonSequence)

global cfN ctN;

codonFrequency=zeros(1,length(codonSequence));

% disp('codon subsequence coding Asn number calculation begins');
Nia=find(ismember(codonSequence,'AAT'));
if (Nia~=0)
    codonFrequency(Nia)=cfN(1,1);
    NumberCodonNia=length(Nia);
else NumberCodonNia=0;
end

Nib=find(ismember(codonSequence,'AAC'));
if (Nib~=0)
    codonFrequency(Nib)=cfN(1,2);
    NumberCodonNib=length(Nib);
else NumberCodonNib=0;
end

Nip=[(Nia)',(Nib)'];
Ni=find(Nip);

if (Ni~=0)
    Nig=Nip(Ni);
    codonSequenceR=(codonSequence)';
    subsequenceN=codonSequenceR((Nig)');
    Nf=codonFrequency(Nig);
    NNn=length(subsequenceN);
    
    Xn=[NumberCodonNia,NumberCodonNib];
    Pn=cfN;
    Yn=mnpdf(Xn,Pn);
    En=Efor(length(ctN),NNn);
    Ynn=Yn/En;
    
else NNn=NaN;
    Ynn=NaN;
    Xn=NaN;
%     subsequenceN=NaN;
%     disp('codons coding Asn are not existed');     %% subsequence for Asn ends
end

Y=Ynn;
end