function [NNz,Xz,Y] = EndAminoAcidH(codonSequence)

global cfZ ctZ;

codonFrequency=zeros(1,length(codonSequence));

% disp('codon subsequence coding End number calculation begins');
Zia=find(ismember(codonSequence,'TGA'));
if (Zia~=0)
    codonFrequency(Zia)=cfZ(1,1);
    NumberCodonZia=length(Zia);
else NumberCodonZia=0;
end

Zib=find(ismember(codonSequence,'TAG'));
if (Zib~=0)
    codonFrequency(Zib)=cfZ(1,2);
    NumberCodonZib=length(Zib);
else NumberCodonZib=0;
end

Zic=find(ismember(codonSequence,'TAA'));
if (Zic~=0)
    codonFrequency(Zic)=cfZ(1,3);
    NumberCodonZic=length(Zic);
else NumberCodonZic=0;
end

Zip=[(Zia)',(Zib)',(Zic)'];
Zi=find(Zip);

if (Zi~=0)
    Zig=Zip(Zi);
    codonSequenceR=(codonSequence)';
    subsequenceZ=codonSequenceR((Zig)');
    Zf=codonFrequency(Zig);
    NNz=length(subsequenceZ);
    
    Xz=[NumberCodonZia,NumberCodonZib,NumberCodonZic];
    Pz=cfZ;
    Yz=mnpdf(Xz,Pz);
    Ez=Efor(length(ctZ),NNz);
    Yzz=Yz/Ez;
    
else NNz=NaN;
    Xz=NaN;
    Yzz=NaN;
%     disp('codons coding End are not existed');
end

Y=Yzz;
end