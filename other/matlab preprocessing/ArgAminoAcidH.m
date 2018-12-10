function [NNr,Xr,Y] = ArgAminoAcidH(codonSequence)


global cfR ctR;

% disp('codon subsequence coding Arg number calculation begins');
% codonFrequency=zeros(1,length(codonSequence));


Ria=find(ismember(codonSequence, 'AGG'));      %%subsequence for Arg begins
if (Ria~=0)
%     codonFrequency(Ria)=cfR(1,1);%% Ria, Rib, Ric, Rid -- indices for codons coding Arg
    NumberCodonRia=length(Ria);
else NumberCodonRia=0;
end


Rib=find(ismember(codonSequence,'AGA'));  %% find relevant codon and assign corresponding probability
if (Rib~=0)
    codonFrequency(Rib)=cfR(1,2);
    NumberCodonRib=length(Rib);
else NumberCodonRib=0;
end


Ric=find(ismember(codonSequence,'CGG'));
if (Ric~=0)
    codonFrequency(Ric)=cfR(1,3);
    NumberCodonRic=length(Ric);
else NumberCodonRic=0;
end


Rid=find(ismember(codonSequence,'CGA'));
if (Rid~=0)
    codonFrequency(Rid)=cfR(1,4);
    NumberCodonRid=length(Rid);
else NumberCodonRid=0;
end

Rie=find(ismember(codonSequence,'CGT'));
if (Rie~=0)
    codonFrequency(Rie)=cfR(1,5);
    NumberCodonRie=length(Rie);
else NumberCodonRie=0;
end

Rif=find(ismember(codonSequence,'CGC'));
if (Rif~=0)
    codonFrequency(Rif)=cfR(1,6);
    NumberCodonRif=length(Rif);
else NumberCodonRif=0;
end
Rip=[(Ria)',(Rib)',(Ric)',(Rid)',(Rie)',(Rif)'];
Ri=find(Rip);

if (Ri~=0)
    Rig=Rip(Ri);   %% Rig -- remove all the zeros
    codonSequenceR=(codonSequence)';
    subsequenceR=codonSequenceR((Rig)'); %%subsequenceG -- codon subsequence coding Gly
%     Rf=codonFrequency(Rig);       %% remove all the zeros within codonFrequency
    NNr=length(subsequenceR);
    
    Xr=[NumberCodonRia,NumberCodonRib,NumberCodonRic,NumberCodonRid,NumberCodonRie,NumberCodonRif];
    Pr=cfR;
    Yr=mnpdf(Xr,Pr);
    Er=Efor(length(ctR),NNr);
    Yrr=Yr/Er;
    
else NNr=NaN;
    Yrr=NaN;
    Xr=NaN;
%     disp('codons coding Arg are not existed');     %% subsequence for Arg ends
end

Y=Yrr;
end