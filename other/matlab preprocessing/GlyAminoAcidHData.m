function [Y,subsequenceG,Xg] = GlyAminoAcidHData(codonSequence)

global cfG ctG;

codonFrequency=zeros(1,length(codonSequence));

%% naming principle : 'p'--primilinary; 'c'-- codon or cell;'i'--index;'f'--frequency;
disp('codon subsequence coding Gly number calculation begins');
Gia=find(ismember(codonSequence, 'GGG'));      %%subsequence for Gly begins
if (Gia~=0)
    codonFrequency(Gia)=cfG(1,1);%% Gia, Gib, Gic, Gid -- indices for codons coding Gly
    NumberCodonGia=length(Gia);
else NumberCodonGia=0;
end


Gib=find(ismember(codonSequence,'GGA'));  %% find relevant codon and assign corresponding probability
if (Gib~=0)
    codonFrequency(Gib)=cfG(1,2);
    NumberCodonGib=length(Gib);
else NumberCodonGib=0;
end


Gic=find(ismember(codonSequence,'GGT'));
if (Gic~=0)
    codonFrequency(Gic)=cfG(1,3);
    NumberCodonGic=length(Gic);
else NumberCodonGic=0;
end


Gid=find(ismember(codonSequence,'GGC'));
if (Gid~=0)
    codonFrequency(Gid)=cfG(1,4);
    NumberCodonGid=length(Gid);
else NumberCodonGid=0;
end

Gip=[(Gia)',(Gib)',(Gic)',(Gid)'];
Gi=find(Gip);

if (Gi~=0)
    Gig=Gip(Gi);   %% Gig -- remove all the zeros among Gia, Gib, Gic, Gid
    codonSequenceR=(codonSequence)';
    subsequenceG=codonSequenceR((Gig)'); %%subsequenceG -- codon subsequence coding Gly
    Gf=codonFrequency(Gig);       %% remove all the zeros within codonFrequency
    disp(subsequenceG);
    NNg=length(subsequenceG);
    
    Xg=[NumberCodonGia,NumberCodonGib,NumberCodonGic,NumberCodonGid] %% if subsequence coding for Gly exists, calculate Yg
    Pg=cfG;          %% cfG is corresponding synonymous fraction values from 'yeast codon table'
    Yg=mnpdf(Xg,Pg);  %% by function 'mnpdf' calculate Yg, which is the probability of 'comination of occurrence'(Xg)
    Eg=Efor(length(ctG),NNg);    %% Eg, maximum probability, call function Efor4
    Ygg=Yg/Eg;  %% Ygg, used for final plotting
    
else NNg=0;%% if no subsequece coding for Gly, assign relevant values of such subsequence as '0', which for further convenient operation
    Xg=NaN;
    subsequenceG=NaN;
    Yg=NaN;
    Ygg=NaN;
    disp('codons coding Gly are not existed');     %% subsequence for Gly ends
end

Y=Ygg


end