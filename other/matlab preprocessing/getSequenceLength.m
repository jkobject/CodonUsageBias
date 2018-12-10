%% genes length longer than 100 for Arg Ser Leu
% for 6syno
% function [NNr1,NNs1,NNl1] = getSequenceLength(pasteCodon)
% for i=1:length(pasteCodon)
% NNr1(i)=ArgAminoAcidH(pasteCodon{1,i});
% NNs1(i)=SerAminoAcidH(pasteCodon{1,i});
% NNl1(i)=LeuAminoAcidH(pasteCodon{1,i});
% end

% for 3 syno
% function NNi = getSequenceLength(pasteCodon)
% for i=1:length(pasteCodon)
% NNi(i)=IleAminoAcidH(pasteCodon{1,i});
% end


% % for 4 syno
% function [NNg,NNe,NNv,NNa,NNt,NNp] = getSequenceLength(pasteCodon)
% for i=1:length(pasteCodon)
% NNg(i)=GlyAminoAcidH(pasteCodon{1,i});
% NNe(i)=GluAminoAcidH(pasteCodon{1,i});
% NNv(i)=ValAminoAcidH(pasteCodon{1,i});
% NNa(i)=AlaAminoAcidH(pasteCodon{1,i});
% NNt(i)=ThrAminoAcidH(pasteCodon{1,i});
% NNp(i)=ProAminoAcidH(pasteCodon{1,i});
% end

% for 2 syno
function [NNd,NNe,NNk,NNn,NNc,NNy,NNf,NNq,NNh] = getSequenceLength(pasteCodon)
for i=1:length(pasteCodon)
[NNd(i),~,~]=AspAminoAcidH(pasteCodon{1,i});
[NNe(i),~,~]=GluAminoAcidH(pasteCodon{1,i});
[NNk(i),~,~]=LysAminoAcidH(pasteCodon{1,i});
[NNn(i),~,~]=AsnAminoAcidH(pasteCodon{1,i});
[NNc(i),~,~]=CysAminoAcidH(pasteCodon{1,i});
[NNy(i),~,~]=TyrAminoAcidH(pasteCodon{1,i});
[NNf(i),~,~]=PheAminoAcidH(pasteCodon{1,i});
[NNq(i),~,~]=GlnAminoAcidH(pasteCodon{1,i});
[NNh(i),~,~]=HisAminoAcidH(pasteCodon{1,i});
end
