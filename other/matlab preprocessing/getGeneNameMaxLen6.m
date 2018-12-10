%% find subsequence maximum length among 20 species 
%% can write a table for clear view

sequenceName={sequenceName11,sequenceName12,sequenceName13,sequenceName14,sequenceName15,...
    sequenceName21,sequenceName22,sequenceName23,...
    sequenceName31,sequenceName32,sequenceName33,sequenceName34,sequenceName35,sequenceName36,sequenceName37...
    sequenceName41,sequenceName42,sequenceName43,sequenceName44,sequenceName45};

PasteCodon={pasteCodon11,pasteCodon12,pasteCodon13,pasteCodon14,pasteCodon15,...
    pasteCodon21,pasteCodon22,pasteCodon23,...
    pasteCodon31,pasteCodon32,pasteCodon33,pasteCodon34,pasteCodon35,pasteCodon36,pasteCodon37,...
    pasteCodon41,pasteCodon42,pasteCodon43,pasteCodon44,pasteCodon45};

SequenceName={'sah6','saeu','saku','Ag','Yl','sp','so','sj','ac','afl','ans','anr','ao','at','afu','ff','fg','fo','fv','tr'};

lenArrR={NNr11,NNr12,NNr13,NNr14,NNr15,NNr21,NNr22,NNr23,NNr31,NNr32,NNr33,NNr34,NNr35,NNr36,NNr37,NNr41,NNr42,NNr43,NNr44,NNr45};
% [lmaxIdR,lmaxR]=lmaxIndex(lenArrR);
[lIdR,lR]=getLabove(lenArrR,200);

for i=1:20
sequenceTempName=sequenceName{i};
lindexTemp=lIdR{i};
pasteCodonTemp=PasteCodon{i};
wholeLenR{i}=length(pasteCodonTemp{lindexTemp});
GeneNameR{i}=sequenceTempName{lindexTemp};
end

lenArrS={NNs11,NNs12,NNs13,NNs14,NNs15,NNs21,NNs22,NNs23,NNs31,NNs32,NNs33,NNs34,NNs35,NNs36,NNs37,NNs41,NNs42,NNs43,NNs44,NNs45};
[lmaxIdS,lmaxS]=lmaxIndex(lenArrS);

for i=1:20
sequenceTempName=sequenceName{i};
lindexTemp=lmaxIdS(i);
pasteCodonTemp=PasteCodon{i};
wholeLenS(i)=length(pasteCodonTemp{lindexTemp});
GeneNameS{i}=sequenceTempName{lindexTemp};
end

lenArrL={NNl11,NNl12,NNl13,NNl14,NNl15,NNl21,NNl22,NNl23,NNl31,NNl32,NNl33,NNl34,NNl35,NNl36,NNl37,NNl41,NNl42,NNl43,NNl44,NNl45};
[lmaxIdL,lmaxL]=lmaxIndex(lenArrL);

for i=1:20
sequenceTempName=sequenceName{i};
lindexTemp=lmaxIdL(i);
pasteCodonTemp=PasteCodon{i};
wholeLenL(i)=length(pasteCodonTemp{lindexTemp});
GeneNameL{i}=sequenceTempName{lindexTemp};
end

tblMaxLen6 = table((lmaxR)',(GeneNameR)',(wholeLenR)',(lmaxS)',(GeneNameS)',(wholeLenS)',(lmaxL)',(GeneNameL)',(wholeLenL)','VariableNames',{'ArgMaxLength' 'ForArgGeneName' 'ForArgGeneWholeLength' 'SerMaxLength' 'ForSerGeneName' 'ForSerGeneWholeLength' 'LeuMaxLength' 'LeuGeneName' 'ForLeuGeneWholeLength'},'RowNames',(SequenceName)');

writetable(tblMaxLen6,'tblMaxLen6.txt','WriteRowNames','true');

%%%find how many subsequences length > 300
%for i=1:length(lenArrR)
% countR(i)=length(find(lenArrR{i}>300));
% end
% 
% for i=1:length(lenArrS)
% countS(i)=length(find(lenArrS{i}>300));
% end
% 
% for i=1:length(lenArrL)
% countL(i)=length(find(lenArrL{i}>300));
% end
% countR
% countS
% countL