%% for 6 synonymous codons find subsequence length above 200 which need sampling method to obtain expected entropy

sequenceName={sequenceName11,sequenceName12,sequenceName13,sequenceName14,sequenceName15,...
    sequenceName21,sequenceName22,sequenceName23,...
    sequenceName31,sequenceName32,sequenceName33,sequenceName34,sequenceName35,sequenceName36,sequenceName37...
    sequenceName41,sequenceName42,sequenceName43,sequenceName44,sequenceName45};

% PasteCodon={pasteCodon11,pasteCodon12,pasteCodon13,pasteCodon14,pasteCodon15,...
%     pasteCodon21,pasteCodon22,pasteCodon23,...
%     pasteCodon31,pasteCodon32,pasteCodon33,pasteCodon34,pasteCodon35,pasteCodon36,pasteCodon37,...
%     pasteCodon41,pasteCodon42,pasteCodon43,pasteCodon44,pasteCodon45};


% SequenceName={'sah6','saeu','saku','Ag','Yl','sp','so','sj','ac','afl','ans','anr','ao','at','afu','ff','fg','fo','fv','tr'};



lenArrR={NNr11,NNr12,NNr13,NNr14,NNr15,NNr21,NNr22,NNr23,NNr31,NNr32,NNr33,NNr34,NNr35,NNr36,NNr37,NNr41,NNr42,NNr43,NNr44,NNr45};
[lIdR,lR]=getLabove(lenArrR,300); %% get interested sublength index and length value. lIdR lR are 'cell' 

fid=fopen('Rover300.txt','a');
fprintf(fid,'%s\n','format:index,subsequence length,gene name, whole length.');
fclose(fid);    

for j=1:20
    
sting1=[num2str(j),'th species begeins for Arg over 300.'];
fid=fopen('Rover300.txt','a');
fprintf(fid,'%s\n',sting1);
fclose(fid);    

lIdtemp=lIdR{j};  %%lIdtemp is double matrix

if ~isempty(lIdtemp)

pasteCodonTemp=PasteCodon{j};
sequenceTempName=sequenceName{j};
ltemp=lR{j};  %%ltemp is double matrix

for i=1:length(lIdtemp)

lIdTemp=lIdtemp(i);
lTemp=ltemp(i);
GeneName=sequenceTempName{lIdTemp};
wholeLength=length(pasteCodonTemp{lIdTemp});

fid=fopen('Rover300.txt','a');
fprintf(fid,'%5u,%4u,%s,%5u\n',lIdTemp,lTemp,GeneName,wholeLength);
fclose(fid);
 
end
clear lIdtemp ltemp sequenceTempName pasteCodonTemp
end
end


lenArrS={NNs11,NNs12,NNs13,NNs14,NNs15,NNs21,NNs22,NNs23,NNs31,NNs32,NNs33,NNs34,NNs35,NNs36,NNs37,NNs41,NNs42,NNs43,NNs44,NNs45};
[lIdS,lS]=getLabove(lenArrS,300);

fid=fopen('Sover300.txt','a');
fprintf(fid,'%s\n','format:index,subsequence length,gene name, whole length.');
fclose(fid);

for j=1:20
    
sting2=[num2str(j),'th species begeins for Ser over 300.'];
fid=fopen('Sover300.txt','a');
fprintf(fid,'%s\n',sting2);
fclose(fid);    

lIdtemp=lIdS{j};

if ~isempty(lIdtemp)

pasteCodonTemp=PasteCodon{j};%% change
sequenceTempName=sequenceName{j};
ltemp=lS{j};

for i=1:length(lIdtemp)

lIdTemp=lIdtemp(i);
lTemp=ltemp(i);
GeneName=sequenceTempName{lIdTemp};
wholeLength=length(pasteCodonTemp{lIdTemp});

fid=fopen('Sover300.txt','a');
fprintf(fid,'%5u,%4u,%s,%5u\n',lIdTemp,lTemp,GeneName,wholeLength);
fclose(fid);
 
end
clear lIdtemp ltemp sequenceTempName pasteCodonTemp
end
end


lenArrL={NNl11,NNl12,NNl13,NNl14,NNl15,NNl21,NNl22,NNl23,NNl31,NNl32,NNl33,NNl34,NNl35,NNl36,NNl37,NNl41,NNl42,NNl43,NNl44,NNl45};
[lIdL,lL]=getLabove(lenArrL,300);

fid=fopen('Lover300.txt','a');
fprintf(fid,'%s\n','format:index,subsequence length,gene name, whole length.');
fclose(fid);

for j=1:20
    
sting3=[num2str(j),'th species begeins for Leu over 300.'];
fid=fopen('Lover300.txt','a');
fprintf(fid,'%s\n',sting3);
fclose(fid);    

lIdtemp=lIdL{j};

if ~isempty(lIdtemp)

pasteCodonTemp=PasteCodon{j};%% change
sequenceTempName=sequenceName{j};
ltemp=lL{j};

for i=1:length(lIdtemp)

lIdTemp=lIdtemp(i);
lTemp=ltemp(i);
GeneName=sequenceTempName{lIdTemp};
wholeLength=length(pasteCodonTemp{lIdTemp});

fid=fopen('Lover300.txt','a');
fprintf(fid,'%5u,%4u,%s,%5u\n',lIdTemp,lTemp,GeneName,wholeLength);
fclose(fid);
 
end
clear lIdtemp ltemp sequenceTempName pasteCodonTemp
end
end