% %% for 4 synonymous codons find subsequence length above 700 which need sampling method to obtain expected entropy
% 
% % sequenceName={sequenceName11,sequenceName12,sequenceName13,sequenceName14,sequenceName15,...
% %     sequenceName21,sequenceName22,sequenceName23,...
% %     sequenceName31,sequenceName32,sequenceName33,sequenceName34,sequenceName35,sequenceName36,sequenceName37...
% %     sequenceName41,sequenceName42,sequenceName43,sequenceName44,sequenceName45};
% 
% % PasteCodon={pasteCodon11,pasteCodon12,pasteCodon13,pasteCodon14,pasteCodon15,...
% %     pasteCodon21,pasteCodon22,pasteCodon23,...
% %     pasteCodon31,pasteCodon32,pasteCodon33,pasteCodon34,pasteCodon35,pasteCodon36,pasteCodon37,...
% %     pasteCodon41,pasteCodon42,pasteCodon43,pasteCodon44,pasteCodon45};
% 
% 
% % SequenceName={'sah6','saeu','saku','Ag','Yl','sp','so','sj','ac','afl','ans','anr','ao','at','afu','ff','fg','fo','fv','tr'};
% 
% %%lenFileName: contain all the sublength for 4 syno among 20 species
% 
% lenFileName={'4sequenceLength11.mat','4sequenceLength12.mat','4sequenceLength13.mat','4sequenceLength14.mat','4sequenceLength15.mat',...
%     '4sequenceLength21.mat','4sequenceLength22.mat','4sequenceLength23.mat',...
%     '4sequenceLength31.mat','4sequenceLength32.mat','4sequenceLength33.mat','4sequenceLength34.mat','4sequenceLength35.mat','4sequenceLength36.mat','4sequenceLength37.mat',...
%     '4sequenceLength41.mat','4sequenceLength42.mat','4sequenceLength43.mat','4sequenceLength44.mat','4sequenceLength45.mat'};
% % 
% % %%pasFileName: all the input codon sequences among 20 species
% pasFileName={'pasteCodon11.mat','pasteCodon12.mat','pasteCodon13.mat','pasteCodon14.mat','pasteCodon15.mat',...
%     'pasteCodon21.mat','pasteCodon22.mat','pasteCodon23.mat',...
%     'pasteCodon31.mat','pasteCodon32.mat','pasteCodon33.mat','pasteCodon34.mat','pasteCodon35.mat','pasteCodon36.mat','pasteCodon37.mat',...
%     'pasteCodon41.mat','pasteCodon42.mat','pasteCodon43.mat','pasteCodon44.mat','pasteCodon45.mat'};
% % 
% % %%seqNaFileName: all the gene names for each input codon sequence
% seqNaFileName={'sequenceName11.mat','sequenceName12.mat','sequenceName13.mat','sequenceName14.mat','sequenceName15.mat',...
%     'sequenceName21.mat','sequenceName22.mat','sequenceName23.mat',...
%     'sequenceName31.mat','sequenceName32.mat','sequenceName33.mat','sequenceName34.mat','sequenceName35.mat','sequenceName36.mat','sequenceName37.mat'...
%     'sequenceName41.mat','sequenceName42.mat','sequenceName43.mat','sequenceName44.mat','sequenceName45.mat'};
% % 
% % %% for memory usage benefit, generate PasteCodon Cell contain all the input codon sequences
% PasteCodon=cell(1,20);
% % 
% for r=1:20
% 
% m=matfile(pasFileName{r});
% 
% v=who(m);
% 
% vn=v{1};
%     
% PasteCodon{r}=m.(vn);
% 
% end
% % 
% % %% for memory usage benefit, generate sequenceName cell contain all the gene names for input codon sequences
% % 
% 
% sequenceName=cell(1,20);
% for r=1:20
% 
% m=matfile(seqNaFileName{r});
% 
% v=who(m);
% 
% vn=v{1};
%     
% sequenceName{r}=m.(vn);
% 
% end
% 
% 
% % %%retieve all the length for 4 syno
% lenArrA=cell(1,20);
% lenArrE=cell(1,20);
% lenArrG=cell(1,20);
% lenArrP=cell(1,20);
% lenArrT=cell(1,20);
% lenArrV=cell(1,20);
% for r=1:20
% 
% m=matfile(lenFileName{r});
% 
% v=who(m);
% 
% va=v{1};
% ve=v{2};
% vg=v{3};
% vp=v{4};
% vt=v{5};
% vv=v{6};
% 
% lenArrA{r}=m.(va);
% lenArrE{r}=m.(ve);
% lenArrG{r}=m.(vg);
% lenArrP{r}=m.(vp);
% lenArrT{r}=m.(vt);
% lenArrV{r}=m.(vv);
% end


[lIdA,lA]=getLabove(lenArrA,300); %% get interested sublength index and length value. lIdA lA are 'cell' 

fid=fopen('Aover300.txt','a');
fprintf(fid,'%s\n','format:index,subsequence length,gene name, whole length.');
fclose(fid);    

for j=1:20
    
sting1=[num2str(j),'th species begeins for Ala over 300.'];
fid=fopen('Aover300.txt','a');
fprintf(fid,'%s\n',sting1);
fclose(fid);    

lIdtemp=lIdA{j};  %%lIdtemp is double matrix

if ~isempty(lIdtemp)

pasteCodonTemp=PasteCodon{j};
sequenceTempName=sequenceName{j};
ltemp=lA{j};  %%ltemp is double matrix

for i=1:length(lIdtemp)

lIdTemp=lIdtemp(i);
lTemp=ltemp(i);
GeneName=sequenceTempName{lIdTemp};
wholeLength=length(pasteCodonTemp{lIdTemp});

fid=fopen('Aover300.txt','a');
fprintf(fid,'%5u,%4u,%s,%5u\n',lIdTemp,lTemp,GeneName,wholeLength);
fclose(fid);
 
end
clear lIdtemp ltemp sequenceTempName pasteCodonTemp
end
end


[lIdE,lE]=getLabove(lenArrE,300);

fid=fopen('Eover300.txt','a');
fprintf(fid,'%s\n','format:index,subsequence length,gene name, whole length.');
fclose(fid);

for j=1:20
    
sting2=[num2str(j),'th species begeins for Glu over 300.'];
fid=fopen('Eover300.txt','a');
fprintf(fid,'%s\n',sting2);
fclose(fid);    

lIdtemp=lIdE{j};

if ~isempty(lIdtemp)

pasteCodonTemp=PasteCodon{j};%% change
sequenceTempName=sequenceName{j};
ltemp=lE{j};

for i=1:length(lIdtemp)

lIdTemp=lIdtemp(i);
lTemp=ltemp(i);
GeneName=sequenceTempName{lIdTemp};
wholeLength=length(pasteCodonTemp{lIdTemp});

fid=fopen('Eover300.txt','a');
fprintf(fid,'%5u,%4u,%s,%5u\n',lIdTemp,lTemp,GeneName,wholeLength);
fclose(fid);
 
end
clear lIdtemp ltemp sequenceTempName pasteCodonTemp
end
end


[lIdG,lG]=getLabove(lenArrG,300);

fid=fopen('Gover300.txt','a');
fprintf(fid,'%s\n','format:index,subsequence length,gene name, whole length.');
fclose(fid);

for j=1:20
    
sting3=[num2str(j),'th species begeins for Gly over 300.'];
fid=fopen('Gover300.txt','a');
fprintf(fid,'%s\n',sting3);
fclose(fid);    

lIdtemp=lIdG{j};

if ~isempty(lIdtemp)

pasteCodonTemp=PasteCodon{j};%% change
sequenceTempName=sequenceName{j};
ltemp=lG{j};

for i=1:length(lIdtemp)

lIdTemp=lIdtemp(i);
lTemp=ltemp(i);
GeneName=sequenceTempName{lIdTemp};
wholeLength=length(pasteCodonTemp{lIdTemp});

fid=fopen('Gover300.txt','a');
fprintf(fid,'%5u,%4u,%s,%5u\n',lIdTemp,lTemp,GeneName,wholeLength);
fclose(fid);
 
end
clear lIdtemp ltemp sequenceTempName pasteCodonTemp
end
end



[lIdP,lP]=getLabove(lenArrP,300);

fid=fopen('Pover300.txt','a');
fprintf(fid,'%s\n','format:index,subsequence length,gene name, whole length.');
fclose(fid);

for j=1:20
    
sting3=[num2str(j),'th species begeins for Pro over 300.'];
fid=fopen('Pover300.txt','a');
fprintf(fid,'%s\n',sting3);
fclose(fid);    

lIdtemp=lIdP{j};

if ~isempty(lIdtemp)

pasteCodonTemp=PasteCodon{j};%% change
sequenceTempName=sequenceName{j};
ltemp=lP{j};

for i=1:length(lIdtemp)

lIdTemp=lIdtemp(i);
lTemp=ltemp(i);
GeneName=sequenceTempName{lIdTemp};
wholeLength=length(pasteCodonTemp{lIdTemp});

fid=fopen('Pover300.txt','a');
fprintf(fid,'%5u,%4u,%s,%5u\n',lIdTemp,lTemp,GeneName,wholeLength);
fclose(fid);
 
end
clear lIdtemp ltemp sequenceTempName pasteCodonTemp
end
end


[lIdT,lT]=getLabove(lenArrT,300);

fid=fopen('Tover300.txt','a');
fprintf(fid,'%s\n','format:index,subsequence length,gene name, whole length.');
fclose(fid);

for j=1:20
    
sting3=[num2str(j),'th species begeins for Thr over 300.'];
fid=fopen('Tover300.txt','a');
fprintf(fid,'%s\n',sting3);
fclose(fid);    

lIdtemp=lIdT{j};

if ~isempty(lIdtemp)

pasteCodonTemp=PasteCodon{j};%% change
sequenceTempName=sequenceName{j};
ltemp=lT{j};

for i=1:length(lIdtemp)

lIdTemp=lIdtemp(i);
lTemp=ltemp(i);
GeneName=sequenceTempName{lIdTemp};
wholeLength=length(pasteCodonTemp{lIdTemp});

fid=fopen('Tover300.txt','a');
fprintf(fid,'%5u,%4u,%s,%5u\n',lIdTemp,lTemp,GeneName,wholeLength);
fclose(fid);
 
end
clear lIdtemp ltemp sequenceTempName pasteCodonTemp
end
end



[lIdV,lV]=getLabove(lenArrV,300);

fid=fopen('Vover300.txt','a');
fprintf(fid,'%s\n','format:index,subsequence length,gene name, whole length.');
fclose(fid);

for j=1:20
    
sting3=[num2str(j),'th species begeins for Val over 300.'];
fid=fopen('Vover300.txt','a');
fprintf(fid,'%s\n',sting3);
fclose(fid);    

lIdtemp=lIdV{j};

if ~isempty(lIdtemp)

pasteCodonTemp=PasteCodon{j};%% change
sequenceTempName=sequenceName{j};
ltemp=lV{j};

for i=1:length(lIdtemp)

lIdTemp=lIdtemp(i);
lTemp=ltemp(i);
GeneName=sequenceTempName{lIdTemp};
wholeLength=length(pasteCodonTemp{lIdTemp});

fid=fopen('Vover300.txt','a');
fprintf(fid,'%5u,%4u,%s,%5u\n',lIdTemp,lTemp,GeneName,wholeLength);
fclose(fid);
 
end
clear lIdtemp ltemp sequenceTempName pasteCodonTemp
end
end

