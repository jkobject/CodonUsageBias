%%%%%%%syno3over400  should use the common initialization as 'syno4over700.m'
lenFileName={'3sequenceLength11.mat','3sequenceLength12.mat','3sequenceLength13.mat','3sequenceLength14.mat','3sequenceLength15.mat',...
    '3sequenceLength21.mat','3sequenceLength22.mat','3sequenceLength23.mat',...
    '3sequenceLength31.mat','3sequenceLength32.mat','3sequenceLength33.mat','3sequenceLength34.mat','3sequenceLength35.mat','3sequenceLength36.mat','3sequenceLength37.mat',...
    '3sequenceLength41.mat','3sequenceLength42.mat','3sequenceLength43.mat','3sequenceLength44.mat','3sequenceLength45.mat'};

lenArrI=cell(1,20);
for r=1:20
m=matfile(lenFileName{r});
v=who(m);
vi=v{1};
lenArrI{r}=m.(vi);
end

[lIdI,lI]=getLabove(lenArrI,400); %% get interested sublength index and length value. lIdA lA are 'cell' 

fid=fopen('Iover400.txt','a');
fprintf(fid,'%s\n','format:index,subsequence length,gene name, whole length.');
fclose(fid);    

for j=1:20
    
sting1=[num2str(j),'th species begeins for Ile over 400.'];
fid=fopen('Iover400.txt','a');
fprintf(fid,'%s\n',sting1);
fclose(fid);    

lIdtemp=lIdI{j};  %%lIdtemp is double matrix

if ~isempty(lIdtemp)

pasteCodonTemp=PasteCodon{j};
sequenceTempName=sequenceName{j};
ltemp=lI{j};  %%ltemp is double matrix

for i=1:length(lIdtemp)

lIdTemp=lIdtemp(i);
lTemp=ltemp(i);
GeneName=sequenceTempName{lIdTemp};
wholeLength=length(pasteCodonTemp{lIdTemp});

fid=fopen('Iover400.txt','a');
fprintf(fid,'%5u,%4u,%s,%5u\n',lIdTemp,lTemp,GeneName,wholeLength);
fclose(fid);
 
end
clear lIdtemp ltemp sequenceTempName pasteCodonTemp
end
end