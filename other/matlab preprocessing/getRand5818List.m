%%% caucious to run this file since it will generate different random gene
%%% list for homology analysis

%%%%%%%generate random list of 5818 genes for homology analysis

randID=randperm(5818);
save('homo5818randID.mat','randID');
fID=fopen('homoList5818Raw.txt');
tt=textscan(fID,'%s','Delimiter',',');
fclose(fID);
homoList5818=tt{1,1};
homo5818Rand=homoList5818(randID);

fID2=fopen('homoList5818Rand.txt','w');
fprintf(fID2,'%s,',homo5818Rand{:});
fclose(fID2);