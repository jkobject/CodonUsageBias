function Homologies = getHomologyList(filename)

fileID = fopen(filename);
tt = textscan(fileID,'%s','delimiter',',');
fclose(fileID);

Homologies = tt{1};

end
% % % % % % % % 5000 homolist process%%%%
% % % fID=fopen('homoList5818Raw.txt');
% % % tt = textscan(fID,'%s','Delimiter','\n');
% % % fclose(fID);
% % % 
% % % homoList5818Raw=tt{1,1};
% % % 
% % % randID=randi(5818,[5818,1]);
% % % 
% % % save('homoList5818randID.mat','randID');
% % % 
% % % homoList5818Rand=homoList5818Raw(randID);
% % % 
% % % fileID=fopen('homoList5818Rand.txt','w');
% % % fprintf(fileID,'%s,',homoList5818Rand{:});
% % % fclose(fileID);
