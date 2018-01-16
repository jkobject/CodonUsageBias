%% get 461 SepeciesName

fileID = fopen('SpeciesName.csv');
tt = textscan(fileID,'%s','delimiter','\n');
fclose(fileID);

SpeciesName = tt{1};

save('SpeciesName.mat','SpeciesName')



