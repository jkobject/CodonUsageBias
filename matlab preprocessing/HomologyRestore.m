filename='HomologySpeciesYOR202W.csv';
fileID = fopen(filename);
tt = textscan(fileID,'%s','delimiter',',');
fclose(fileID);

xx = tt{1};

save 'HomologySpeciesYOR202W.mat' xx


filename='HomologyGeneYOR202W.csv';
fileID = fopen(filename);
tt = textscan(fileID,'%s','delimiter',',');
fclose(fileID);

xx = tt{1};

save 'HomologyGeneYOR202W.mat' xx
