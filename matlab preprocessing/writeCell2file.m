fileID=fopen('speciesName.csv','w');
fprintf(fileID,'%s\n',speciesName{1:end,1});
fclose(fileID);