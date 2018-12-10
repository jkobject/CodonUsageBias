%%read from a file to creat input codon sequences

function SequenceName = getSequenceName(filename)

%%read the sequence from file 'inputSequence'

fileID = fopen(filename);
tt = textscan(fileID,'%s','delimiter',';');
fclose(fileID);


%%seperate the sequence name & sequence itself
xx = tt{1};
l = length(xx);
SequenceName = xx(1:2:l,:);

% save 'SequenceName6.mat' SequenceName6
end