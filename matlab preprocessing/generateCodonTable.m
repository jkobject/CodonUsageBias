% generate codon usage table for different species

fileID=fopen('sah6FF.fa')
strSah6=fscanf(fileID,'%c');
fclose(fileID);
s='ATGATCACTCTGCCCCGTATGATCACTCTGCCCCGT'
CodonFreq = codonbias(s,'GeneticCode',1)