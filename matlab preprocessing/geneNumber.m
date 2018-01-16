function NNg= geneNumber(codonSequence)  %%calculate gene number in a particular sequence

disp('A new gene number calculation begins');

codonString=strjoin(codonSequence); %% convert CELL array to a normal string 

geneAindice=find(codonString,'A');
geneAnumber=length(geneAindice)

geneTindice=find(codonString,'T');
geneTnumber=length(geneTindice)

geneCindice=find(codonString,'C');
geneCnumber=length(geneCindice)

geneGindice=find(codonString,'G');
geneGnumber=length(geneGindice)

disp('The current gene number calculation ends');

end

