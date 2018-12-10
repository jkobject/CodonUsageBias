%% this method doesn't work well

strut=webread('http://rest.ensemblgenomes.org/homology/id/YGL123W?content-type=application/json;compara=fungi;sequence=cdna;type=orthologues;format=condensed');

sourceSpecies=strut.data.id; %% yeast species
sourceSq=getgenbank(sourceSpecies,'SequenceOnly',1);
fileID=fopen('homologyGene.csv','a');
fprintf(fileID,'%s %s\n','>',sourceSpecies);
fprintf(fileID,'%s\n',sourceSq);
fclose(fileID);

for i=1:141    
    tempSpecies=strut.data.homologies(i).target.species; %% retrieve homology species
    try
        geneID=strut.data.homologies(i).target.id;
        tempSq=getgenbank(geneID,'SequenceOnly',1);
        fileID=fopen('homologyGene.csv','a');
        fprintf(fileID,'%s %s %s %s\n','>',tempSpecies,'>',geneID);
        fprintf(fileID,'%s\n',tempSq);
        fclose(fileID);
    catch
        %         fileID=fopen('homologyGene.csv','a');
        %         fprintf(fildID,'%s\n','>no relevant gene found');
        %         fclose(fileID);
    end
end