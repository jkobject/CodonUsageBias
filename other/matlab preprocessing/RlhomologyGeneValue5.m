

%%%%replace random sequence and compare with originial 
setSynonymousCodonTable;
global ctD ctI ctG ctR

homoList={'YJL200C','YGR204W','YAL012W','YIL116W','YML027W','YDR201W','YDR293C'};   %%homology sourse
%%%%%here, consider Asp, Ile, Gly,Arg

for homoc=1:7 %%homology counter
    
    stringName=homoList{homoc};
    
    
    fileID5=fopen([stringName,'RlhomologyAsp.txt'],'a'); %%%%%Rl: means Replaced
    fprintf(fileID5,'%s,%s,%s,%s,%s\n','species','gene','entropyValue','entropyLocation','length');
    
    fileID10=fopen([stringName,'RlhomologyIle.txt'],'a');
    fprintf(fileID10,'%s,%s,%s,%s,%s\n','species','gene','entropyValue','entropyLocation','length');
    
    fileID11=fopen([stringName,'RlhomologyGly.txt'],'a');
    fprintf(fileID11,'%s,%s,%s,%s,%s\n','species','gene','entropyValue','entropyLocation','length');
    
    fileID16=fopen([stringName,'RlhomologyArg.txt'],'a');
    fprintf(fileID16,'%s,%s,%s,%s,%s\n','species','gene','entropyValue','entropyLocation','length');
    
    
    fileID19=fopen([stringName,'RlhomologyError.txt'],'a');
    
    [homoSpecies,homoGene,l] = getHomoInfor(stringName);%%%call function 'getHomoInfor' to retrieve homology information
    
    for j=1:l
        try
            speciesName=homoSpecies{j};
            fileName1=sprintf('%s',[speciesName,'.fa']);%% find corresponding species file
            fileName2=sprintf('%s',[speciesName,'Gene.fa']);%% only gene names in such file
            geneNamesWhole=getHomologyList(fileName2); %% list all the gene names in the interesed species
            geneName=homoGene{j};
            
            id=find(contains(geneNamesWhole,geneName));
            pasteCodonWhole=getCodonSequence(fileName1); %%%%% find interested homology
            pasteCodon=pasteCodonWhole{id};
            
            %%%%%Asp%%%%%%
            [NNd,~,~]=AspAminoAcidH(pasteCodon);
            pasteD = replaceE(NNd,ctD);  %%%%%%replace with random sequence: pasteD
            [~,~,Drp]=AspAminoAcidH(pasteD);
            Dr=-log(Drp)/NNd;                  %%%%%%Dr: replaced D value
            if NNd<=800
                Plocation = entropy4Location(2,NNd,Dr);
            else
                Plocation=NaN;
            end
            
            fprintf(fileID5,'%s,%s,%f,%f,%i\n',speciesName,geneName,Dr,Plocation,NNd);
            
            %%%%%Ile%%%%%%
            [NNi,~,~]=IleAminoAcidH(pasteCodon);
            pasteI=replaceE(NNi,ctI);
            [~,~,Irp]=IleAminoAcidH(pasteI);
            Ir=-log(Irp)/NNi;
            if NNi<=400
                Plocation = entropy4Location(3,NNi,Ir);
            else
                Plocation=NaN;
            end
            
            fprintf(fileID10,'%s,%s,%f,%f,%i\n',speciesName,geneName,Ir,Plocation,NNi);
            
            %%%Gly%%%%%%%
            [NNg,~,~]=GlyAminoAcidH(pasteCodon);
            pasteG=replaceE(NNg,ctG);
            [~,~,Grp]=GlyAminoAcidH(pasteG);
            Gr=-log(Grp)/NNg;
            if NNg<=300
                Plocation = entropy4Location(4,NNg,Gr);
            else
                Plocation=NaN;
            end
            
            fprintf(fileID11,'%s,%s,%f,%f,%i\n',speciesName,geneName,Gr,Plocation,NNg);

            %%Arg%%%%%%%
            [NNr,~,~]=ArgAminoAcidH(pasteCodon);
            pasteR=replaceE(NNr,ctR);
            [~,~,Rrp]=ArgAminoAcidH(pasteR);
            Rr=-log(Rrp)/NNr;
            if NNr<=400
                Plocation = entropy4Location(6,NNr,Rr);
            else
                Plocation=NaN;
            end
            
            fprintf(fileID16,'%s,%s,%f,%f,%i\n',speciesName,geneName,Rr,Plocation,NNr);

            
        catch
            fprintf(fileID19,'%s,%s,%s,%s\n',[speciesName,":",geneName,':error!']);
        end
        clear pasteCodon
    end
    
    fclose(fileID5);
    fclose(fileID10);
    fclose(fileID11);
    fclose(fileID16);
    fclose(fileID19);
    
end

