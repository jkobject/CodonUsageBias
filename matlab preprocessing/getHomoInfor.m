function [homoSpecies,homoGene,l] = getHomoInfor(stringName)

weblink=['http://rest.ensemblgenomes.org/homology/id/',stringName,'?compara=fungi&content-type=application/json&sequence=cdna&type=orthologues&format=condensed'];
option=weboptions('Timeout',120);
try
strut=webread(weblink,option);

l=length(strut.data.homologies);
homoSpecies=cell(l,1);
homoGene=cell(l,1);

for fdID=1:l    %%%fdID:filed id
    homoSpecies{fdID}=strut.data.homologies(fdID).species;
    homoGene{fdID}=strut.data.homologies(fdID).id;
end

catch         %%% avoid computing pausing 
    fid=fopen('webreadFail.txt','a');  %%%%store webread failure species for further recalculate
    fprintf(fid,'%s,',stringName); %%%webreadFail.txt store homo read failure items   
    fclose(fid);
    homoSpecies={};  %%%%if webread fail return empty value.otherwise main body error
    homoGene={};
    l=0;
    
end

end