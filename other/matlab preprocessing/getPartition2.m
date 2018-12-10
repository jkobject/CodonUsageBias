function [] = getPartition2(Syno,sublength) %% Syno: synonymous codons number; Sublength: length of sequence

parfor Sublength=401:sublength  %% initialize start length and stop length

P=zeros(1,Syno);
P(1,:)=1/Syno;

filename2=['partition2o',num2str(Sublength),'.txt'];
    for i=0:Sublength        
        j=Sublength-i;
        mnvect=[i,j];    %% mnvect: all the possible partitions (order matters here) 
        p=mnpdf(mnvect,P); %% for each combination, calculate multinomial distribution density
        fileID = fopen(filename2,'a');
        fprintf(fileID,'%i,%i;%d;',mnvect,p);
        fclose(fileID);  
        end
    end
disp([num2str(Sublength),'done']);
end
end
 