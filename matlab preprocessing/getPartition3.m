function [] = getPartition3(Syno) %% Syno: synonymous codons number; Sublength: length of sequence

for Sublength=[462 459 429 443 441 431 451 453 447]  %% initialize start length and stop length

P=zeros(1,Syno);
P(1,:)=1/Syno;

filename3=['partition3o',num2str(Sublength),'.txt'];
 for i=0:Sublength    %%i, j, k represents amounts of 3 synonymous codons
    for j=0:Sublength    
        k=Sublength-i-j;
        if k>=0
        mnvect=[i,j,k];    %% mnvect: all the possible partitions (order matters here) 
        p=mnpdf(mnvect,P); %% for each combination, calculate multinomial distribution density
        fileID = fopen(filename3,'a');
        fprintf(fileID,'%d,',p);   %% choose right data type can save space and memory
        fclose(fileID);
        end
    end
 end
disp([num2str(Sublength),' p3done']);
end
end
