function [] = getPartition4(Syno,sublength) %% Syno: synonymous codons number; Sublength: length of sequence

for Sublength=301:sublength

P=zeros(1,Syno);
P(1,:)=1/Syno;

filename4=['partition4o',num2str(Sublength),'p.txt'];

for i=0:Sublength    %%i, j, k represents amounts of 4 synonymous codons
    for j=0:Sublength  
        for k=0:Sublength
        r=Sublength-i-j-k;
        if r>=0
        mnvect=[i,j,k,r];%% mnvect: all the possible partitions (order matters here) 
        p=mnpdf(mnvect,P); %% for each combination, calculate multinomial distribution density        
        fileID = fopen(filename4,'a');
        fprintf(fileID,'%d,',p);   %% choose right data type can save space and memory
        fclose(fileID);

        end
        end
    end
end
disp([num2str(Sublength),'done']);
end
end

