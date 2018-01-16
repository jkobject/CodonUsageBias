function [] = getPartition6(Syno,sublength) %% Syno: synonymous codons number; Sublength: length of sequence

for Sublength=10:sublength  %% initialize start length and stop length

P=zeros(1,Syno);
P(1,:)=1/Syno;

filename6=['partition6o',num2str(Sublength),'ptest.txt'];
fileID = fopen(filename6,'w');

for i=0:Sublength    %%i, j, k represents amounts of 6 synonymous codons
    for j=0:Sublength  
        for k=0:Sublength
            for r=0:Sublength
                for s=0:Sublength
                    t=Sublength-i-j-k-r-s;
                    if t>=0
                    mnvect=[i,j,k,r,s,t];    %% mnvect: all the possible partitions (order matters here) 
                    p=mnpdf(mnvect,P);                     
                    fprintf(fileID,'%d,',p);                  
                    end
                end
            end
        end
    end
end
fclose(fileID); 
disp([num2str(Sublength),' p6done']);
end
end


