function [] = EntropyMean(Syno,sublength) %% Syno: synonymous codons number; Sublength: length of sequence

for Sublength=31:sublength  %% initialize start length and stop length

P=zeros(1,Syno);
P(1,:)=1/Syno;
pmax=Efor(Syno,Sublength);
p=zeros(Sublength^(Syno-2),1); %% guarantee enough storing space,matlab has limit for variable size
n=1;
i=0;
j=0;
k=0;
r=0;
s=0;
t=0;

% filename6=sprintf('Partition6F%u.txt',Sublength);%%save partitions
% filename6p=sprintf('pFpartition6%u.txt',Sublength); %%save mnpdf value for each partition
switch Syno
%     case 2
%     for i=0:Sublength        
%         j=Sublength-i;
%         if j>=0
%         mnvect=[i,j];    %% mnvect: all the possible partitions (order matters here) 
%         p(1,n)=mnpdf(mnvect,P); %% for each combination, calculate multinomial distribution density
%         n=n+1;                  %% maybe 0 in vector p, but doesn't matter since will delet later.
%         end
%     end

        
    
%     case 3
% for i=0:Sublength    %%i, j, k represents amounts of 3 synonymous codons
%     for j=0:Sublength    
%         k=Sublength-i-j;
%         if k>=0
%         mnvect=[i,j,k];    %% mnvect: all the possible partitions (order matters here) 
%         p(1,n)=mnpdf(mnvect,P); %% for each combination, calculate multinomial distribution density
%         n=n+1;                  %% maybe 0 in vector p, but doesn't matter since will delet later.
%         end
%     end
% end

case 4
for i=0:Sublength    %%i, j, k represents amounts of 4 synonymous codons
    for j=0:Sublength  
        for k=0:Sublength
        r=Sublength-i-j-k;
        if r>=0
        mnvect=[i,j,k,r];%% mnvect: all the possible partitions (order matters here) 
        p(1,n)=mnpdf(mnvect,P); %% for each combination, calculate multinomial distribution density
        n=n+1;                  %% maybe 0 in vector p, but doesn't matter since will delet later.
        end
        end
    end
end

% case 6
% for i=0:Sublength    %%i, j, k represents amounts of 6 synonymous codons
%     for j=0:Sublength  
%         for k=0:Sublength
%             for r=0:Sublength
%                 for s=0:Sublength
%                     t=Sublength-i-j-k-r-s;
%                     if t>=0
%                     mnvect=[i,j,k,r,s,t];    %% mnvect: all the possible partitions (order matters here) 
%                     p=mnpdf(mnvect,P);                     
%                     dlmwrite(filename6,mnvect,'-append','precision','%3u','delimiter',' ');
%                     dlmwrite(filename6p,p,'-append') %% for each combination, calculate multinomial distribution density
%                     n=n+1;                  %% maybe 0 in vector p, but doesn't matter since will delet later.
%                     end
%                 end
%             end
%         end
%     end
% end

end

% SumOfP(i)=sum(pf);
% X=dlmread(filename6p);
% pf=X(X~=0);
pf=p(p~=0);
% AveEntropy=sum(pf.*(log(pf/pmax))); %%average entropy; 
save 'pf4o31.mat' pf
% fileID = fopen('AveEntropy3f1001t1500.txt','a');
% fprintf(fileID,'%f\n',AveEntropy);   %% choose right data type can save space and memory
% fclose(fileID);
% figure

% histogram(pf,100);
% xlabel('(mnpdf)p');
% ylabel('frequency');
% title(['histogram of length ',num2str(Sublength),'; Synonymous codon ',num2str(Syno)]);

end
end
% x=1:length(AveEntropyNoPm);plot(x,AveEntropyNoPm,'o')