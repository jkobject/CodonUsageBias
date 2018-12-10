%%check the partition codes by comparison to binomial distribution

Sublength=100;
Syno=6;

P=zeros(1,Syno);
P(1,:)=1/Syno;
pmax=Efor(Syno,Sublength);
p=zeros(1,Sublength^(Syno-2)); %% guarantee enough storing space 
n=1;
i=0;
j=0;
k=0;
r=0;
s=0;
t=0;


%for i=0:Sublength    %%i, j, k represents amounts of 3 synonymous codons
for i=0:Sublength    %%i, j, k represents amounts of 3 synonymous codons
    for j=0:Sublength  
        for k=0:Sublength
            for r=0:Sublength
                for s=0:Sublength
                    t=Sublength-i-j-k-r-s;
                    if t>=0
                    mnvect=[i,j,k,r,s,t];    %% mnvect: all the possible partitions (order matters here) 
                    p(1,n)=mnpdf(mnvect,P); %% for each combination, calculate multinomial distribution density
                    n=n+1;                  %% maybe 0 in vector p, but doesn't matter since will delet later.
                    end
                end
            end
        end
    end
end

