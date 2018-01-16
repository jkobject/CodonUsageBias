function Evalue = Efor(cLeng,subSleng) %% for synonymous codons equally likely show up, to assign occurrence as Domonique suggested

P=zeros(1,cLeng);  %%cLeng: number of synonymous codons for one amino acid

P(:,:)=1/cLeng;

i=mod(subSleng,cLeng);  %%subSleng: the length of subsequence

X=zeros(1,cLeng);


for j=1:(i+1)
    X(j)=ceil(subSleng/cLeng);
end


for k=j:cLeng
    X(k)=floor(subSleng/cLeng);
end

Evalue=mnpdf(X,P); %% calculate maximum probability of such equally likely case
end