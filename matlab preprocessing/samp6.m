function [] = samp6(syno,NN,itr)

P=zeros(1,syno);
P(:,:)=1/syno;
R=mnrnd(NN,P,itr);


for i=1:itr
    mnp(i)=mnpdf(R(i,:),P);
end

fileName=['partition6o',num2str(NN),'p.csv'];
csvwrite(fileName,mnp);

end

