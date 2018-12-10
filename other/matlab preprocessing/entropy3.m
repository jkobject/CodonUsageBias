clc
clear


LL =[10 15 20 25 100];

for j=1:5
    L = LL(j);
    for i=1:100000
    D=rand(1,L);
    vec=zeros(1,L);
    vec(D<0.2) = 0.2;
    vec(D>0.2 & D< 0.31) = 0.31-0.2;
    vec(D>0.31 & D< 0.98) = 0.98-0.31;
    vec( D> 0.98) = 1-0.98;
    S(i)=sum(vec.*log(vec));
    end
    subplot(5, 1, j)
    hist(S)
    title(sprintf('length = %.d',L));

end







