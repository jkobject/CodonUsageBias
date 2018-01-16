calculate sample mean

filename1 = 'p1Sample6o500.mat';
m1 = matfile(filename1);
sp61=m1.p1Sample6;


filename2 = 'p2Sample6o500.mat';
m2 = matfile(filename2);
sp62=m2.p2Sample6;

filename3 = 'p3Sample6o500.mat';
m3 = matfile(filename3);
sp63=m3.p3Sample6;



filename4 = 'p4Sample6o500.mat';
m4 = matfile(filename4);
sp64=m4.p4Sample6;


filename5 = 'p5Sample6o500.mat';
m5 = matfile(filename5);
sp65=m5.p5Sample6;

filename6 = 'p6Sample6o500.mat';
m6 = matfile(filename6);
sp66=m6.p6Sample6;

sp67=p44Sample6{500};
length(sp67)


mean1=sum(sp61(100:1000000))/(1000000-100);
mean2=sum(sp62(100:1000000))/(1000000-100);
mean3=sum(sp63(100:1000000))/(1000000-100);
mean4=sum(sp64(100:1000000))/(1000000-100);
mean5=sum(sp65(100:1000000))/(1000000-100);
mean6=sum(sp66(100:1000000))/(1000000-100);
mean7=sum(sp67(100:1000000))/(1000000-100);

A=[mean1 mean2 mean3 mean4 mean5 mean6 mean7]

Am=sum(A)/6

Ea=std(A)
