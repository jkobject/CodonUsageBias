%% fro 6 synonymous codons: compare true value of length 200 and sample values
pmax=Efor(6,200);
Pmax=log(pmax);

sp21=sp2001{1,200};
sp22=sp2002{1,200};
sp23=sp2003{1,200};
sp24=sp2004{1,200};
sp25=sp2005{1,200};
sp26=sp2006{1,200};

sp200=[sp21(100000:2000000),sp22(100000:2000000),sp23(100000:3000000),sp24(100000:1000000),sp25(100000:2000000),sp26(100000:3000000)];
save 'SampF200.mat' sp200; %% sp200 original mnpdf log value

PMAX=zeros(1,12400006);
PMAX(1:12400006)=Pmax;

SP200=sp200-PMAX; %% for each partition, normalization by maximum value;

save 'SampF200f.mat' SP200 %% SP200 normalized by maximum pmax

ave1=mean(SP200(10000:110000))
ave2=mean(SP200(20000:220000))
ave3=mean(SP200(30000:330000))
ave4=mean(SP200(40000:4040000))  %%400,0000 is enough
ave5=mean(SP200(50000:5050000))
ave6=mean(SP200(60000:6060000))
ave7=mean(SP200)
compV=sum(SP200)/length(SP200);
true=-2.4573
x=[100000,200000,300000,4000000,5000000,6000000,12000000];
y=[ave1,ave2,ave3,ave4,ave5,ave6,ave7];
plot(x,y);
xlabel('sample size');
ylabel('sample entropy');


% hold on
% x=0;y=-2.4573;
% plot(x,y,'o')
% hold off