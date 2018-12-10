%% hypothesis test

u0=-2.4573;
samp1=SP200(4000000:10000000);
samp2=SP200(4100000:10100000);
samp3=SP200(5000000:11000000);
samp4=SP200(10:12000010);
samp={samp1,samp2,samp3,samp4};

for i=1:length(samp)
n(i)=length(samp{i});
smean(i)=mean(samp{i});
svar(i)=var(samp{i});
t(i)=(smean(i)-u0)/(sqrt(svar(i))/sqrt(n(i)));
end

h=ttest(samp4)


