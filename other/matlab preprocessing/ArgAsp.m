% % % plots for Asp Arg in saeu, afl:
fileName1='AspComparison6.mat';
m1=matfile(fileName1);
D1=m1.D;
Dr1=m1.Dr;


fileName2='AspComparison10.mat';
m2=matfile(fileName2);
D2=m2.D;
Dr2=m2.Dr;

fileName3='ArgComparison2.mat';
m3=matfile(fileName3);
R1=m3.R;
Rr1=m3.Rr;


fileName4='ArgComparison10.mat';
m4=matfile(fileName4);
R2=m4.R;
Rr2=m4.Rr;

% % % Asp for Sp
Ddis1=(D1(~isnan(D1)))./(Dr1(~isnan(Dr1)));
n1bins=floor(length(D1)/15);

figure
histogram(Ddis1,n1bins);
xlabel('ratio between real entropy and expected entropy');
ylabel('frequency');
title('Distribution: Schizzosaccharomyces pombe--Asp(2sysnonymous)');

%%Asp in afl
Ddis2=(D2(~isnan(D2)))./(Dr2(~isnan(Dr2)));
n2bins=floor(length(D2)/15);

figure 
histogram(Ddis2,n2bins);
xlabel('ratio between real entropy and expected entropy');
ylabel('frequency');
title('Distribution: Aspergillus flavus--Asp(2sysnonymous)');

%%histogram Arg in saeu
Rdis1=(R1(~isnan(R1)))./(Rr1(~isnan(Rr1)));
n3bins=floor(length(R1)/15);

figure
histogram(Rdis1,n3bins);
xlabel('ratio between real entropy and expected entropy');
ylabel('frequency');
title('Distribution: Saccharomyces eubayanus--Arg(6sysnonymous)');

%% Arg in afl
Rdis2=(R2(~isnan(R2)))./(Rr2(~isnan(Rr2)));
n4bins=floor(length(R2)/15);

figure 
histogram(Rdis2,n4bins);
xlim([0 25]);
xlabel('ratio between real entropy and expected entropy');
ylabel('frequency');
title('Distribution: Aspergillus flavus--Arg(6sysnonymous)');
