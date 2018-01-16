%% function of this file: save average entropy for 6 synonymous codons for equal replacement

% N=1000000;  %% data process from mnpdf to log(pi/pmax) 
% for i=101:200
% currentSamp=p11Sample6{1,i};
% pmax=Efor(6,i);
% PMAX=zeros(1,(N-1));
% PMAX(1:(N-1))=log(pmax);
% sampN=currentSamp(2:N);
% AveEntropy6f101t200(i)=sum(sampN-PMAX)/(N-1);
% end
% save 'AveEntropy6f101t200.mat' AveEntropy6f101t200

% N=1000000;  %% data process from mnpdf to log(pi/pmax)
% for i=201:300
% currentSamp=p22Sample6{1,i};
% pmax=Efor(6,i);
% PMAX=zeros(1,(N-1));
% PMAX(1:(N-1))=log(pmax);
% sampN=currentSamp(2:N);
% AveEntropy6f201t300(i)=sum(sampN-PMAX)/(N-1);
% end
% save 'AveEntropy6f201t300.mat' AveEntropy6f201t300

% N=1000000;  %% data process from mnpdf to log(pi/pmax)
% for i=301:400
% currentSamp=p33Sample6{1,i};
% pmax=Efor(6,i);
% PMAX=zeros(1,(N-1));
% PMAX(1:(N-1))=log(pmax);
% sampN=currentSamp(2:N);
% AveEntropy6f301t400(i)=sum(sampN-PMAX)/(N-1);
% end
% save 'AveEntropy6f301t400.mat' AveEntropy6f301t400

% N=1000000;  %% data process from mnpdf to log(pi/pmax)
% for i=401:500
% currentSamp=p44Sample6{1,i};
% pmax=Efor(6,i);
% PMAX=zeros(1,(N-1));
% PMAX(1:(N-1))=log(pmax);
% sampN=currentSamp(2:N);
% AveEntropy6f401t500(i)=sum(sampN-PMAX)/(N-1);
% end
% save 'AveEntropy6f401t500.mat' AveEntropy6f401t500

% N=1000000;  %% data process from mnpdf to log(pi/pmax)
% for i=501:600
% currentSamp=p55Sample6{1,i};
% pmax=Efor(6,i);
% PMAX=zeros(1,(N-1));
% PMAX(1:(N-1))=log(pmax);
% sampN=currentSamp(2:N);
% AveEntropy6f501t600(i)=sum(sampN-PMAX)/(N-1);
% end
% save 'AveEntropy6f501t600.mat' AveEntropy6f501t600

AveEntropy6f=zeros(1,600);
AveEntropy6f(1:55)=AveEntropy6(1:55);
AveEntropy6f(56:100)=AveEntropy6f56t100(56:100);
AveEntropy6f(101:200)=AveEntropy6f101t200(101:200);
AveEntropy6f(201:300)=AveEntropy6f201t300(201:300);
AveEntropy6f(301:400)=AveEntropy6f301t400(301:400);
AveEntropy6f(401:500)=AveEntropy6f401t500(401:500);
AveEntropy6f(501:600)=AveEntropy6f501t600(501:600);

save 'AveEntropy6f.mat' AveEntropy6f

