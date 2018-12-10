% function [Svalue,freq,interestL,aveFreq] = getSFforT(XX,aa,Pbound)  %%%% get x y datasets for temperature fit
function [Svalue,freq,aveFreq] = getSFforT(XX,aa,interestL)

% lowP=Pbound(1);
% upP=Pbound(2);

aaID=find(ismember(XX{1},aa));  %%%extract the same aa
Sublength1=XX{2}(aaID,1);
Pi1=XX{3}(aaID,1);
DegPi1=XX{4}(aaID,1);
Pmax1=XX{5}(aaID,1);
% DegPmax1=XX{6}(aaID,1);
Frequency1=XX{7}(aaID,1);

%%%%histogram(Sublength,'BinWidth',1);
%%%%xlabel('length');
%%%%ylabel('fequency');
%%%%titile(['length distribution of ',aaList{aaCount},' in ',speciesName{tP}]);

% [N1,edges1]=histcounts(Sublength1,'BinWidth',1,'Normalization','cdf');
% [N2,edges2]=histcounts(Sublength1,'BinWidth',1);
% [~,interestID]=max(N2(find(N1>lowP,1):find(N1>upP,1)));%%%%choose lengths whose accumulative percentile: 60%-90%, then among them find the most largest one
% edges2f=edges2((find(N1>lowP)+1):(find(N1>upP,1)+1));
% interestL=edges2f(interestID(1));%%%% find which sublength has the optimal sample size to fit

SubLid=find(Sublength1==interestL);%%%%% certain sublength ID
Pi2=Pi1(SubLid,1);
DegPi2=DegPi1(SubLid,1);
Pmax2=Pmax1(SubLid,1);
% DegPmax2=DegPmax1(SubLid,1);
Frequency2=Frequency1(SubLid,1);

[Pi3,uniID,~]=unique(Pi2);
DegPi3=DegPi2(uniID);
Pmax3=Pmax2(uniID);
Frequency3=Frequency2(uniID);
% DegPmax3=DegPmax2(uniID);
Svalue=log(Pmax3./(Pi3.*DegPi3));


freq=log(Frequency3/sum(Frequency3));
aveFreq=sum(Frequency3)/length(Svalue);

end