% %%%%don't forget to divide length
syno=6;
NN=75;
itr=1000;

% P=zeros(1,syno);
% P(:,:)=1/syno;
% R=mnrnd(NN,P,itr);
% for i=1:itr
%     mnp(i)=mnpdf(R(i,:),P);
% end
% figure
% plot(mnp,'.');

% figure
% histogram(mnp);
% % 
pEx=dlmread('partition6o75.txt');
filename=['partition6o',num2str(NN),'p.csv'];
pSa=csvread(filename);
% % % figure
% % % histogram(pf1);
pmax=Efor(6,75);
s1=sum(log(1./pEx))/length(pEx);
s2=sum(pSa.*log(1./pSa));
% s3=sum(mnp.*log(1./mnp));


syno=6;
NN=50;
itr=100000;
pmax=Efor(syno,NN);

P=zeros(1,syno);
P(:,:)=1/syno;
R=mnrnd(NN,P,itr);

filename=['partition6o',num2str(NN),'p.csv'];
pf=csvread(filename);
[NHIST1,edges1]=histcounts(pf,'BinWidth',0.00000001);
x1=log(1./(edges1(1:(end-1))));
y1=edges1(1:(end-1)).*(NHIST1/sum(NHIST1));
figure
plot(x1,y1,'.');
title('exhaustive method');
xlabel('log(pmax/pi)');
ylabel('probability');

for i=1:itr
    mnp(i)=mnpdf(R(i,:),P);
end

[NHIST2,edges2]=histcounts(mnp,'BinWidth',0.00000001);
x2=log(1./edges2(1:(end-1)));
y2=edges2(1:(end-1)).*(NHIST2/sum(NHIST2));
figure
plot(x2,y2,'.');
title('sampling method (size 100000)');
xlabel('log(pmax/pi)');
ylabel('probability');

y1done=y1(2:end);
y2done=y2(2:end);
sumId=find((y1done.*y2done)~=0);
KL=sum(y2done(sumId).*log(y2done(sumId)./y1done(sumId)))