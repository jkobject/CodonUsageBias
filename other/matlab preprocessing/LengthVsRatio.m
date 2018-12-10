%%subsequence length Vs raito, to analyze length influence
x14=MidFq1{1,4}{:,1};
y14=MidFq1{1,4}{:,2};

% N=cell(1,140);
for i=1:140
sqind{i}=find(Y<=exp(0.5-i/2)&Y>exp(-i/2))
N{1,i}=NN(sqind{i});
Navg(i)=round(sum(N{1,i})/length(N{1,i}));
end

plot(Navg,'o')

xlabel('corresponding bin');
ylabel('subsequence length');

xlim([0,70]);

txt={'0','-10','-20','-30','-40','-50','-60','-70'};

set(gca,'xTickLabel',txt); 

title('Val, Sah6,subsequence average length distribution against bin');



% sqind1=find(Y<=exp(0)&Y>exp(-2));
% N{1,1}=NN(sqind1);
% Navg1=round(sum(N{1,1})/length(N{1,1}));
% 
% sqind2=find(Y<=exp(-2)&Y>exp(-4));
% N{1,2}=NN(sqind2);
% Navg2=round(sum(N{1,2})/length(N{1,2}));
% 
% sqind3=find(Y<=exp(-4)&Y>exp(-6));
% N{1,3}=NN(sqind3);
% Navg3=round(sum(N{1,3})/length(N{1,3}));
% 
% sqind4=find(Y<=exp(-6)&Y>exp(-8));
% N{1,4}=NN(sqind4);
% Navg4=round(sum(N{1,4})/length(N{1,4}));
% 
% sqind5=find(Y<=exp(-8)&Y>exp(-10));
% N{1,5}=NN(sqind5);
% Navg5=round(sum(N{1,5})/length(N{1,5}));
% 
% sqind6=find(Y<=exp(-10)&Y>exp(-30));
% N{1,6}=NN(sqind6);
% Navg6=round(sum(N{1,6})/length(N{1,6}));
% 
% sqind7=find(Y<=exp(-30)&Y>exp(-50));
% N{1,7}=NN(sqind7);
% Navg7=round(sum(N{1,7})/length(N{1,7}));
% 
% 
% sqind8=find(Y<=exp(-50)&Y>exp(-70));
% N{1,8}=NN(sqind8);
% Navg8=round(sum(N{1,8})/length(N{1,8}));

% sqind9=find(Y<=exp(-70)&Y>exp(-100));
% N{1,9}=NN(sqind9);

% for i=1:8
% figure
% Xinterval={(-2:0),'(-4,-2]','(-6,-4]','(-8,-6]','(-10,-8]','(-30,-10]','(-50,-30]','(-70,-50]'};
% plot(N{1,i},'o');
% hold on
% end
% hold off
% 
% % xlabel('indices');
% % ylabel('subsequence length');
% % lgd=legend('(-2,0]','(-4,-2]','(-6,-4]','(-8,-6]','(-10,-8]','(-30,-10]','(-50,-30]','(-70,-50]');
% % title('Val, Sah6,subsequence length distribution corresponding to log(pr/pm) value shown as legend');
% 
% 
figure
plot(x14,y14,'o');
% xlim([-10,0]);
xlabel('midpoint of bin');
ylabel('log of frequency');
title('Val, Sah6,(Pr/Pm) probability ratio distribution before replacement');
