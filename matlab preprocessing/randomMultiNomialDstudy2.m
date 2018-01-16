% nl=[4,8,100,200,400];
% for i=1:length(nl)
% [midG,NgbinG]=randomMultiNomialDstudy(4,nl(i),10000)
% 
% % figure
% plot(midG,NgbinG,'o'); 
% xlabel('midpoint of bin');
% ylabel('frequency');
%     
% hold on
% end
% 
% legend('length 4','length 8','leng 100','leng 400','Location','NorthWest');
% title('Occurrences N=10000 for each length, category m=4');


[midG,NgbinG]=randomMultiNomialDstudy(4,1,1000)

plot(midG,NgbinG,'o'); 
xlabel('midpoint of bin');
ylabel('frequency');
hold on
title('Occurrences N=100000 for lengthNl=2, category m=4');
