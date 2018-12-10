function [] = PlotPmax(Syno,sublength) 
pmax=zeros(1,94);
for Sublength=106:sublength
pmax(Sublength)=Efor(Syno,Sublength);
end
x=106:200;

plot(log(pmax(106:200)),'o')
xlabel('sequence length');
ylabel('log of maximum pi')
title('study about maximum probability for 4 synonymous codons of different length')

end



% figure;   %% plot average entropy
% 
% x=(106:200);
% 
% plot(x,AveEntropy4(106:200),'o')
% 
% title('multinomial distribution density analysis (Syno=4)');
% xlabel('sequence length');
% ylabel('sum(pi.*log(pi/pmax))');



% figure    %% plot sum(logp) or plot sum(p*logp)
% j=1;
% for i=106:200
% % pt(j)=sum(pf4f106t200{i}.*log(pf4f106t200{i}));
% %pt(j)=sum(log(pf4f106t200{i}));
% j=j+1;
% end
% 
% x=(106:200);
% 
% plot(x,pt,'o')
% 
% x=(100:200);
% title('multinomial distribution density analysis (Syno=4)');
% xlabel('sequence length');
% ylabel('sum(log(pi))');