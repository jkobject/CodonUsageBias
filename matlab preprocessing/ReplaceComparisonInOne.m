%% species sah6 pr/pm probability ratio comparison

NameAminoAcid={'Gly','Asp','Glu','Val','Ala','Arg','Ser','Lys','Asn','Ile','Thr','Cys','Tyr','Leu','Phe','Gln','His','Pro'};

for i=1:18
    figure
    plot(xb{1,i},yb{1,i},'o');
    hold on
    plot(xa{1,i},ya{1,i},'*');
    hold on
    plot(xaT{1,i},yaT{1,i},'^');
   

xlabel('midpoint of bin');
  
ylabel('log of frequency');

legend('before','afterEqualReplacement','afterTableReplacement','Location','Northwest');

title([NameAminoAcid(i),'sj log(Pr/Pm) distribution comparison before&after replacement']);
     hold off
end