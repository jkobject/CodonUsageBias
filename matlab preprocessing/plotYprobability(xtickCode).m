function [] = plotYprobability(Y) %% plot one Y
    
figure

plot(Y,'O')

xlabel('amino acid');

ylabel('probability');

xlim([0,22]);

set(gca,'xTick',1:1:21);  

xtxt={'G','D','E','V','A','R','S','K','N','M','I','T','W','Z','C','Y','L','F','Q','H','P'};

set(gca,'xTickLabel',xtxt); %% set xTickLabel for each amino acid

grid on

shg %% make the current figure visible

end