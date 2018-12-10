clc 
clear


pasteCodon=getCodonSequence; %% call function getCodonSequence to create input 

sequenceName={'YAL042W','YAL039C','YAL044C'};

for i=1:length(pasteCodon)
    
Y=AminoAcidNumberSubsequence(pasteCodon{1,i});

X(i,:)=Y;

end

plot(X','o');%%box plot 

xlabel('all mRNAs');

ylabel('probability');

xlim([0,22]);

%%ylim([-180,0]);

set(gca,'xTick',1:1:21);  

xtxt={'G','D','E','V','A','R','S','K','N','M','I','T','W','Z','C','Y','L','F','Q','H','P'};

set(gca,'xTickLabel',xtxt); %% set xTickLabel for each amino acid

shg %% make the current figure visible

grid off



