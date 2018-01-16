setSynonymousCodonTable

global ctD cfD
 
for i=1:length(pasteCodon1) %% note: pasteCodon is column vector
    
disp([num2str(i), 'codon sequence calculation begins']);
    
[NNd,~,D(i,:)]=AspAminoAcidH(pasteCodon1{1,i}); %%get Ygg for each mRNA

if isnan(NNd) 
    continue;
else   

Dr(i,:)=AveEntropy2(NNd);

end 
end

figure    
    
plot(log(D),Dr,'.');

 hold on
 
x=(-35):0;

y=x;

plot(x,y);

hold off

set(gca,'xaxislocation','bottom','yaxislocation','left','xdir','reverse','ydir','reverse')

xlabel('entropy of original sequence');
 
ylabel('average entropy of sequence after equal codon replacement');

title(['Species--sah6: ','Amino Acid--Asp']);

% hold on;
% x=(-40:0); 
% y=x;
% plot(x,y);
% hold off;    