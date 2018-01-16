setSynonymousCodonTable

global ctG cfG
 
for i=1:length(pasteCodon1) %% note: pasteCodon is column vector
    
disp([num2str(i), 'codon sequence calculation begins']);
    
[NNg,~,G(i,:)]=GlyAminoAcidH(pasteCodon1{1,i}); %%get Ygg for each mRNA

if isnan(NNg) 
    continue;
else   
Go(i,:)=log(G(i,:))/NNg;
Gr(i,:)=AveEntropy4(NNg)/NNg;

end 
end

SelWeakId=find((Go-Gr)>0);

figure    
    
plot(Go,Gr,'.');

xlim([-1.4,0]);

hold on
 
x=-0.25:0.1:0;

y=x;

plot(x,y);

hold off

set(gca,'xaxislocation','bottom','yaxislocation','left','xdir','reverse','ydir','reverse')

xlabel('per codon entropy of original sequence');
 
ylabel('per codon average entropy of sequence after equal codon replacement');

title(['Species--sah6: ','Amino Acid--Gly(4 synonymous)']);


