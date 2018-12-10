setSynonymousCodonTable

global ctR cfR
 
for i=1:length(pasteCodon1) %% note: pasteCodon is column vector
    
disp([num2str(i), 'codon sequence calculation begins']);
    
[NNr,~,R(i,:)]=ArgAminoAcidH(pasteCodon1{1,i}); %%get Ygg for each mRNA

if isnan(NNr) 
    continue;
else   
Ro(i,:)=log(R(i,:))/NNr;
Rr(i,:)=AveEntropy6f(NNr)/NNr;

end 
end

SelWeakId=find((Ro-Rr)>0);

figure    
    
plot(Ro,Rr,'.');

hold on
 
x=(-0.3):0.1:0;

y=x;

plot(x,y);

hold off

set(gca,'xaxislocation','bottom','yaxislocation','left','xdir','reverse','ydir','reverse')

xlabel('per codon entropy of original sequence');
 
ylabel('per codon average entropy of sequence after equal codon replacement');

title(['Species--sah6: ','Amino Acid--Arg(6 synonymous)']);

