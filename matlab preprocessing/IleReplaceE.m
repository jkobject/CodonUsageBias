setSynonymousCodonTable

global ctI cfI
 
for i=1:length(pasteCodon1) %% note: pasteCodon is column vector
    
disp([num2str(i), 'codon sequence calculation begins']);
    
[NNi,~,I(i,:)]=IleAminoAcidH(pasteCodon1{1,i}); %%get Ygg for each mRNA

if isnan(NNi) 
    continue;
else   

Ir(i,:)=AveEntropy3(NNi);

end 
end

figure    
    
plot(log(I),Ir,'.');

% hold on
% 
% x=(-40):0;
% 
% y=x;
% 
% plot(x,y);
% 
% hold off

set(gca,'xaxislocation','bottom','yaxislocation','left','xdir','reverse','ydir','reverse')

xlabel('entropy of original sequence');
 
ylabel('entropy of sequence after equal codon replacement');

title(['Species--sah6: ','Amino Acid--Ile']);
% 
% hold on;
% x=(-40:0); 
% y=x;
% plot(x,y);
% hold off;    