function []= histY1(n,q)  %% n: weights length; q: sequence length vector

for i=1:length(q)
    
l(:,1:100000)=q(i); 

    arrayfun(@(x) entropyD(n,x),l); %% for a certain length calculate its entropy
    
    subplot(length(q), 1, i);

    hist(arrayfun(@(x) entropyD(n,x),l)); %% for a certain length draw its histogram
     
    xlim ([-0.1 0]);
 
     xlabel('entropy');
     ylabel('frequency');
     legend(sprintf('length = %.d',q(i)));

end
end

