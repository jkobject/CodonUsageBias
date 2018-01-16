
for l=75:100

pout=sampEntropyMean6(6,l,10000000);

filename=['Samp6o',num2str(l),'p.mat'];

save(filename,'pout')

disp([num2str(l),'done']);

end


