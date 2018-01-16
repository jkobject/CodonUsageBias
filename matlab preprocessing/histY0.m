clc
clear


l(:,1:100000)=25;

arrayfun(@entropyY0,l);

hist(arrayfun(@entropyY0,l));

title(sprintf('length = %.d',25));

xlabel('entropy');

legend('frequency');

