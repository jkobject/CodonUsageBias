load('KL20pre.mat');
Y1=XXp;
load('KL20.mat')
Y2=XXp;
AAlist={'Asp','Asn','Cys','Gln','Glu','His','Lys','Phe','Tyr','Ile','Ala','Gly','Pro','Thr','Val','Arg','Ser','Leu'};


for i=1:18
figure
plot(Y1(:,i),Y2(:,i),'.');
xlabel('Selection Pressure by Previous Method');
ylabel('Selection Pressure by KL divergence');
title(['amino acid: ',AAlist{i}]);
end