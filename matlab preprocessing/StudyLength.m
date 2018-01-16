  %% Glybe into length influence on results


setSynonymousCodonTable;

i=1; 
for i=1:length(pasteCodon1)
[NN1(i,1),X1{i,1},Y1(i,1)]=GlyAminoAcidH(pasteCodon1{1,i});
end

Y1(isnan(Y1))=0;
NN1(isnan(NN1))=0;
NNp1=[Y1 NN1];
NNs1=sortrows(NNp1);
NNrp1=NNs1(:,2);
NNr1= NNrp1(end:-1:1);
NNc1p=cumsum(NNr1);
NNc1=NNc1p./NNc1p(length(NNc1p));

figure
plot(NNc1,'o')
xlim([0,5500]);
xlabel('indices');
ylabel('cummulative length');
title('sah6: Gly Cummulative Subsequence Length Distribution');
 
% pasteCodon1=getCodonSequence('sah6.csv'); %% call function getCodonSequence to creat input 

i=1; 
for i=1:length(pasteCodon2)
[NN2(i,1),X2{i,1},Y2(i,1)]=GlyAminoAcidH(pasteCodon2{1,i});
end

Y2(isnan(Y2))=0;
NN2(isnan(NN2))=0;
NNp2=[Y2 NN2];
NNs2=sortrows(NNp2);
NNrp2=NNs2(:,2);
NNr2= NNrp2(end:-1:1);
NNc2p=cumsum(NNr2);
NNc2=NNc2p./NNc2p(length(NNc2p));

figure
plot(NNc2,'o');
xlim([0,5500]);
xlabel('indices');
ylabel('cummulative length');
title('saeu: Gly Cummulative Subsequence Length Distribution');

i=1; 
for i=1:length(pasteCodon3)
[NN3(i,1),X3{i,1},Y3(i,1)]=GlyAminoAcidH(pasteCodon3{1,i});
end

Y3(isnan(Y3))=0;
NN3(isnan(NN3))=0;
NNp3=[Y3 NN3];
NNs3=sortrows(NNp3);
NNrp3=NNs3(:,2);
NNr3= NNrp3(end:-1:1);
NNc3p=cumsum(NNr3);
NNc3=NNc3p./NNc3p(length(NNc3p));

figure
plot(NNc3,'o');
xlim([0,5500]);
xlabel('indices');
ylabel('cummulative length');
title('saku: Gly Cummulative Subsequence Length Distribution');

i=1; 
for i=1:length(pasteCodon4)
[NN4(i,1),X4{i,1},Y4(i,1)]=GlyAminoAcidH(pasteCodon4{1,i});
end

Y4(isnan(Y4))=0;
NN4(isnan(NN4))=0;
NNp4=[Y4 NN4];
NNs4=sortrows(NNp4);
NNrp4=NNs4(:,2);
NNr4= NNrp4(end:-1:1);
NNc4p=cumsum(NNr4);
NNc4=NNc4p./NNc4p(length(NNc4p));

figure
plot(NNc4,'o');
xlim([0,5500]);
xlabel('indices');
ylabel('cummulative length');
title('sp: Gly Cummulative Subsequence Length Distribution');

i=1; 
for i=1:length(pasteCodon5)
[NN5(i,1),X5{i,1},Y5(i,1)]=GlyAminoAcidH(pasteCodon5{1,i});
end

Y5(isnan(Y5))=0;
NN5(isnan(NN5))=0;
NNp5=[Y5 NN5];
NNs5=sortrows(NNp5);
NNrp5=NNs5(:,2);
NNr5= NNrp5(end:-1:1);
NNc5p=cumsum(NNr5);
NNc5=NNc5p./NNc5p(length(NNc5p));

figure
plot(NNc5,'o');
xlim([0,5500]);
xlabel('indices');
ylabel('cummulative length');
title('so: Gly Cummulative Subsequence Length Distribution');

for i=1:length(pasteCodon6)
[NN6(i,1),X6{i,1},Y6(i,1)]=GlyAminoAcidH(pasteCodon6{1,i});
end
% Xg{i,:},

% tbl=table(D,NNd,'VariableNames',{'RatioGlyue' 'length'});
% Tbl = sortrows(tbl,{'RatioGlyue'},'descend');
% sort(NNd,'Ascend');
% figure
% 
% h=histogram(NN);
% 
% xlabel('subsequence length');
%  
% ylabel('frequency');
% 
% title('sah6:Gly Histogram');
% 
% [Nbin,edges]=histcounts(NN);
% 
% k=1;
% mid=zeros(1,(length(edges)-1)); %%get (mid) midpoints vector
% for k=1:(length(edges)-1)
% %      if Nbin(k)==0    %%discard outliers
% %         mid(k)=NaN;
% %         Nbin=NaN;
% %     else
%         mid(k)=(edges(k)+edges(k+1))/2;
% %     end
% end

Y6(isnan(Y6))=0;
NN6(isnan(NN6))=0;
NNp6=[Y6 NN6];
NNs6=sortrows(NNp6);
NNrp6=NNs6(:,2);
NNr6= NNrp6(end:-1:1);
NNc6p=cumsum(NNr6);
NNc6=NNc6p./NNc6p(length(NNc6p));

figure
plot(NNc6,'o');
xlim([0,5500]);
xlabel('indices');
ylabel('cummulative length');
title('sj: Gly Cummulative Subsequence Length Distribution');
% [R,P]=corrcoef(mid,Nbin)
