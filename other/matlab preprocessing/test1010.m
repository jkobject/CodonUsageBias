setSynonymousCodonTable;
order={'11','12','13','14','15','21','22','23','31','32','33','34','35','36','37','41','42','43','44','45'};

for o=1:1

fileName1=(['pasteCodon',order{o},'.mat']);
m1=matfile(fileName1);
v1p=who(m1);
v1=v1p{1};
pasteCodon=m1.(v1);


fileName6=(['6sequenceLength',order{o},'.mat']);
m6=matfile(fileName6);
v6p=who(m6);
v61=v6p{1};
v62=v6p{2};
v63=v6p{3};
NNl=m6.(v61);
NNr=m6.(v62);
NNs=m6.(v63);

%%Arg%%%%%
%%%%%%%%%%%%%theoretical reference calculation
SumHist=zeros(1,1);
for i=1:length(NNr)
    NN=NNr(i);%
    if NN<=400
        pref=getPref(NN,6); %% remember to change syno value;  read all the Entropy value for possible configurations--pref
        SumHist=sumHist(SumHist,pref,6,NN);%% SumHist resore sum value
    end
end
[SumHistF,edgesF,ThoeryPercent]=getTheoryPlot(SumHist);

%%%%%%%%%%%%real entropy distribution calculation
for i=1:length(pasteCodon) %% note: pasteCodon is column vector
    [NNrs,~,Rp]=ArgAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA
    if NNrs<=400
        R(i)=-log(Rp)/NNrs;
    else
        R(i)=NaN;
    end
end
[NHISTreal,EdgesReal]=getRealPlot(R);

figure
plot(edgesF,SumHistF,'r.');
hold on
plot(EdgesReal,NHISTreal,'g.');
plot(ThoeryPercent,0,'b*');
hold off
title('distribution of miexed lengths for Arg in Sah6');
xlabel('entropy value');
ylabel('probability');
legend('theoretical','real','99%percentil');
end