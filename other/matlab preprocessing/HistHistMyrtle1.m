setSynonymousCodonTable;
order={'11','12','13','14','15','21','22','23','31','32','33','34','35','36','37','41','42','43','44','45'};

for o=1:4

fileName1=(['pasteCodon',order{o},'.mat']);
m1=matfile(fileName1);
v1p=who(m1);
v1=v1p{1};
pasteCodon=m1.(v1);

fileName2=(['2sequenceLength',order{o},'.mat']);
m2=matfile(fileName2);
v2p=who(m2);
v21=v2p{1};
v22=v2p{2};
v23=v2p{3};
v24=v2p{4};
v25=v2p{5};
v26=v2p{6};
v27=v2p{7};
v28=v2p{8};
v29=v2p{9};
NNc=m2.(v21);
NNd=m2.(v22);
NNe=m2.(v23);
NNf=m2.(v24);
NNh=m2.(v25);
NNk=m2.(v26);
NNn=m2.(v27);
NNq=m2.(v28);
NNy=m2.(v29);

fileName3=(['3sequenceLength',order{o},'.mat']);
m3=matfile(fileName3);
v3p=who(m3);
v3=v3p{1};
NNi=m3.(v3);

fileName4=(['4sequenceLength',order{o},'.mat']);
m4=matfile(fileName4);
v4p=who(m4);
v41=v4p{1};
v43=v4p{3};
v44=v4p{4};
v45=v4p{5};
v46=v4p{6};
NNa=m4.(v41);
NNg=m4.(v43);
NNp=m4.(v44);
NNt=m4.(v45);
NNv=m4.(v46);

fileName6=(['6sequenceLength',order{o},'.mat']);
m6=matfile(fileName6);
v6p=who(m6);
v61=v6p{1};
v62=v6p{2};
v63=v6p{3};
NNl=m6.(v61);
NNr=m6.(v62);
NNs=m6.(v63);

%%% Gly %%%
%%%%%%%%%%%theoretical reference calculation%%%
SumHist=zeros(1,1);
for i=1:length(NNg)
    NN=NNg(i);
    if NN<=300
    pref=getPref(NN,4); %%read all the Entropy value for possible configurations--pref 
    SumHist=sumHist(SumHist,pref,4,NN);%% SumHist resore sum value
    end
end
[SumHistF,edgesF,ThoeryPercent]=getTheoryPlot(SumHist);
% 
% %%%%%%%%%%%real entropy distribution calculation
for i=1:length(pasteCodon) %% note: pasteCodon is column vector
[NNgs,~,Gp]=GlyAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA
if NNgs<=300
G(i)=-log(Gp)/NNgs;
else
G(i)=NaN;
end
end
[NHISTreal,EdgesReal]=getRealPlot(G);

sFile1=(['GlyAccurate',order{o},'.mat']);
save(sFile1,'SumHistF','edgesF','NHISTreal','EdgesReal');
%%%generate random sequences for comparison
% % % global ctG
% % % for i=1:length(pasteCodon) %% note: pasteCodon is column vector
% % % [NNg,~,~]=GlyAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA
% % % if ~(isnan(NNg)||NNg>300)
% % % InGd=replaceE(NNg,ctG);
% % % [~,~,Gp]=GlyAminoAcidH(InGd); %%get Ygg for each mRNA
% % % G(i)=-log(Gp)/NNg;
% % % else
% % % G(i)=NaN;
% % % end
% % % end
% % % [NHISTrealRandom,EdgesRealRandom]=getRealPlot(G);
% % % 
% % % plot(EdgesRealRandom,NHISTrealRandom,'y.');

%%%Val%%%%%
SumHist=zeros(1,1);
for i=1:length(NNv)
    NN=NNv(i);
    if NN<=300
    pref=getPref(NN,4); %%read all the Entropy value for possible configurations--pref 
    SumHist=sumHist(SumHist,pref,4,NN);%% SumHist resore sum value
    end
end
[SumHistF,edgesF,ThoeryPercent]=getTheoryPlot(SumHist);

%%%%%%%%%%%real entropy distribution calculation
for i=1:length(pasteCodon) %% note: pasteCodon is column vector
[NNvs,~,Vp]=ValAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA
if NNvs<=300
V(i)=-log(Vp)/NNvs;
else
V(i)=NaN;
end
end
[NHISTreal,EdgesReal]=getRealPlot(V);

sFile2=(['ValAccurate',order{o},'.mat']);
save(sFile2,'SumHistF','edgesF','NHISTreal','EdgesReal');

%%%Ala%%%%%%
SumHist=zeros(1,1);
for i=1:length(NNa)
    NN=NNa(i);
    if NN<=300
    pref=getPref(NN,4); %%read all the Entropy Alaue for possible configurations--pref 
    SumHist=sumHist(SumHist,pref,4,NN);%% SumHist resore sum Alaue
    end
end
[SumHistF,edgesF,ThoeryPercent]=getTheoryPlot(SumHist);

%%%%%%%%%%%real entropy distribution calculation
for i=1:length(pasteCodon) %% note: pasteCodon is column vector
[NNas,~,Ap]=AlaAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA
if NNas<=300
A(i)=-log(Ap)/NNas;
else
A(i)=NaN;
end
end
[NHISTreal,EdgesReal]=getRealPlot(A);

sFile3=(['AlaAccurate',order{o},'.mat']);
save(sFile3,'SumHistF','edgesF','NHISTreal','EdgesReal');

%%%%Thr%%%%%
SumHist=zeros(1,1);
for i=1:length(NNt)
    NN=NNt(i);
    if NN<=300
    pref=getPref(NN,4); %%read all the Entropy value for possible configurations--pref 
    SumHist=sumHist(SumHist,pref,4,NN);%% SumHist resore sum value
    end
end
[SumHistF,edgesF,ThoeryPercent]=getTheoryPlot(SumHist);

%%%%%%%%%%%real entropy distribution calculation
for i=1:length(pasteCodon) %% note: pasteCodon is column vector
[NNts,~,Tp]=ThrAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA
if NNts<=300
T(i)=-log(Tp)/NNts;
else
T(i)=NaN;
end
end
[NHISTreal,EdgesReal]=getRealPlot(T);

sFile4=(['ThrAccurate',order{o},'.mat']);
save(sFile4,'SumHistF','edgesF','NHISTreal','EdgesReal');

%%%Pro%%%
SumHist=zeros(1,1);
for i=1:length(NNp)
    NN=NNp(i);
    if NN<=300
    pref=getPref(NN,4); %%read all the Entropy value for possible configurations--pref 
    SumHist=sumHist(SumHist,pref,4,NN);%% SumHist resore sum value
    end
end
[SumHistF,edgesF,ThoeryPercent]=getTheoryPlot(SumHist);

%%%%%%%%%%%real entropy distribution calculation
for i=1:length(pasteCodon) %% note: pasteCodon is column vector
[NNps,~,Pp]=ProAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA
if NNps<=300
P(i)=-log(Pp)/NNps;
else
P(i)=NaN;
end
end
[NHISTreal,EdgesReal]=getRealPlot(P);

sFile5=(['ProAccurate',order{o},'.mat']);
save(sFile5,'SumHistF','edgesF','NHISTreal','EdgesReal');


%%%Ile%%%%%
%%%%%%%%%%%%theoretical reference calculation
SumHist=zeros(1,1);
for i=1:length(NNi)
    try
    NN=NNi(i);
    if NN<=500
    pref=getPref(NN,3); %% remember to change syno value;  read all the Entropy value for possible configurations--pref 
    SumHist=sumHist(SumHist,pref,3,NN);%% SumHist resore sum value
    end
    catch
        disp(['Ile',num2str(NN),'error']);
    end
end
[SumHistF,edgesF,ThoeryPercent]=getTheoryPlot(SumHist);

%%%%%%%real entropy distribution calculation
for i=1:length(pasteCodon) %% note: pasteCodon is column vector  
[NNis,~,Ip]=IleAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA
if NNis<=500
I(i)=-log(Ip)/NNis;
else
I(i)=NaN;
end
end
[NHISTreal,EdgesReal]=getRealPlot(I);

sFile6=(['IleAccurate',order{o},'.mat']);
save(sFile6,'SumHistF','edgesF','NHISTreal','EdgesReal');

%%%%Asp%%%%%
%%%%%%%%%%%%theoretical reference calculation
SumHist=zeros(1,1);
for i=1:length(NNd)
     NN=NNd(i);
     if NN<=800
     pref=getPref(NN,2); %% function 'getPfre', read all the Entropy value for possible configurations--pref 
     SumHist=sumHist(SumHist,pref,2,NN);%% function 'sumHist' resore sum value
     end
end
[SumHistF,edgesF,ThoeryPercent]=getTheoryPlot(SumHist); %%function 'getTheoryPlot' calculate data for plot

%%%%%%%%%%%real entropy distribution calculation
for i=1:length(pasteCodon) %% note: pasteCodon is column vector
[NNds,~,Dp]=AspAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA
if NNds<=800
D(i)=-log(Dp)/NNds;
else
D(i)=NaN;
end
end
[NHISTreal,EdgesReal]=getRealPlot(D);  %% function 'getRealPlot' calculate data for plot

sFile7=(['AspAccurate',order{o},'.mat']);
save(sFile7,'SumHistF','edgesF','NHISTreal','EdgesReal');

%%%%Lys%%%%%
SumHist=zeros(1,1);
for i=1:length(NNk)
     NN=NNk(i);
     if NN<=800
     pref=getPref(NN,2); %% function 'getPfre', read all the Entropy value for possible configurations--pref 
     SumHist=sumHist(SumHist,pref,2,NN);%% function 'sumHist' resore sum value
     end
end
[SumHistF,edgesF,ThoeryPercent]=getTheoryPlot(SumHist); %%function 'getTheoryPlot' calculate data for plot

%%%%%%%%%%real entropy distribution calculation
for i=1:length(pasteCodon) %% note: pasteCodon is column vector
[NNks,~,Kp]=LysAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA
if NNks<=800
K(i)=-log(Kp)/NNks;
else
K(i)=NaN;
end
end
[NHISTreal,EdgesReal]=getRealPlot(K);  %% function 'getRealPlot' calculate data for plot

sFile8=(['LysAccurate',order{o},'.mat']);
save(sFile8,'SumHistF','edgesF','NHISTreal','EdgesReal');

%%Asn%%%
SumHist=zeros(1,1);
for i=1:length(NNn)
     NN=NNn(i);
     if NN<=800
     pref=getPref(NN,2); %% function 'getPfre', read all the Entropy value for possible configurations--pref 
     SumHist=sumHist(SumHist,pref,2,NN);%% function 'sumHist' resore sum value
     end
end
[SumHistF,edgesF,ThoeryPercent]=getTheoryPlot(SumHist); %%function 'getTheoryPlot' calculate data for plot

%%%%%%%%%%real entropy distribution calculation
for i=1:length(pasteCodon) %% note: pasteCodon is column vector
[NNns,~,Np]=AsnAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA
if NNns<=800
N(i)=-log(Np)/NNns;
else
N(i)=NaN;
end
end
[NHISTreal,EdgesReal]=getRealPlot(N);  %% function 'getRealPlot' calculate data for plot

sFile9=(['AsnAccurate',order{o},'.mat']);
save(sFile9,'SumHistF','edgesF','NHISTreal','EdgesReal');

%%%Cys%%%%%
SumHist=zeros(1,1);
for i=1:length(NNc)
     NN=NNc(i);
     if NN<=800
     pref=getPref(NN,2); %% function 'getPfre', read all the Entropy value for possible configurations--pref 
     SumHist=sumHist(SumHist,pref,2,NN);%% function 'sumHist' resore sum value
     end
end
[SumHistF,edgesF,ThoeryPercent]=getTheoryPlot(SumHist); %%function 'getTheoryPlot' calculate data for plot

%%%%%%%%%%real entropy distribution calculation
for i=1:length(pasteCodon) %% note: pasteCodon is column vector
[NNcs,~,Cp]=CysAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA
if NNcs<=800
C(i)=-log(Cp)/NNcs;
else
C(i)=NaN;
end
end
[NHISTreal,EdgesReal]=getRealPlot(C);  %% function 'getRealPlot' calculate data for plot

sFile10=(['CysAccurate',order{o},'.mat']);
save(sFile10,'SumHistF','edgesF','NHISTreal','EdgesReal');

%%Tyr%%%
SumHist=zeros(1,1);
for i=1:length(NNy)
     NN=NNy(i);
     if NN<=800
     pref=getPref(NN,2); %% function 'getPfre', read all the Entropy value for possible configurations--pref 
     SumHist=sumHist(SumHist,pref,2,NN);%% function 'sumHist' resore sum value
     end
end
[SumHistF,edgesF,ThoeryPercent]=getTheoryPlot(SumHist); %%function 'getTheoryPlot' calculate data for plot

%%%%%%%%%%real entropy distribution calculation
for i=1:length(pasteCodon) %% note: pasteCodon is column vector
[NNys,~,Np]=TyrAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA
if NNys<=800
Y(i)=-log(Np)/NNys;
else
Y(i)=NaN;
end
end
[NHISTreal,EdgesReal]=getRealPlot(Y);  %% function 'getRealPlot' calculate data for plot

sFile11=(['TyrAccurate',order{o},'.mat']);
save(sFile11,'SumHistF','edgesF','NHISTreal','EdgesReal');

%%%Phe%%%%%
SumHist=zeros(1,1);
for i=1:length(NNf)
     NN=NNf(i);
     if NN<=800
     pref=getPref(NN,2); %% function 'getPfre', read all the Entropy value for possible configurations--pref 
     SumHist=sumHist(SumHist,pref,2,NN);%% function 'sumHist' resore sum value
     end
end
[SumHistF,edgesF,ThoeryPercent]=getTheoryPlot(SumHist); %%function 'getTheoryPlot' calculate data for plot

%%%%%%%%%%real entropy distribution calculation
for i=1:length(pasteCodon) %% note: pasteCodon is column vector
[NNfs,~,Fp]=PheAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA
if NNfs<=800
F(i)=-log(Fp)/NNfs;
else
F(i)=NaN;
end
end
[NHISTreal,EdgesReal]=getRealPlot(F);  %% function 'getRealPlot' calculate data for plot

sFile12=(['PheAccurate',order{o},'.mat']);
save(sFile12,'SumHistF','edgesF','NHISTreal','EdgesReal');

%%Gln%%%%
SumHist=zeros(1,1);
for i=1:length(NNq)
     NN=NNq(i);
     if NN<=800
     pref=getPref(NN,2); %% function 'getPfre', read all the Entropy value for possible configurations--pref 
     SumHist=sumHist(SumHist,pref,2,NN);%% function 'sumHist' resore sum value
     end
end
[SumHistF,edgesF,ThoeryPercent]=getTheoryPlot(SumHist); %%function 'getTheoryPlot' calculate data for plot

%%%%%%%%%%real entropy distribution calculation
for i=1:length(pasteCodon) %% note: pasteCodon is column vector
[NNqs,~,Qp]=GlnAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA
if NNqs<=800
Q(i)=-log(Qp)/NNqs;
else
Q(i)=NaN;
end
end
[NHISTreal,EdgesReal]=getRealPlot(Q);  %% function 'getRealPlot' calculate data for plot

sFile13=(['GlnAccurate',order{o},'.mat']);
save(sFile13,'SumHistF','edgesF','NHISTreal','EdgesReal');

%%His%%%%
SumHist=zeros(1,1);
for i=1:length(NNh)
     NN=NNh(i);
     if NN<=800
     pref=getPref(NN,2); %% function 'getPfre', read all the Entropy value for possible configurations--pref 
     SumHist=sumHist(SumHist,pref,2,NN);%% function 'sumHist' resore sum value
     end
end
[SumHistF,edgesF,ThoeryPercent]=getTheoryPlot(SumHist); %%function 'getTheoryPlot' calculate data for plot

%%%%%%%%%%real entropy distribution calculation
for i=1:length(pasteCodon) %% note: pasteCodon is column vector
[NNhs,~,Hp]=HisAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA
if NNhs<=800
H(i)=-log(Hp)/NNhs;
else
H(i)=NaN;
end
end
[NHISTreal,EdgesReal]=getRealPlot(H);  %% function 'getRealPlot' calculate data for plot

sFile14=(['HisAccurate',order{o},'.mat']);
save(sFile14,'SumHistF','edgesF','NHISTreal','EdgesReal');

%%%%Glu%%%%
SumHist=zeros(1,1);
for i=1:length(NNe)
     NN=NNe(i);
     if NN<=800
     pref=getPref(NN,2); %% function 'getPfre', read all the Entropy value for possible configurations--pref 
     SumHist=sumHist(SumHist,pref,2,NN);%% function 'sumGlut' resore sum value
     end
end
[SumHistF,edgesF,ThoeryPercent]=getTheoryPlot(SumHist); %%function 'getTheoryPlot' calculate data for plot

%%%%%%%%%%real entropy distribution calculation
for i=1:length(pasteCodon) %% note: pasteCodon is column vector
[NNes,~,Ep]=GluAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA
if NNes<=800
E(i)=-log(Ep)/NNes;
else
E(i)=NaN;
end
end
[NHISTreal,EdgesReal]=getRealPlot(E);  %% function 'getRealPlot' calculate data for plot

sFile15=(['GluAccurate',order{o},'.mat']);
save(sFile15,'SumHistF','edgesF','NHISTreal','EdgesReal');


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

sFile16=(['ArgAccurate',order{o},'.mat']);
save(sFile16,'SumHistF','edgesF','NHISTreal','EdgesReal');


%%Ser%%%%%
SumHist=zeros(1,1);
for i=1:length(NNs)
    NN=NNs(i);
    if NN<=400
        pref=getPref(NN,6); %% remember to change syno value;  read all the Entropy value for possible configurations--pref
        SumHist=sumHist(SumHist,pref,6,NN);%% SumHist resore sum value
    end
end
[SumHistF,edgesF,ThoeryPercent]=getTheoryPlot(SumHist);

%%%%%%%%%%%%real entropy distribution calculation
for i=1:length(pasteCodon) %% note: pasteCodon is column vector
    [NNss,~,Sp]=SerAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA
    if NNss<=400
        S(i)=-log(Sp)/NNss;
    else
        S(i)=NaN;
    end
end
[NHISTreal,EdgesReal]=getRealPlot(S);

sFile17=(['SerAccurate',order{o},'.mat']);
save(sFile17,'SumHistF','edgesF','NHISTreal','EdgesReal');

%%%Leu%%%
SumHist=zeros(1,1);
for i=1:length(NNl)
    NN=NNl(i);%
    if NN<=400
        pref=getPref(NN,6); %% remember to change syno value;  read all the Entropy value for possible configurations--pref
        SumHist=sumHist(SumHist,pref,6,NN);%% SumHist resore sum value
    end
end
[SumHistF,edgesF,ThoeryPercent]=getTheoryPlot(SumHist);

%%%%%%%%%%%%real entropy distribution calculation
for i=1:length(pasteCodon) %% note: pasteCodon is column vector
    [NNls,~,Lp]=LeuAminoAcidH(pasteCodon{1,i}); %%get Ygg for each mRNA
    if NNls<=400
        L(i)=-log(Lp)/NNls;
    else
        L(i)=NaN;
    end
end
[NHISTreal,EdgesReal]=getRealPlot(L);

sFile18=(['LeuAccurate',order{o},'.mat']);
save(sFile18,'SumHistF','edgesF','NHISTreal','EdgesReal');

end
