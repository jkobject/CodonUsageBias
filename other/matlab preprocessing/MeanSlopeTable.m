%% regression for slope and r2


NumSynonymous={'4','2','4','4','4','6','6','2','2','3','4','2','2','6','2','2','2','4'};

NameSynonymous={'Gly(4synonymous)','Asp(2synonymous)','Glu(4synonymous)','Val(4synonymous)','Ala(4synonymous)','Arg(6synonymous)','Ser(6synonymous)','Lys(2synonymous)','Asn(2synonymous)','Ile(3synonymous)','Thr(4synonymous)','Cys(2synonymous)','Tyr(2synonymous)','Leu(6synonymous)','Phe(2synonymous)','Gln(2synonymous)','His(2synonymous)','Pro(4synonymous)'};

NameAminoAcid={'Gly','Asp','Glu','Val','Ala','Arg','Ser','Lys','Asn','Ile','Thr','Cys','Tyr','Leu','Phe','Gln','His','Pro'};

for i=1:18
    
figure

x1= MidFq1{1,i}{:,1};
y1= MidFq1{1,i}{:,2};
x2= MidFq2{1,i}{:,1};
y2= MidFq2{1,i}{:,2};
x3= MidFq3{1,i}{:,1};
y3= MidFq3{1,i}{:,2};
x4= MidFq4{1,i}{:,1};
y4= MidFq4{1,i}{:,2};
x5= MidFq5{1,i}{:,1};
y5= MidFq5{1,i}{:,2};
x6= MidFq6{1,i}{:,1};
y6= MidFq6{1,i}{:,2};

% pl=plot(x1,y1,x2,y2,x3,y3,x4,y4,x5,y5,x6,y6);
% 
% NameArray={'LineStyle','Marker','Color'};
% 
% ValueArray={'none','o',[1,0,0];'none','o',[1,0.3,0];'none','o',[1,0.5,0];'none','o',[0,0,1];'none','o',[0,0.3,1];'none','o',[0,0.5,1]};
% 
% set(pl,NameArray,ValueArray);
% 
% legend('Saccharomyces arboricola_h_6','Saccharomyces eubayanus','Saccharomyces kudriavzevii_ifo_1802','Schizosaccharomyces pombe','Schizosaccharomyces octosporus','Schizosaccharomyces japonicus','Location','Northwest');
% 
% xlabel('midpoint of bin');
%   
% ylabel('log of frequency');
% 
% title([NameAminoAcid(i),'bin plot']);
mean1=sum(x1.*exp(y1))./sum(exp(y1));
slopeOneC1 = polyfit(x1,y1,1); 
slopeOne1=slopeOneC1(1);
% filin1=polyval(slopeOneC1,x1);
mdl1 = fitlm(x1,y1);
figure
plotResiduals(mdl1,'fitted');
title(['Sah6',NameAminoAcid{1}]);
meanAll1(1,i)=mean1;
RsqAll1(1,i)=mdl1.Rsquared.Ordinary;
SlopeAll1(1,i)=slopeOne1;

mean2=sum(x2.*exp(y2))./sum(exp(y2));
slopeOneC2 = polyfit(x2,y2,1); 
slopeOne2=slopeOneC2(1);
% filin2=polyval(slopeOneC2,x1);
mdl2 = fitlm(x2,y2);
figure
plotResiduals(mdl2,'fitted');
title(['Saeu',NameAminoAcid{2}]);
meanAll2(1,i)=mean2;
RsqAll2(1,i)=mdl2.Rsquared.Ordinary;
SlopeAll2(1,i)=slopeOne2;

mean3=sum(x3.*exp(y3))./sum(exp(y3));
slopeOneC3 = polyfit(x3,y3,1); 
slopeOne3=slopeOneC3(1);
mdl3 = fitlm(x3,y3);
figure
plotResiduals(mdl3,'fitted');
title(['Saku',NameAminoAcid{3}]);
meanAll3(1,i)=mean3;
RsqAll3(1,i)=mdl3.Rsquared.Ordinary;
SlopeAll3(1,i)=slopeOne3;

mean4=sum(x4.*exp(y4))./sum(exp(y4));
slopeOneC4 = polyfit(x4,y4,1); 
slopeOne4=slopeOneC4(1);
mdl4 = fitlm(x4,y4);
figure
plotResiduals(mdl4,'fitted');
title(['Sp',NameAminoAcid{4}]);
meanAll4(1,i)=mean4;
RsqAll4(1,i)=mdl4.Rsquared.Ordinary;
SlopeAll4(1,i)=slopeOne4;

mean5=sum(x5.*exp(y5))./sum(exp(y5));
slopeOneC5 = polyfit(x5,y5,1); 
slopeOne5=slopeOneC5(1);
mdl5 = fitlm(x5,y5);
figure
plotResiduals(mdl5,'fitted');
title(['So',NameAminoAcid{5}]);
meanAll5(1,i)=mean5;
RsqAll5(1,i)=mdl5.Rsquared.Ordinary;
SlopeAll5(1,i)=slopeOne5;

mean6=sum(x6.*exp(y6))./sum(exp(y6));
slopeOneC6 = polyfit(x6,y6,1); 
slopeOne6=slopeOneC6(1);
mdl6 = fitlm(x6,y6);
figure
plotResiduals(mdl6,'fitted');
title(['Sj',NameAminoAcid{6}]);
meanAll6(1,i)=mean6;
RsqAll6(1,i)=mdl6.Rsquared.Ordinary;
SlopeAll6(1,i)=slopeOne6;

end

tbl = table((NumSynonymous)',(SlopeAll1)',(RsqAll1)',(meanAll1)',(SlopeAll2)',(RsqAll2)',(meanAll2)',(SlopeAll3)',(RsqAll3)',(meanAll3)',(SlopeAll4)',(RsqAll4)',(meanAll4)',(SlopeAll5)',(RsqAll5)',(meanAll5)',(SlopeAll6)',(RsqAll6)',(meanAll6)','VariableNames',{'NumSynonymous' 'Slope1' 'Rsq1' 'mean1' 'Slope2' 'Rsq2' 'mean2' 'Slope3' 'Rsq3' 'mean3' 'Slope4' 'Rsq4' 'mean4' 'Slope5' 'Rsq5' 'mean5' 'Slope6' 'Rsq6' 'mean6'},'RowNames',(NameAminoAcid)');
     
Tbl = sortrows(tbl,{'NumSynonymous'},'ascend')

setHeading(Tbl,'After Replacement');
