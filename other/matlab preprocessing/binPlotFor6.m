% clc  %% bin plot for 6 species
% clear

NumSynonymous={'4','2','4','4','4','6','6','2','2','3','4','2','2','6','2','2','2','4'};

NameSynonymous={'Gly(4synonymous)','Asp(2synonymous)','Glu(4synonymous)','Val(4synonymous)','Ala(4synonymous)','Arg(6synonymous)','Ser(6synonymous)','Lys(2synonymous)','Asn(2synonymous)','Ile(3synonymous)','Thr(4synonymous)','End(3synonymous)','Cys(2synonymous)','Tyr(2synonymous)','Leu(6synonymous)','Phe(2synonymous)','Gln(2synonymous)','His(2synonymous)','Pro(4synonymous)'};

NameAminoAcid={'Gly','Asp','Glu','Val','Ala','Arg','Ser','Lys','Asn','Ile','Thr','Cys','Tyr','Leu','Phe','Gln','His','Pro'};


% setSynonymousCodonTable;
% 
% 
% % pasteCodon1=getCodonSequence('sah6.csv'); %% call function getCodonSequence to creat input 
% 
% MidFqT1=AllHMidF(pasteCodon1);
% 
% % pasteCodon2=getCodonSequence('saeu.csv');
% 
% MidFqT2=AllHMidF(pasteCodon2);
% 
% % pasteCodon3=getCodonSequence('saku.csv');
% 
% MidFqT3=AllHMidF(pasteCodon3);
% 
% % pasteCodon4=getCodonSequence('sp.csv');
% 
% MidFqT4=AllHMidF(pasteCodon4);
% 
% % pasteCodon5=getCodonSequence('so.csv');
% 
% MidFqT5=AllHMidF(pasteCodon5);
% 
% % pasteCodon6=getCodonSequence('sj.csv');
% 
% MidFqT6=AllHMidF(pasteCodon6);


for i=1:18
    
figure

x1= MidFqT1{1,i}{:,1};
y1= MidFqT1{1,i}{:,2};
x2= MidFqT2{1,i}{:,1};
y2= MidFqT2{1,i}{:,2};
x3= MidFqT3{1,i}{:,1};
y3= MidFqT3{1,i}{:,2};
x4= MidFqT4{1,i}{:,1};
y4= MidFqT4{1,i}{:,2};
x5= MidFqT5{1,i}{:,1};
y5= MidFqT5{1,i}{:,2};
x6= MidFqT6{1,i}{:,1};
y6= MidFqT6{1,i}{:,2};

pl=plot(x1,y1,x2,y2,x3,y3,x4,y4,x5,y5,x6,y6);


NameArray={'LineStyle','Marker','Color'};

ValueArray={'none','o',[1,0,0];'none','o',[1,0.3,0];'none','o',[1,0.6,0];'none','o',[0,0,1];'none','o',[0,0.3,1];'none','o',[0,0.6,1]};

set(pl,NameArray,ValueArray);

legend('Saccharomyces arboricola_h_6','Saccharomyces eubayanus','Saccharomyces kudriavzevii_ifo_1802','Schizosaccharomyces pombe','Schizosaccharomyces octosporus','Schizosaccharomyces japonicus','Location','Northwest');

xlabel('midpoint of bin');
  
ylabel('log of frequency');

title([NameAminoAcid(i),'(Pr/Pm) probability ratio distribution after replacement according to Codon Table']);

end