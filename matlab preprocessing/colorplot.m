pl=plot(x1,y1,x2,y2,x3,y3,x4,y4,x5,y5,x6,y6);
NameArray={'LineStyle','Marker','Color'};
ValueArray={'none','o',[1,0,0];'none','o',[1,0.3,0];'none','o',[1,0.6,0];'none','o',[0,0,1];'none','o',[0,0.3,1];'none','o',[0,0.6,1]};
set(pl,NameArray,ValueArray);
legend('Sccharomyces arboricola_h_6','Saccharomyces eubayanus','Saccharomyces kudriavzevii_ifo_1802','Schizosaccharomyces pombe','Schizosaccharomyces octosporus','Schizosaccharomyces japonicus','Location','Northwest');