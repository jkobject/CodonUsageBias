setSynonymousCodonTable;

PasteCodon={pasteCodon11,pasteCodon12,pasteCodon13,pasteCodon14,pasteCodon15,...
    pasteCodon21,pasteCodon22,pasteCodon23,...
    pasteCodon31,pasteCodon32,pasteCodon33,pasteCodon34,pasteCodon35,pasteCodon36,pasteCodon37,...
    pasteCodon41,pasteCodon42,pasteCodon43,pasteCodon44,pasteCodon45};
% PasteCodon={pasteCodon36,pasteCodon43,pasteCodon44};

SpeciesName={'sah6','saeu','saku','Ag','Yl','sp','so','sj','ac','afl','ans','anr','ao','at','afu','ff','fg','fo','fv','tr'};

% SpeciesName={'at','fo','fv'};

for i=1:length(PasteCodon)
    
[Sp(i),Sid{i}]=SelectionPressure(PasteCodon{i});

% figure    
%     
% plot(log(G),Gr,'.');

% hold on
% 
%x=(-0.9):0;
% 
% y=x;
% 
% plot(x,y);
% 
% hold off

% set(gca,'xaxislocation','bottom','yaxislocation','left','xdir','reverse','ydir','reverse')

xlabel('per codon entropy cost of observed sequence');
 
ylabel('per codon average entropy cost of sequence after equal codon replacement');

title([SpeciesName{i},'Amino Acid--Ala(4 synonymous)']);

% disp(['Selection pressure:',num2str(Sp(i))]);

end

tblAla = table((SpeciesName)',(Sp)','VariableNames',{'SpeciesName' 'SelectionPressure'})
% setHeading(tbl, 'Gly Slection Pressure Analysis');
writetable(tblAla);
