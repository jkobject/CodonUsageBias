%%compare with codon table, get the sequence indice of the input sequence
%%for further probability assignment.

function si = getSequenceIndice(pasteCodon)
codonTable={'GGG','GGA','GGT','GGC','GAG','GAA','GAT','GAC',...
    'GTG','GTA','GTT','GTC','GCG','GCA','GCT','GCC',...
    'AGG','AGA','AGT','AGC','AAG','AAA','AAT','AAC',...
    'ATG','ATA','ATT','ATC','ACG','ACA','ACT','ACC',...
    'TGG','TGA','TGT','TGC','TAG','TAA','TAT','TAC',...
    'TTG','TTA','TTT','TTC','TCG','TCA','TCT','TCC',...
    'CGG','CGA','CGT','CGC','CAG','CAA','CAT','CAC',...
    'CTG','CTA','CTT','CTC','CCG','CCA','CCT','CCC'};

tableLength=length(codonTable);

%%manipulate input codon sequence
codonStr=pasteCodon;
strLength=length(codonStr);

for i=1:strLength
    for j=1:tableLength
        if (strcmp(codonStr(i),codonTable(j))); %%find corresponding frequency for each codon
        %%if codonStr{i}==codonTable{j}
            si(i)=j;
           break;
        end
    end
end
end