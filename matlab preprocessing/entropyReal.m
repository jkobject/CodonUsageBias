%% engtropyReal can calculate entropy of one real input sequence, codon table is from SGD.
function S = entropyReal(pasteCodon)

%%codon talbe
codonTable={'GGG','GGA','GGT','GGC','GAG','GAA','GAT','GAC',...
    'GTG','GTA','GTT','GTC','GCG','GCA','GCT','GCC',...
    'AGG','AGA','AGT','AGC','AAG','AAA','AAT','AAC',...
    'ATG','ATA','ATT','ATC','ACG','ACA','ACT','ACC',...
    'TGG','TGA','TGT','TGC','TAG','TAA','TAT','TAC',...
    'TTG','TTA','TTT','TTC','TCG','TCA','TCT','TCC',...
    'CGG','CGA','CGT','CGC','CAG','CAA','CAT','CAC',...
    'CTG','CTA','CTT','CTC','CCG','CCA','CCT','CCC'};
codonFraction=[0.12,0.23,0.45,0.20,0.30,0.70,0.65,0.35,...
    0.20,0.22,0.39,0.20,0.11,0.30,0.37,0.22,...
    0.21,0.47,0.16,0.11,0.42,0.58,0.60,0.40,...
    1.00,0.28,0.46,0.26,0.14,0.31,0.34,0.21,...
    1.00,0.30,0.62,0.38,0.23,0.47,0.57,0.43,...
    0.28,0.28,0.59,0.41,0.10,0.21,0.26,0.16,...
    0.04,0.07,0.14,0.06,0.32,0.68,0.64,0.36,...
    0.11,0.14,0.13,0.06,0.12,0.41,0.31,0.16];
tableLength=length(codonTable);

%%manipulate input codon sequence
inCodonstr=pasteCodon;
codonStr=strsplit(inCodonstr);
strLength=length(codonStr);
codonFrequency=zeros(1,strLength);

for i=1:strLength
    for j=1:tableLength
        if (strcmp(codonStr(i),codonTable(j))); %%find corresponding frequency for each codon
        %%if codonStr{i}==codonTable{j}
           codonFrequency(i)=codonFraction(j);
           break;
        end
    end
end

fprintf('length=');
disp(strLength);

S = sum(codonFrequency .* log(codonFrequency));%% calculate entropy for input codon sequence 
fprintf('entropy=');
disp(S);

end








