function setSynonymousCodonTable

global cfG cfGT ctG cfD cfDT ctD cfE cfET ctE cfV cfVT ctV cfA cfAT ctA cfR cfRT ctR cfS cfST ctS cfK cfKT ctK cfN cfNT ctN cfI cfIT ctI cfT cfTT ctT cfZ cfZT ctZ cfC cfCT ctC cfY cfYT  ctY cfL cfLT ctL cfF cfFT ctF cfQ cfQT ctQ cfH cfHT ctH cfP cfPT ctP;

% prompt = 'Do you want equal probability for each synonymous codon? If so type Y';
% str=input(prompt,'s');

% if isempty(str)
    
%     cfGT=[0.12,0.23,0.45,0.20];  %% this codon table is for sac.. vita.
%     cfDT=[0.65,0.35];
%     cfET=[0.30,0.70];
%     cfVT=[0.20,0.21,0.39,0.20];
%     cfAT=[0.11,0.30,0.37,0.22];
%     cfRT=[0.22,0.47,0.04,0.07,0.14,0.06];%%change 0.21
%     cfST=[0.16,0.11,0.10,0.21,0.26,0.16];
%     cfKT=[0.42,0.58];
%     cfNT=[0.60,0.40];
% %     cfMT=[1.00];
%     cfIT=[0.28,0.46,0.26];
%     cfTT=[0.14,0.31,0.34,0.21];
% %     cfWT=[1.00];
%     cfZT=[0.30,0.23,0.47];
%     cfCT=[0.62,0.38];
%     cfYT=[0.57,0.43];
%     cfLT=[0.28,0.28,0.11,0.14,0.13,0.06];
%     cfFT=[0.59,0.41];
%     cfQT=[0.32,0.68];
%     cfHT=[0.64,0.36];
%     cfPT=[0.12,0.41,0.31,0.16]; 

% else
    cfG=[1/4,1/4,1/4,1/4];
    cfD=[1/2,1/2];
    cfE=[1/2,1/2];
    cfV=[1/4,1/4,1/4,1/4];
    cfA=[1/4,1/4,1/4,1/4];
    cfR=[1/6,1/6,1/6,1/6,1/6,1/6];
    cfS=[1/6,1/6,1/6,1/6,1/6,1/6];
    cfK=[1/2,1/2];
    cfN=[1/2,1/2];
%     cfM=[1.00];
    cfI=[1/3,1/3,1/3];
    cfT=[1/4,1/4,1/4,1/4];
%     cfW=[1.00];
    cfZ=[1/3,1/3,1/3];
    cfC=[1/2,1/2];
    cfY=[1/2,1/2];
    cfL=[1/6,1/6,1/6,1/6,1/6,1/6];
    cfF=[1/2,1/2];
    cfQ=[1/2,1/2];
    cfH=[1/2,1/2];
    cfP=[1/4,1/4,1/4,1/4];
    
% end
%%preprocessed to 1 long sequence which connecte all the genes
% [cfGT,cfDT,cfET,cfVT,cfAT,cfRT,cfST,cfKT,cfNT,cfIT,cfTT,cfZT,cfCT,cfYT,cfLT,cfFT,cfQT,cfHT,cfPT]=CodonFrequency('sah6FF.csv');


ctG={'GGG','GGA','GGT','GGC'}; %%Gly
ctD={'GAT','GAC'};  %%Asp
ctE={'GAG','GAA'};  %%Glu
ctV={'GTG','GTA','GTT','GTC'}; %%Val
ctA={'GCG','GCA','GCT','GCC'}; %%Ala
ctR={'AGG','AGA','CGG','CGA','CGT','CGC'};   %%Arg
ctS={'AGT','AGC','TCG','TCA','TCT','TCC'};   %%Ser
ctK={'AAG','AAA'};   %%Lys
ctN={'AAT','AAC'};   %%Asn
% ctM={'ATG'};  %%Met
ctI={'ATA','ATT','ATC'}; %%Ile
ctT={'ACG','ACA','ACT','ACC'};  %%Thr
% ctW={'TGG'};   %%Trp
ctZ={'TGA','TAG','TAA'};  %%End
ctC={'TGT','TGC'};  %%Cys
ctY={'TAT','TAC'};   %%Tyr
ctL={'TTG','TTA','CTG','CTA','CTT','CTC'};  %%Leu
ctF={'TTT','TTC'}; %%Phe
ctQ={'CAG','CAA'};  %%Gln
ctH={'CAT','CAC'};  %%His
ctP={'CCG','CCA','CCT','CCC'};  %%Pro

end