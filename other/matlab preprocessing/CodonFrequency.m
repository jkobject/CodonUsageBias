function [cfGT,cfDT,cfET,cfVT,cfAT,cfRT,cfST,cfKT,cfNT,cfIT,cfTT,cfZT,cfCT,cfYT,cfLT,cfFT,cfQT,cfHT,cfPT]=CodonFrequency(filename)
SeqNT=fileread(filename);
cc=codoncount(SeqNT);
cfGT=[cc.GGG,cc.GGA,cc.GGT,cc.GGC]./(cc.GGG+cc.GGA+cc.GGT+cc.GGC);
cfDT=[cc.GAT,cc.GAC]./(cc.GAT+cc.GAC);
cfET=[cc.GAG,cc.GAA]./(cc.GAG+cc.GAA);
cfVT=[cc.GTG,cc.GTA,cc.GTT,cc.GTC]./(cc.GTG+cc.GTA+cc.GTT+cc.GTC);
cfAT=[cc.GCG,cc.GCA,cc.GCT,cc.GCC]./(cc.GCG+cc.GCA+cc.GCT+cc.GCC);
cfRT=[cc.AGG,cc.AGA,cc.CGG,cc.CGA,cc.CGT,cc.CGC]./(cc.AGG+cc.AGA+cc.CGG+cc.CGA+cc.CGT+cc.CGC);
cfST=[cc.AGT,cc.AGC,cc.TCG,cc.TCA,cc.TCT,cc.TCC]./(cc.AGT+cc.AGC+cc.TCG+cc.TCA+cc.TCT+cc.TCC);
cfKT=[cc.AAG,cc.AAA]./(cc.AAG+cc.AAA);
cfNT=[cc.AAT,cc.AAC]./(cc.AAT+cc.AAC);
cfIT=[cc.ATA,cc.ATT,cc.ATC]./(cc.ATA+cc.ATT+cc.ATC);
cfTT=[cc.ACG,cc.ACA,cc.ACT,cc.ACC]./(cc.ACG+cc.ACA+cc.ACT+cc.ACC);
cfZT=[cc.TGA,cc.TAG,cc.TAA]./(cc.TGA+cc.TAG+cc.TAA);
cfCT=[cc.TGT,cc.TGC]./(cc.TGT+cc.TGC);
cfYT=[cc.TAT,cc.TAC]./(cc.TAT+cc.TAC);
cfLT=[cc.TTG,cc.TTA,cc.CTG,cc.CTA,cc.CTT,cc.CTC]./(cc.TTG+cc.TTA+cc.CTG+cc.CTA+cc.CTT+cc.CTC);
cfFT=[cc.TTT,cc.TTC]./(cc.TTT+cc.TTC);
cfQT=[cc.CAG,cc.CAA]./(cc.CAG+cc.CAA);
cfHT=[cc.CAT,cc.CAC]./(cc.CAT+cc.CAC);
cfPT=[cc.CCG,cc.CCA,cc.CCT,cc.CCC]./(cc.CCG+cc.CCA+cc.CCT+cc.CCC);
end















