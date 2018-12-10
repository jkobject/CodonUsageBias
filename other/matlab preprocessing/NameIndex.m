run on myrtle command
nohup matlab -nodisplay -nojvm < savAveEntropy101t200.m > output file
scp yd46@myrtle.kent.ac.uk:/home/cug/yd46/MATLAB/mydatatest.mat /Users/yd46/Documents/MATLAB/ 

Partitions of length 100 for 6 synonymous codons: 96560647, running about 3 hours.

1.
http://downloads.yeastgenome.org/unpublished_data/codon/ysc.orf.cod ??codon table

http://downloads.yeastgenome.org/sequence/S288C_reference/orf_dna/

SPECIES NAME( http://fungi.ensembl.org/info/website/ftp/index.html) 

Saccharomyces arboricola_h_6

Saccharomyces eubayanus

Saccharomyces kudriavzevii

Schizosaccharomyces pombe

Schizosaccharomyces octosporus

Schizosaccharomyces japonicus

2. Amino Acid 

NameSynonymous={'Gly(4synonymous)G','Asp(2synonymous)D','Glu(4synonymous)E','Val(4synonymous)V',...
    'Ala(4synonymous)A','Arg(6synonymous)R','Ser(6synonymous)S','Lys(2synonymous)K',...
    'Asn(2synonymous)N','Ile(3synonymous)I','Thr(4synonymous)T','End(3synonymous)Z',...
    'Cys(2synonymous)C','Tyr(2synonymous)Y','Leu(6synonymous)L','Phe(2synonymous)F',...
    'Gln(2synonymous)Q','His(2synonymous)H','Pro(4synonymous)P'};
short={G,D,E,V,A,R,S,K,N,I,T,Z,C,Y,L,F,Q,H,P};

3. Rawdata File name
input sequence: pasteCodon1 2 3 4 5 6
subsequence for each AA: SqG1--SqG6 
p ratio midpoint for one species: MidFq1{1,1}--MidFq1{1,19}
                                  MidFq2{1,1}--MidFq2{1,19}
                                  for example each G in 6 species MidFq1{1,1}--MidFq6{1,1}
                                      
4. important potencial record file
all21AAPrange: 21 Amino Acid of SV species, probability ratio range 
MeanSlopeTable: slopes among 6 species
midsah6,midsaeu,midsaku,midsp,midso,midsj:before replacement
midafterR:after replacement with equal mutate probability

5. terms for writing 
probability densities not probabilities
sampling distribution of mean :http://stattrek.com/sampling/sampling-distribution.aspx
normal approximation

monte carlo method
master equation