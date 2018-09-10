# STORY

we load the data from ensembl on 'ensembl2' workspace for SC and fungi with 6692 homologies saved under 'full_16august'
we retrieve this list of species

## remarks about the run

- the load takes about 1h30 to run 
- memory leakage still
- There seems to be islands of genes with a lot of recorded SNPs
- saving the full data seems to take approx 5mn and require around 20gb of RAM for 6700 homologies and 550 species, it takes 3.5Gb uncompress and 0.3Gb compresssed
- it takes 10mn and 16Gb of RAM to load the full homology matrix into the homoset object

- all the runs have been made on a c5d.xlarge Amazon Cloud compute instance with 150Gb of storage, Intel Xeon Scalable 3,0 GHz turbo boost 3,5 GHz materialy accelerated hypervisor  of 4 vCPU, >1gb internet connection, 8Gb of RAM and additional 50Gb of slow Swap memory


the number of genes retrieved is 1730558
and the final number of valid homologies is 5807

## clustering homosets

plotting all the homologies we have :

![png](output_19_0.png)

homologies per species
![png](output_21_0.png)

the final clustering of the homologies into 4 with scores of 
the quality of the clustering is: [silhouette_score,calinski_harabaz_score]
    0.435136981254
    6397.88903542

![png](output_25_2.png)


![png](output_25_3.png)


#### using PCA

it is a bit puzzling to see that there is no homology that is shared by all species, thought this should be the case for the set of genes essential to life. 

We postulate that this might show the alignment errors which only looks at similarity between DNA sequences and not at the final function of protein coding genes.

It might be interesting to have a look at what factors influence the accuracy of such algorithm, to be able to account for it
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4267518/
More data (for example a mix of DNA sequence alignment, protein alignment, RNA seq and taxonomical information) might help greatly improve the accuracy of such software. 

Values seem to be placed as one would think with a high density clusters ( maybe of entropy not driven by any major effects) and then a spreading of entropy from many different directions (amino acids) 
is there amino acids with less spread in their CUB? the high values seems to be randomly distributed but we can notice 2 things: 
- there is a limit to which one amino acid entropy is decorelated to another (this triangle shape)
- there seems to be specific entropy direction than other but it could be noise.

#### using tSNE
This plot seems pretty magical, expecially compared to what PCA outputs, However knowing the flaws and advantages of tSNE: https://distill.pub/2016/misread-tsne/, here what might come to mind. Here is really highlighted any relationship in the datapoint accross the 18 dimensions. 
we can see that this plot creates very distinctive shapes: circle. Taken as is, this would mean that among the genes there is groups where the entropy values are not the same but the amount of entropy the totality of the genome has is mostly the same. this means, the entropy over the 18 amino acid is distributed such that genes are pushed to have a total entropy of each distribution of their synonimous codon around the same point. althought being an interesting point, this clustering around a particular radius could be tSNE amplifying a noise in the distribution, especially when an important number of points center on such circular distribution, which is effectivelly the case with an important number of low entropy genes.

What is more intriguing however is that it is not one but 4 distinct circle that appear, as if there would be 4 distinct pairs of amino acids for where entropy tends to be correlated when the entropy is high. This might make more sense given the fact that in amino acids with 4+ codons, the entropy is easily higher than elsewhere especially for ARG, LEU, SER where there is 6 such codons. 
trying to have some form of coloring of the plot might help understanding what may drive such correlation 

overhall the why of such circular shapes remains an open question which would need to be further explored.

---------
The tSNE and PCA all_plot are the same for working and all homoset.
Meaning that there is no latent process driving the distribution in the highly shared homologies.

One can see such plot as a visual representation of non linear correlation amongst the variables. Here whatever is the values and even while changing the data, we can see that there is 4 axis correlated to their center in some ways. what drives it is unknown and should be explored further.

-----------

*see plots*

when extracting a working homology (homology with high number of species)
we are left with 3192 comprising 1443683 genes

-----------

### looking at one group of homologies

homo_namelist[810:816] what we can see here is that there is some interesting things. first, even if some homologies can be spread out, some of them are grouped pretty well hinting on some factors related to translation. Many things may be able to explain that and it is not so much the mean entropy or length, more complex factors thus have to explain that.

homo_namelist[810:816] Another look here shows that one of the homologies (the only secreted one) has an extremelly high CUB, there may be an explanation as it is an highly expressed excreted protein. this might lead to such an important information gain and might be further tested in the following plots

------------
On this project there is many things to see and to investigate. Several hypothesis can be make on the many different plots and one would need the functions of different homologous genes and what the species are 



----------------

### contiunued first exploration
we can see here that there is a strongcorrelation of the entropy to the number of synonimous codons, however, it is not enough to be able to classify it (for example PRO(14,4 codons) is less informative than LYS(12, 2 codons)

Moreover, you can see it has been decided not to normalize the vectors. this is because we are never comparing the values of the vector to any other values for classifying or clustering etc.. 
We are only using the vectors. and we might want to compare the different values within the vectors, different vectors among homologies or different vectors among different hmologies etc. rescaling the vectors within themselves would not change a thing in comparing the different values since we would still have the same disrepency between different dimensions and we would not be able to compare different vectors, rescaling all the vectors from within an homology would prevent us from comparing different vectors together.. etc. not rescaling here seems to be the better option (for any other variable in the regression step we preprocess and normalize all of them.

homology averages : [  5.94535188,   7.55111827,   2.21075474,   2.11802318,
             0.74706038,   1.99656174,   3.07157284,   6.80222758,
             1.17308691,   5.8847418 ,  12.63998386,   4.18494766,
             1.99909645,   3.98521084,   6.06791147,   4.28587105,
             1.6369395 ,   6.39921517])
codons = {
    \'ALA': ['GCA', 'GCC', 'GCG', 'GCT'],  # GC
    'ARG': ['CGA', 'CGC', 'CGG', 'CGT', 'AGG', 'AGA'],  # CG
    'ASN': ['AAC', 'AAT'],  # ASN - LYS
    'ASP': ['GAT', 'GAC'],  # ASP - GLU
    'CYS': ['TGT', 'TGC'],  # ~
    'GLN': ['CAA', 'CAG'],  # GLN - HIS
    'GLU': ['GAG', 'GAA'],
    \'GLY': ['GGT', 'GGG', 'GGA', 'GGC'],  # GG
    'HIS': ['CAT', 'CAC'],
    'ILE': ['ATC', 'ATA', 'ATT'],  # ILE -- MET
    'LEU': ['TTA', 'TTG', 'CTC', 'CTT', 'CTG', 'CTA'],  # CT
    'LYS': ['AAA', 'AAG'],
    # 'MET': ['ATG'],
    'PHE': ['TTT', 'TTC'],  # PHE - LEU
    'PRO': ['CCT', 'CCG', 'CCA', 'CCC'],  # CC
    'SER': ['TCT', 'TCG', 'TCA', 'TCC', 'AGC', 'AGT'],  # TC
    # 'STOP': ['TAG', 'TGA', 'TAA'],
    'THR': ['ACC', 'ACA', 'ACG', 'ACT'],  # AC
    # 'TRP': ['TGG'],
    'TYR': ['TAT', 'TAC'],
    'VAL': ['GTA', 'GTC', 'GTG', 'GTT']}  # GT

homology variances [  7.65164553,   8.51687922,   3.58228369,   3.8571531 ,
             1.16062417,   3.38140755,   5.18205945,   8.31587079,
             1.85351603,   6.52986744,  13.64420273,   6.37602431,
             3.04094619,   5.04495174,   6.78347536,   5.14261895,
             2.53622075,   7.12209085]))

same for working
[  5.97628942,   7.70757895,   2.22471043,   2.13199829,
             0.74290188,   2.01099741,   3.17823769,   6.85060898,
             1.19385793,   5.96524715,  12.82902082,   4.46611862,
             2.00879056,   3.99377676,   6.03202785,   4.25497258,
             1.65132861,   6.49788342]

[  7.67603982,   8.76717122,   3.57895547,   3.92763532,
             1.15799442,   3.40407144,   5.348315  ,   8.26130029,
             1.88114428,   6.61396437,  13.98781425,   6.66422731,
             3.02849258,   5.05068292,   6.70552172,   5.01830508,
             2.55054707,   7.22239885]

We can see that the values are not the same amongst them (even if low, the size of the set >10^6 is enough to give signicance to such differences) Wich hint toward different factors on the highly shared homologies

the max values are 

array([ 315.11608917,  327.41540819,  160.0414408 ,  487.89411216,
             52.94063673,  167.46844147,  179.16199446,  355.80320172,
             83.20523309,  174.17504316,  703.43954027,  163.36772378,
            135.52066175,  407.63142791,  317.63055617,  263.0384299 ,
             92.81086287,  349.44444376])

and are conserved onto the working homoset


----------
Here is the values of a sampled homology 

homology: YMR186W
    size: 550
    fullmax: 95.030388428
    var: [  95.27018358   64.0558798    28.77263303   39.97858985    0.60478574
       26.22099706  380.76446341   39.49440197    3.86260229   78.66340124
      290.91663999  257.0192746    36.23927396   40.34022531  122.87893448
      100.83216627   12.61892141   70.24443896]
    mean: [ 18.2098458   16.90886703  10.07872097   5.33678356   0.48330023
       6.84888559  28.51495954  15.52612261   2.98357867  22.34445412
      39.59356623  32.93171482   8.94501112  10.16085624  22.35686529
      19.95389706   5.72739935  19.77422222]
    metrics: {}
    nans mean: 0.154545454545
    lengths mean/var: 663.372727273, 151.73919366
    doublon sum: 54
    GCcount mean/var: 0.511269728852, 0.0720579496983
    KaKs_Scores: None
    isrecent: False
    ishighpreserved: False
    ref: None
    refprot: YMR186W
    refgene: YMR186W
    ecai mean/var: 0.49489271782, 0.0974550502757
    cai mean/var: 0.589823682579, 0.0972128372504
    protein_abundance: 218787
    weight: 80899
    conservation: 5.13216885007
    mRNA_abundance: 6.4
    cys_elements: 0
    is_secreted: True
    decay_rate: 0.0007
    tot_volume: 95785.0
    mean_hydrophobicity: -0.6
    glucose_cost: 673.77
    synthesis_steps: 3483
    isoelectricpoint: 6.085
    othercods: 18.9054545455

Even if the entropy has to and is sometimes correlated with the length of the vectors and the number of codons, there is much stronger correlation to other factors as we can see.

----------

the clustering of the homologies was always pretty efficient (except when high variances amongst the homologies (freq / entropy)) finding always below averages clusters (between 2-10)

*see plots*
![png](output_29_1.png)


## Cluster Stats

What we have here is quite interesting: 

first there seems to be some homologies (0-620)  that are entirely comprised of primary cluster (meaning that all the CUB vectors are gathering around a way smaller  nsphere than anywhere else. There seems to be  (864,764,989,1107,1179,1227,1257,1527) and much more at higher values that seems to be mostly comprised of secondary cluster. this means that there is, within these homologies, small clusters of species which CUBs are very similar. and finally there seems to be groups of homologies (5430-5800) that are entirely composed of outliers these are homologies where the variance between gene entropy values is very high and it may also be that for some of these, the size of the homology is much lower, increasing the probability of not finding any cluster (to a lesser extent, there is also such homologies in the working homoset (comprised only of highly shared homologies)

second, for the species there seems to be a less widespreed difference but the variation are still very much continuous. for a group of the first 4, most of their genes are in a very small amount of principal cluster and it appears also (even more with the smaller n-sphare radius value) on the working-homo one. the percentage of secondary cluster varies a bit for species that have much more presence in primary cluster and a group of species (514-519) are almost only in primary groups.


Here we can see that amongst the outlier species we have only species in very few homologies except for __hepatospora_eriocheir__ which is present in ~700 and still present such CUB behavior.

Having a look at one species:
    name: xylona_heveae_tc161
    ----
    genome_size: 24337763
    ----
    metadata: [True, False, False, False, False]
    ----
    average_entropy: [ 22.4642271   26.60688608  10.99742797   3.02795858   0.31169614
      12.12520353  14.76008424   9.9544623    2.05144308  15.22045262
      41.84146647  18.65835059   2.5560994   13.15196309  17.60777759
      11.44584344   9.31768831  14.03726015]
    ----
    average_size: 783.3
    ----
    meanGChomo: 0.288298791225
    ----
    tot_homologies: 10
    ----
    meanecai: None
    ----

## plotting homologies

There is a lot to say here. first we can see that ecai and similarity scores correlate more with entropy than anything else and correlate similarly with each other. It means that mean entropy retrieves evolutionary drift information of different genes. 
Moreover, we can see that some very similar cluster in their entropy and shape, reveal that two different genes are encoding very similar proteins (see their size, etc). It is also almost exactly the same cluster. However their ecai is not the same for both (as if compared with two different genes)
It can mean two things. There is a problem with ensembl's database and it shows the wrong homologous genes to a reference one. Or we have homologous genes that have not been found and instead placed into different clusters. Or we are seeing two genes that produce very similar proteins and are extremelly related to one another (as they have drifted and evolved together) this is very well picked up by entropy and homologous comparisons. It shows even more when plotting the homologies together and seeing that they overlaps perfeclty and that the length of the genes are the same as well. YBL027W;YBR084C-A and YOL120C;YNL301C and YER131W;YGL189C

Moreover, these homologies being exactly the same, might hint on what were the density clusters in the all_homoplot, earlier. 

Also we can also see that some other homologies are completely different in their CUB but lead to similar species being into similar clusters which very much shows a similar drift of the genome toward a species CUB for some genes.

It is also intersting to see that the correlation often differ amongst homologies, even if similarity scores and ecai stay somewhat the same, length, nans ad GCcount sometimes explain quite a few of the differences in codon usage bias we see. YGL189C;YER131W;YPR102C


'YER131W'

similarity_scores: SpearmanrResult(correlation=0.43487331185717687, pvalue=2.6003821222295293e-25)
    nans: SpearmanrResult(correlation=-0.27068816195697354, pvalue=3.7721552777484254e-10)
    ecai: SpearmanrResult(correlation=0.43033486715175473, pvalue=9.1365917867182727e-25)
    cai: SpearmanrResult(correlation=-0.0789140703800634, pvalue=0.072731142201338167)
    lenmat: SpearmanrResult(correlation=0.026909476813379948, pvalue=0.54114636079113398)
    GCcount: SpearmanrResult(correlation=0.1025444230874635, pvalue=0.019575315835816339)

YGL189C

similarity_scores: SpearmanrResult(correlation=0.42720146613820714, pvalue=1.9436464141625602e-24)
    nans: SpearmanrResult(correlation=-0.27136737481256884, pvalue=3.2660076243783963e-10)
    ecai: SpearmanrResult(correlation=0.42240622026574848, pvalue=7.1041628456677805e-24)
    cai: SpearmanrResult(correlation=-0.068786189053128929, pvalue=0.11755165598522778)
    lenmat: SpearmanrResult(correlation=0.02285945138870555, pvalue=0.60335183722590546)
    GCcount: SpearmanrResult(correlation=0.092770755313180533, pvalue=0.034607094161626988)

 silhouette: 0.388873696407
    cal_hara: 218.980391519
    cluster_phylodistance: [0.8504936969340643, 0.8063694420020829, 1.009453945168029]

Correlation between similarity and ecai:
 SpearmanrResult(correlation=0.47852460467219776, pvalue=4.655130219170088e-31)

we also retrieve the distances:

![png](output_15_0.png)

they are in concordances with the species we have but not very precise when we span further away. this is because we don't have more data. we would need a bigger tree to have sufficient statistics for the 500 species we have here. the other ones would support information on the real distances between the 500


### homologies outliers

very high variances, similarities that do not relate to the change in CUB at all, inverse correlation to many things but correlation to GCcount. hinting on a different kind of bias, showing a really strong one however. There is strong correlation to length and GCcount here 
The similarities are really low. Meaning that there is a really strong evolutionary constraint for them. It would be wise for this group, to extract it and compute other clustering hyperparameters.

Similarly here we have instances of almost exactly the same homologies (twin homologies)
It also seems that some of them are completly different than the rest. YKL182W;YNR016C
For them, the similarity score is much higher.
FOr all the other, it seems that the reference gene has also a very highly entropy value 

YKL182W

similarity_scores: SpearmanrResult(correlation=0.037859815466143781, pvalue=0.37551565432676903)
    nans: SpearmanrResult(correlation=-0.1669341654415212, pvalue=8.3654587165791271e-05)
    ecai: SpearmanrResult(correlation=-0.013795450543254805, pvalue=0.74683816587670093)
    cai: SpearmanrResult(correlation=-0.18451131554398742, pvalue=1.3320448316152529e-05)
    lenmat: SpearmanrResult(correlation=0.49269384700572383, pvalue=5.6149276120470033e-35)
    GCcount: SpearmanrResult(correlation=0.24454986571941831, pvalue=6.224591628537765e-09)

YNR016C

similarity_scores: SpearmanrResult(correlation=0.15296426892106091, pvalue=0.00056934708193189247)
    nans: SpearmanrResult(correlation=-0.17110272856249881, pvalue=0.00011328878575171294)
    ecai: SpearmanrResult(correlation=-0.15894624695793266, pvalue=0.00034061809294703635)
    cai: SpearmanrResult(correlation=-0.3503784401533328, pvalue=5.284184358566423e-16)
    lenmat: SpearmanrResult(correlation=0.1803169711666025, pvalue=4.6716433281590285e-05)
    GCcount: SpearmanrResult(correlation=0.38903938597346588, pvalue=1.1677370420646218e-19)

mean ecai: 0.44298859627
    mean cai: 0.608844955289


It is very interesting to see in primary homo, how some homologies have a huge entropy for one or two dimension only, meanin that the entropy is only driven by one amino acid

YPL249C-A


similarity_scores: SpearmanrResult(correlation=0.15668458797295806, pvalue=0.00067900733673671384)
    nans: SpearmanrResult(correlation=0.28586330406273419, pvalue=3.1180931951463311e-10)
    ecai: SpearmanrResult(correlation=0.35815742616977347, pvalue=1.3990479861887904e-15)
    cai: SpearmanrResult(correlation=-0.059886576258249898, pvalue=0.19640866101880053)
    lenmat: SpearmanrResult(correlation=0.0051118698202435491, pvalue=0.91227208771407242)
    GCcount: SpearmanrResult(correlation=0.1392534124048882, pvalue=0.0025618902813443235)

avg similarity Scores
    0.679533539952
    ------------------------------------
    mean ecai: 0.754447187918
    mean cai: 0.551622420676


YDL130W

similarity_scores: SpearmanrResult(correlation=0.32963398880492711, pvalue=2.2941910765484818e-14)
    nans: SpearmanrResult(correlation=-0.050120172803386744, pvalue=0.25902776911849618)
    ecai: SpearmanrResult(correlation=0.02066270840375984, pvalue=0.64187740024601836)
    cai: SpearmanrResult(correlation=-0.084657244755652339, pvalue=0.05630257313438284)
    lenmat: SpearmanrResult(correlation=-0.45581051665020356, pvalue=1.7659741848818821e-27)
    GCcount: SpearmanrResult(correlation=0.21769726581613691, pvalue=7.0858251668679864e-07)

We alleviate the high variance problem for freq & entropy by computing different hyperparameters for different groups of homologies with more or less variances (overhall the change is not that important (ranging from 5 to 15) for entropy) here.

#### example for lowest ones

![png](output_66_1.png)


## final computations

We can very much see by the shape that mostly one direction is present where all homology go. Moreover, the highest homology species, have de facto the highest variance

*amongst the things that do not work: distance to tRNA UB, full/mean GC diff, tRNA number, phylodistance*


### regression on species

num_genes: SpearmanrResult(correlation=-0.20392120936213901, pvalue=2.635831327646218e-06)
    genome_size: SpearmanrResult(correlation=-0.14812823660067045, pvalue=0.00068636971631414411)
    fullGCcount: SpearmanrResult(correlation=0.26672770045972816, pvalue=5.9432859937195404e-10)
    varGCcount: SpearmanrResult(correlation=0.18647347951225177, pvalue=1.8039311033672066e-05)
    phylogenetic distances: SpearmanrResult(correlation=-0.10299262826299072, pvalue=0.018586274765915109)
    the R^2 score is of: -0.104337545306
    -------------------------------
    average_size: -0.703058659725
    num_genes: -1.75739703204
    genome_size: 0.459321538785
    fullGCcount: 2.29165929882
    varGCcount: -0.763439912152
    tot_homologies: 0.238266815121
    isplant_symbiotic: 0.465421263163
    isanimal_pathogen: -0.307053477627



ow, it seems that all highly conserved homologies are mostly present on one part of the plot and the other and highly preserved ones on the other part. this plot possess a really high amount of information. It seems they both come from the same place when they have a low entropy but the higher the entropy the higher the divergence in how this high entropy is represented. which mean that the entropy distribution amongst the more highly conserved genes and the more recent ones is very different.

the similarity score also shows some interesting pattern. there really seems to be some entropy related with the GC content and some other with how recent is the gene (i.e. other factors) and other with the length

There again seems to be cluster of secreted proteins and a relation of cys reg elements to lengthy/ high variance/ high entropy homologies

the hydrophobicity also seems to predict some information about the relative position and thus entropy values !

Evolutionary CAI also seems to be highly predictive of this shape. we can see that low entropy homologies and high entropy/high variance/conserved homologies are at the opposite of the spectrum with a degraded between the two (which is not the case of the similarity values)

*we really required both entropy and a value, each allowed me to show and discover different relationships*

*we basically used homology to understand a notion of gene ages etc..*


**see plots**

-----------

Here the future seems also grim for the algorithm, however some hope exists, we have found that the lasso was finally not that the interaction may not be linear and may be important amongst the variable as the score of the MLP was consistently better that the one of the lasso and that varying the available parameters was, for some set of paramters, increasing the R^2 score. we may perform a better regression by carefully choosing the set of parameters to give to the lasso


similarity_scores =-0.081261407786012513, pvalue=5.5983630443298877e-10
nans =-0.3555959460371706, pvalue=1.1441764655930061e-172
lenmat =0.80583957969973008, pvalue=0.0
GCcount =0.22619613607892744, pvalue=2.8877807326984728e-68
weight =0.72817624858641561, pvalue=0.0
protein_abundance =0.18651697015748744, pvalue=1.2999207161464825e-46
mRNA_abundance =0.050599719595574855, pvalue=0.00011454073869885887
decay_rate =0.14022360231537928, pvalue=6.8717537057815878e-27
is_secreted =0.090758808355216727, pvalue=4.2391813217099995e-12
cys_elements =0.12857575316090689, pvalue=7.8244087477063335e-23
tot_volume =0.77385917280841998, pvalue=0.0
mean_hydrophobicity =0.11672375791852342, pvalue=4.5226765754239889e-19
glucose_cost =0.77214059050694916, pvalue=0.0
synthesis_steps =0.77459143311298995, pvalue=0.0
isoelectricpoint =-0.026971626011661889, pvalue=0.039853086935354678
meanecai =-0.3734691190980633, pvalue=1.2365551211497936e-191
meancai =-0.1372190049285606, pvalue=8.2741500233278835e-26
conservation =0.080801201912223494, pvalue=6.9991869867232414e-10


##entropy location
then we compute the entropy location 
it takes 2H for all $ 2\*10^6 $ genes
with possible partition function sizes of 10^20

It is very interesting. we can see here that there is multiple clusters, independantly of the homology or the length. the 9 clearly defined clusters are representing some other notion. slightly related with the species but not enough. We can postulate that the clusters define some form of shared constraint neither related to homology nor to species. 


the averages are
homology averages : [ 0.94144441  0.89892266  0.76866767  0.76787882  0.59010412  0.76331989
      0.77705589  0.94065878  0.71435589  0.93837663  0.92280362  0.78456446
      0.76631687  0.91948455  0.9116723   0.93393733  0.75076303  0.94312841]

var : [ 0.2237579 ,  0.28926838,  0.31572851,  0.31721181,  0.37081365,
            0.31820945,  0.31420915,  0.22549738,  0.33437126,  0.21387471,
            0.24471889,  0.31170976,  0.31744569,  0.26116283,  0.26529497,
            0.23693527,  0.32200822,  0.2209886 ]

Right of the bat we can see that some amino acids are very litlle pushed to an unusual entorpy value. This could be explained by some reasons regarding the general role of this amino acid. But it relatively vary between species. So there is some other factors.

THis is interesting that this homology has little to do with anything and it should not be a secondary cluster However, we can see that tSNE shows a non linear secondary cluster.
For the seconda one however, there is this notion of clusters that is even more present. Clusters show much more using $$ A_{value}$$
For most of them, there seems to be an equal correlation of the value to nans, cai length and GC whether it is high or low, which also hints on the idea that it is not correlation to the measure in itself.

mRNA abundance seems also to be linked with the secondary cluster presence/outlier cluster presence

Looking at the plot of many homologies, we kinda see different small groups and bigger groups within them. It resembles a lot a small world graph...
If one would overlay previous knowledge about the CUB it seems that different values cluster show the different drivers of the CUB, it can be about the general protein function (with a cluster of one homology). It can be related to the species environment with small clusters and or usage of the protein. 

------
we also find, looking at the plot, how the entropy cost allows to cluster the datapoints much better. however, it appears that it does not show simple correlation amongst variables.  the very good cluster found, do not even relate to anything of phylogenetic distances. some other latent variables might be at play here. we can see that the entropy cost contains much more scrambled information about indirect factors and latent variables but also remove some not important information. this can be shown by the way the NN is able to infer it from different values


###secondary
YMR186W

imilarity_scores: SpearmanrResult(correlation=-0.15984702808231441, pvalue=0.00016691732866527833)
    nans: SpearmanrResult(correlation=-0.10446836580972899, pvalue=0.014240595035458457)
    ecai: SpearmanrResult(correlation=0.021057431026088324, pvalue=0.62217395694803812)
    cai: SpearmanrResult(correlation=0.041378414733043593, pvalue=0.33273634999554524)
    lenmat: SpearmanrResult(correlation=-0.059038206659138968, pvalue=0.16677993138434313)
    GCcount: SpearmanrResult(correlation=-0.06612591732312488, pvalue=0.12139189160360152)


silhouette: 0.0920315361753
    cal_hara: 28.3819986672
    cluster_phylodistance: [1.0064131983336058, 1.1363601011083024, 1.003642183120282, 
     0.9544638941877908, 1.0158370600816644]

Here it tends to be mostly the lengths and number of NaN values that seem to be corelated with the entropy values. If they are correlated to the first cluster, they don't explain all the other ones.

Compared to other metrics, it seems to be for now the best one for clustering and finding clusters etc..

### outliers

YBR089C-A

similarity_scores: SpearmanrResult(correlation=-0.056817674870850376, pvalue=0.20836000240005659)
    nans: SpearmanrResult(correlation=-0.39803122016203385, pvalue=3.9754782355307412e-20)
    ecai: SpearmanrResult(correlation=0.16328502680745435, pvalue=0.00027576475430247596)
    cai: SpearmanrResult(correlation=-0.089018302220878379, pvalue=0.048447932653799514)
    lenmat: SpearmanrResult(correlation=0.43207324673915537, pvalue=8.5192368836702005e-24)
    GCcount: SpearmanrResult(correlation=0.037584279363295073, pvalue=0.40550318325433032)

 silhouette: 0.11192530708
    cal_hara: 23.206754635

We can really see that there is differences amongst homologies with a lot of secondary clusters and others here. Moreover we can see in the full plot of mulitple homologies that an


###primary
YLR142W

similarity_scores: SpearmanrResult(correlation=0.0371730540226008, pvalue=0.28592102682190379)
    nans: SpearmanrResult(correlation=-0.10419323606742718, pvalue=0.0027161042246365901)
    ecai: SpearmanrResult(correlation=-0.054798878661751413, pvalue=0.11555007679193391)
    cai: SpearmanrResult(correlation=-0.16789221953558919, pvalue=1.2195465121695145e-06)
    lenmat: SpearmanrResult(correlation=0.1112294842727367, pvalue=0.0013654399651707328)
    GCcount: SpearmanrResult(correlation=0.17078327176356722, pvalue=7.9239682671512234e-07)

YMR121C

similarity_scores: SpearmanrResult(correlation=0.045255382799603602, pvalue=0.31011147342294748)
    nans: SpearmanrResult(correlation=-0.080488084645829738, pvalue=0.070731805083924348)
    ecai: SpearmanrResult(correlation=-0.015135724467302483, pvalue=0.7343764607186829)
    cai: SpearmanrResult(correlation=-0.051478028866851711, pvalue=0.24820353133591458)
    lenmat: SpearmanrResult(correlation=0.14009461817731003, pvalue=0.0015994284596261669)
    GCcount: SpearmanrResult(correlation=0.075505049891718148, pvalue=0.090079126599269893)

silhouette: 0.223301491637
    cal_hara: 13.2400501282



### final entloc


The species with the highest variance seems to be all surrounding the one with less variance in their CUB, species satoella_complicata seems to be far outside the regular cluster. In PCA there really seems to be an elongated shape meaning that one dimension contains way more variance than any other. 
Plant symbiotic species seem to cluster more toward low entropy, low variance values.
Only one part of the plot seems to have brown rot bacteria but this is of low statistical significance.


it really seems to be difficult to regress on that low number of datapoint, we achieve a low score, which may be the maximum given the probably low correlation.
We can see that the average size correlates with the CUB which makes sense given the measure we use. Moreover the GCcount seems to be really important for the regressor algorithm. finally, we have a small inverse correlation of the CUB to the total number of homology (maybe that the highly shared/conserved genes tend to be ones which codon usage bias is little influenced by exterior factors. 

average_size: SpearmanrResult(correlation=0.2441321746988882, pvalue=1.6034758413378369e-08)
    num_genes: SpearmanrResult(correlation=0.025754952030612056, pvalue=0.55712420802658313)
    genome_size: SpearmanrResult(correlation=0.02465102393132142, pvalue=0.57415346354289198)
    tot_homologies: SpearmanrResult(correlation=-0.15643674596871962, pvalue=0.00033360229590838011)

score of Lasso regression: 0.2555102083452766
highest possible

-------------
What a strangely well defined shape. we can see a cluster of high mean eCAI homologies: a form of preservation which is also a place where the mean entropy is the lowest.

Homologies with a high cost/length  also seem to cluster together around one point of high entropy (but not the highest). 

a correlated thing is the number of cys regulatory elements to the length and cost of the gene. 

secreted homologies seem to sometimes be more clustered than non secreted ones. 

High mRNA abundances seems to be mostly clustered on the outside of the plot. we can also see just how many. 

Some low entropy homologies also seem to be driven by a lot of NaN values. 

similarity score seems to be influenced by the length (which is not trivial) And there also seems to be some clusters of higher similarity and the one of low entorpy highlited earlier. 

Moreover, it is the one that has been left out as outlier by DBSCAN. 

Most homologies seems to have high entropy and the others seem only to be a few (may be a lot of outliers due to all homoset

we can see an important cluster in PCA displaying the influence of hydrophobicity on the homologies.

* it is interesting that despite looking at averages over many different genes of different species, such relation exist despite looking mostly at the metadata of only SC*

#### DOING working
The shape is now way more rounded than before. and some of the clusters we had spotted earlier have disapeared. we have mostly the same correlation as before




---------
The results are pretty grim for this regression algorithm, there is traces of hope however, some important correlation have been found and we can see an almost linear increase of the r^2 score with the increase of datapoint.

    similarity_scores = -0.072817143336051104, pvalue=2.7742966040990751e-08
    nans = -0.67673563479785559, pvalue=0.0
    lenmat = 0.61317143434654797, pvalue=0.0
    GCcount = 0.1042831576169921, pvalue=1.6286476139846346e-15
    weight = 0.59488263810025221, pvalue=0.0
    protein_abundance = 0.056439099440282671, pvalue=1.6821875511060589e-05
    mRNA_abundance = -0.086933272838433981, pvalue=3.2280434415417848e-11
    decay_rate = 0.1567485306934171, pvalue=2.9071670461463889e-33
    is_secreted = 0.0013553018355325524, pvalue=0.91775916604867347
    cys_elements = 0.077968371827372554, pvalue=2.6936563465393485e-09
    tot_volume = 0.61404859129722777, pvalue=0.0
    mean_hydrophobicity = 0.12100588754115996, pvalue=2.1802017159482846e-20
    glucose_cost = 0.61905446084924243, pvalue=0.0
    synthesis_steps = 0.62413262866369279, pvalue=0.0
    isoelectricpoint = 0.020678172263525068, pvalue=0.11512218340330567
    meanecai = -0.31898066567312289, pvalue=1.7902743204999836e-137
    meancai = -0.11200359394823546, pvalue=1.1267798927103194e-17
    conservation = 0.091562523649363861, pvalue=2.7373927220031136e-12

the R^2 score is of: 0.632627415777

Moreover, using only a subset of values that are correlated with the entropy and only using the average CUB, a 2 layer (2\*11 neurons) MLP is able to regress on the values with a very good accuracy (0.6 < r^2 < 0.76) this is because of nans values (removing them get's you back to a r^2 of (0.2 - 0.35)

*with have seen that, for regression, the entropy loc seems to be better than regular entropy. Moreover, the NN is also way better than the lasso. removing uncorrelated parameters is also very important*

we might be able to have a better look at the influence of the nans by plotting the error of regression against the nan values with and without the nan values. to see how they influence the regression

## getting the frequencies:

we load as usual.

Silhouette: The best value is 1 and the worst value is -1. Values near 0 indicate overlapping clusters. Negative values generally indicate that a sample has been assigned to the wrong cluster, as a different cluster is more similar.

Harabaz: The score is defined as ratio between the within-cluster dispersion and the between-cluster dispersion

AIC and BIC are bayesian information criterions. one for the parameters that we use ( to remove some) and the other for the model that we use to compare it and choose the best one


-----

Well what we had here for the homologies: 
we can see that there is some very interesting plot where only a few lines are describing hte CUB of some homologies and the rest is clustering around the center. what can explain that ?

However, not only one or two CUB variance can explain that as shown in the mean variances of this homology
for example in   ['YDL224C', 'YNL315C', 'YDL183C', 'YCL039W', 'YER140W', 'YOR386W'] using PCA


MEAN 

[ 0.21637936  0.32080878  0.18636864  0.27644322  0.14272234  0.24840703
      0.11710233  0.16635646  0.12342153  0.20199031  0.61637185  0.38362815
      0.4761617   0.5238383   0.42164512  0.57835488  0.44717273  0.55282727
      0.56648396  0.43351604  0.29441138  0.13370346  0.21963885  0.35224631
      0.43861137  0.56138863  0.51523569  0.13789745  0.34686686  0.08945905
      0.19547382  0.26863177  0.1617783   0.20551245  0.07914461  0.34062367
      0.65937633  0.37764455  0.62235545  0.26229166  0.1973185   0.25675042
      0.28363942  0.18351299  0.17407077  0.14152359  0.21027929  0.17828906
      0.1123243   0.32952318  0.2206005   0.21471063  0.23516569  0.38095949
      0.61904051  0.11534297  0.37655913  0.24581345  0.26228444]

VAR [ 0.1313701 ,  0.15693788,  0.12874605,  0.14560639,  0.11901046,
            0.18630218,  0.10431298,  0.13361226,  0.10702646,  0.20899046,
            0.22540011,  0.22540011,  0.20512176,  0.20512176,  0.29797907,
            0.29797907,  0.2396791 ,  0.2396791 ,  0.22084201,  0.22084201,
            0.17983354,  0.11111848,  0.13170709,  0.1866077 ,  0.26102624,
            0.26102624,  0.20913771,  0.1333468 ,  0.16361986,  0.12796655,
            0.12882443,  0.16626248,  0.0938507 ,  0.13063002,  0.06612581,
            0.21331124,  0.21331124,  0.21221855,  0.21221855,  0.14775439,
            0.14652639,  0.18250794,  0.18185308,  0.11566318,  0.11482386,
            0.10416015,  0.12572926,  0.10775521,  0.08494582,  0.17407936,
            0.14023019,  0.15193334,  0.14537714,  0.23516215,  0.23516215,
            0.10099699,  0.18069381,  0.13775152,  0.15272984])

## homoplots

e can see for this one, how much the cai is linked with the similarity , the GC content and the frequency of codon usage

there is also a strong link of ecai to CUF. we can see that there is very strong links to 
 we can see some different homologies that have much more variances and thus could not be clustered together. For them, the GC content explains entirely the 2D plot that we have
 
 We can strongly see that most of the time, the manifold of data that is retrieved by tSNE is exactly the one of the GC content or eCAI or CAI etc..

*Here is correlation is against the Endres Schrindelin distance of the CUF of each gene to the reference gene*

*it seems, thoughout this pipeline, using frequency, that it does not merely show as much information that the other ones. this can be explained because the pipeline was meant thinking about the types of general information that entropy could give: it both validate the effiency of entropy as a measure and validate also this pipeline and this project in its goals, intuitions etc*

'YEL022W'

    similarity_scores =-0.32206253140336, pvalue=3.0647714617907084e-13
    ecai =-0.89114450221958186, pvalue=6.4115801870425671e-169
    cai =-0.88864101235781234, pvalue=1.1694596527175149e-166
    lenmat =0.26560276040697112, pvalue=2.5272798170047021e-09
    GCcount =0.88227074022528273, pvalue=3.8477272520583511e-161

    silhouette: 0.33600633148
    cal_hara: 6.76887799758
    cluster_phylodistance: [0.9032958624032591, 1.0108045214673116]



ecai = -0.94665512849427791, pvalue=5.8145543274871832e-254
cai = -0.93119502023189893, pvalue=1.513910210748799e-226
lenmat = 0.12333150654137384, pvalue=0.0051101342158700819
GCcount = 0.93139784774252099, pvalue=7.302251468645603e-227

homology: YLR024C
------------------------------------
silhouette: 0.339863129829
cal_hara: 3.02683639626
cluster_phylodistance: [0.8684993690891437, 1.00296683546272]



We can really see clusters here. they certainly mostly represent subtypes of species that are highly related. this makes that every metadata creates clusters which are in fact already highly related species 

Again here, we see full well the difference between the two groups of homologies. GC content, CAI and eCAI are high determinant of this groups.
Here a contrario to the other measures, we can see a cluster from all genes with high mRNA abundances. It can also be seen using PCA.
Moreover, high similarity homologies seemed to be all on the outside of the second cluster. working homoset is basically just one of them (the first one)


## compute 3DGD:

we found very different values:

computing with the homologous similarity led to interesting results.