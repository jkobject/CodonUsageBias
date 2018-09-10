# Remarks about the code and what I did there


## utils

- because of hyper optimised C function it might is faster to run some unoptimized function relying heavily on C libraries (pandas, numpy, scipy...) than to try to optimize one's algorithm in python (needs to be tested with %%timeit)

- the multivariate normal approximation to the multinomial distribution does not seem to work if I provide a set of probabilities that are the same for each draw dimension because the matrix is not invertible since it is not SDP. SO the function might only work given non-symetrical probailities

- computejerem is faster than any other partition function and stays faster than even the randomdraw until reaching approx 30 000 000 possibilities. it uses the repetition in the multinomial computation to bypass the multinomial of scipy and directly do dynamic programming on it

- trying to compute the complexity of the regular algorithm with for loops of for loops, I did not managed to do it for nbcod = 6 (see paper document on the subject)

- computepartition_sorted_full is basically the fullest and simplest computation, it will compute all possible values in a sorted manner, one by one with groups of for-loops (almost like yun's algorithm)

-computepartition_without_permutation is supposed to be much more efficient as it does not compute the values that will produce the same output value ( any permutation of the values within the random vector) and compute just the number of permutation of a given vector to get the actual output

- the random draw uses heavily dynamic programming and coding tricks to speed itself up. one of the trick is to 

- overall, the compute partition functions have allowed me to explore HPC and programmation issue when doing such highly repetitive computation on python as well as the many tricks to improve performance and speed. A great exploration of computational complexity

-getloc is really about recoding what someone wrote and enabled me to really fully understand the code of Yun Deng, it uses also heavily dynamic programming


- because of the way I do my save, basically storing everything into a dict and then a json, I have quite a few problem (typing issues, memory needed (approx 5times more than the final output file..))

- even using python's ability to display foerever bigger and smaller number, using numpy may force some very very small values 10^-100 to be zero,for example when drawing from the multinomial

- I use many times, the lnexp trick to have faster and more reliable computations, however this non linearity trick could have also been used in jerem-compute (and assessed the optimization) and in other context such as when displaying the color gradients for different values to better highlight the differences)

- the architecture used makes it pretty straightforward to use and allows for many developments however it might be complicated to perform some things and is only one amongst many possible I have thought about (see papers on it)

- I use simple grid search on different occasion to continue on the process of not infering anything about the data

## other

- the saving function is really not the most efficient and shows how complicated it might be to save some set of data using python. (I had thought about serializing the data (like pickle) but this is not possible on python as it is on other languages). I am using gzip, json and csv formats for the storage

- I thought about using hdfs to do some large scale computation. This would allow my code to be used on more complex data sets than fungi (like vertebrate or plants for example. However the time cost of implementing it would have been really high)

- tSNE, with its requirement to look at the full data distribution, is very memory ineficient, I used 3 different version of it, all using barn's huts SGD and a parallelized one seemed to be more than 10x faster than the regular scikit one (sota)

## homoset

- I found very usefull to use heatmaps and datashaders to plot point clouds that are too big to be displayed with regular plots, using real cloud of points.
heatmaps with binning can be very usefull as well

- I added random homologies to see if the plot would make the same kind of sense and it did not.

- I often use displaying where the computation is is really usefull to take decisions about the computation and the code being performed. 

## homology

- I was abl to plot 4 dimensions using 3d plot+ size. works ok.

AIC and BIC measures or silhouette/CH scores are not really usefull in the context I am in. I don't really know what a good cluster is and what a good amount of parameters are, and this might change accross homologies 

## pycub

parallelizing to a high extent (up to 16 processes when doing parallel loading of the data.) It is very efficient and increases the speed almost linearly as long as there is internet speed.

- using 3 different ways to load datas from internet 

- the data is almost always not formated, very hard to process. need to use techniques to solve problems of namespaces and references..

- using R within python with Rpy2, 

- 4 big functions to finally compares species and homologies together with all the data, using 2 distinct techniques, a visual one with plottings and interactive color gradients and another with a regression. both together should be enough to investigate the explanation power of different variables to the codon usage bias.

- always being carefull for both to use different methods one with a more simple/consistent/used approach and another more  to investigate 

- it has been very difficult to reproduce diament work as it has just been explained orally throught a few number of line without any availability of pseudo code or other in the supplementary material or the supplementary documents, very hard to find the data as well as it was not lined and some of the cited work which should have corresponded to the data did not linked the data they had either. I have used 4 different versions for the code and runned it for a week before seeing similar results Moreover, each dataset is different so hard to generalize code. used masked arrays to allow fast computation over arrays.

- I have found a search method that is very pythonic and related to what I was looking for.

#assumptions
- the codons are supposed to be geometrically randomly distributed,

- CAI can work accross species

- for most of the metadata on the underlying protein, the values are roughly the same for each proteins in the homology. for the costs, volume etc. the simplification is to the point of only adding the underlying values of each amino acids for the sequence making the protein.

- assuming that dividing perc_id to perc_pos will give me the similarity of a given sequence to the reference one

- Assuming that databases/Tobias/Diament data is accurate  

- I have assumed that I can basically replace the value representing unknown 

- I have assumed that the pseudo phylogenetic distances I create from the phylogenetic dendogram tree are somewhat related to an actual concept of evolutionary distance.

- I am assuming that the amount of tRNA genomic position even if uninformative about the actual tRNA is accurate in itself. 

#statistics about the code

- the code is about 6500 lines of pure code 
- it is 90% python, 5% javascript, 3% bash, 2% R and it contains 13 additional python packages and 2 additional R packages

- using the dataset from NCBI, ensembl, from labs : Von der Harr, Diament, 

- using json/csv for data storage
- using gzip for compression
- using aws for computing
- using conda and pip for managing packages in python and homebrew, yum for managing other packages,
- using pipreqs for maanging requirements file, git/github for the file versionning, with gitignore.io, doxypypy + doxygen for the documentation,
github pages for the website
-using stackoverflow, python, developpez, quora, github, biostar, ...

- development pipeline, linux ubuntu + macOS, Iterm, sublime + many packages, rsub + remote jupyter for on-server development. jupyter notebook.

- using:
 os, json, zipfile, shutil, ftplib, gzip, copy, requests, sklearn, urllib2, joblib, functools32, rpy2, ete2, pandas, numpy, scipy, bokeh
matplotlib, Bio, pdb, glob
math, holoviews, kmodes, MulticoreTSNE, collections, mpl_toolkits

8 functions over 25 are parallel

# Diament plots 
2.5million contacts retrieved on dist3D
500Mb matrices

endresdistance
DIAMENT2 (0.038297713917720327,
 0.086844486492188941, CUB ENT
 0.71438053837945337,
 4.0047867171581658e-312, CUF
 0.84295299789908806,
 0.0) SYNCUF

euclidean distance

DiAMENT2(-0.13840943669212902,
 5.1022318884020248e-10,
 0.59963739504768032,
 1.386413533191727e-195,
 0.75864680930450634,
 0.0)

0.44650853756252862
 (-0.2685301534522273, CUB HOMOLOGY
 2.236504987814927e-34,

 correlatio between endres and euclidean (correlation=0.96924282826462216, pvalue=0.0) better with 98%=cub > cuf > entropy = 94%

## intriguing stuff

- no shared homologies accross each species: 
- some strange forming in plot all homologies with tsne
- some homologies with really clear cut clusters
- some homologies with secondary clusters bigger than primary cluster 
- Diament plot
- one species with most of its values as outliers
- very high level of entropy for som


ENTLOC

ENTLOC - HOMO
similarity_scores: SpearmanrResult(correlation=-0.072817143336051104, pvalue=2.7742966040990751e-08)
nans: SpearmanrResult(correlation=-0.67673563479785559, pvalue=0.0)
lenmat: SpearmanrResult(correlation=0.61317143434654797, pvalue=0.0)
GCcount: SpearmanrResult(correlation=0.1042831576169921, pvalue=1.6286476139846346e-15)
weight: SpearmanrResult(correlation=0.59488263810025221, pvalue=0.0)
protein_abundance: SpearmanrResult(correlation=0.056439099440282671, pvalue=1.6821875511060589e-05)
mRNA_abundance: SpearmanrResult(correlation=-0.086933272838433981, pvalue=3.2280434415417848e-11)
decay_rate: SpearmanrResult(correlation=0.1567485306934171, pvalue=2.9071670461463889e-33)
is_secreted: SpearmanrResult(correlation=0.0013553018355325524, pvalue=0.91775916604867347)
cys_elements: SpearmanrResult(correlation=0.077968371827372554, pvalue=2.6936563465393485e-09)
tot_volume: SpearmanrResult(correlation=0.61404859129722777, pvalue=0.0)
mean_hydrophobicity: SpearmanrResult(correlation=0.12100588754115996, pvalue=2.1802017159482846e-20)
glucose_cost: SpearmanrResult(correlation=0.61905446084924243, pvalue=0.0)
synthesis_steps: SpearmanrResult(correlation=0.62413262866369279, pvalue=0.0)
isoelectricpoint: SpearmanrResult(correlation=0.020678172263525068, pvalue=0.11512218340330567)
meanecai: SpearmanrResult(correlation=-0.31898066567312289, pvalue=1.7902743204999836e-137)
meancai: SpearmanrResult(correlation=-0.11200359394823546, pvalue=1.1267798927103194e-17)
conservation: SpearmanrResult(correlation=0.091562523649363861, pvalue=2.7373927220031136e-12)


ENTLOC - SPECE
average_size: 0.2441321746988882, pvalue=1.6034758413378369e-08
num_genes: SpearmanrResult(correlation=0.025754952030612056, pvalue=0.55712420802658313)
genome_size: SpearmanrResult(correlation=0.02465102393132142, pvalue=0.57415346354289198)
tot_homologies: SpearmanrResult(correlation=-0.15643674596871962, pvalue=0.00033360229590838011)

HOMOPLOT ENTLOC
nans: SpearmanrResult(correlation=-0.080488084645829738, pvalue=0.070731805083924348)
ecai: SpearmanrResult(correlation=-0.015135724467302483, pvalue=0.7343764607186829)
cai: SpearmanrResult(correlation=-0.051478028866851711, pvalue=0.24820353133591458)
lenmat: SpearmanrResult(correlation=0.14009461817731003, pvalue=0.0015994284596261669)
GCcount: SpearmanrResult(correlation=0.075505049891718148, pvalue=0.090079126599269893)

similarity_scores: SpearmanrResult(correlation=0.37483760198015076, pvalue=7.7121359818286199e-25)
nans: SpearmanrResult(correlation=-0.10806248167061397, pvalue=0.0041511485637422685)
ecai: SpearmanrResult(correlation=0.010218256542274859, pvalue=0.78696010067632338)
cai: SpearmanrResult(correlation=-0.017079558552661085, pvalue=0.65144603963385839)
lenmat: SpearmanrResult(correlation=-0.15684365422195573, pvalue=2.9920382844967417e-05)
GCcount: SpearmanrResult(correlation=0.044066508556495902, pvalue=0.24359812284101062)


similarity_scores: SpearmanrResult(correlation=-0.056817674870850376, pvalue=0.20836000240005659)
nans: SpearmanrResult(correlation=-0.39803122016203385, pvalue=3.9754782355307412e-20)
ecai: SpearmanrResult(correlation=0.16328502680745435, pvalue=0.00027576475430247596)
cai: SpearmanrResult(correlation=-0.089018302220878379, pvalue=0.048447932653799514)
lenmat: SpearmanrResult(correlation=0.43207324673915537, pvalue=8.5192368836702005e-24)
GCcount: SpearmanrResult(correlation=0.037584279363295073, pvalue=0.40550318325433032)

similarity_scores: SpearmanrResult(correlation=-0.056817674870850376, pvalue=0.20836000240005659)
nans: SpearmanrResult(correlation=-0.39803122016203385, pvalue=3.9754782355307412e-20)
ecai: SpearmanrResult(correlation=0.16328502680745435, pvalue=0.00027576475430247596)
cai: SpearmanrResult(correlation=-0.089018302220878379, pvalue=0.048447932653799514)
lenmat: SpearmanrResult(correlation=0.43207324673915537, pvalue=8.5192368836702005e-24)
GCcount: SpearmanrResult(correlation=0.037584279363295073, pvalue=0.40550318325433032)

MEAN

$5.945 \pm 2.766$ & $0.945 \pm 0.223$ & $0.212\pm 0.121, 0.324\pm 0.157, 0.188\pm 0.129, 0.274\pm 0.145$ \\
$7.551 \pm 2.918$ & $0.907 \pm 0.289$ & $0.144\pm 0.119, 0.254\pm 0.185, 0.118\pm 0.104, 0.169\pm 0.135, 0.121\pm 0.105, 0.192\pm 0.202, $\\
$2.210 \pm 1.892$ & $0.747 \pm 0.315$ & $0.625 \pm 0.224, 0.374 \pm 0.224$\\
$2.118 \pm 1.963$ & $0.753 \pm 0.317$ & $0.470 \pm 0.203, 0.529 \pm 0.203$\\
$0.747 \pm 1.077$ & $0.532 \pm 0.370$ & $0.415 \pm 0.298, 0.584 \pm 0.298$\\
$1.996 \pm 1.838$ & $0.722 \pm 0.318$ & $0.439 \pm 0.238, 0.560 \pm 0.238$\\
$3.071 \pm 2.276$ & $0.750 \pm 0.314$ & $0.574 \pm 0.219, 0.425 \pm 0.219$\\
$6.802 \pm 2.883$ & $0.948 \pm 0.225$ & $ $
$1.173 \pm 1.361$ & $0.701 \pm 0.334$ & $ $
$5.884 \pm 2.555$ & $0.932 \pm 0.213$ & $ $
$12.639 \pm 3.693$ & $0.936 \pm 0.244$ & $ $
$4.184 \pm 2.525$ & $0.791 \pm 0.311$ & $ $
$1.999 \pm 1.743$ & $0.772 \pm 0.317$ & $ $
$3.985 \pm 2.246$ & $0.924 \pm 0.261$ & $ $
$6.067 \pm 2.604$ & $0.927 \pm 0.265$ & $ $
$4.285 \pm 2.267$ & $0.940 \pm 0.236$ & $ $
$1.636 \pm 1.592$ & $0.686 \pm 0.322$ & $ $
$6.399 \pm 2.668$ & $0.958 \pm 0.220$ & $ $

  ,  ,  ,  ,  ,
            ,  ,  ,  ,  ,
            ,  ,  ,  ,  ,
            ,  ,   




ENT

ENT HOMOLOGIES
similarity_scores: SpearmanrResult(correlation=-0.081261407786012513, pvalue=5.5983630443298877e-10)
nans: SpearmanrResult(correlation=-0.3555959460371706, pvalue=1.1441764655930061e-172)
lenmat: SpearmanrResult(correlation=0.80583957969973008, pvalue=0.0)
GCcount: SpearmanrResult(correlation=0.22619613607892744, pvalue=2.8877807326984728e-68)
weight: SpearmanrResult(correlation=0.72817624858641561, pvalue=0.0)
protein_abundance: SpearmanrResult(correlation=0.18651697015748744, pvalue=1.2999207161464825e-46)
mRNA_abundance: SpearmanrResult(correlation=0.050599719595574855, pvalue=0.00011454073869885887)
decay_rate: SpearmanrResult(correlation=0.14022360231537928, pvalue=6.8717537057815878e-27)
is_secreted: SpearmanrResult(correlation=0.090758808355216727, pvalue=4.2391813217099995e-12)
cys_elements: SpearmanrResult(correlation=0.12857575316090689, pvalue=7.8244087477063335e-23)
tot_volume: SpearmanrResult(correlation=0.77385917280841998, pvalue=0.0)
mean_hydrophobicity: SpearmanrResult(correlation=0.11672375791852342, pvalue=4.5226765754239889e-19)
glucose_cost: SpearmanrResult(correlation=0.77214059050694916, pvalue=0.0)
synthesis_steps: SpearmanrResult(correlation=0.77459143311298995, pvalue=0.0)
isoelectricpoint: SpearmanrResult(correlation=-0.026971626011661889, pvalue=0.039853086935354678)
meanecai: SpearmanrResult(correlation=-0.3734691190980633, pvalue=1.2365551211497936e-191)
meancai: SpearmanrResult(correlation=-0.1372190049285606, pvalue=8.2741500233278835e-26)
conservation: SpearmanrResult(correlation=0.080801201912223494, pvalue=6.9991869867232414e-10)


ENT SPECIES

average_size: SpearmanrResult(correlation=0.027770541559695276, pvalue=0.52667911319693816)
num_genes: SpearmanrResult(correlation=-0.20392120936213901, pvalue=2.635831327646218e-06)
genome_size: SpearmanrResult(correlation=-0.14812823660067045, pvalue=0.00068636971631414411)
fullGCcount: SpearmanrResult(correlation=0.26672770045972816, pvalue=5.9432859937195404e-10)
varGCcount: SpearmanrResult(correlation=0.18647347951225177, pvalue=1.8039311033672066e-05)
tot_homologies: SpearmanrResult(correlation=-0.016224109108492327, pvalue=0.71152173456985102)
phylogenetic distances: SpearmanrResult(correlation=-0.10299262826299072, pvalue=0.018586274765915109)


ENT HOMOPLOTS

nans: SpearmanrResult(correlation=-0.19590289969118865, pvalue=1.540731325834782e-05)
ecai: SpearmanrResult(correlation=0.35135518228446844, pvalue=2.155877607469682e-15)
cai: SpearmanrResult(correlation=-0.30619021569389493, pvalue=7.0738882403011127e-12)
lenmat: SpearmanrResult(correlation=0.10524978740516901, pvalue=0.021093649469112005)
GCcount: SpearmanrResult(correlation=0.40033350586522315, pvalue=6.6791437489435558e-20)

similarity_scores: SpearmanrResult(correlation=0.32963398880492711, pvalue=2.2941910765484818e-14)
nans: SpearmanrResult(correlation=-0.050120172803386744, pvalue=0.25902776911849618)
ecai: SpearmanrResult(correlation=0.02066270840375984, pvalue=0.64187740024601836)
cai: SpearmanrResult(correlation=-0.084657244755652339, pvalue=0.05630257313438284)
lenmat: SpearmanrResult(correlation=-0.45581051665020356, pvalue=1.7659741848818821e-27)
GCcount: SpearmanrResult(correlation=0.21769726581613691, pvalue=7.0858251668679864e-07)

similarity_scores: SpearmanrResult(correlation=0.15668458797295806, pvalue=0.00067900733673671384)
nans: SpearmanrResult(correlation=0.28586330406273419, pvalue=3.1180931951463311e-10)
ecai: SpearmanrResult(correlation=0.35815742616977347, pvalue=1.3990479861887904e-15)
cai: SpearmanrResult(correlation=-0.059886576258249898, pvalue=0.19640866101880053)
lenmat: SpearmanrResult(correlation=0.0051118698202435491, pvalue=0.91227208771407242)
GCcount: SpearmanrResult(correlation=0.1392534124048882, pvalue=0.0025618902813443235)

similarity_scores: SpearmanrResult(correlation=0.15296426892106091, pvalue=0.00056934708193189247)
nans: SpearmanrResult(correlation=-0.17110272856249881, pvalue=0.00011328878575171294)
ecai: SpearmanrResult(correlation=-0.15894624695793266, pvalue=0.00034061809294703635)
cai: SpearmanrResult(correlation=-0.3503784401533328, pvalue=5.284184358566423e-16)
lenmat: SpearmanrResult(correlation=0.1803169711666025, pvalue=4.6716433281590285e-05)
GCcount: SpearmanrResult(correlation=0.38903938597346588, pvalue=1.1677370420646218e-19)

similarity_scores: SpearmanrResult(correlation=0.43487331185717687, pvalue=2.6003821222295293e-25)
nans: SpearmanrResult(correlation=-0.27068816195697354, pvalue=3.7721552777484254e-10)
ecai: SpearmanrResult(correlation=0.43033486715175473, pvalue=9.1365917867182727e-25)
cai: SpearmanrResult(correlation=-0.0789140703800634, pvalue=0.072731142201338167)
lenmat: SpearmanrResult(correlation=0.026909476813379948, pvalue=0.54114636079113398)
GCcount: SpearmanrResult(correlation=0.1025444230874635, pvalue=0.019575315835816339)

similarity_scores: SpearmanrResult(correlation=0.44126966286818109, pvalue=9.3715047790922373e-28)
nans: SpearmanrResult(correlation=-0.26189260032207634, pvalue=3.9955072037803424e-10)
ecai: SpearmanrResult(correlation=0.16599490070839287, pvalue=8.7840155071909685e-05)
cai: SpearmanrResult(correlation=-0.12527002782781424, pvalue=0.0031694191946236841)
lenmat: SpearmanrResult(correlation=0.028970437037223521, pvalue=0.4965876805756344)
GCcount: SpearmanrResult(correlation=0.16289944019868202, pvalue=0.00011918414037122116)

MAX

[ 315.11608917,  327.41540819,  160.0414408 ,  487.89411216,
         52.94063673,  167.46844147,  179.16199446,  355.80320172,
         83.20523309,  174.17504316,  703.43954027,  163.36772378,
        135.52066175,  407.63142791,  317.63055617,  263.0384299 ,
         92.81086287,  349.44444376])

MEAN 
[  ,   ,   ,   ,
         ,   ,   ,   ,
         ,    ,  ,   ,
         ,   ,   ,   ,
          ,   ])

VAR

[  7.65164553,   8.51687922,   3.58228369,   3.8571531 ,
         1.16062417,   3.38140755,   5.18205945,   8.31587079,
         1.85351603,   6.52986744,  13.64420273,   6.37602431,
         3.04094619,   5.04495174,   6.78347536,   5.14261895,
         2.53622075,   7.12209085]





FREQ

SPECIES

average_size: SpearmanrResult(correlation=0.1014259282649264, pvalue=0.020463725722187187)
num_genes: SpearmanrResult(correlation=-0.16797255799919525, pvalue=0.00011527965531646965)
genome_size: SpearmanrResult(correlation=-0.041720350497750823, pvalue=0.34143823651282212)
tot_homologies: SpearmanrResult(correlation=-0.62331068834553016, pvalue=1.6136566554619786e-57)




HOMOPLOT

similarity_scores: SpearmanrResult(correlation=-0.29360862840299046, pvalue=1.7693945350444917e-11)
nans: SpearmanrResult(correlation=0.047406344872149661, pvalue=0.28813377148145641)
ecai: SpearmanrResult(correlation=-0.85158605786146635, pvalue=6.5235917202319086e-143)
cai: SpearmanrResult(correlation=-0.94311459661580255, pvalue=3.1261784365209552e-242)
lenmat: SpearmanrResult(correlation=0.2034507556660296, pvalue=4.1373692151291406e-06)
GCcount: SpearmanrResult(correlation=0.94977125336212787, pvalue=1.9888246307598616e-255)

similarity_scores: SpearmanrResult(correlation=-0.12832706625541768, pvalue=0.0078585998710924056)
nans: SpearmanrResult(correlation=0.20111973844421688, pvalue=2.7707361049821816e-05)
ecai: SpearmanrResult(correlation=-0.59910819942715443, pvalue=4.7692635428902217e-43)
cai: SpearmanrResult(correlation=-0.63888612323809668, pvalue=1.8079958777708946e-50)
lenmat: SpearmanrResult(correlation=-0.42055694436941071, pvalue=8.9854302603987605e-20)
GCcount: SpearmanrResult(correlation=0.53530088243316443, pvalue=4.1947927536458557e-33)

similarity_scores: SpearmanrResult(correlation=-0.028244575046229594, pvalue=0.59531530432540936)
nans: SpearmanrResult(correlation=0.017922305177279493, pvalue=0.73612053778495834)
ecai: SpearmanrResult(correlation=-0.6390267674459581, pvalue=2.9830988217421072e-42)
cai: SpearmanrResult(correlation=-0.51202114899306939, pvalue=3.4929686071265405e-25)
lenmat: SpearmanrResult(correlation=-0.21624117601482046, pvalue=3.8808021129908748e-05)
GCcount: SpearmanrResult(correlation=0.51582667724233799, pvalue=1.3534639816636899e-25)

nans: SpearmanrResult(correlation=-0.043110406049348099, pvalue=0.34243736212629583)
ecai: SpearmanrResult(correlation=-0.8378660882737945, pvalue=1.336030500534991e-129)
cai: SpearmanrResult(correlation=-0.81904834384192315, pvalue=4.1328456095225179e-119)
lenmat: SpearmanrResult(correlation=-0.21268897585395116, pvalue=2.1812031434509663e-06)
GCcount: SpearmanrResult(correlation=0.74766999004906232, pvalue=2.8728014364658284e-88)

similarity_scores: SpearmanrResult(correlation=-0.19842552415786888, pvalue=1.2144098879373411e-05)
nans: SpearmanrResult(correlation=0.12363035934495531, pvalue=0.0067463317674184868)
ecai: SpearmanrResult(correlation=-0.80216057733610335, pvalue=6.8398817697200166e-109)
cai: SpearmanrResult(correlation=-0.81635835124569833, pvalue=8.4491517255428707e-116)
lenmat: SpearmanrResult(correlation=-0.1026970974754615, pvalue=0.024594746783326182)
GCcount: SpearmanrResult(correlation=0.80729545911710099, pvalue=2.5300046835900362e-111)



MEAN

[ 0.21253173,  0.32414871,  0.18897989,  0.27433966,  0.14459283,
        0.25434144,  0.11827204,  0.16901116,  0.12129166,  0.19249087,
        ,  ,  ,,  ,  ,
        0.29363671,  0.13169634,  0.21718323,  0.35748372,  0.42929203,
        0.57070797,  0.52437031,  0.13219167,  0.34343802,  0.08415495,
        0.19262564,  0.27633943,  0.16215862,  0.20782163,  0.07689973,
        0.32978873,  0.67021127,  0.36915724,  0.63084276,  0.26198427,
        0.19828697,  0.24977401,  0.28995476,  0.18238289,  0.17646958,
        0.13850166,  0.21347218,  0.17962491,  0.10954877,  0.33435679,
        0.21729486,  0.2172119 ,  0.23113645,  0.37296075,  0.62703925,
        0.11210939,  0.38441595,  0.24511579,  0.25835886]

VAR 

[  ,  ,  ,  ,  ,
        ,  ,  ,  ,  ,
                 0.17953086,  0.11056564,  0.13106786,  0.18632887,  0.26021176,
        0.26021176,  0.20706937,  0.12986139,  0.16349044,  0.12353504,
        0.12888701,  0.16655317,  0.09480435,  0.13062479,  0.06488233,
        0.2108924 ,  0.2108924 ,  0.21167717,  0.21167717,  0.1484071 ,
        0.14718296,  0.18022325,  0.18390385,  0.1168238 ,  0.11685122,
        0.10289038,  0.12800847,  0.10842605,  0.08409559,  0.17581536,
        0.13988518,  0.15354668,  0.14469953,  0.23506393,  0.23506393,
        0.09853923,  0.17982803,  0.13750459,  0.15121574]


TRNA EXAMPLE
galerina_marginata_cbs_339_88
{u'ALA': {u'CGA': 6, u'CGC': 4, u'CGG': 2, u'CGU': 4},
 u'ARG': {u'GCA': 5, u'GCC': 3, u'GCG': 2, u'GCU': 4, u'UCC': 3, u'UCU': 3},
 u'ASN': {u'UUA': 2, u'UUG': 5},
 u'ASP': {u'CUA': 2, u'CUG': 7},
 u'CYS': {u'ACA': 2, u'ACG': 4},
 u'GLN': {u'GUC': 4, u'GUU': 5},
 u'GLU': {u'CUC': 6, u'CUU': 5},
 u'GLY': {u'CCA': 2, u'CCC': 4, u'CCG': 7, u'CCU': 5},
 u'HIS': {u'GUA': 2, u'GUG': 5},
 u'ILE': {u'UAA': 16, u'UAG': 2, u'UAU': 3},
 u'LEU': {u'AAC': 4, u'AAU': 4, u'GAA': 6, u'GAC': 4, u'GAG': 2, u'GAU': 3},
 u'LYS': {u'UUC': 6, u'UUU': 5},
 u'PHE': {u'AAA': 2, u'AAG': 5},
 u'PRO': {u'GGA': 5, u'GGC': 3, u'GGG': 2, u'GGU': 3},
 u'SER': {u'AGA': 5, u'AGC': 6, u'AGG': 2, u'AGU': 10, u'UCA': 2, u'UCG': 5},
 u'THR': {u'UGA': 6, u'UGC': 3, u'UGG': 2, u'UGU': 3},
 u'TYR': {u'AUA': 3, u'AUG': 5},
 u'VAL': {u'CAA': 6, u'CAC': 4, u'CAG': 2, u'CAU': 3},
 u'datapoints': 292,
 u'num': 127}

 SPECIES

 {0: u'_candida_glabrata',
     1: u'_candida_intermedia',
     2: u'_candida_tenuis_atcc_10573',
     3: u'absidia_glauca',
     4: u'acidomyces_richmondensis_bfw',
     5: u'acremonium_chrysogenum_atcc_11550',
     6: u'agaricus_bisporus_var_burnettii_jb137_s8',
     7: u'allomyces_macrogynus_atcc_38327',
     8: u'alternaria_alternata',
     9: u'amanita_muscaria_koide_bx008',
     10: u'amphiamblys_sp_wsbs2006',
     11: u'anthracocystis_flocculosa_pf_1',
     12: u'arthrobotrys_oligospora_atcc_24927',
     13: u'arthroderma_otae_cbs_113480',
     14: u'aschersonia_aleyrodis_rcef_2490',
     15: u'ascochyta_rabiei',
     16: u'ascoidea_rubescens_dsm_1968',
     17: u'ascosphaera_apis_arsef_7405',
     18: u'aspergillus_aculeatus_atcc_16872',
     19: u'aspergillus_bombycis',
     20: u'aspergillus_brasiliensis_cbs_101740',
     21: u'aspergillus_calidoustus',
     22: u'aspergillus_carbonarius_item_5010',
     23: u'aspergillus_clavatus',
     24: u'aspergillus_cristatus',
     25: u'aspergillus_flavus',
     26: u'aspergillus_fumigatus',
     27: u'aspergillus_glaucus_cbs_516_65',
     28: u'aspergillus_lentulus',
     29: u'aspergillus_luchuensis',
     30: u'aspergillus_niger_atcc_1015',
     31: u'aspergillus_nomius_nrrl_13137',
     32: u'aspergillus_ochraceoroseus',
     33: u'aspergillus_oryzae',
     34: u'aspergillus_parasiticus_su_1',
     35: u'aspergillus_rambellii',
     36: u'aspergillus_ruber_cbs_135680',
     37: u'aspergillus_sydowii_cbs_593_65',
     38: u'aspergillus_terreus',
     39: u'aspergillus_tubingensis_cbs_134_48',
     40: u'aspergillus_udagawae',
     41: u'aspergillus_versicolor_cbs_583_65',
     42: u'aspergillus_wentii_dto_134e9',
     43: u'aureobasidium_melanogenum_cbs_110374',
     44: u'aureobasidium_namibiae_cbs_147_97',
     45: u'aureobasidium_pullulans_exf_150',
     46: u'aureobasidium_subglaciale_exf_2481',
     47: u'babjeviella_inositovora_nrrl_y_12698',
     48: u'batrachochytrium_salamandrivorans',
     49: u'baudoinia_panamericana_uamh_10762',
     50: u'beauveria_bassiana',
     51: u'bipolaris_maydis_c5',
     52: u'bipolaris_oryzae_atcc_44560',
     53: u'bipolaris_sorokiniana_nd90pr',
     54: u'bipolaris_victoriae_fi3',
     55: u'bipolaris_zeicola_26_r_13',
     56: u'blastomyces_dermatitidis_er_3',
     57: u'blastomyces_gilchristii_slh14081',
     58: u'blastomyces_percursus',
     59: u'blumeria_graminis',
     60: u'botryobasidium_botryosum_fd_172_ss1',
     61: u'botrytis_cinerea',
     62: u'brettanomyces_bruxellensis_awri1499',
     63: u'byssochlamys_spectabilis_no_5',
     64: u'calocera_cornea_hhb12733',
     65: u'calocera_viscosa_tufc12733',
     66: u'candida_albicans_sc5314',
     67: u'candida_arabinofermentans_nrrl_yb_2248',
     68: u'candida_dubliniensis_cd36',
     69: u'candida_maltosa_xu316',
     70: u'candida_orthopsilosis_co_90_125',
     71: u'candida_parapsilosis_cdc317',
     72: u'candida_tanzawaensis_nrrl_y_17324',
     73: u'capronia_coronata_cbs_617_96',
     74: u'capronia_epimyces_cbs_606_96',
     75: u'ceraceosorus_bombacis',
     76: u'ceratocystis_platani',
     77: u'chaetomium_globosum_cbs_148_51',
     78: u'chaetomium_thermophilum_var_thermophilum_dsm_1495',
     79: u'choanephora_cucurbitarum',
     80: u'cladophialophora_bantiana_cbs_173_52',
     81: u'cladophialophora_carrionii',
     82: u'cladophialophora_immunda',
     83: u'cladophialophora_psammophila_cbs_110553',
     84: u'cladophialophora_yegresii_cbs_114405',
     85: u'claviceps_purpurea_20_1',
     86: u'coccidioides_immitis_h538_4',
     87: u'coccidioides_posadasii_str_silveira',
     88: u'colletotrichum_fioriniae_pj7',
     89: u'colletotrichum_gloeosporioides',
     90: u'colletotrichum_graminicola',
     91: u'colletotrichum_incanum',
     92: u'colletotrichum_nymphaeae_sa_01',
     93: u'colletotrichum_orbiculare',
     94: u'colletotrichum_orchidophilum',
     95: u'colletotrichum_salicis',
     96: u'colletotrichum_simmondsii',
     97: u'colletotrichum_sublineola',
     98: u'colletotrichum_tofieldiae',
     99: u'conidiobolus_coronatus_nrrl_28638',
     100: u'coniochaeta_ligniaria_nrrl_30616',
     101: u'coniosporium_apollinis_cbs_100218',
     102: u'coprinopsis_cinerea_okayama7_130',
     103: u'cordyceps_brongniartii_rcef_3172',
     104: u'cordyceps_confragosa_rcef_1005',
     105: u'cordyceps_militaris_cm01',
     106: u'cryptococcus_amylolentus_cbs_6039',
     107: u'cryptococcus_depauperatus_cbs_7841',
     108: u'cryptococcus_gattii_ca1873',
     109: u'cryptococcus_gattii_vgii_r265',
     110: u'cryptococcus_gattii_vgiv_ind107',
     111: u'cryptococcus_gattii_wm276',
     112: u'cryptococcus_neoformans',
     113: u'cutaneotrichosporon_oleaginosus',
     114: u'cyberlindnera_fabianii',
     115: u'cyberlindnera_jadinii',
     116: u'cylindrobasidium_torrendii_fp15055_ss_10',
     117: u'cyphellophora_europaea_cbs_101466',
     118: u'dacryopinax_primogenitus',
     119: u'dactylellina_haptotyla_cbs_200_50',
     120: u'daedalea_quercina_l_15889',
     121: u'debaryomyces_fabryi',
     122: u'debaryomyces_hansenii_cbs767',
     123: u'diaporthe_ampelina',
     124: u'diaporthe_helianthi',
     125: u'diplodia_corticola',
     126: u'diplodia_seriata',
     127: u'dothistroma_septosporum',
     128: u'drechmeria_coniospora',
     129: u'drechslerella_stenobrocha_248',
     130: u'emergomyces_pasteuriana_ep9510',
     131: u'emmonsia_crescens_uamh_3008',
     132: u'emmonsia_parva_uamh_139',
     133: u'emmonsia_sp_cac_2015a',
     134: u'endocarpon_pusillum_z07020',
     135: u'eremothecium_cymbalariae_dbvpg_7215',
     136: u'eremothecium_gossypii_fdag1',
     137: u'eremothecium_sinecaudum',
     138: u'erysiphe_necator',
     139: u'escovopsis_weberi',
     140: u'eutypa_lata_ucrel1',
     141: u'exidia_glandulosa_hhb12029',
     142: u'exophiala_aquamarina_cbs_119918',
     143: u'exophiala_dermatitidis_nih_ut8656',
     144: u'exophiala_mesophila',
     145: u'exophiala_oligosperma',
     146: u'exophiala_sideris',
     147: u'exophiala_spinifera',
     148: u'exophiala_xenobiotica',
     149: u'fibroporia_radiculosa',
     150: u'fibulorhizoctonia_sp_cbs_109695',
     151: u'fistulina_hepatica_atcc_64428',
     152: u'fomitopsis_pinicola_fp_58527_ss1',
     153: u'fonsecaea_erecta',
     154: u'fonsecaea_multimorphosa_cbs_102226',
     155: u'fonsecaea_nubica',
     156: u'fonsecaea_pedrosoi_cbs_271_37',
     157: u'fusarium_culmorum',
     158: u'fusarium_fujikuroi',
     159: u'fusarium_graminearum',
     160: u'fusarium_langsethiae',
     161: u'fusarium_mangiferae',
     162: u'fusarium_oxysporum',
     163: u'fusarium_poae',
     164: u'fusarium_proliferatum_et1',
     165: u'fusarium_pseudograminearum',
     166: u'fusarium_solani',
     167: u'fusarium_verticillioides',
     168: u'gaeumannomyces_graminis',
     169: u'galerina_marginata_cbs_339_88',
     170: u'gelatoporia_subvermispora_b',
     171: u'geotrichum_candidum',
     172: u'gloeophyllum_trabeum_atcc_11539',
     173: u'gonapodya_prolifera_jel478',
     174: u'grosmannia_clavigera_kw1407',
     175: u'hanseniaspora_guilliermondii',
     176: u'hanseniaspora_opuntiae',
     177: u'hanseniaspora_osmophila',
     178: u'hanseniaspora_valbyensis_nrrl_y_1626',
     179: u'hebeloma_cylindrosporum_h7',
     180: u'heterobasidion_irregulare_tc_32_1',
     181: u'hirsutella_minnesotensis_3608',
     182: u'histoplasma_capsulatum_nam1',
     183: u'hypholoma_sublateritium_fd_334_ss_4',
     184: u'hyphopichia_burtonii_nrrl_y_1933',
     185: u'hypsizygus_marmoreus',
     186: u'jaapia_argillacea_mucl_33604',
     187: u'kalmanozyma_brasiliensis_ghg001',
     188: u'kazachstania_africana_cbs_2517',
     189: u'kazachstania_naganishii_cbs_8797',
     190: u'kluyveromyces_lactis',
     191: u'kluyveromyces_marxianus_dmku3_1042',
     192: u'komagataella_pastoris',
     193: u'kuraishia_capsulata_cbs_1993',
     194: u'kwoniella_bestiolae_cbs_10118',
     195: u'kwoniella_dejecticola_cbs_10117',
     196: u'kwoniella_heveanensis_bcc8398',
     197: u'kwoniella_mangroviensis_cbs_10435',
     198: u'kwoniella_pini_cbs_10737',
     199: u'laccaria_amethystina_laam_08_1',
     200: u'laccaria_bicolor_s238n_h82',
     201: u'lachancea_dasiensis_cbs_10888',
     202: u'lachancea_fermentati',
     203: u'lachancea_lanzarotensis',
     204: u'lachancea_meyersii_cbs_8951',
     205: u'lachancea_mirantina',
     206: u'lachancea_nothofagi_cbs_11611',
     207: u'lachancea_sp_cbs_6924',
     208: u'lachancea_thermotolerans_cbs_6340',
     209: u'laetiporus_sulphureus_93_53',
     210: u'leptosphaeria_maculans',
     211: u'lichtheimia_ramosa',
     212: u'lipomyces_starkeyi_nrrl_y_11557',
     213: u'lodderomyces_elongisporus_nrrl_yb_4239',
     214: u'macrophomina_phaseolina_ms6',
     215: u'madurella_mycetomatis',
     216: u'magnaporthe_oryzae',
     217: u'magnaporthe_poae',
     218: u'malassezia_pachydermatis',
     219: u'malassezia_sympodialis_atcc_42132',
     220: u'marssonina_brunnea_f_sp_multigermtubi_mb_m1',
     221: u'melampsora_laricipopulina',
     222: u'metarhizium_acridum_cqma_102',
     223: u'metarhizium_album_arsef_1941',
     224: u'metarhizium_anisopliae_arsef_23',
     225: u'metarhizium_anisopliae_brip_53293',
     226: u'metarhizium_brunneum_arsef_3297',
     227: u'metarhizium_guizhouense_arsef_977',
     228: u'metarhizium_majus_arsef_297',
     229: u'metarhizium_rileyi_rcef_4871',
     230: u'metschnikowia_bicuspidata_var_bicuspidata_nrrl_yb_4993',
     231: u'meyerozyma_guilliermondii_atcc_6260',
     232: u'microbotryum_violaceum',
     233: u'microdochium_bolleyi',
     234: u'millerozyma_farinosa_cbs_7064',
     235: u'mixia_osmundae_iam_14324',
     236: u'moesziomyces_antarcticus',
     237: u'moniliophthora_perniciosa_fa553',
     238: u'moniliophthora_roreri_mca_2997',
     239: u'mortierella_elongata_ag_77',
     240: u'mucor_ambiguus',
     241: u'mucor_circinelloides_f_circinelloides_1006phl',
     242: u'mycosphaerella_eumusae',
     243: u'nadsonia_fulvescens_var_elongata_dsm_6958',
     244: u'nannizzia_gypsea_cbs_118893',
     245: u'naumovozyma_castellii_cbs_4309',
     246: u'naumovozyma_dairenensis_cbs_421',
     247: u'neofusicoccum_parvum_ucrnp2',
     248: u'neolecta_irregularis_dah_3',
     249: u'neolentinus_lepideus_hhb14362_ss_1',
     250: u'neonectria_ditissima',
     251: u'neosartorya_fischeri',
     252: u'neurospora_crassa',
     253: u'neurospora_tetrasperma_fgsc_2509',
     254: u'ogataea_parapolymorpha_dl_1',
     255: u'ogataea_polymorpha',
     256: u'oidiodendron_maius_zn',
     257: u'ophiocordyceps_sinensis_co18',
     258: u'ophiocordyceps_unilateralis',
     259: u'ophiostoma_piceae_uamh_11346',
     260: u'paracoccidioides_brasiliensis_pb18',
     261: u'paracoccidioides_sp_lutzii_pb01',
     262: u'paraphaeosphaeria_sporulosa',
     263: u'parasitella_parasitica',
     264: u'paxillus_rubicundulus_ve08_2h10',
     265: u'penicilliopsis_zonata_cbs_506_65',
     266: u'penicillium_antarcticum',
     267: u'penicillium_arizonense',
     268: u'penicillium_brasilianum',
     269: u'penicillium_camemberti_fm_013',
     270: u'penicillium_chrysogenum',
     271: u'penicillium_coprophilum',
     272: u'penicillium_decumbens',
     273: u'penicillium_digitatum_phi26',
     274: u'penicillium_expansum',
     275: u'penicillium_flavigenum',
     276: u'penicillium_freii',
     277: u'penicillium_griseofulvum',
     278: u'penicillium_italicum',
     279: u'penicillium_nalgiovense',
     280: u'penicillium_nordicum',
     281: u'penicillium_oxalicum_114_2',
     282: u'penicillium_polonicum',
     283: u'penicillium_roqueforti_fm164',
     284: u'penicillium_rubens_wisconsin_54_1255',
     285: u'penicillium_solitum',
     286: u'penicillium_steckii',
     287: u'penicillium_subrubescens',
     288: u'penicillium_vulpinum',
     289: u'peniophora_sp_cont',
     290: u'pestalotiopsis_fici_w106_1',
     291: u'phaeoacremonium_minimum_ucrpa7',
     292: u'phaeosphaeria_nodorum',
     293: u'phanerochaete_carnosa_hhb_10118_sp',
     294: u'phialocephala_scopiformis',
     295: u'phialocephala_subalpina',
     296: u'phialophora_americana',
     297: u'phialophora_attae',
     298: u'phlebia_centrifuga',
     299: u'phlebiopsis_gigantea_11061_1_cr5_6',
     300: u'phycomyces_blakesleeanus_nrrl_1555_',
     301: u'pichia_membranifaciens_nrrl_y_2026',
     302: u'piloderma_croceum_f_1598',
     303: u'pisolithus_microcarpus_441',
     304: u'pisolithus_tinctorius_marx_270',
     305: u'pleurotus_ostreatus_pc15',
     306: u'pneumocystis_carinii_b80',
     307: u'pneumocystis_jirovecii_ru7',
     308: u'pneumocystis_murina_b123',
     309: u'pochonia_chlamydosporia_170',
     310: u'podospora_anserina_s_mat_',
     311: u'pseudocercospora_fijiensis_cirad86',
     312: u'pseudocercospora_musae',
     313: u'pseudogymnoascus_destructans_20631_21',
     314: u'pseudogymnoascus_sp_03vt05',
     315: u'pseudogymnoascus_sp_05ny08',
     316: u'pseudogymnoascus_sp_23342_1_i1',
     317: u'pseudogymnoascus_sp_24mn13',
     318: u'pseudogymnoascus_sp_vkm_f_103',
     319: u'pseudogymnoascus_sp_vkm_f_3557',
     320: u'pseudogymnoascus_sp_vkm_f_3775',
     321: u'pseudogymnoascus_sp_vkm_f_3808',
     322: u'pseudogymnoascus_sp_vkm_f_4246',
     323: u'pseudogymnoascus_sp_vkm_f_4281_fw_2241_',
     324: u'pseudogymnoascus_sp_vkm_f_4513_fw_928_',
     325: u'pseudogymnoascus_sp_vkm_f_4514_fw_929_',
     326: u'pseudogymnoascus_sp_vkm_f_4515_fw_2607_',
     327: u'pseudogymnoascus_sp_vkm_f_4516_fw_969_',
     328: u'pseudogymnoascus_sp_vkm_f_4517_fw_2822_',
     329: u'pseudogymnoascus_sp_vkm_f_4518_fw_2643_',
     330: u'pseudogymnoascus_sp_vkm_f_4519_fw_2642_',
     331: u'pseudogymnoascus_sp_vkm_f_4520_fw_2644_',
     332: u'pseudogymnoascus_sp_wsf_3629',
     333: u'pseudogymnoascus_verrucosus',
     334: u'pseudozyma_hubeiensis_sy62',
     335: u'puccinia_graminis',
     336: u'puccinia_sorghi',
     337: u'puccinia_striiformis',
     338: u'puccinia_triticina',
     339: u'purpureocillium_lilacinum',
     340: u'pyrenochaeta_sp_ds3say3a',
     341: u'pyrenophora_teres',
     342: u'pyrenophora_triticirepentis',
     343: u'pyronema_omphalodes_cbs_100304',
     344: u'rachicladosporium_antarcticum',
     345: u'rasamsonia_emersonii_cbs_393_64',
     346: u'rhinocladiella_mackenziei_cbs_650_93',
     347: u'rhizoctonia_solani_ag_1_ib',
     348: u'rhizophagus_irregularis_daom_181602_gca_000439145',
     349: u'rhizopogon_vesiculosus',
     350: u'rhizopogon_vinicolor_am_or11_026',
     351: u'rhizopus_delemar_ra_99_880',
     352: u'rhizopus_microsporus',
     353: u'rhodotorula_graminis_wp1',
     354: u'rhodotorula_sp_jg_1b',
     355: u'rhodotorula_toruloides_np11',
     356: u'rhynchosporium_agropyri',
     357: u'rhynchosporium_commune',
     358: u'rhynchosporium_secalis',
     359: u'rozella_allomycis_csf55',
     360: u'saccharomyces_arboricola_h_6',
     361: u'saccharomyces_cerevisiae_x_saccharomyces_kudriavzevii_vin7',
     362: u'saccharomyces_eubayanus',
     363: u'saccharomyces_kudriavzevii_ifo_1802',
     364: u'saccharomyces_sp_boulardii_',
     365: u'saccharomycetaceae_sp_ashbya_aceri_',
     366: u'scedosporium_apiospermum',
     367: u'scheffersomyces_stipitis_cbs_6054',
     368: u'schizophyllum_commune_h4_8',
     369: u'schizopora_paradoxa',
     370: u'schizosaccharomyces_cryophilus',
     371: u'schizosaccharomyces_japonicus',
     372: u'schizosaccharomyces_octosporus',
     373: u'schizosaccharomyces_pombe',
     374: u'scleroderma_citrinum_foug_a',
     375: u'sclerotinia_borealis_f_4128',
     376: u'serendipita_indica_dsm_11827',
     377: u'serendipita_vermifera_maff_305830',
     378: u'serpula_lacrymans_var_lacrymans_s7_3',
     379: u'setosphaeria_turcica_et28a',
     380: u'sistotremastrum_niveocremeum_hhb9708',
     381: u'sistotremastrum_suecicum_hhb10207_ss_3',
     382: u'sordaria_macrospora',
     383: u'spathaspora_passalidarum_nrrl_y_27907',
     384: u'sphaerobolus_stellatus_ss14',
     385: u'sphaerulina_musiva_so2202',
     386: u'spizellomyces_punctatus_daom_br117',
     387: u'sporidiobolus_salmonicolor',
     388: u'sporisorium_reilianum',
     389: u'sporisorium_scitamineum',
     390: u'sporothrix_brasiliensis_5110',
     391: u'sporothrix_insectorum_rcef_264',
     392: u'sporothrix_schenckii_atcc_58251',
     393: u'stachybotrys_chartarum_ibt_40288',
     394: u'stachybotrys_chlorohalonata_ibt_40285',
     395: u'stagonospora_sp_src1lsm3a',
     396: u'stemphylium_lycopersici',
     397: u'suillus_luteus_uh_slu_lm8_n1',
     398: u'talaromyces_islandicus',
     399: u'talaromyces_marneffei_atcc_18224',
     400: u'talaromyces_stipitatus_atcc_10500',
     401: u'termitomyces_sp_j132',
     402: u'tetrapisispora_blattae_cbs_6284',
     403: u'tetrapisispora_phaffii_cbs_4417',
     404: u'thermothelomyces_thermophila_atcc_42464',
     405: u'thielavia_terrestris_nrrl_8126',
     406: u'thielaviopsis_punctulata',
     407: u'tilletia_caries',
     408: u'tilletia_controversa',
     409: u'tilletia_indica',
     410: u'tilletia_walkeri',
     411: u'tilletiaria_anomala_ubc_951',
     412: u'tolypocladium_ophioglossoides_cbs_100239',
     413: u'torrubiella_hemipterigena',
     414: u'tortispora_caseinolytica_nrrl_y_17796',
     415: u'torulaspora_delbrueckii',
     416: u'trametes_pubescens',
     417: u'tremella_mesenterica_dsm_1558',
     418: u'trichoderma_atroviride_imi_206040',
     419: u'trichoderma_gamsii',
     420: u'trichoderma_guizhouense',
     421: u'trichoderma_harzianum',
     422: u'trichoderma_reesei',
     423: u'trichoderma_virens',
     424: u'trichophyton_benhamiae_cbs_112371',
     425: u'trichophyton_equinum_cbs_127_97',
     426: u'trichophyton_interdigitale_mr816',
     427: u'trichophyton_rubrum_cbs_118892',
     428: u'trichophyton_soudanense_cbs_452_61',
     429: u'trichophyton_tonsurans_cbs_112818',
     430: u'trichophyton_verrucosum_hki_0517',
     431: u'trichophyton_violaceum',
     432: u'trichosporon_asahii_var_asahii_cbs_8904',
     433: u'tsuchiyaea_wingfieldii_cbs_7118',
     434: u'tuber_melanosporum',
     435: u'tulasnella_calospora_mut_4182',
     436: u'umbilicaria_pustulata',
     437: u'uncinocarpus_reesii_1704',
     438: u'ustilaginoidea_virens',
     439: u'ustilago_bromivora',
     440: u'ustilago_hordei',
     441: u'ustilago_maydis',
     442: u'valsa_mali_var_pyri',
     443: u'vanderwaltozyma_polyspora_dsm_70294',
     444: u'verruconis_gallopava',
     445: u'verticillium_alfalfae_vams_102',
     446: u'verticillium_dahliae',
     447: u'verticillium_longisporum',
     448: u'wallemia_ichthyophaga_exf_994',
     449: u'wallemia_mellicola_cbs_633_66',
     450: u'wickerhamomyces_anomalus_nrrl_y_366_8',
     451: u'xanthophyllomyces_dendrorhous',
     452: u'xylona_heveae_tc161',
     453: u'yarrowia_lipolytica',
     454: u'zygosaccharomyces_bailii_isa1307',
     455: u'zygosaccharomyces_parabailii',
     456: u'zymoseptoria_brevis',
     457: u'zymoseptoria_tritici',
     458: u'_candida_auris',
     459: u'aspergillus_nidulans',
     460: u'candida_tropicalis_mya_3404',
     461: u'clavispora_lusitaniae_atcc_42720',
     462: u'colletotrichum_chlorophyti',
     463: u'colletotrichum_higginsianum',
     464: u'dichomitus_squalens_lyad_421_ss1',
     465: u'fonsecaea_monophora',
     466: u'glarea_lozoyensis_74030',
     467: u'hanseniaspora_uvarum_dsm_2768',
     468: u'hydnomerulius_pinastri_md_312',
     469: u'lentinula_edodes',
     470: u'lichtheimia_corymbifera_jmrc_fsu_9682',
     471: u'paxillus_involutus_atcc_200175',
     472: u'phaeomoniella_chlamydospora',
     473: u'pichia_kudriavzevii',
     474: u'saitoella_complicata_nrrl_y_17804',
     475: u'stereum_hirsutum_fp_91666_ss1',
     476: u'wickerhamomyces_ciferrii',
     477: u'zygosaccharomyces_rouxii',
     478: u'anncaliia_algerae_pra339',
     479: u'edhazardia_aedis_usnm_41457',
     480: u'encephalitozoon_cuniculi_gb_m1',
     481: u'encephalitozoon_hellem_atcc_50504',
     482: u'encephalitozoon_intestinalis_atcc_50506',
     483: u'encephalitozoon_romaleae_sj_2008',
     484: u'grifola_frondosa',
     485: u'isaria_fumosorosea_arsef_2679',
     486: u'komagataella_phaffii_gs115',
     487: u'mitosporidium_daphniae',
     488: u'moesziomyces_aphidis_dsm_70725',
     489: u'nematocida_displodere',
     490: u'nematocida_parisii_ertm3',
     491: u'nematocida_sp_1_ertm2',
     492: u'nematocida_sp_ertm5',
     493: u'nosema_apis_brl_01',
     494: u'nosema_bombycis_cq1',
     495: u'nosema_ceranae_brl01',
     496: u'ordospora_colligata_oc4',
     497: u'pseudoloma_neurophilia',
     498: u'sclerotinia_sclerotiorum',
     499: u'spraguea_lophii_42_110',
     500: u'sugiyamaella_lignohabitans',
     501: u'trachipleistophora_hominis',
     502: u'trametes_cinnabarina',
     503: u'vavraia_culicis_subsp_floridensis',
     504: u'vittaforma_corneae_atcc_50505',
     505: u'gymnopus_luxurians_fd_317_m1',
     506: u'plicaturopsis_crispa_fd_325_ss_3',
     507: u'postia_placenta_mad_698_r',
     508: u'mortierella_verticillata_nrrl_6337',
     509: u'punctularia_strigosozonata_hhb_11173_ss5',
     510: u'aspergillus_ustus',
     511: u'coniophora_puteana_rwd_64_598_ss2',
     512: u'batrachochytrium_dendrobatidis_jel423',
     513: u'trametes_versicolor_fp_101664_ss1',
     514: u'enterocytozoon_bieneusi_h348',
     515: u'enterocytozoon_hepatopenaei',
     516: u'enterospora_canceri',
     517: u'hepatospora_eriocheir',
     518: u'talaromyces_cellulolyticus',
     519: u'fomitiporia_mediterranea_mf3_22',
     520: u'leucoagaricus_sp_symc_cos',
     521: u'pachysolen_tannophilus_nrrl_y_2460'}


