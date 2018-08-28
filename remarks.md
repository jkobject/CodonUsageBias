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

 correlatio between endres and euclidean (correlation=0.96924282826462216, pvalue=0.0) better with 98%=cub > cuf > entropy = 94%

