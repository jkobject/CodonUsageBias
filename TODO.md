# CODE

## to test
 - download datasets


## todo 

--> set documentation for each functions and objects __done__
--> ask Yun for the phylogenetic distance if she has it and the function to compute the phylogenetic distance matrix from a pseudo half life of genes with the gamma function. __done__
>I have thought to create a pseudo evolution distance matrix based on phylogenetic tree classification by way of calculating ‘half-life’ (using gamma distribution) of real amino acid or gene which I read on book 'molecular evolution'


- plot for each species, how much its genes are outliers, how much are belonging to a secondary cluster and how much are belonging to the principal cluster.
--> create a long stacked bar plot with these values

- compare amongst the working homoset homologies, the clusters together, by what species they contains by creating a new vector of species presence in each cluster and plotting the similarity matrix of each of those vectors. 
     --> look at similarity matrix within genes of 1 species and get some statistics over this matrix. full rank, etc ? 
	---> look at similarity matrix between species of their genes pertaining to same cluster (cluster matrix species A === cluster Matrix of species B for all B and A) then sum of the bools and divide by total number of homologies
	or by total number of common genes (or just cosine similarity of their has homo matrix.)
    --> create a compare function in homoset of homologies clusters similarity matrix and ordering. 

- find if we are close to ancestry tree, here we need to represent a comparison of the closeness in a phylogenetic tree to a cluster of species
	--> can create an interactive plot showing a gradient of color intensity on the points computed on the phylogenetic distance to the point currently hovered. 
or	--> given a grouping of phylogenetic tree, what cluster is the most similar to it
	--> compute the F1 score (amount of species that should be in the cluster and are not, that shouldnot be but are in the cluster ...) regularized with respect to the size of each cluster, for each homologies.
	--> set the eps values by finding clusters that best relate to ancestry tree
	 	--> find for 20 random homologies, what is the eps value that increase the F1 score (regularized by some value) the most
then 	--> how does it behaves on the full set 

-> use gaussian clustering and look at the variance of the kernels (requested by dominique to maybe have some ideas of variance
as it is not well displayed by eps) add this as another information when plotting and make gaussian clustering work well, then work using the values found by eps. 

other idea:
	--> given the distance between each two data points (weighted) averaged by each homologies. find the difference with the pseudo phylogenetic distance, how much does it explains this distance. what are the weights that better helps explaining it. what are those weights exactly for each homology and what can they tell us (if they can) on the importance and role of each homology. is it similar for each two species? (try for a group of distant ones and of closer ones)

- ensembl might think that it means they are not related and are not part of the same homology as they are not using similarity of the transcribed protein function but only similarity of the DNA strings (this is also normal since two individuals could have gained the same functions from two different events) 
--> so, one should look at the mean variance in CUB value and mean range for each homology (add this when plotting)

- use a mean sequence values and compare it to all sequences values for one species and one homology
--> variance in data (for full species) can be explained by a value in tRNA copy number.

--> Write two new functions one that computes codon frequency from ensembl and another that computes
entropy location from entropy in homoset

-->Use entropy, entropy location, codon frequency, (normalized and unormalized ones.), use random values as well (random values can be created from the same thing as the partition functions for entropy location. 
--> average over sequences/species/homologies

---> download or find more meta datasets and see what can be done with them http://ensemblgenomes.org/info/data/protein_features and http://ensemblgenomes.org/info/data/repeat_features ( protein features and repeat features)
	 
- look at similarity distance in sequences and compare it to the CUB values I have whether similarity 
(is entropy)

---> 


---- STOP end of JULY and do this ------------
 - set everything in order and working + documentation(lateX, html) + website (3d)
 - create a nice explaining notebook (1d)
 - extract figures, create tabs, diagrams (2d)
 - have each other doument ready ( logbook, readme, explanations, ...) (1d)
 - write a 60 page thesis (10d)
 - format it in lateX (1d)


# INFO 

## ideas 

- see if there is a similarity in the cluster position of each of their genes, compare it and find relationships
- compare the similarities between the species for each types of related genes between the species
- find a relationship between 
- plot entropy location values similarity between genes of one particular individual
- look for philogenetic similarity using ncbi taxonomy information/website or the information for Tobias

* To find Metadata 
	* from ensembl :  http://ensemblgenomes.org/info/data/protein_features and http://ensemblgenomes.org/info/data/repeat_features ( protein features and repeat features)
	* look for other websites 


- look for the trna genomic content using the metadata and then use it as a discriminatory factor (following paper by Mario dos Reis et. al.)

- evaluate the distance between clusters. 
- find the cluster number with algorithm X-means (gaussian distribution and extract the inter-cluster distance)
- prepare a pipeline for intra species homologies

- create a classification NN P(species | genename )

 -  I should look at species that are not represented in the dataset.

##  other things to do 

- the species are ordered according to philogenetic/taxonomic similarities
send an email to tobias for the name of the comparer species.
ask for the philogenetic list amongs the 461 species
the species present here are all the same ones. 

- find how to replace Nans as they convey info and couldnot be replaced by 0 or 1 (maybe 0.5 ?)
`http://scikit-learn.org/stable/modules/preprocessing.html`

- look at the equations on the thesis. 

- WE can look for classification methods such as :

CLUSTER

 * K-means 
 * Gaussian Mixture (clustering) http://scikit-learn.org/stable/modules/mixture.html#gaussian-mixture

 * Affinity propagation (cl)
>don't forget performance testing of clustering algo http://scikit-learn.org/stable/modules/clustering.html#clustering-performance-evaluation

 * DB-scan 

CLUSTERER-CLASSIFIER

 * self organizing maps > 1000 neurons
NOVELTY/OUTLIER DETECTOR http://scikit-learn.org/stable/modules/outlier_detection.html

FEATURE EXTRACTOR

 * Var autoencoder
 * graphical models http://scikit-learn.org/stable/modules/neural_networks_unsupervised.html#graphical-model-and-parametrization
 
DIM REDUCTION

 * PCA - ICA
 * tsn-e

OTHER

 * NMF (given this gene values,  to which species it belongs| given this gene cluster, to which species it belongs | given this species what could be the value if it had this homology) http://scikit-learn.org/stable/modules/decomposition.html#non-negative-matrix-factorization-nmf-or-nnmf
 


## Pipeline of the matlab code
so you have one species and you look for homologies to a subset of 5818 genes of this species. 
th
main = homologyvalue()
then you gethomoinformation() to get the information about the homology to 461 species
then you look at preinstalled genotype of 461 species and get this particular gene 


- read the paper
- cancerous cells should have a particular kind of codon usage bias (more than regular cells)

-  Entropy location : you compute the entropy value for each possible configuration, it makes a distribution to wich you compare the measured entropy value. with a integration.

- there is a whole genome part as well 
KL method : for the comparison of the whole genome per species (only exons !) we sum up all the distributions and normalize them ( by the length of the gene) and we compute the KL divergence to the distribution made by the measured distribution made by the plotting of all the exons
Expected entropy: (see notebook)


- focus first on the classification 
- be really carefull about what you are goin to present on your thesis. 


## communication management

__Dominique Chu__

>Really do the same plot with the 5000 species to look for structure
Do it with entropyValue and Entropy\*Length as well


__Yun__

- what to do with similar homologies in first 500 and 600-1000?


__Tobias__

- what to do with double homologies ( homology in the same species)


__Alex Freitas__

- density based clustering (and ensemble)
- look for way to assess clustering quality (inter intra cluster)
- regularize by the amount of clusters (pareto)

> the gene clustering (of i-one gene ) can be seen as a mini phylogenetic tree
> once I have the meta data how much can each features or groups of features can explain the clustering that I found. 


### question

- how are you finding your homologies I feel as though they are not the same for every species. How can I compare them in this way ? 

> you should not look at the homologies like that. the genes can actually be very different and the way they find homologies is by comparing how it behaves or if they know it is an homologies or the RNA or DNA list...

- should I normalize the data using the lenghth to compare them ?

>no

- how do I read your csv ?

> it is ordered as a tree from similarity alignment matching ( look for similarity alignment phylogenetic matching algorithm)
BUT for Yun's it is according to the names (family) taxonomical data + similarity alignement matching (just getting it from itol.embl.de)

- where can I get the philo tree that you have ? 

>to recreate it from the CSV file you need to send it to the itol.embl.de website

- unicel or not ?

>not

 - pathogene or not (plant/animal) ?

> some yes 

- what is the difference between Nc value and entropy location ? 



### To say


problem on what I am going to write and do as a master thesis
	python package with full pipeline (versatile)
	how much can I say about yun's work ( she was stressed)

homologues not related to their functions but their common ancestors
the way 

- How do i find the full exome of eahc species I have ?

### require from Tobias

 - phylogenetic tree
 - temperature
 - replication speed
 - frequency of use of the gene
 - nocive species
 - gene size
 - type of species
 - more ...

## what we said

- there is a problem with the homologies try to look at yun's data YFL039C else maybe a problem of the database or a problem of the sequencing. look at it find ideas to debug this (infer this). it is due to the moment where I do not take into account the duplicates and the files where less species are written due to 
lack of some amino acids, so you need to rewrite your function for processing Yun's data as well as look at what is outputed is 
consistent with what should naturally happening: an has homo matrix that is almost full. 
_did it and it worked._ 

I think how to give a better explanation for why cluster genes into 2 or 3 clusters, we can use artificial sequences (replaced according to codon usage table) to test, if the principle cluster genes have the same performance as artificial ones, which means these genes have the same trend as the whole genome wide.  

-> where to make them ?


## work pipeline 

Phase (1) – Business understanding
Understanding the business objectives and requirements, and converting them into a data mining problem

Phase (2) – Data understanding
Understand the data to be mined, consider data quality issues

Phase (3) – Data preparation 
Data cleaning, attribute and record selection, data transformation, etc.; dependent on the data mining algorithm(s) to be used in (4)

Phase (4) – Modelling / data mining phase
Choose one or more data mining algorithms to be applied to the data, adjust its(their) parameters, build a model of the data
Includes an initial objective, data-driven validation

Phase (5) – Model evaluation
Interpret and validate the models from a business perspective 

Phase (6) – Deployment of the model
Deploy the data mining results in the form of a plan, monitoring and maintaining the plan deployment

