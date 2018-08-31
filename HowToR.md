# How to :
## create a phylogenetic distance matrix
	- get the _ncbi taxonomic id_ (from ensembl for example)
	- create a phylogenetic tree from a tree creating algorithm
	- on R write:
	```R
	 treeText <- readLines(tree.phy)
	        treeText <- paste0(treeText, collapse="")
	        library(treeio)
	        tree <- read.tree(text = treeText) ## load tree 
	    distMat <- cophenetic(tree) ## generate dist matrix
	```
[source](https://www.biostars.org/p/312148/)
[could have used](http://etetoolkit.org/documentation/ete-compare/)