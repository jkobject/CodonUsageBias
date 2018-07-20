""""
Created by Jeremie KALFON
Date : 21 FEV 2018
University of Kent, ECE paris
jkobject.com

"""


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#import holoviews as hv
#from holoviews.operation.datashader import datashade, dynspread, shade
from joblib import Parallel, delayed

import utils
import homology as h

from sklearn.cluster import KMeans
from sklearn.cluster import MiniBatchKMeans
from sklearn.metrics.pairwise import cosine_similarity
#from kmodes.kmodes import KModes
from sklearn import manifold as man
from sklearn import metrics
from sklearn.decomposition import PCA
import tsne
from scipy import random
from scipy import sparse

import collections

import pdb


class HomoSet(collections.MutableMapping):
    """HomoSet is the object containing evrey homology as a dictionnary according to thie rhomology code
        from Homoset you can do much computation that requires the set of homologies

        Object where we store an homology group basically where we do our entire
        Computation from.

                params:
                ------
                hashomo_matrix : a numpy boolean array (species, homologies) that store the matrix of gene presence in species
                homo_matrix : a numpy array float (homologies*species(inthehomologies),aminoacids)
                                    but containing the codon entropy vectors instead
                homodict : dictionnary of homology object containing codon usage per species
                    for the gene coresponding to this homology
                homo_namelist : list of all the homology names
                species_namelist: set of all the species names
                clusters: a list of clusters for the homoset clustering can be of size of nb of species
                    or of nb of homologies
                homogroupnb: number of group for the homology clustering
                red_homomatrix: numpy array (homologies*species(inthehomologies),x*y)the reduced 2D version
                    of the homology matrix for the full homology matrix
                wasclusterized: if the homologies have been clustered or not. usefull for processing requiring clusters
                homocluster: np array homologies*species of int of cluster number

    """

    hashomo_matrix = None
    homo_matrix = None
    homo_matrixnames = None
    fulleng = None
    homodict = {}
    homo_namelist = []
    species_namelist = set([])
    clusters = []
    homogroupnb = 2
    red_homomatrix = None
    wasclusterized = False
    homo_clusters = None
    datatype = ''
    averagehomo_matrix = None

    def __init__(self, **kwargs):
        """
        will initialize the object with the different values you might have from another project
        use the data dictionnary to add any type of data
        """
        data = kwargs.get("data", None)
        if data is not None:
            self.homo_matrix = np.asarray(data["homo_matrix"]) if data.get(
                "homo_matrix", None) is not None else None
            self.homo_namelist = data.get("homo_namelist", [])
            self.species_namelist = set(data.get("species_namelist", []))
            self.homogroupnb = data.get("homogroupnb", 2)
            self.clusters = data.get("clusters", [])
            self.datatype = data.get("datatype", '')
            self.phylo_distance = pd.read_json(data["phylo_distance"], orient='split') if data.get(
                "phylo_distance", None) is not None else None
            self.red_homomatrix = np.asarray(data["red_homomatrix"]) if data.get(
                "red_homomatrix", None) is not None else None
            utils.speciestable = data.get("speciestable", {})
            utils.phylo_distances = data.get("phylo_distances", None)
            utils.meandist = utils.phylo_distances.sum().sum() / (len(
                utils.phylo_distances)**2 - len(utils.phylo_distances)) if utils.phylo_distances is not None else None
            self.hashomo_matrix = np.asarray(data["hashomo_matrix"]) if data.get(
                "hashomo_matrix", None) is not None else None
            self.fulleng = np.asarray(data["fulleng"]) if data.get(
                "fulleng", None) is not None else None
            self.wasclusterized = data.get('wasclusterized', False)
            self.homo_clusters = np.asarray(data["homo_clusters"]) if data.get(
                "homo_clusters", None) is not None else None
            self.averagehomo_matrix = np.asarray(data["averagehomo_matrix"]) if data.get(
                "averagehomo_matrix", None) is not None else None
            tp = {}
            for key, val in data["homodict"].iteritems():
                tp.update({key: h.homology(data=val)})
            self.homodict = tp
        else:
            self.homo_matrix = kwargs.get("homo_matrix", None)
            self.homo_namelist = kwargs.get("homo_namelist", [])
            self.species_namelist = kwargs.get("species_namelist", set([]))
            self.homogroupnb = kwargs.get("homogroupnb", 2)
            self.clusters = kwargs.get("clusters", [])
            self.datatype = kwargs.get("datatype", '')
            self.phylo_distance = kwargs.get("phylo_distance", None)
            self.red_homomatrix = kwargs.get("red_homomatrix", None)
            utils.speciestable = kwargs.get("speciestable", {})
            utils.phylo_distances = kwargs.get("phylo_distances", None)
            utils.meandist = utils.phylo_distances.sum().sum() / (len(
                utils.phylo_distances)**2 - len(utils.phylo_distances)) if utils.phylo_distances is not None else None
            self.hashomo_matrix = kwargs.get("hashomo_matrix", None)
            self.fulleng = kwargs.get("fulleng", None)
            self.wasclusterized = kwargs.get('wasclusterized', False)
            self.homo_clusters = kwargs.get("homo_clusters", None)
            self.averagehomo_matrix = kwargs.get("averagehomo_matrix", None)
        # TODO: -- add loads and save everywhere
    # magic methods https://rszalski.github.io/magicmethods/

    def __getitem__(self, key):
        if type(key) is str:
            return self.homodict[key]
        if type(key) is int:
            return self.homodict[self.homo_namelist[key]]

    def __setitem__(self, key, value):
        self.homodict[key] = value

    def __delitem__(self, key):
        del self.homodict[key]

    def __iter__(self):
        return iter(self.homodict)

    def iteritems(self):
        return self.homodict.iteritems()

    def __len__(self):
        return len(self.homodict)

    def update(self, val):
        self.homodict.update(val)


##############################################################

    def plot_all(self, With='tsne', perplexity=60, interactive=False, iteration=300):
        """
        will plot all the homologies in the full_homo_matrix (and compute it)
        (sometimes around 800 000 datapoints) to look at any kind of relationships
        as there is too much datapoints, the plots are density ones.

        Params:
        -------
        With: flag the dim reduction algorithm to use (tsne: need >16gigs of RAM,now use another version
        of tsne for large datasets )(PCA: works well)(lsta/hessian:untested)
            perplexity,iteration: ints of basic tsne hyperparams
        interactive: bool if true should use bokeh else matplot lib

        Returns:
        -------
        the desired plot
        """
        if self.homo_matrix is None:
            self.loadfullhomo()

        if self.red_homomatrix is None:
            red = None
            if With == 'tsne':
                # TODO: to debug
                # red = man.TSNE(n_components=2, perplexity=perplexity, n_iter=iteration, verbose=3)
                # .fit_transform(self.homo_matrix)
                red = tsne.tsne(self.homo_matrix, no_dims=2, initial_dims=self.homo_matrix.shape[1], perplexity=30.0)
            if With == 'PCA':
                red = PCA(n_components=2).fit_transform(self.homo_matrix)
            if With == 'ltsa' or With == 'hessian':
                red = man.LocallyLinearEmbedding(self.homo_matrix.shape[0] / 100, 2,
                                                 eigen_solver='auto', method=With).fit_transform(self.homo_matrix)
            self.red_homomatrix = red
        if not interactive:
            fig = plt.figure(figsize=(40, 40))
            ax = fig.add_subplot(111)
            ax.scatter(self.red_homomatrix[:, 0], self.red_homomatrix[:, 1])
            plt.show()
        else:
            """
            #hv.notebook_extension('bokeh')
            dynspread.max_px = 200
            dynspread.threshold = 0.5
            shade.cmap = "#30a2da"  # to match HV Bokeh default
            points = hv.Points(self.red_homomatrix, label="all " + str(self.red_homomatrix.shape[0]) + " homologies")

            def heatmap(coords, bins=10, offset=20.0, transform=lambda d, m: d, label=None):
                ""
                Given a set of coordinates, bins them into a 2d histogram grid
                of the specified size, and optionally transforms the counts
                and/or compresses them into a visible range starting at a
                specified offset between 0 and 1.0.
                ""
                hist, xs, ys = np.histogram2d(coords.T[0], coords.T[1], bins=bins)
                counts = hist[:, ::-1].T
                transformed = transform(counts, counts != 0)
                span = transformed.max() - transformed.min()
                compressed = np.where(counts != 0, offset + (1.0 - offset) * transformed / span, 0)
                args = dict(label=label) if label else {}
                return hv.Image(compressed, bounds=(xs[-1], ys[-1], xs[1], ys[1]), **args)
            print 'please write %output size=200" before calling this function'
            return heatmap(self.red_homomatrix, 100)(style=dict(cmap="fire")) + datashade(points)
            """

    def loadfullhomo(self):
        """
        function to concatenate all the homologies in one big array(practicle for certain computations)
        """
        self.homo_matrix = self.homodict[self.homo_namelist[0]].full
        self.homo_matrixnames = self.homodict[self.homo_namelist[0]].names
        self.fulleng = self.homodict[self.homo_namelist[0]].lenmat
        for x, homo in enumerate(self.homo_namelist[1:]):
            try:
                self.homo_matrix = np.vstack((self.homo_matrix, self.homodict[homo].full))
                self.homo_matrixnames.extend(self.homodict[homo].names)
                self.fulleng = np.vstack((self.fulleng, self.homodict[homo].lenmat))
            except ValueError:
                print x, self.homodict[homo].full.shape
                pdb.set_trace()
            print '{0}%\r'.format((x * 100) / len(self.homo_namelist)),
        self.homo_matrixnames = np.asarray(self.homo_matrixnames)

    def loadhashomo(self):
        """
        function to compute the matrix of bool saying wether species X has a gene or more in homology Y
        """
        hashomo = np.zeros((len(self.homo_namelist),
                            len(self.species_namelist)), dtype=np.bool)
        for i, homo in enumerate(self.homo_namelist):
            ind = [na for y, na in enumerate(self.homodict[homo].names) if not self.homodict[homo].doub[y]]
            # get the indeces from each species name
            hashomo[i][ind] = True
        self.hashomo_matrix = hashomo

    def size(self):
        if self.has_homomatrix is None:
            self.loadhashomo()
        return np.count_nonzero(self.hashomo_matrix)

    def add_random_homology(self):
        """
        a function to populate a homology with random values
        """
        # TODO: to test
        full = np.random.rand(len(self.species_namelist), 18)
        nans = np.zeros(len(self.species_namelist))
        names = []
        for val in len(self.species_namelist):
            names.append("species_" + str(int(random.rand() * 3000)))
        similarity_scores = np.random.rand(len(self.species_namelist))
        proteinids = []
        for val in len(self.species_namelist):
            proteinids.append("proteinids_" + str(int(random.rand() * 6)))
        Kaks_Scores = np.random.rand(len(self.species_namelist))
        lenmat = np.random.rand(len(self.species_namelist), 18) * 100
        lenmat = lenmat.astype(int)
        homo = h.homology(full=full, nans=nans, names=names, similarity_scores=similarity_scores,
                          proteinids=proteinids, Kaks_Scores=Kaks_Scores, lenmat=lenmat)
        self.homodict.update({'random_' + str(int(random.rand() * 1000)): homo})

    def preprocessing(self, withtaxons=False):
        """
        will compute the full list of names, find doublons, and set the names to ints instead
        of strings. called after loading from ensembl and associate namelist in each homologies to a number
        in utils.speciestable

        Returns:
        ------
        taxons,species : list, if the names contains an additional list of taxon ids.
                        and the corresponding species in the same order
                        only if with taxons
        """
        species_namelist = set([])
        taxons = []
        species = []  # we need to have another species list for the ordering
        # to be kept as it should (a set is ordered to be accessed in log time)
        utils.speciestable = {}
        if self.homodict is not None:
            i = 0
            helper = {}
            for key, val in self.homodict.iteritems():
                if withtaxons:
                    doub = np.zeros(len(val.names[0]), dtype=np.bool)
                    names = []
                    temp = ''
                    for j, name in enumerate(val.names[0]):
                        if name == temp:
                            doub[j] = True
                        else:
                            temp = name
                            old_l = len(species_namelist)
                            species_namelist.add(name)
                            if len(species_namelist) != old_l:
                                species.append(val.names[0][j])
                                taxons.append(val.names[1][j])
                                utils.speciestable.update({i: name})
                                helper.update({name: i})
                                i += 1
                        names.append(helper[name])
                    val.names = names
                    val.doub = doub
                else:
                    np.zeros(len(val.names), dtype=np.bool)
                    names = []
                    temp = ''
                    for j, name in enumerate(val.names):
                        if name == temp:
                            doub[j] = True
                        else:
                            temp = name
                            old_l = len(species_namelist)
                            species_namelist.add(name)
                            if len(species_namelist) != old_l:
                                utils.speciestable.update({i: name})
                                helper.update({name: i})
                                i += 1
                        names.append(helper[name])
                    val.names = names
                    val.doub = doub
            self.species_namelist = species_namelist
            if withtaxons:
                return taxons, species

    def compute_entropyloc(self, using='normal'):
        """
        called if need entropy location and used ensembl data. you can always compute entropy location
        from entropy data.

        Will be much faster than doing it directly when calling
        ensembl's data as it computes the partition function
        only one for each lengths
        """
        # TODO: totest
        if self.datatype is 'entropy':
            if self.homo_matrix is None:
                self.homomatrix = self.homodict[self.homo_namelist[0]].full
                self.fulleng = self.homodict[self.homo_namelist[0]].lenmat
                for x, homo in enumerate(self.homo_namelist[1:]):
                    try:
                        self.homo_matrix = np.vstack((self.homo_matrix, self.homodict[homo].full))
                        self.fulleng = np.vstack((self.fulleng, self.homodict[homo].lenmat))
                    except ValueError:
                        print x, self.homodict[homo].full.shape
                        pdb.set_trace()
                    print '{0}%\r'.format((x * 50) / len(self.homo_namelist)),
            else:
                self.fulleng = self.homodict[self.homo_namelist[0]].lenmat
                for x, homo in enumerate(self.homo_namelist[1:]):
                    try:
                        self.fulleng = np.vstack((self.fulleng, self.homodict[homo].lenmat))
                    except ValueError:
                        print x, self.homodict[homo].full.shape
                    print '{0}%\r'.format((x * 50) / len(self.homo_namelist)),
            self.homomatrix = utils.getloc(self.homo_matrix, self.fulleng, using=using)
            pos = 0
            for x, homo in enumerate(self.homo_namelist):
                posi = pos + len(self.homodict[homo].full) + 1
                self.homodict[homo].full = self.homo_matrix[pos:posi]
                pos = posi
                print '{0}%\r'.format(50 + ((x * 50) / len(self.homo_namelist))),

    def remove(self, species):
        """
        remove this list of species from the homologies
        """
        # TODO: to finish coding
        for key in self.homodict.keys():
            self.homodict[key].remove(species)

    def plot_homoperspecies(self):
        """
        will plot the number of homology per spcies in this homology group
        """
        sumed = np.sum(self.hashomo_matrix, axis=1)
        plt.figure(figsize=(40, 10))
        plt.title('number of homologies per species')
        plt.bar(range(len(sumed)), sumed)
        print "you can always look at a particular range of species with 'homo_namelist' "

    def plot_speciesperhomo(self):
        """
        will plot the number of species per homologies in this homology group
        """
        sumed = np.sum(self.hashomo_matrix, axis=0)
        plt.figure(figsize=(40, 10))
        plt.title('number of species per homologies')
        plt.bar(range(len(sumed)), sumed)
        print "you can always look at a particular range of species with 'homo_namelist' "

    def order_from_matrix(self, clustering='kmeans', byspecie=False, order=True,
                          plot_ordering=True, homogroupnb=2, findnb=False,
                          reducedim=False):
        """
        Compute an homology group :
        from matrix computation using the homo_matrix
        (or from network computation in homologize_from_network)

        Can be computed many times and will updata homoset with the most recent homoset found
        if homoset exists, it will save it.

        Params:
        -------
        clustering: flags to 'kmeans', 'spectral', 'xmeans' to use different sk-learn algorithms

        plot: flags to true for the function to output ploting
        of the affinity matrix with and without the
        clusters

        homogroupnb: nb of groups you want to extract

        """
        # TODO: find more tests and visualisations to perform
        if reducedim:
            # TODO: try to reduce the dimensionality here before clustering
            print "tocode yet"
        if findnb is True:
            homogroupnb = 2
        clusts = []
        while True:
            if clustering == "kmeans":
                # The parallel version of K-Means is broken on OS X
                # when numpy uses the Accelerate Framework.
                # This is expected behavior: Accelerate can be called after a fork
                # but you need to execv the subprocess with the Python binary
                # (which multiprocessing does not do under posix)
                alg = KMeans(n_clusters=homogroupnb, n_jobs=-1)
            elif clustering == "fast":
                # TODO: totest
                alg = MiniBatchKMeans(n_clusters=homogroupnb)

            # elif clustering == "kmodes":
                # https://github.com/nicodv/kmodes/blob/master/kmodes/kmodes.py
                # alg = KModes(n_clusters=homogroupnb,
                #            init='Huang', n_init=2, verbose=1)

            else:
                print "you entered a wrong clustering algorithm"
                return False

            self.hashomo_matrix = self.hashomo_matrix.T if byspecie else self.hashomo_matrix
            # BYSPECIES is a special case where one would like to keep all homologies and
            # remove species instead. we would advise the opposite though
            alg.fit(self.hashomo_matrix)
            clust = alg.labels_
            metricA = metrics.silhouette_score(self.hashomo_matrix, clust, metric='euclidean')
            metricB = metrics.calinski_harabaz_score(self.hashomo_matrix, clust)
            if findnb is True:
                clusts[homogroupnb] = clust
                if homogroupnb < 10:
                    homogroupnb += 1
                    print "for " + homogroupnb + " cluster, the metric is " + metricA, metricB
                else:
                    print "you just select the clust by doing PyCUB.homo.clusters =\
                    output['groupnbyouwant']"
                    print "you an run the 'orderfromclust' also to have the plottings and orderings"
                    return clusts
            else:
                break
        print "the quality of the clustering is: "
        print metricA
        print metricB
        if order:
            self.hashomo_matrix = self.orderfromclust(homogroupnb, clust, byspecie=byspecie,
                                                      plot_ordering=plot_ordering)
        self.hashomo_matrix = self.hashomo_matrix.T if byspecie else self.hashomo_matrix
        self.homogroupnb = homogroupnb

    def orderfromclust(self, homogroupnb, clust, byspecie, plot_ordering):
        """
        creates an ordering of every elements (names, homologies according to the found clusters)
        from an ordered cluster of species or of homologies
        self.hashomomatrix should reflect this orientation as well
        """
        orderedhas = np.zeros(self.hashomo_matrix.shape)
        ltemp = [0] * len(self.homo_namelist) if not byspecie else [0] * len(self.species_namelist)
        self.clusters = [0] * len(clust)
        begin = 0
        # reorder all matrices
        species_namelist = list(self.species_namelist)
        for i in range(homogroupnb):
            ind = np.argwhere(clust == i)[:, 0]
            orderedhas[begin:begin + len(ind)] = self.hashomo_matrix[ind]
            # the list as well
            sublist = [self.homo_namelist[e] for e in ind] if not byspecie else [species_namelist[e] for e in ind]
            ltemp[begin:begin + len(ind)] = sublist
            self.clusters[begin:begin + len(ind)] = clust[ind]
            begin += len(ind)

        if byspecie:
            self.species_namelist = set(ltemp)
        else:
            self.homo_namelist = ltemp
        if plot_ordering:
            self._plot_clust(self.hashomo_matrix, orderedhas)
        return orderedhas

    def get_clusterstats(self, by_proportion=True, sort=True):
        """
        will find the number of cluster i per homologies and per species

        plot for each species, how much its genes are outliers, how much are belonging
        to a secondary cluster and how much are belonging to the principal cluster.
        --> create a long stacked bar plot with these values
        Params:
        ------
        by_proportion: bool to compute the prop or the raw numbers
        sort: to sort by 'species','homologies' or None
        """
        # TODO: totest
        homoclusters = np.array((len(self.homodict), 3))
        specluster = np.array((len(self.species_namelist), 3))
        if self.wasclusterized:
            if self.stats is None:
                for j, val in enumerate(self.homodict.iteritems()):
                    x = 0
                    hprop = np.zeros(3)
                    for i in range(len(val.clusters)):
                        if val.clusters[i] == 0:  # primary clusters
                            hprop[0] += 1
                            specluster[val.names[x]][0] += 1
                        elif val.clusters[i] == -1:  # outliers clusters
                            hprop[1] += 1
                            specluster[val.names[x]][1] += 1
                        else:   # secondary clusters
                            hprop[2] += 1
                            specluster[val.names[x]][2] += 1
                        if val.doub[i + 1] is False:
                            x += 1
                    homoclusters[j] = np.divide(hprop, i) if by_proportion else hprop
                self.stats = {'homologies': homoclusters, 'species': specluster}
                if sort:
                    ind = self.stats['homologies'].argsort(axis=0)
                    self.stats['homologies'][:] = self.stats['homologies'][ind]
                    if self.hashomo_matrix is not None:
                        self.hashomo_matrix[:] = self.hashomo_matrix[ind]
                    if self.homo_namelist is not None:
                        self.homo_namelist[:] = self.homo_namelist[ind]
                    if self.homo_matrix is not None:
                        self.homo_matrix[:] = self.homo_matrix[ind]
                    if self.red_homomatrix is not None:
                        self.red_homomatrix[:] = self.red_homomatrix[ind]
                    if self.homo_matrixnames is not None:
                        self.homo_matrixnames[:] = self.homo_matrixnames[ind]
                    if self.clusters is not None:
                        self.clusters[:] = self.clusters[ind]
                    self.stats["species"] = self.stats["species"].tolist()
                    self.stats["species"].sort()
            self._barplot()
        else:
            print "you need to find clusters first"

    def compare_clusters(self, cubdistance_matrix=True, plot=True, size=40):
        """
        for each clusters in homologies, will compare them with a similarity matrix
        and a distance matrix

        compare amongst the working homoset homologies, the clusters together,
        by what species they contains by creating a new vector of species presence
        in each cluster and plotting the similarity matrix of each of those vectors.
        --> create a compare function in homoset of
        homologies clusters similarity matrix and ordering.

        basically the distance should be nan if it has not the species,
        -1 if outlier to other and 1 if one cluster to another and zeros if the same to the same
        """
        # TODO: to test
        se = np.zeros(len(self.species_namelist), dtype=int)
        j = 0
        simimatrix = np.zeros((len(self.homodict), len(self.homodict)))
        for _, homo in self.homodict.iteritems():
            a = set(homo.names)
            se[homo.names] = homo.clusters[:]
            comp = np.zeros(len(self.species_namelist), dtype=int)
            i = 0
            for _, homocmp in self.homodict.iteritems():
                similar = len(a.intersection(homocmp.names))
                comp[homo.names] = homocmp.clusters[:]
                simimatrix[j, i] = ((se == comp).sum() - (len(self.species_namelist) - similar)) / similar
                i += 1
            j += 1
        if plot:
            plt.figure(figsize=(size, 200))
            plt.title('similarity amongst clusters per homologies')
            plt.imshow(simimatrix)
            plt.savefig("utils/clustersimilarity.pdf")
            plt.show()
        self.cluster_similarity = simimatrix.copy()
        if cubdistance_matrix:
            se = np.zeros((len(self.species_namelist), 18))
            j = 0
            simimatrix = np.zeros((len(self.homodict), len(self.homodict)))
            for _, homo in self.homodict.iteritems():
                a = set(homo.names)
                se[homo.names] = homo.full[:]
                comp = np.zeros((len(self.species_namelist), 18))
                i = 0
                for _, homocmp in self.homodict.iteritems():
                    similar = len(a.intersection(homocmp.names))
                    comp[homo.names] = homocmp.full[:]
                    simimatrix[j, i] = np.NaN if similar == 0 else np.sqrt(np.einsum('ij,ij->i',
                                                                                     se, comp)).sum() / similar
                    i += 1
                j += 1
            if plot:
                plt.figure(figsize=(size, 200))
                plt.title('similarity amongst clusters per homologies')
                plt.imshow(simimatrix)
                plt.savefig("utils/CUBsimilarity.pdf")
                plt.show()
        self.cub_similarity = simimatrix

    def find_clusters(self, clustering='dbscan', homogroupnb=None,
                      assess=True, eps=0.8, best_eps=True, trainingset=30):
        """
        Finds, for each homologies in the working homoset, groups that are part of compact clusters
        it will be using gaussian mixture clustering or DBSCAN and order them according
        to the density of each cluster (we are interested in the densest ones) and assess
        the quality using 3 criterion:BIC, number of outliers,

        also: - find if we are close to ancestry tree,
        here we need to represent a comparison of the closeness
        in a phylogenetic tree to a cluster of species
        --> given a grouping of phylogenetic tree, what cluster is the most similar to it


        Params:
        -----
        clustering: method (DBSCAN, gaussian mixture)
        homogroupnb: the number of groups can be a number or else will look for the better number of
        cluster according to assessments.
        assess: plot or not the assessments
        """
        # TODO: totest
        if best_eps:
            eps = self.findbest_eps(trainingset)
        for val in self.homo_namelist:
            print val
            self.homodict[val].clusterize_(clustering=clustering, eps=eps,
                                           homogroupnb=homogroupnb, assess=assess)
            print "-----------------------------"
        self.wasclusterized = True
        if best_eps:
            score = 0
            for _, homo in self.homodict.iteritems():
                score += homo.metrics["avg_phylodistance"][2:].mean() - (homo.metrics["avg_phylodistance"][0] - 1)
            print "-----------------TOT------------------"
            print score

    def findbest_eps(self, trainingset=400, clustering="dbscan", homogroupnb=None,
                     epoch=20, ranges=[0.2, 0.9], size=40):
        """
        will find the best eps hyperparameter (the one that minimizes the evolutionary distance
        within its clusters)

        Params:
        ------
        """
        score = 0
        scores = []
        best_score = 10000000
        best_eps = 0
        val = ranges[0]
        if utils.phylo_distances is None:
            raise LookupError("you need to compute the phylogenetic distances first to have some form of labels")
        for _ in range(epoch):
            i = 0
            for _, homo in self.homodict.iteritems():
                if i > trainingset:
                    break
                homo.clusterize_(clustering=clustering,
                                 eps=val, homogroupnb=homogroupnb, assess=True)
                score += homo.metrics["avg_phylodistance"][2:].mean() - (homo.metrics["avg_phylodistance"][0] - 1)
                i += 1
            # get the best eps of all
            # compute its similarity on the full set of homologies
            val += (ranges[1] - ranges[0]) / epoch
            scores.append(score)
            if score < best_score:
                best_score = score
                best_eps = val
        plt.figure(figsize=(size, 200))
        plt.title('scores for best clustering hyperparameter')
        plt.plot(scores)
        plt.savefig("utils/scores_" + clustering + "_" + str(epoch) + "epochs.pdf")
        plt.show()
        return best_eps

    def similarity_perspecies(self):
        """
        """
    # TODO have a function to compute similarity matrix for each genes of a species. and put it in
    # espece.py

    def get_effectonvalues(self, values=['doublons,nans']):
        """
        find if having doublons affects the values by comparing
        the mean position of doublons (if higher than average) and plotting them

        Here 0.5 means 0 effect, 1 and 0 are positive and negative effects
        """
        # TODO: totest
        doubpositions = []
        nanpositions = []
        keys = []
        for key, val in self.homodict.iteritems():
            posdict = val.getdoubpos()
            # we might encounter a triplon which is counted as two doublons here. should not impact
            # the answer to the fundamental question here.
            if posdict:
                keys.append(key)
                for subkey, subval in posdict.iteritems():
                    doubpositions.append(subval)
            posdict = val.getnanpos()
            # we might encounter a triplon which is counted as two doublons here. should not impact
            # the answer to the fundamental question here.
            if posdict:
                keys.append(key)
                for subkey, subval in posdict.iteritems():
                    nanpositions.append(subval)
        # TODO: have a plot with a cursor to choose one of the dimensions
        print '----------------------------------'
        print "the effect of doublons of the entropy is of: "
        print np.array(doubpositions).mean(0)
        print "the effect of nans of the entropy is of: "
        print np.array(nanpositions).mean(0)


# ###################################################################

    def _barplor(self):
        """
        called by statistics function to plot a barplot of
        the proportion of different cluster per homologies
        and per species
        """
        # TODO: totest
        dims = dict(kdims='prop', vdims='homologies')
        primary = hv.Area(self.stats[0], label='primary', **dims)
        outlier = hv.Area(self.stats[1], label='outlier', **dims)
        secondary = hv.Area(self.stats[2], label='secondary', **dims)
        # could add xaxis = self.homoname_list
        overlay = (primary * outlier * secondary).options('Area', fill_alpha=0.5, name='testdata')
        hv.Area.stack(overlay).relabel("cluster info: primary_clust +outlier +secondary_clust")

    def _dictify(self):
        """
        Used by the saving function. transform the object into a dictionary that can be
        json serializable
        """
        dictihomo = {}
        homodict = self.homodict.iteritems()
        for key, val in homodict:
            dictihomo.update({key: val._dictify()})
        return {"hashomo_matrix": self.hashomo_matrix.tolist() if self.hashomo_matrix is not None else None,
                "homo_matrix": self.homo_matrix.tolist() if self.homo_matrix is not None else None,
                "clusters": self.clusters,
                "homo_namelist": self.homo_namelist,
                "homodict": dictihomo,
                "species_namelist": list(self.species_namelist),
                "speciestable": utils.speciestable,
                "homogroupnb": self.homogroupnb,
                "datatype": self.datatype,
                "phylo_distance": self.phylo_distance.to_json(orient='split') if utils.phylo_distance is not None else None,
                "red_homomatrix": self.red_homomatrix.tolist() if self.red_homomatrix is not None else None,
                "fulleng": self.fulleng.tolist() if self.fulleng is not None else None,
                "wasclusterized": self.wasclusterized,
                "homo_clusters": self.homo_clusters.tolist() if self.homo_clusters is not None else None,
                "averagehomo_matrix": self.averagehomo_matrix.tolist() if self.averagehomo_matrix is not None else None}

    def plot(self, invert=False, size=40, interactive=False):
        """
        plot the has homo matrix

        Params:
        -------
        invert: bool to invert the ordering of the matrix
        size: size of the matrix

        """
        # TODO: code the intereactive version where you can hover
        # over particular lines and interact with them (seeing what homologies they are/ what species)
        plt.figure(figsize=(size, 200))
        plt.title('the regular matrix')
        plt.imshow(self.hashomo_matrix.T if invert else self.hashomo_matrix)
        plt.savefig("utils/hashomomatrix.pdf")
        plt.show()

    def _plot_clust(self, mat, orderedmat):
        """
        will plot the correlation matrix of the has_homomatrix before and after ordering
        allowing one to show its effect
        """
        # the regular matrix
        plt.figure(figsize=(40, 40))
        plt.title('the regular matrix')
        plt.imshow(mat.T)
        plt.show()
        # the affinity matrix
        mat_sparse = sparse.csr_matrix(mat)
        similarities = cosine_similarity(mat_sparse)
        plt.figure(figsize=(40, 40))
        plt.title('the affinity matrix')
        plt.imshow(similarities)
        plt.show()

        # the ordered matrix
        plt.figure(figsize=(40, 40))
        plt.title('the ordered matrix')
        plt.imshow(orderedmat.T)
        plt.show()

        # affinity of the ordered matrix
        mat_sparseo = sparse.csr_matrix(orderedmat)
        similaritieso = cosine_similarity(mat_sparseo)
        plt.figure(figsize=(40, 40))
        plt.title('the affinity of the ordered matrix')
        plt.imshow(similaritieso)
        plt.show()
