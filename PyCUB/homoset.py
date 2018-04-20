""""
Created by Jeremie KALFON
Date : 21 FEV 2018
University of Kent, ECE paris
jkobject.com

"""


import numpy as np
import matplotlib.pyplot as plt

import utils
import homology as h

from sklearn.cluster import KMeans
from sklearn.cluster import MiniBatchKMeans
from sklearn.metrics.pairwise import cosine_similarity
from kmodes.kmodes import KModes
from sklearn import manifold as man
from sklearn import metrics


from scipy import sparse

import collections

import pdb


class HomoSet(collections.MutableMapping):
    """docstring for HomoSet

                Object where we store an homology group basically where we do our entire
                Computation from.

                params:
                ------
                hashomo_matrix : a numpy boolean array (species, homologies) that store the matrix of gene presence in species
                homo_matrix : a numpy array (homologies*species(inthehomologies),aminoacids)
                                    but containing the codon entropy vectors instead
                homodict = dictionnary of homology object containing codon usage per species 
                    for the gene coresponding to this homology
                homo_namelist : list of all the homology names
                species_namelist: list of all the species names
                clusters: a list of clusters for the homology clustering can be of size of nb of species
                    or of nb of homologies
                homogroupnb: number of group for the homology clustering
                red_homomatrix: numpy array (homologies*species(inthehomologies),x*y)the reduced 2D version 
                    of the homology matrix for the full homology matrix
    """

    hashomo_matrix = None
    homo_matrix = None
    homodict = None
    homo_namelist = []
    species_namelist = []
    clusters = []
    homogroupnb = 2
    red_homomatrix = None

    def __init__(self, data=None):
        """
        will..
        """
        if data is not None:
            self.homo_matrix = np.asarray(data["homo_matrix"]) if not (
                data["homo_matrix"] is None) else None
            self.homo_namelist = data["homo_namelist"]
            self.species_namelist = data["species_namelist"]
            self.homogroupnb = data["homogroupnb"]
            self.clusters = data["clusters"]
            self.red_homomatrix = data.get("red_homomatrix", None)
            utils.speciestable = data["speciestable"]
            self.hashomo_matrix = np.asarray(data["hashomo_matrix"]) if not (
                data["hashomo_matrix"] is None) else None
            tp = {}
            for key, val in data["homodict"].iteritems():
                tp.update({key: h.homology(data=val)})
            self.homodict = tp
    # magic methods https://rszalski.github.io/magicmethods/

    def __getitem__(self, key):
        return self.homodict[key]

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

    def plot_all(self, perplexity=60, interactive=False, iteration=300):
        """
        will..
        """
        # here alg == 'tsne':
        if self.homo_matrix is None:
            self.loadfullhomo()

        if self.red_homomatrix is None:
            red = man.TSNE(n_components=2, perplexity=perplexity, n_iter=iteration, verbose=3).fit_transform(self.homo_matrix)
            self.red_homomatrix = red
        if not interactive:
            fig = plt.figure(figsize=(40, 40))
            ax = fig.add_subplot(111)
            ax.scatter(self.red_homomatrix[:, 0], self.red_homomatrix[:, 1])
            plt.show()
        else:
            print " not yet coded"
            # TODO: code this

    def loadfullhomo(self):
        self.homo_matrix = self.homodict[self.homo_namelist[0]].full
        for x, homo in enumerate(self.homo_namelist[1:]):
            self.homo_matrix = np.vstack((self.homo_matrix, self.homodict[homo].full))
            print '{0}%\r'.format((x * 100) / len(self.homo_namelist)),

    def loadhashomo(self, isbyhomo):
        # TODO: to test
        hashomo = np.zeros((len(self.homo_namelist),
                            len(self.species_namelist)), dtype=np.bool)
        if isbyhomo:
            hashomo = hashomo.T
        inde = {}
        for i, n in enumerate(self.species_namelist):
            inde.update({n: i})
        for homo in self.homo_namelist:
            ind = [inde[utils.speciestable[na]] for na in self.homodict[homo].names]
            # try to get the indeces from each species name
            hashomo[ind] = True

    def find_clusters(self, clustering='dbscan', homogroupnb=None, assess=True):
        """
        Finds, for each homologies in the working homoset, groups that are part of compact clusters
        it will be using gaussian mixture clustering or DBSCAN and order them according
        to the density of each cluster (we are interested in the densest ones) and assess
        the quality using 3 criterion:BIC, number of outliers, .

        Params:
        -----
        clustering: method (DBSCAN, gaussian mixture)
        homogroupnb: the number of groups can be a number or else will look for the better number of
        cluster according to assessments.
        plotassess: plot or not the assessments
        """
        for val in self.homo_namelist:
            self.homodict[val].clusterize_(clustering=clustering, homogroupnb=homogroupnb, assess=assess)

    def remove(self, species):
        """
        remove this list of species from the homologies
        """
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

    def order_from_matrix(self, clustering='kmeans', plot_ordering=True, homogroupnb=2, findnb=False,
                          reducedim=False, orderspecie=False):
        """
        Compute an homology group :
        from matrix computation using the homo_matrix
        (or from network computation in homologize_from_network)

        Can be computed many times and will updata homoset with the most recent homoset found
        if homoset exists, it will save it.

        Params:
        -------
        clustering: flags to 'kmeans', 'spectral', 'xmeans' to use different sk-learn algorithms

        plot: flags to true for the function to output ploting of the affinity matrix with and without the
        clusters

        homogroupnb: nb of groups you want to extract

        """
        # TODO: test all the new clusterings
        # TODO: find more tests and visualisations to perform
        if reducedim:
            # TODO: try to reduce the dimensionality here before clustering
            print "tocode yet"
        if findnb is True:
            homogroupnb = 2
        clusts = []
        while True:
            if clustering == "kmeans":
                # The parallel version of K-Means is broken on OS X when numpy uses the Accelerate Framework.
                # This is expected behavior: Accelerate can be called after a fork
                # but you need to execv the subprocess with the Python binary
                # (which multiprocessing does not do under posix)
                alg = KMeans(n_clusters=homogroupnb, n_jobs=-1)
            elif clustering == "fast":
                alg = MiniBatchKMeans(n_clusters=homogroupnb)

            elif clustering == "kmodes":
                # https://github.com/nicodv/kmodes/blob/master/kmodes/kmodes.py
                alg = KModes(n_clusters=homogroupnb,
                             init='Huang', n_init=2, verbose=1)

            else:
                print "you entered a wrong clustering algorithm"
                return False

            if not orderspecie:
                # special case where one would like to keep all homologies and
                # remove species instead. we would advise the opposite though
                self.hashomo_matrix = self.hashomo_matrix.T

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
                    print "you just select the clust by doing PyCUB.homo.clusters = output['groupnbyouwant']"
                    print "you an run the 'orderfromclust' also to have the plottings and orderings"
                    return clusts
            else:
                break
        print "the quality of the clustering is: "
        print metricA
        print metricB
        if plot_ordering:
            print "plotting the ordering... might take time"
            self.orderfromclust(homogroupnb, clust, fromspecies=orderspecie)
        self.homogroupnb = homogroupnb
        return True

    def orderfromclust(self, homogroupnb, clust, fromspecies):
        orderedhas = np.zeros(self.hashomo_matrix.shape)
        orderedfull = np.zeros(self.hashomo_matrix.shape) if not (
            self.homo_matrix is None) else None
        ltemp = [0] * len(self.homo_namelist) if not fromspecies else [0] * len(self.species_namelist)
        self.clusters = [0] * len(clust)
        begin = 0
        # reorder all matrices
        for i in range(homogroupnb):
            ind = np.argwhere(clust == i)[:, 0]
            orderedhas[begin:begin + len(ind)] = self.hashomo_matrix[ind]
            if self.homo_matrix is not None:
                orderedfull[begin:begin + len(ind)] = self.homo_matrix[ind]
            # the list as well
            sublist = [self.homo_namelist[e] for e in ind] if not fromspecies else [self.species_namelist[e] for e in ind]
            ltemp[begin:begin + len(ind)] = sublist
            self.clusters[begin:begin + len(ind)] = clust[ind]
            begin += len(ind)

        self._plot_clust(self.hashomo_matrix, orderedhas)
        if fromspecies:
            self.species_namelist = ltemp
        else:
            self.homo_namelist = ltemp
        self.hashomo_matrix = orderedhas
        self.homo_matrix = orderedfull


###################################################################

    def _dictify(self):
        """

        Used by the saving function. transform the object into a dictionary that can be
        json serializable

        """
        dictihomo = {}
        for key, val in self.homodict.iteritems():
            dictihomo.update({key: val._dictify()})
        return {"hashomo_matrix": self.hashomo_matrix.tolist() if self.hashomo_matrix is not None else None,
                "homo_matrix": self.homo_matrix.tolist() if self.homo_matrix is not None else None,
                "clusters": self.clusters,
                "homo_namelist": self.homo_namelist,
                "homodict": dictihomo,
                "species_namelist": self.species_namelist,
                "speciestable": utils.speciestable,
                "homogroupnb": self.homogroupnb}

    def _plot_clust(self, mat, orderedmat):
        # the regular matrix
        plt.figure(figsize=(40, 40))
        plt.title('the regular matrix')
        plt.imshow(mat)
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
        plt.imshow(orderedmat)
        plt.show()

        # affinity of the ordered matrix
        mat_sparse = sparse.csr_matrix(orderedmat)
        similarities = cosine_similarity(mat_sparse)
        plt.figure(figsize=(40, 40))
        plt.title('the affinity of the ordered matrix')
        plt.imshow(similarities)
        plt.show()
