""""
Created by Jeremie KALFON
Date : 21 FEV 2018
University of Kent, ECE paris
jkobject.com

"""


import json
import os
import glob

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import espece as spe
import utils
import homology as h

from sklearn.cluster import SpectralClustering
from sklearn.cluster import KMeans
from sklearn.cluster import MiniBatchKMeans
from sklearn.metrics.pairwise import cosine_similarity
from kmodes.kmodes import KModes

from scipy import sparse


class HomoSet(object):
    """docstring for HomoSet

                Object where we store an homology group basically where we do our entire
                Computation from.

                params:
                ------
                hashomo_matrix : a numpy boolean array that store the matrix of gene presence in species
                homo_matrix : a numpy array similar to hashomo_matrix
                                    but containing the codon entropy vectors instead
                homodict = dictionnary of dataframes of codon usage per species from homology names
                homo_namelist : list of all the homology names
    """

    hashomo_matrix = None
    homo_matrix = None
    homodict = None
    homo_namelist = []
    species_namelist = []
    clusters = []
    homogroupnb = 2

    def __init__(self, data=None):
        """
        will..
        """
        if data not is None:
            self.homo_matrix = np.asarray(data["homo_matrix"]) if not (data["homo_matrix"] is None) else None
            self.homo_namelist = data["homo_namelist"]
            self.species_namelist = data["species_namelist"]
            self.homogroupnb = data["homogroupnb"]
            self.clusters = data["clusters"]
            self.hashomo_matrix = np.asarray(data["hashomo_matrix"]) if not (data["hashomo_matrix"] is None) else None
            for key, val in data["homodict"].iteritems():
                self.homodict.update({key: h.homology(data=val)})

    def plot_all(self):
        """
        will..
        """

    def plot_homo_per_species(self):
        """
        will plot the number of homology per spcies 
        """

        sumed = np.sum(self.hashomo_matrix, axis=1)
        plt.figure(figsize=(40, 10))
        plt.title('number of homologies per species')
        plt.bar(range(len(sumed)), sumed)
        print "you can always look at a particular range of species with 'homo_namelist' "

    def gene_plot(self, species):
        """
        will..
        """
        pass

    def order_from_matrix(self, clustering='kmeans', plot_ordering=True, homogroupnb=2):
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

        def plot(self, homogroupnb, clust):
            orderedhas = np.zeros(self.hashomo_matrix.shape)
            orderedfull = np.zeros(self.hashomo_matrix.shape) if not (self.homo_matrix is None) else None
            ltemp = [0] * len(self.species_namelist)
            self.clusters = [0] * len(clust)
            begin = 0
            # reorder all matrices
            for i in range(homogroupnb):
                a = np.argwhere(clust == i)[:, 0]
                orderedhas[begin:begin + len(a)] = self.hashomo_matrix[a]
                if self.homo_matrix not is None:
                    orderedfull[begin:begin + len(a)] = self.homo_matrix[a]
                # the list as well
                c = [self.species_namelist[e] for e in a]
                ltemp[begin:begin + len(a)] = c
                self.clusters[begin:begin + len(a)] = clust[a]
                begin += len(a)
            self.species_namelist = ltemp
            self.hashomo_matrix = orderedhas
            self.homo_matrix = orderedfull
            _plot_clust(self.hashomo_matrix, orderedhas)

        if clustering == "spectral":
            alg = SpectralClustering(n_clusters=homogroupnb, n_jobs=-1)

        elif clustering == "kmeans":
            # The parallel version of K-Means is broken on OS X when numpy uses the Accelerate Framework.
            # This is expected behavior: Accelerate can be called after a fork
            # but you need to execv the subprocess with the Python binary
            # (which multiprocessing does not do under posix)
            alg = KMeans(n_clusters=homogroupnb, n_jobs=-1)

        elif clustering == "fast":
            alg = MiniBatchKMeans(n_clusters=homogroupnb)

        elif clustering == "kmodes":
            alg = KModes(n_clusters=homogroupnb, init='Huang', n_init=3, verbose=1)

        else:
            print "you entered a wrong clustering algorithm"
            return False

        alg.fit(self.hashomo_matrix)
        clust = alg.labels_
        #---sklearn.metrics.silhouette_score - --
        if plot_ordering:
            print "plotting the ordering... might take time"
            plot(self, homogroupnb, clust)
        self.homogroupnb = homogroupnb
        self.clusters = clust.tolist()
        return True

    def _dictify(self):
        """
        Used by the saving function. transform the object into a dictionary that can be 
        json serializable
        """
        dictihomo = {}
        for key, val in self.homodict.iteritems():
            dictihomo.update({key: val._dictify()})
        return {"hashomo_matrix": self.hashomo_matrix.to_dict() if self.hashomo_matrix is not None else None,
                "homo_matrix": self.homo_matrix.tolist() if self.homo_matrix is not None else None,
                "clusters": self.clusters,
                "homo_namelist": self.homo_namelist,
                "homodict": dictihomo,
                "species_namelist": self.species_namelist,
                "homogroupnb": self.homogroupnb}

    def _plot_clust(mat, orderedmat)
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
