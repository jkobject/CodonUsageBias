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
from scipy import sparse


class HomoSet(object):
    """docstring for HomoSet

                Object where we store an homology group basically where we do our entire
                Computation from.

                params:
                ------
                has_homo_matrix : a numpy boolean array that store the matrix of gene presence in species
                full_homo_matrix : a numpy array similar to has_homo_matrix
                                    but containing the codon entropy vectors instead
                homodict = dictionnary of dataframes of codon usage per species from homology names
                homo_namelist : list of all the homology names
    """

    has_homo_matrix = False
    full_homo_matrix = False
    homodict = {}
    homo_namelist = []
    clusters = []
    homogroupnb = 2

    def __init__(self, data=False):
        """
        will..
        """
        if data:
            self.full_homo_matrix = np.asarray(data["full_homo_matrix"]) if not (type(data["full_homo_matrix"]) is bool) else False
            self.homo_namelist = data["homo_namelist"]
            self.homogroupnb = data["homogroupnb"]
            self.clusters = data["clusters"]
            self.has_homo_matrix = pd.DataFrame.from_dict(data["has_homo_matrix"])
            for key, val in data["homodict"].iteritems():
                self.homodict.update({key: h.homology(data=val)})

    def plot_all():
        """
        will..
        """

    def plot_homo_per_species(self):
        """
        will plot the number of homology per spcies 
        """

        sumed = np.sum(self.has_homo_matrix, axis=1)
        plt.figure(figsize=(40, 10))
        plt.title('number of homologies per species')
        plt.bar(range(len(sumed)), sumed)
        print "you can always look at a particular range of species with 'homo_namelist' "

    def gene_plot(species):
        """
        will..
        """
        pass

    def order_from_matrix(self, clustering='kmeans', plot=True, homogroupnb=2):
        """
        Compute an homology group :
        from matrix computation using the full_homo_matrix
        (or from network computation in homologize_from_network)

        Can be computed many times and will updata homoset with the most recent homoset found
        if homoset exists, it will save it.

        Params:
        -------
        clustering: flags to 'kmeans', 'spectral', ... to use different sk-learn algorithms

        plot: flags to true for the function to output ploting of the affinity matrix with and without the
        clusters

        homogroupnb: nb of groups you want to extract

        """

        def compute(homo_namelist, mat, homogroupnb, clust, names):
            orderedmat = np.zeros(mat.shape)
            ltemp = [0] * len(names)
            begin = 0
            for i in range(homogroupnb):
                a = np.argwhere(clust == i)[:, 0]
                orderedmat[begin:begin + len(a)] = mat[a]
                c = [names[i] for i in a]
                ltemp[begin:begin + len(a)] = c
                begin += len(a)
            names = ltemp
            if plot:
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
            return pd.DataFrame(orderedmat, names, homo_namelist)

        mat = self.has_homo_matrix.as_matrix()
        names = self.has_homo_matrix.index.tolist()
        if clustering == "spectral":
            alg = SpectralClustering(n_clusters=homogroupnb, n_jobs=-1)

        elif clustering == "kmeans":
            alg = KMeans(n_clusters=homogroupnb, n_jobs=-1)

        elif clustering == "fast":
            alg = MiniBatchKMeans(n_clusters=homogroupnb)

        else:
            print "you entered a wrong clustering algorithm"
            return False

        alg.fit(mat)
        clust = alg.labels_
        self.has_homo_matrix = compute(self.homo_namelist, mat, homogroupnb, clust, names)
        self.homogroupnb = homogroupnb
        self.clusters = clust.tolist()
        return True

    def get_homoset(self, size, clustering=True, clusternb=2, max_clique=False, per_specie=True,
                    per_gene=True, name=[]):
        """
        To use once you have clustered homology groups, else takes everything. 


        Params:
        ------
        number: set the number of the group you want to get need to be between 1 and homogroupnb

        size :  a value which will determine the size (max size of the retrived homology)
        50 means the first 50 in the sim
        if similarity alignment

        simi_align: flag for computing the matrix on similarity alignement

        max_clique: flag for computing on the max clique 

        Return:
        ------
        a HomoSet object (see homoset.py)
        """

    def _clusterize_kmeans():
        """
        will..
        """
        pass

    def _clusterize_gaussmixture():
        """
        will..
        """

    def _clusterize_affinityprop():
        """
        will..
        """

    def _clusterize_DB_scan():
        """
        will..
        """

    def _assess_clust():
        """
        will..
        """

    def _dictify(self):
        """
        will..
        """
        dictihomo = {}
        for key, val in self.homodict.iteritems():
            dictihomo.update({key: val._dictify()})
        return {"has_homo_matrix": self.has_homo_matrix.to_dict() if not (type(self.has_homo_matrix) is bool) else False,
                "full_homo_matrix": self.full_homo_matrix.tolist() if not (type(self.full_homo_matrix) is bool) else False,
                "clusters": self.clusters,
                "homo_namelist": self.homo_namelist,
                "homodict": dictihomo,
                "homogroupnb": self.homogroupnb}
