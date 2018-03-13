""""
Created by Jeremie KALFON 
Date : 21 FEV 2018
University of Kent, ECE paris
jkobject.com

"""


import pandas as pd
import json


class HomoSet(object):
    """docstring for HomoSet

                Object where we store an homology group basically where we do our entire 
                Computation from. 

                params: 
                ------
                reduced = PD.DF of genes *species  containing 2D vectors ( dimensionality reduction using T-SNE)
                clusters =  matrix of genes * species containing  cluster index ( note the same between genes)
                full = PD.DF of genes *species  containing the entire 18D vectors of their entropy location
                homodict = dictionnary of dataframes of codon usage per species from homology names
    """

    reduced = None
    clusters = None
    full = None
    homodict = {}

    def __init__(self, data=False):
        """
        will..
        """
        if data:
            reduced = pd.from_dict(data["reduced"])
            clusters = np.array(data["clusters"])
            full = pd.from_dict(data["full"])
            for key, val in data["homodict"].iteritems():
                homodict.update({key: pd.from_dict(val)})

    def plot_all():
        """
        will..
        """

    def gene_plot(species):
        """
        will..
        """
        pass

    def save(filename):
        """
        will save the object as a json string 
        """

    def clusterize_kmeans():
        """
        will..
        """
        pass

    def clusterize_gaussmixture():
        """
        will..
        """

    def clusterize_affinityprop():
        """
        will..
        """

    def clusterize_DB_scan():
        """
        will..
        """

    def assess_clust():
        """
        will..
        """

    def dictify():
        """
        will..
        """
        dictihomo = {}
        for key, val in self.homodict.iteritems():
            dictihomo.update({key: val.to_dict()})
        return {"reduced": reduced.to_dict(),
                "clusters": clusters.tolist(),
                "full": full.to_dict(),
                "homodict": dictihomo}
