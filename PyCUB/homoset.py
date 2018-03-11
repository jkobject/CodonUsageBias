""""
Created by Jeremie KALFON 
Date : 21 FEV 2018
University of Kent, ECE paris
jkobject.com

"""


import pandas as pd


class HomoSet(object):
    """docstring for HomoSet

                Object where we store an homology group basically where we do our entire 
                Computation from. 

                params: 
                ------
                reduced = PD.DF of genes *species  containing 2D vectors ( dimensionality reduction using T-SNE)
                clusters =  matrix of genes * species containing  cluster index ( note the same between genes)
                full = PD.DF of genes *species  containing the entire 18D vectors of their entropy location
    """

    reduced = None
    clusters = None
    full = None
    homodict = {}

    def __init__(self):
        """
        will..
        """

    def plot_all():
        """
        will..
        """

    def gene_plot(species):
        """
        will..
        """
        pass

    def save(filename,):
        """
        will..
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

    def function():
        """
        will..
        """
        pass
