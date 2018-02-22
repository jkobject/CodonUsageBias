""""
Created by Jérémie KALFON 
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
    clusters = [][]
    full = None
    homodict = {}

    def __init__(self):

    def plot_all():

    def gene_plot(species):
        pass

    def save(filename,):

    def clusterize():
        pass

    def function():
        pass
