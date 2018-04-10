""""
Created by Jeremie KALFON
Date : 21 FEV 2018
University of Kent, ECE paris
jkobject.com

"""
import pandas as pd
import glob
import numpy as np


class homology(object):
    """in homology we store an homology with its dataframe, 
    it reduced reduced df with dim reduction and its clusters
    it is supposed to be store in a homoogy dictionary 

    Params:
    ------
    full : df of one homology with entropy value vector per species
    reduced :  df of one homology dimensionality reduced
    clusters : list of cluster labels
    """

    reduced = False
    clusters = False
    full = False

    def __init__(self, data=False, full=False):
        """
        will intialize an instance of the object and be used for the loading mechanism
        """
        if not (type(data) is bool):
            self.reduced = pd.DataFrame.from_dict(data["reduced"]) if not (type(data["reduced"]) is bool) else False
            self.clusters = data["clusters"]
            self.full = pd.DataFrame.from_dict(data["full"]) if not (type(data["full"]) is bool) else False

        elif not (type(full) is bool):
            self.full = full

    def clusterize_(self, clustering='gaussian', plot=False, homogroupnb=4):
        """
        will clusterize the homology using gaussian mixture clustering or DBSCAN and order them according 
        to the density of each cluster (we are interested in the dense ones) and assess the quality using 3 criterion:
        BIC, , . 
        """
        if clustering == 'gaussian':

        elif clustering == 'dbscan':

        return clusters

    def assess_clust():
        """
        Assesses the clusterization quality with different metrics (BIC, )
        """

    def reduce_dim(self, alg='tsne', n=2, perplexity=40):
        """
        reduce the dimensionality of your gene dataset to a defined dimension 
        using the t-SNE algorithm

        :param gene:  a matrix of gene codon usage per species
                  n:  the desired dimension
        perplexity :  an optional value when you know about tsne 

        :return tsned: the reduced dataset
        """
        if alg == 'tsne':
            red = man.TSNE(n_components=n, perplexity=perplexity).fit_transform(self.scaled)
        elif alg == 'pca':
            break
        self.reduced = red
        return red

    def plot():
        """
        will plot the object's data using bokeh to create nice interactive plots (either on jupyter notebook or as
         html-javascripts pages)
        in these plots of your gene dataset you can have a look at your previously 
        dimensionality reduced gene dataset and hover over points to have a look at the species or display 
        colors and centroids according to the clusterize function's output

        contains an option for simple plots.

        :param gene:  a matrix of gene codon usage per species reduced into 2 dimensions
               species : a list of all the species in your gene dataset (usually one per row)
               getimage: (optional) flag to true if you want it to ouput an html file 
               labels, centroids : (optional) the labels and centroids returned
                by the clusterize function

        :return p: your figure as a bokeh object.
        """
        if centroids.any() and labels.any():
            #colormap = [[rand(256), rand(256), rand(256)] for _ in range(100)]
            colormap = ["#1abc9c", "#3498db", "#2ecc71", "#9b59b6", '#34495e', '#f1c40f',
                        '#e67e22', '#e74c3c', '#7f8c8d', '#f39c12']
            colors = [colormap[x] for x in labels]
        else:
            colors = '#1abc9c'

        source = ColumnDataSource(
            data=dict(
                x=tsnedgene[:, 0],
                y=tsnedgene[:, 1],
                color=colors,
                label=["species : %s" % (x_) for x_ in species]

            )
        )
        hover = HoverTool(tooltips=[
            ("label", "@label"),
        ])
        p = figure(title="T-sne of homologous gene X for each species",
                   tools=[hover, BoxZoomTool(), WheelZoomTool(), SaveTool(), ResetTool()])
        p.circle('x', 'y', size=10, source=source)
        show(p)
        return p

    def _dictify(self):
        """
        Used by the saving function. transform the object into a dictionary that can be 
        json serializable
        """
        return {"reduced": self.reduced.to_dict() if not (type(self.reduced) is bool) else False,
                "clusters": self.clusters,
                "full": self.full.to_dict() if not (type(self.full) is bool) else False}
