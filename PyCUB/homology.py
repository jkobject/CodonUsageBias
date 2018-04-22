""""
Created by Jeremie KALFON
Date : 21 FEV 2018
University of Kent, ECE paris
jkobject.com

"""

import numpy as np

from sklearn import manifold as man
from sklearn import cluster, mixture
from sklearn import metrics

import pdb
import utils

from bokeh.plotting import *
from bokeh.models import *
import matplotlib.pyplot as plt


class homology(object):
    """in homology we store an homology with its dataframe,
    it reduced reduced df with dim reduction and its clusters
    it is supposed to be store in a homoogy dictionary

    Params:
    ------
    names : list of int corresponding to names
    full : np array (species, amino) of one homology with entropy value vector per species
    reduced :  np array (species, x*y) of one homology dimensionality reduced
    clusters : list of cluster val for each species
    """
    names = None
    reduced = None
    clusters = None
    full = None
    centroids = None
    metrics = {}
    nans = None
    lenmat = None

    def __init__(self, data=None, full=None, names=None, nans=None, lenmat=None):
        """
        will intialize an instance of the object and be used for the loading mechanism
        """
        if not (data is None):
            self.reduced = np.asarray(data["reduced"]) if not (data["reduced"] is None) else None
            self.clusters = data["clusters"]
            self.full = np.asarray(data["full"]) if not (data["full"] is None) else None
            self.names = data["names"]
            self.centroids = data["centroids"]
            self.metrics = data["metrics"]
            self.nans = np.asarray(data["nans"]) if not (data["nans"] is None) else None
            self.lenmat = np.asarray(data["lenmat"]) if not (data["lenmat"] is None) else None
        elif not (full is None):
            self.full = full
            self.names = names
            self.clusters = None
            self.nans = nans
            self.lenmat = lenmat

    def remove(self, species):
        """
        removes the list of species from this homology if it exists there
        """
        names = [utils.speciestable[str(na)] for na in self.names]
        mask = np.ones(len(self.reduced.shape[0]), dtype=bool)
        for spe in species:
            for i, na in enumerates(names):
                if na == spe:
                    self.names.pop(i)
                    mask[i] = False
                    if self.cluster is not None:
                        self.clusters.pop(i)
        self.reduced = self.reduced[mask, :]
        self.full = self.full[mask, :]

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
            red = man.TSNE(n_components=n, perplexity=perplexity).fit_transform(self.full)
        # elif alg == 'pca':

        self.reduced = red
        return red

    def plot(self, per=40, interactive=False):
        if self.clusters is not None:
            # colormap = [[rand(256), rand(256), rand(256)] for _ in range(100)]
            colormap = ["#1abc9c", "#3498db", "#2ecc71", "#9b59b6", '#34495e', '#f1c40f',
                        '#e67e22', '#e74c3c', '#7f8c8d', '#f39c12']
            colors = [colormap[x] for x in self.clusters]
        else:
            colors = '#1abc9c'
        if self.reduced is None:
            self.reduce_dim(perplexity=per)

        if interactive:
            # TODO: debug this part
            # TODO: show the clusterisation in bokeh with centroids
            print " if you are on a notebook you should write 'from bokeh.io import output_notebook'"
            if self.clusters is None:
                colors = [colors] * len(self.names)
            source = ColumnDataSource(data=dict(x=self.reduced[:, 0], y=self.reduced[:, 1],
                                                label=["species : %s" % utils.speciestable[str(x__)] for x__ in self.names], color=colors))
            output_notebook()
            hover = HoverTool(tooltips=[
                ("label", "@label"),
            ])
            p = figure(title="T-sne of homologous gene X for each species",
                       tools=[hover, BoxZoomTool(), WheelZoomTool(), SaveTool(), ResetTool()])
            p.circle(x='x', y='y', source=source, color='color')

            show(p)
            return p
        else:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.scatter(self.reduced[:, 0], self.reduced[:, 1], c=colors)
            plt.show()

    def clusterize_(self, clustering='gaussian', homogroupnb=None, assess=True):
        """
        will clusterize the homology using gaussian mixture clustering or DBSCAN and order them according
        to the density of each cluster (we are interested in the dense ones) and assess the quality using 3 criterion:
        BIC, , .
        """
        # TODO: find more tests and visualisations to perform
        if clustering == 'gaussian':
            # http://scikit-learn.org/stable/modules/generated/sklearn.mixture.GaussianMixture.html
            alg = mixture.GaussianMixture(n_components=homogroupnb, n_init=2, init_params='random')
            alg.fit(self.full)
            if assess:
                aic = alg.aic(self.full)
                bic = alg.bic(self.full)
                self.metrics.update({'aic': aic, 'bic': bic})
                print "the BIC scores for the GMM is"
                print aic
                print "the AIC scores for the GMM is"
                print bic
            self.clusters = alg.predict(self.full)

        elif clustering == 'dbscan' and homogroupnb is not None:
            # http://scikit-learn.org/stable/modules/generated/sklearn.cluster.DBSCAN.html
            val = np.around(len(self.names) / (2 * (homogroupnb + 1)))
            print val
            alg = cluster.DBSCAN(eps=0.8, min_samples=7,
                                 algorithm='auto', n_jobs=-1)
            self.clusters = alg.fit_predict(self.full)
            n_clusters_ = len(set(self.clusters)) - (1 if -1 in self.clusters else 0)
            print "Estimated number of clusters using DBscan: " + str(n_clusters_)
        elif clustering == 'jkmeans':
            # TODO:  create your own clustering algorithm with gaussian kernels that explains as much as possible the data and get a threshold and set all others as outliers
            """
                here we want to find the smallest group of gaussian that explains
                the highest number of data point with the smallest variance possible
                as
            """
            print "tocode"
        if assess:
            try:
                silhouette = metrics.silhouette_score(self.full, self.clusters)
                cal_hara = metrics.calinski_harabaz_score(self.full, self.clusters)
            except ValueError:
                silhouette = 0
                cal_hara = 0
            self.metrics.update({'silhouette': silhouette,
                                 'cal_hara': cal_hara})
        return self.clusters

    def _dictify(self):
        """
        Used by the saving function. transform the object into a dictionary that can be
        json serializable
        """
        return {"reduced": self.reduced.tolist() if not (self.reduced is None) else None,
                "clusters": self.clusters,
                "full": self.full.tolist() if not (self.full is None) else None,
                "names": self.names,
                "centroids": self.centroids,
                "metrics": self.metrics,
                "nans": self.nans.tolist() if self.nans is not None else None,
                "lenmat": self.lenmat.tolist() if self.lenmat is not None else None}
