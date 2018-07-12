""""
Created by Jeremie KALFON
Date : 21 FEV 2018
University of Kent, ECE paris
jkobject.com

"""

import numpy as np

from sklearn import manifold as man
from sklearn.decomposition import PCA
from sklearn import cluster, mixture
from sklearn import metrics
from mpl_toolkits.mplot3d import Axes3D
import utils

from bokeh.plotting import *
from bokeh.models import *
import matplotlib.pyplot as plt

import pdb


class homology(object):
    """in homology we store an homology with its dataframe,
    it reduced reduced df with dim reduction and its clusters
    it is supposed to be store in a homoogy dictionary

    Params:
    ------
    names : list of int corresponding to names in utils.speciestable
    full : np array of float (species, amino) of one homology with entropy value vector per species
    reduced :  np array float (species, x*y) of one homology dimensionality reduced
    clusters : list[int] of cluster val for each species
    centroids: the position in aminoacid# Dimension of each centroids (number of centroid == number of cluster)
    metrics: a dict of metrics names and values for the cluster of this homology
    nans: numpy array of bool wether or not this position is a nan
    lenmat: a numpy array species amino of int of length of each amino acids (number)
    doub: nupy array of [bool] of wether or not this position is a doublon ( a copy number of the gene for the species)
    reduced_algo: dimensionality reduction algorithm used on this homology
    """
    names = None
    reduced = None
    clusters = None
    full = None
    centroids = None
    metrics = {}
    nans = None
    lenmat = None
    doub = None
    reduced_algo = None
    KaKs_Scores = None
    similarity_scores = None
    proteinids = []
    # TODO: compute the mean variance in CUB value and mean range for each homology (add this when plotting)
    # TODO: compute a distance matrix between the species
    # Replace long params by **kwargs (dict for the data)

    def __init__(self, data=None, full=None, names=None, nans=None, lenmat=None, doub=None,
                 similarity_scores=None, KaKs_Scores=None, proteinids=None, clusters=None):
        """
        will intialize an instance of the object and be used for the loading mechanism
        """
        if data is not None:
            self.reduced = np.asarray(data["reduced"]) if not (data["reduced"] is None) else None
            self.clusters = data["clusters"]
            self.full = np.asarray(data["full"]) if not (data["full"] is None) else None
            self.names = data["names"]
            self.centroids = data["centroids"]
            self.metrics = data["metrics"]
            self.nans = np.asarray(data["nans"]) if not (data["nans"] is None) else None
            self.lenmat = np.asarray(data["lenmat"]) if not (data["lenmat"] is None) else None
            self.doub = np.asarray(data["doub"]) if not (data["doub"] is None) else None
            self.KaKs_Scores = np.asarray(data(["KaKs_Scores"])) if (data["KaKs_Scores"])\
                is not None else None
            self.similarity_scores = np.asarray(data(["similarity_scores"])) if\
                (data["similarity_scores"]) is not None else None
            self.proteinids = data(["proteinids"])
        if full is not None:
            self.full = full
        if names is not None:
            self.names = names
        if clusters is not None:
            self.clusters = None
        if nans is not None:
            self.nans = nans
        if lenmat is not None:
            self.lenmat = lenmat
        if doub is not None:
            self.doub = doub
        if KaKs_Scores is not None:
            self.KaKs_Scores = KaKs_Scores
        if proteinids is not None:
            self.proteinids = proteinids
        if similarity_scores is not None:
            self.similarity_scores = similarity_scores

    def remove(self, species):
        """
        removes the list of species from this homology if it exists there
        """
        names = [utils.speciestable[na] for na in self.names]
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

    def nb_unique_species(self):
        """
        compute the number of unique species in this homologies
        (basically count the number of doub)
        """
        return self.doub.shape[0] - np.count_nonzero(self.doub)

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
        elif alg == 'pca':
            red = PCA(n_components=n).fit_transform(self.full)

        self.reduced = red
        self.reduced_algo = alg
        return red

    def plot(self, reducer="tsne", per=40, interactive=False, D=2, size=40, shownans=True):
        """
        plot the dimensionality reduced datapoint of the homology matrix

        colors represents the different clusters

        Params:
        ------
        reducer: flag the algorithm to reduce the matrix of 18D to D D
            per: int he perplexity when using tsne algorithm
        size: int the size of the plot
        interactive: (recommended) wether to use bokeh or not to represend the data
                    will show name of the species when hovering over datapoint
                    will show a gradient of evolutionary distance of the species of
                    the datapoint currently hovered to each other species
            D: int the goal dimension (2-3-4) when using matplotlib
        """
        # TODO: use the ancestry distance from homoset (put it in utils)
        # and add the possibility to plot it with colors
        if self.clusters is not None:
            # colormap = [[rand(256), rand(256), rand(256)] for _ in range(100)]
            colormap = ["#1abc9c", "#3498db", "#2ecc71", "#9b59b6", '#34495e', '#f1c40f',
                        '#e67e22', '#e74c3c', '#7f8c8d', '#f39c12']
            colors = [colormap[x] for x in self.clusters]
        else:
            colors = '#1abc9c'
        if shownans:
            for i, val in enumerate(self.nans):
                if val:
                    colormap[i] = "#000000"
        if (self.reduced is None or self.reduced_algo != reducer):
            self.reduce_dim(alg=reducer, n=D, perplexity=per)
        elif self.reduced.shape[1] == D:
            self.reduce_dim(alg=reducer, n=D, perplexity=per)

        if interactive:
            # TODO: to test
            print " if you are on a notebook you should write 'from bokeh.io import output_notebook'"
            if self.clusters is None:
                colors = [colors] * len(self.names)
            source = ColumnDataSource(data=dict(x=self.reduced[:, 0], y=self.reduced[:, 1],
                                                label=["species : %s" % utils.speciestable[x__] for x__ in self.names],
                                                doub=self.doub if self.doub is not None else False,
                                                clusters=self.clusters if self.clusters is not None else False,
                                                similarity_scores=self.similarity_scores if self.similarity_scores is not None else False,
                                                KaKs_Scores=self.KaKs_Scores if self.KaKs_Scores is not None else False,
                                                nans=self.nans if self.nans is not None else False,
                                                color=colors))
            output_notebook()

            radio_button_group = widgets.RadioButtonGroup(
                labels=["show Cluster", "show Doublon", "show KaKs_Scores", "show similarity_scores", "Show PhyloDistance"], active=0)
            callback = CustomJS(args=dict(source=source, colors=colormap), code="""
            // JavaScript code goes here

            // the model that triggered the callback is cb_obj:
            var b = cb_obj.value;
            // models passed as args are automagically available
            var data = source.data;
            var len = data.color.length;
            var col = [];
            if(b === "show Doublon" && data.doub){
                var j = 0;
                for(i=0;i<len;++i){
                    if(doub[i]){
                        col[i-1] = colors[j];
                        col.push(colors[j]);
                    }else{col.push("#999999");}
                }
            }
            if(b === "show Cluster"){
                var j = 0;
                othercol = '#1abc9c'
                if(data.clusters){
                    for(i=0;i<len;++i){col.push(colors[data.clusters[i]]);}
                }else{
                    for(i=0;i<len;++i){col.push(othercol);}
                }
            }
            if(b === "show KaKs_Scores" && data.KaKs_Scores){
                var arr = Array(len);
                while(i--) arr[i] = data.KaKs_Scores[i];
                arr.sort();
                var max = arr[len-1];
                //console.log("max Ka/Ks score here is" + max);
                for(i=0;i<len;++i){
                    // normalizing
                    col.push(rgb(34,256*(data.KaKs_Scores[i]/max),76))
                }
            }
            if(b === "show similarity_scores" && data.similarity_scores){
                for(i=0;i<len;++i){
                    col.push(rgb(34,76,256*data.similarity_scores[i]))
                }
            }
            data.color = col;
            """)
            radio_button_group.on_click(callback)
            show(radio_button_group)

            hover = HoverTool(tooltips=[("label", "@label")])
            p = figure(title="T-sne of homologous gene X for each species",
                       tools=[hover, BoxZoomTool(), WheelZoomTool(), SaveTool(), ResetTool()])
            p.circle(x='x', y='y', source=source, color='color', size=10)
            if self.centroids is not None:
                p.square(self.centroids[0], self.centroids[0], size=12, color="olive", alpha=0.3)
            show(p)

            print "------------------------------------"
            for key, val in self.metrics.iteritems():
                print key + ': ' + str(val)
            return p
        else:
            fig = plt.figure(figsize=(40, size))
            if D == 2:
                ax = fig.add_subplot(111)
                ax.scatter(self.reduced[:, 0], self.reduced[:, 1], c=colors)
            if D == 3:
                ax = Axes3D(fig)
                ax.scatter(self.reduced[:, 0], self.reduced[:, 1], self.reduced[:, 2], c=colors)
            if D == 4:
                # TODO: totest
                ax = fig.add_subplot(111, projection='3d')
                ax.scatter(self.reduced[:, 0], self.reduced[:, 1], self.reduced[:, 2], s=self.reduced[:, 3] * 60, c=colors)
            else:
                print "please choose a D between 2 and 4"
                return
            plt.show()
            print "------------------------------------"
            for key, val in self.metrics.iteritems():
                print key + ': ' + str(val)

    def order(self, withtaxons=False):
        """
        order the names by numerical increasing order
        and sorts every representations as well according to this ordering
        """
        names = self.names[0] if withtaxons else self.names
        indices = sorted(range(len(names)), key=lambda k: names[k])
        names.sort()
        if self.full is not None:
            self.full[:] = self.full[indices]
        if self.reduced is not None:
            self.reduced[:] = self.reduced[indices]
        if self.clusters is not None:
            self.clusters = [self.clusters[i] for i in indices]
        if self.lenmat is not None:
            self.lenmat[:] = self.lenmat[indices]
        if withtaxons:
            self.names[0] = names
            for j, i in enumerate(indices):
                self.names[1][j] = self.names[1][i]
        else:
            self.names = names

    def clusterize_(self, clustering='gaussian', eps=0.8, homogroupnb=None, assess=True):
        """
        will clusterize the homology using gaussian mixture clustering or DBSCAN and order them according
        to the density of each cluster (we are interested in the dense ones) and assess the quality using 3 criterion:
        BIC, AIC ,silhouette, cal_hara .
        """
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
            self.clusters = alg.predict(self.full).tolist()

        elif clustering == 'dbscan' and homogroupnb is not None:
            # http://scikit-learn.org/stable/modules/generated/sklearn.cluster.DBSCAN.html
            val = np.around(len(self.names) / (2 * (homogroupnb + 1)))
            print val
            alg = cluster.DBSCAN(eps=eps, min_samples=7,
                                 algorithm='auto', n_jobs=-1)
            self.clusters = alg.fit_predict(self.full).tolist()
            n_clusters_ = len(set(self.clusters)) - (1 if -1 in self.clusters else 0)
            # TODO: use gaussian clustering and look at the variance of the kernels
            # (requested by dominique to maybe have some ideas of variance
            # as it is not well displayed by eps) add this as another information
            # when plotting and make gaussian clustering work well, then work using the values found by eps.

            print "Estimated number of clusters using DBscan: " + str(n_clusters_)
        elif clustering == 'jkmeans':
            # TODO:  create your own clustering algorithm with gaussian kernels
            # that explains as much as possible the data and get a threshold and set all others as outliers
            """
                here we want to find the smallest group of gaussian that explains
                the highest number of data point with the smallest variance possible
                as
            """
            print "tocode"
        if assess:
            # TODO: add F1 score if possible

            try:
                silhouette = metrics.silhouette_score(self.full, self.clusters).item()
                cal_hara = metrics.calinski_harabaz_score(self.full, self.clusters).item()
            except ValueError:
                silhouette = 0
                cal_hara = 0
            self.metrics = {}
            print 'silhouette_score ' + str(silhouette)
            print 'cal_hara ' + str(cal_hara)
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
                "KaKs_Scores": self.KaKs_Scores.tolist() if self.KaKs_Scores is not None else None,
                "similarity_scores": self.similarity_scores.tolist() if self.similarity_scores is not None else None,
                "proteinids": self.proteinids,
                "nans": self.nans.tolist() if self.nans is not None else None,
                "lenmat": self.lenmat.tolist() if self.lenmat is not None else None,
                "doub": self.doub.tolist() if self.doub is not None else None}
