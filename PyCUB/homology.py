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
    full : np array of float (species, amino) of one homology
    with entropy value vector per species
    reduced :  np array float (species, x*y) of one homology dimensionality reduced
    clusters : list[int] of cluster val for each species
    centroids: the position in aminoacid# Dimension of each centroids
    (number of centroid == number of cluster)
    metrics: a dict of metrics names and values for the cluster of this homology
    nans: numpy array of bool wether or not this position is a nan
    lenmat: a numpy array species amino of int of length of each amino acids (number)
    doub: nupy array of [bool] of wether or not this position is a doublon
    ( a copy number of the gene for the species)
    reduced_algo: dimensionality reduction algorithm used on this homology
    """

    full = None
    var = None
    mean = None

    clusters = None
    centroids = None
    metrics = {}

    nans = None
    nanspos = None
    lenmat = None
    names = None
    doub = None
    doublonpos = None
    GCcount = None

    reduced = None
    reduced_algo = None

    KaKs_Scores = None
    similarity_scores = None
    proteinids = []
    isrecent = None
    ishighpreserved = None

    # TODO: compute the mean variance in CUB value and mean range
    # for each homology (add this when plotting)
    # TODO: compute a distance matrix between the species

    def __init__(self, **kwargs):
        """
        will intialize an instance of the object and be used for the loading mechanism
        """
        # TODO: Replace long params by **kwargs (dict for the data)
        data = kwargs.get("data", None)
        if data is not None:
            self.full = np.asarray(data.get("full")) if data.get("full", None) is not None else None
            self.var = np.asarray(data.get("var")) if data.get("var", None) is not None else None
            self.mean = np.asarray(data.get("mean")) if data.get("mean", None) is not None else None
            self.clusters = np.asarray(data.get("clusters")) if data.get("clusters", None) is not None else None
            self.centroids = np.asarray(data.get("centroids")) if data.get("centroids", None) is not None else None
            self.metrics = data.get("metrics", {})
            self.nans = np.asarray(data.get("nans")) if data.get("nans", None) is not None else None
            self.nanspos = np.asarray(data.get("nanspos")) if data.get("nanspos", None) is not None else None
            self.lenmat = np.asarray(data.get("lenmat")) if data.get("lenmat", None) is not None else None
            self.names = data.get("names", None)
            self.doub = np.asarray(data.get("doub")) if data.get("doub", None) is not None else None
            self.doublonpos = np.asarray(data.get("doublonpos")) if data.get("doublonpos", None) is not None else None
            self.GCcount = np.asarray(data.get("GCcount")) if data.get("GCcount", None) is not None else None
            self.reduced = np.asarray(data.get("reduced")) if data.get("reduced", None) is not None else None
            self.reduced_algo = data.get("reduced_algo", None)
            self.KaKs_Scores = data.get("KaKs_Scores", None)
            self.similarity_scores = data.get("similarity_scores", None)
            self.proteinids = data.get("proteinids", [])
            self.isrecent = data.get("isrecent", None)
            self.ishighpreserved = data.get("ishighpreserved", None)
        else:
            self.full = kwargs.get("full", None)
            self.var = kwargs.get("var", None)
            self.mean = kwargs.get("mean", None)
            self.clusters = kwargs.get("clusters", None)
            self.centroids = kwargs.get("centroids", None)
            self.metrics = kwargs.get("metrics", {})
            self.nans = kwargs.get("nans", None)
            self.nanspos = kwargs.get("nanspos", None)
            self.lenmat = kwargs.get("lenmat", None)
            self.names = kwargs.get("names", None)
            self.doub = kwargs.get("doub", None)
            self.doublonpos = kwargs.get("doublonpos", None)
            self.GCcount = kwargs.get("GCcount", None)
            self.reduced = kwargs.get("reduced", None)
            self.reduced_algo = kwargs.get("reduced_algo", None)
            self.KaKs_Scores = kwargs.get("KaKs_Scores", None)
            self.similarity_scores = kwargs.get("similarity_scores", None)
            self.proteinids = kwargs.get("proteinids", [])
            self.isrecent = kwargs.get("isrecent", None)
            self.ishighpreserved = kwargs.get("ishighpreserved", None)

    def remove(self, species):
        """
        removes the list of species from this homology if it exists there
        """
        names = [utils.speciestable[na] for na in self.names] if type(species[0]) == str else self.names
        mask = np.ones(len(self.names), dtype=bool)
        pdb.set_trace()
        for i, na in enumerate(names):
            if na in species:
                self.names.pop(i)
                mask[i] = False
                if self.clusters is not None:
                    self.clusters.pop(i)
        self.reduced = self.reduced[mask, :] if self.reduced is not None else None
        self.full = self.full[mask, :] if self.full is not None else None

    def nb_unique_species(self):
        """
        compute the number of unique species in this homologies
        (basically count the number of doub)
        """
        return self.doub.shape[0] - np.count_nonzero(self.doub)

    def getdoubpos(self):
        pos = {}
        if self.doub is not None:
            for i, val in enumerate(self.doub):
                if val:
                    mean = full[i - 1:i + 1].mean(0)
                    pos.update({utils.speciestable[self.names[i]]:
                                np.divide(np.count_nonzero(mean > self.full, 0), self.nb_unique_species())})
        return pos

    def getnanpos(self):
        pos = {}
        if self.nans is not None:
            for i, val in enumerate(self.nans):
                if val:
                    pos.update({utils.speciestable[self.names[i]]:
                                np.divide(np.count_nonzero(val > self.full, 0), self.nb_unique_species())})
        return pos

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

    def compute_halflife(self):
        """
        retrieve the half life of the mRNAs according to a db

        """
        # TODO to code
        pass

    def compute_averages(self):
        """
        """
        self.mean = self.full.mean(0)
        self.var = self.full.var(0)

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
            source = ColumnDataSource(
                data=dict(x=self.reduced[:, 0], y=self.reduced[:, 1],
                          label=["species : %s" % utils.speciestable[x__] for x__ in self.names],
                          doub=self.doub if self.doub is not None else False,
                          clusters=self.clusters if self.clusters is not None else False,
                          similarity_scores=self.similarity_scores if
                          self.similarity_scores is not None else False,
                          KaKs_Scores=self.KaKs_Scores if self.KaKs_Scores is not None else False,
                          nans=self.nans if self.nans is not None else False,
                          color=colors))
            output_notebook()

            radio_button_group = widgets.RadioButtonGroup(
                labels=["show Cluster", "show Doublon", "show KaKs_Scores",
                        "show similarity_scores", "Show PhyloDistance"], active=0)
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
            print "------------------------------------"
            print self.full.mean(axis=0)
            print self.full.var(axis=0)
            print self.KaKs_Scores.mean()
            print self.similarity_scores.mean()
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
                ax.scatter(self.reduced[:, 0], self.reduced[:, 1],
                           self.reduced[:, 2], s=self.reduced[:, 3] * 60, c=colors)
            else:
                print "please choose a D between 2 and 4"
                return
            plt.show()

            print "------------------------------------"
            for key, val in self.metrics.iteritems():
                print key + ': ' + str(val)
            print "------------------------------------"
            print self.full.mean(axis=0)
            print self.full.var(axis=0)
            print self.KaKs_Scores.mean()
            print self.similarity_scores.mean()

    def clusterize_(self, clustering='gaussian', eps=0.8, homogroupnb=None, assess=True):
        """
        will clusterize the homology using gaussian mixture clustering or DBSCAN and order
        them according
        to the density of each cluster (we are interested in the dense ones)
        and assess the quality using 3 criterion:
        BIC, AIC ,silhouette, cal_hara, phylodistance.
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
            # when plotting and make gaussian clustering work well,
            # then work using the values found by eps.

            print "Estimated number of clusters using DBscan: " + str(n_clusters_)
        if assess:
            if utils.phylo_distances is not None:
                avg_phylodistance = []
                spe = [utils.speciestable[j] for j in self.names]
                div = utils.phylo_distances[spe].loc[spe].sum().sum() / (len(spe)**2 - len(spe))
                for i in range(-1, n_clusters_ - 1):
                    # Here the first value is for the outliers and the second for the unclusterized
                    # data points
                    species = [utils.speciestable[j] for j in self.names[np.argwhere(self.clusters == i)]]
                    avg_phylodistance.append((utils.phylo_distances[species].
                                              loc[species].sum().sum() / (len(species)**2 - len(species))) / div)
                # print 'average phylo distance of the clusters: '+ str(avg_phylodistance)
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
                                 'cal_hara': cal_hara,
                                 'avg_phylodistance': np.array(avg_phylodistance)})
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
                "KaKs_Scores": self.KaKs_Scores,
                "similarity_scores": self.similarity_scores,
                "proteinids": self.proteinids,
                "nans": self.nans.tolist() if self.nans is not None else None,
                "doub": self.doub.tolist() if self.doub is not None else None,
                "var": self.var.tolist() if self.var is not None else None,
                "mean": self.mean.tolist() if self.mean is not None else None,
                "nanspos": self.nanspos.tolist() if self.nanspos is not None else None,
                "lenmat": self.lenmat.tolist() if self.lenmat is not None else None,
                "doublonpos": self.doublonpos.tolist() if self.doublonpos is not None else None,
                "GCcount": self.GCcount.tolist() if self.GCcount is not None else None,
                "reduced_algo": self.reduced_algo,
                "isrecent": self.isrecent,
                "ishighpreserved": self.ishighpreserved}
