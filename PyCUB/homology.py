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
from bokeh.io import save, show
from bokeh.models import *
from bokeh.layouts import column
import matplotlib.pyplot as plt
from scipy.spatial.distance import euclidean

import pdb


class homology(object):
    """in homology we store an homology with all its related data,

    it reduced matrix with dim reduction and its clusters for example
    it is supposed to be store in a dictionary of homologies
    an homology is a set of genes from different species related by a common ancester
    gene and generally a common function.
    the unique metadatas are generally from the reference species/genome

    Args:
        names: list of int corresponding to names in utils.speciestable
        full: np.array[float] (species, amino) of one homology
            with entropy value vector per species
        reduced:  np.array[float] (species, x*y) of one homology dimensionality reduced
        clusters: list[int] of cluster val for each species
        centroids: np.array[float] the position in aminoacid# Dimension of each centroids
            (number of centroid == number of cluster), dimension is reduced when plotted
        metrics: a dict of metrics names and values for the cluster of this homology
        nans: np.array[bool] wether or not this position is a nan
        lenmat: a np.array[int] species amino of int of length of each amino acids (number)
        doub: np.array[bool] of wether or not this position is a doublon
            ( a copy number of the gene for the species)
        reduced_algo: str dimensionality reduction algorithm used on this homology
        var: np.array of the variance in the CUB value for each datapoint
        mean: np.array[float] of the mean CUB for each datapoint
        homocode: str the code of the homology
        GCcount: np.array[float] GC bias for each datapoint
        KaKs_Scores: np.array[float] a form of similarity score between genes/datapoints
        similarity_scores: np.array[float] similarity score between genes/datapoints
        proteinids: list[str] the protein names for each datapoints
        isrecent: float a proxy for wether or not this homology has appeated recently
        ishighpreserved: bool a proxy for wether or not this homology is conserved throughout evolution
        geneids: list[str] the name of the genes
        ref: np.array[float] the reference CUB value of cerevisiae gene
        refprot: str the reference protein name of cerevisiae gene
        refgene: str the reference gene name of cerevisiae gene
        ecai: np.array[float] the codon adaptation index of each gene
        meanecai: float the mean ecai of the ecai
        protein_abundance: float the average abundance of the protein enoded by this gene in cerevisiae cells
        weight: int the molecular weight of this protein
        mRNA_abundance: float the average abundance of the messenger RNA of this coding gene in cerevisiae cells
        cys_elements: int the number of cys regulatory elements known for cerevisiae on this gene
            (the value is only referenced for secreted genes)
        is_secreted: bool wether or not this protein is secreted out of the cell
        decay_rate: float the half life in minute of this protein
        tot_volume: int a proxy of the molecular volume of this protein
        mean_hydrophobicity: a very distant proxy of the hydrophobicity of this protein
        glucose_cost: the glucose cost of creating the amino acids required to build up this protein
        synthesis_steps: the number of steps required by the cell to build up the amino acids of this protein
        isoelectricpoint: float, a proxy of the Pi of this protein.
        conservation: float, the total conservation of each amino of the corresponding protein
        othercods: float, the average number of codons other than the ones of the 18 amino acids we are
            looking at per species on the homology

    """

    full = None  # np array len(species)xCUBD
    var = None  # np array len(species)
    mean = None  # np array float

    clusters = None  # list[int]
    centroids = None  # np array Xx2/3/4
    metrics = {}
    homocode = 'None'
    nans = None  # array len(species)
    lenmat = None  # array len(species)x 18
    names = None  # list[int]
    doub = None  # array len(species)
    GCcount = None  # array len(species)

    reduced = None  # array len(species)x2/3/4
    reduced_algo = None

    KaKs_Scores = None  # array len(species)
    similarity_scores = None  # array len(species)
    proteinids = []  # list[str]
    isrecent = None  # float
    ishighpreserved = None  # bool
    geneids = None  # list[str]
    ref = None  # array CUBD
    refprot = None  # str
    refgene = None  # str
    ecai = None  # array len(species)
    meanecai = None  # float
    cai = None  # array len(species)
    meancai = None  # float

    protein_abundance = 0.  # float
    weight = 0  # int
    conservation = None  # float
    mRNA_abundance = 0.  # float
    cys_elements = 0  # int
    is_secreted = False  # bool
    decay_rate = 0.  # float
    tot_volume = None  # float
    mean_hydrophobicity = None  # float
    glucose_cost = None  # float
    synthesis_steps = None  # float
    isoelectricpoint = None  # float
    othercods = None  # float

    def __init__(self, **kwargs):
        """
        can intialize the file from kwargs as a raw dictionnary for json format (output of dictify) or
        from regular args.
        """
        data = kwargs.get("data", None)
        if data is not None:
            self.full = np.asarray(data.get("full")) if data.get("full", None) is not None else None
            self.var = np.asarray(data.get("var")) if data.get("var", None) is not None else None
            self.mean = np.asarray(data.get("mean")) if data.get("mean", None) is not None else None
            self.clusters = data.get("clusters", None)
            self.centroids = np.asarray(data.get("centroids")) if data.get("centroids", None) is not None else None
            self.metrics = data.get("metrics", {})
            self.nans = np.asarray(data.get("nans")) if data.get("nans", None) is not None else None
            self.lenmat = np.asarray(data.get("lenmat")) if data.get("lenmat", None) is not None else None
            self.names = data.get("names", None)
            self.doub = np.asarray(data.get("doub")) if data.get("doub", None) is not None else None
            self.GCcount = np.asarray(data.get("GCcount")) if data.get("GCcount", None) is not None else None
            self.reduced = np.asarray(data.get("reduced")) if data.get("reduced", None) is not None else None
            self.reduced_algo = data.get("reduced_algo", None)
            self.KaKs_Scores = np.asarray(data.get("KaKs_Scores", None)) if data.get("KaKs_Scores", None) is not None else None
            self.homocode = str(data.get("homocode", 'None'))
            self.similarity_scores = np.asarray(data.get("similarity_scores", None))\
                if data.get("similarity_scores", None) is not None else None
            self.proteinids = [str(i) for i in data.get("proteinids", [])]
            self.isrecent = data.get("isrecent", None)
            self.ishighpreserved = data.get("ishighpreserved", None)
            self.geneids = [str(i) for i in data.get("geneids", None)]
            self.conservation = data.get("conservation", None)
            self.cai = np.asarray(data["cai"]) if data.get("cai") else None
            self.meancai = data.get("meancai", None)
            self.ecai = np.asarray(data["ecai"]) if data.get("ecai", False) else None
            self.meanecai = data.get("meanecai", None)
            self.protein_abundance = data.get("protein_abundance", 0.)
            self.weight = data.get("weight", 0)
            self.mRNA_abundance = data.get("mRNA_abundance", 0.)
            self.cys_elements = data.get("cys_elements", 0)
            self.is_secreted = data.get("is_secreted", False)
            self.decay_rate = data.get("decay_rate", 0.)
            self.tot_volume = data.get("tot_volume", None)
            self.mean_hydrophobicity = data.get("mean_hydrophobicity", None)
            self.glucose_cost = data.get("glucose_cost", None)
            self.synthesis_steps = data.get("synthesis_steps", None)
            self.isoelectricpoint = data.get("isoelectricpoint", None)
            self.ref = np.asarray(data.get("ref", None)) if data.get("ref", None) is not None else None
            self.refprot = str(data.get("refprot", None))
            self.refgene = str(data.get("refgene", None))
            self.othercods = data.get("othercods", None)
        else:
            self.full = kwargs.get("full", None)
            self.var = kwargs.get("var", None)
            self.mean = kwargs.get("mean", None)
            self.clusters = kwargs.get("clusters", None)
            self.centroids = kwargs.get("centroids", None)
            self.metrics = kwargs.get("metrics", {})
            self.nans = kwargs.get("nans", None)
            self.lenmat = kwargs.get("lenmat", None)
            self.names = kwargs.get("names", None)
            self.doub = kwargs.get("doub", None)
            self.homocode = kwargs.get("homocode", 'None')
            self.GCcount = kwargs.get("GCcount", None)
            self.reduced = kwargs.get("reduced", None)
            self.reduced_algo = kwargs.get("reduced_algo", None)
            self.KaKs_Scores = kwargs.get("KaKs_Scores", None)
            self.similarity_scores = kwargs.get("similarity_scores", None)
            self.proteinids = kwargs.get("proteinids", [])
            self.isrecent = kwargs.get("isrecent", None)
            self.ishighpreserved = kwargs.get("ishighpreserved", None)
            self.geneids = kwargs.get("geneids", None)
            self.ecai = kwargs.get("ecai", None)
            self.meanecai = kwargs.get("meanecai", None)
            self.cai = kwargs.get("cai", None)
            self.meancai = kwargs.get("meancai", None)
            self.mRNA_abundance = kwargs.get("mRNA_abundance", 0.)
            self.protein_abundance = kwargs.get("mRNA_abundance", 0.)
            self.weight = kwargs.get("weight", 0)
            self.cys_elements = kwargs.get("cys_elements", 0)
            self.is_secreted = kwargs.get("is_secreted", False)
            self.decay_rate = kwargs.get("decay_rate", 0.)
            self.tot_volume = kwargs.get("tot_volume", None)
            self.mean_hydrophobicity = kwargs.get("mean_hydrophobicity", None)
            self.glucose_cost = kwargs.get("glucose_cost", None)
            self.synthesis_steps = kwargs.get("synthesis_steps", None)
            self.isoelectricpoint = kwargs.get("isoelectricpoint", None)
            self.ref = kwargs.get("ref", None)
            self.refprot = kwargs.get("refprot", None)
            self.refgene = kwargs.get("refgene", None)
            self.othercods = kwargs.get("othercods", None)
            self.conservation = kwargs.get("conservation", None)

    def __str__(self):
        speciestable = dict(utils.speciestable)
        returnd = '\nhomology: ' + self.homocode +\
            '\nsize: ' + str(len(self.full)) +\
            '\nfullmax: ' + str(self.full.max()) +\
            '\nvar: ' + str(self.var) +\
            '\nmean: ' + str(self.mean) +\
            '\nmetrics: ' + str(self.metrics) +\
            '\nnans mean: ' + str(self.nans.mean()) +\
            '\nlengths mean/var: ' + str(self.lenmat.sum(1).mean()) + ', ' + str(self.lenmat.sum(1).var()**(0.5)) +\
            '\ndoublon sum: ' + str(self.doub.sum()) +\
            '\nGCcount mean/var: ' + str(self.GCcount.mean()) + ', ' + str(self.GCcount.var()**(0.5)) +\
            '\nKaKs_Scores: ' + str(self.KaKs_Scores) +\
            '\nisrecent: ' + str(self.isrecent) +\
            '\nishighpreserved: ' + str(self.ishighpreserved) +\
            '\nref: ' + str(self.ref) +\
            '\nrefprot: ' + str(self.refprot) +\
            '\nrefgene: ' + str(self.refgene) +\
            '\necai mean/var: ' + str(self.ecai.mean()) + ', ' + str(self.ecai.var()**(0.5)) +\
            '\ncai mean/var: ' + str(self.cai.mean()) + ', ' + str(self.cai.var()**(0.5)) +\
            '\nprotein_abundance: ' + str(self.protein_abundance) +\
            '\nweight: ' + str(self.weight) +\
            '\nconservation: ' + str(self.conservation) +\
            '\nmRNA_abundance: ' + str(self.mRNA_abundance) +\
            '\ncys_elements: ' + str(self.cys_elements) +\
            '\nis_secreted: ' + str(self.is_secreted) +\
            '\ndecay_rate: ' + str(self.decay_rate) +\
            '\ntot_volume: ' + str(self.tot_volume) +\
            '\nmean_hydrophobicity: ' + str(self.mean_hydrophobicity) +\
            '\nglucose_cost: ' + str(self.glucose_cost) +\
            '\nsynthesis_steps: ' + str(self.synthesis_steps) +\
            '\nisoelectricpoint: ' + str(self.isoelectricpoint) +\
            '\nothercods: ' + str(self.othercods)
        return returnd

# CpCpGpApApTpApTpApTpTpCpCpGpApApTpApTpApTpTpCpCpGpApApTpApTpApTpTpCpCpGpApApTpApTpApTpTpTpTpCpCpGpApApTpApTpApTpTp
# GbGbCbTbTbAbTbAbTbAbAbGbGbCbTbTbAbTbAbTbAbAbGbGbCbTbTbAbTbAbTbAbAbGbGbCbTbTbAbTbAbTbAbAbAbGbGbCbTbTbAbTbAbTbAbAbAb

    def remove(self, species):
        """
        removes the list of species from this homology if it exists there

        Args:
            species: list[str] of species to remove

        """
        speciestable = dict(utils.speciestable)
        tradnames = [speciestable[na] for na in self.names] if type(species[0]) == unicode or \
            type(species[0]) == unicode else self.names
        names = list(self.names)
        cluster = []
        mask = np.zeros(len(self.names), dtype=bool)
        self.names = []
        for i, na in enumerate(tradnames):
            if na not in species:
                self.names.append(names[i])
                mask[i] = True
                if self.clusters is not None:
                    cluster.append(self.clusters[i])
        self.cluster = np.array(cluster) if self.clusters is not None else None
        self.reduced = self.reduced[mask, :] if self.reduced is not None else None
        self.full = self.full[mask, :] if self.full is not None else None

    def nb_unique_species(self):
        """
        compute the number of unique species in this homologies

        (basically count the number of doub)

        Args:
            None

        Returns:
            The number of unique species in this homology
        """
        return self.doub.shape[0] - np.count_nonzero(self.doub)

    def order(self, withtaxons=False):
        """
        order the names by numerical increasing order

        and sorts every representations as well according to this ordering

        Args:
            withtaxons: bool to true if there is taxonomic data (present before preprocessing)
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
        if self.mean is not None:
            self.mean[:] = self.mean[indices]
        if self.nans is not None:
            self.nans[:] = self.nans[indices]
        if self.doub is not None:
            self.doub[:] = self.doub[indices]
        if self.var is not None:
            self.var[:] = self.var[indices]
        if self.GCcount is not None:
            self.GCcount[:] = self.GCcount[indices]
        if self.KaKs_Scores is not None:
            if len(self.KaKs_Scores) < len(self.full):
                self.KaKs_Scores = None
            else:
                self.KaKs_Scores[:] = self.KaKs_Scores[indices]
        if self.similarity_scores is not None:
            if len(self.similarity_scores) < len(self.full):
                self.similarity_scores = None
            else:
                self.similarity_scores[:] = self.similarity_scores[indices]
        if self.proteinids is not None:
            self.proteinids = [self.proteinids[i] for i in indices]
        if self.geneids is not None:
            self.geneids = [self.geneids[i] for i in indices]
        if self.ecai is not None:
            self.ecai[:] = self.ecai[indices]
        if self.cai is not None:
            self.cai[:] = self.cai[indices]
        if withtaxons:
            self.names[0] = names
            for j, i in enumerate(indices):
                self.names[1][j] = self.names[1][i]
        else:
            self.names = names

    def compute_averages(self):
        """
        Computes the mean, var and mean of the homology

        Args:
            None
        """
        self.mean = self.full.mean(0)
        self.var = self.full.var(0)**(0.5)
        if self.ecai is not None:
            self.meanecai = self.ecai.mean()
        if self.cai is not None:
            self.meancai = self.cai.mean()

    def reduce_dim(self, alg='tsne', n=2, perplexity=40):
        """
        reduce the dimensionality of your gene dataset to a defined dimension using the t-SNE algorithm

        Args:
            alg: a matrix of gene codon usage per species
            n: the desired dimension
            perplexity: an optional value when you know about tsne

        Returns:
            The reduced dataset
        """
        full = self.full if self.centroids is None else np.vstack((self.full, self.centroids))

        if alg == 'tsne':
            red = man.TSNE(n_components=n, perplexity=perplexity).fit_transform(full)
        elif alg == 'pca':
            red = PCA(n_components=n).fit_transform(full)
        else:
            raise AttributeError("wrong algo my friend choose pca or tsne")
        if self.centroids is None:
            self.reduced = red
        else:
            self.reduced = red[:len(self.full)]
            self.centroids = red[len(self.full):]
        self.reduced_algo = alg
        return red

    def plot(self, reducer="tsne", per=40, interactive=True, D=2, size=20):
        """
        plot the dimensionality reduced datapoint of the homology matrix

        colors represents the different clusters you can set the interactivity to gain
        deeper knowledge of the dataset. It can dim reduce the data and will also save
        the figure in utils/save. Moreover, it will print some interesting data about the
        homology

        Args:
            reducer: str flag the algorithm to reduce the matrix of utils.cubD to D D
            per: int he perplexity when using tsne algorithm
            size: int the size of the plot
            interactive: (recommended) wether to use bokeh or not to represend the data
                        will show name of the species when hovering over datapoint
                        will show a gradient of evolutionary distance of the species of
                        the datapoint currently hovered to each other species
            D: int the goal dimension (2-3-4) when using matplotlib

        Returns:
            The plot interactive or not and some informations (as prints)

        Raises:
            UnboundLocalError: "you need to have dimensionality reduced to 2D to have a right plotting of the centroids"
        """
        colormap = list(utils.colormap)
        spetable = dict(utils.speciestable)
        if self.clusters is not None:
            # colormap = [[rand(256), rand(256), rand(256)] for _ in range(100)]
            colors = [colormap[x + 1] for x in self.clusters]
        else:
            colors = [colormap[0]] * len(self.names)
        if self.reduced is None or self.reduced_algo != reducer:
            self.reduce_dim(alg=reducer, n=D, perplexity=per)
        elif self.reduced.shape[1] != D:
            self.reduce_dim(alg=reducer, n=D, perplexity=per)
        if interactive:
            print " if you are on a notebook you should write 'from bokeh.io import output_notebook'"
            data = dict(x=self.reduced[:, 0], y=self.reduced[:, 1],
                        species=[str(spetable[n]) for n in self.names],
                        meanentropy=["%.2f" % self.full[i].mean() for i in range(len(self.names))],
                        color=colors)
            # TODO: show the ref gene
            if self.doub is not None:
                data.update({'doub': self.doub})
            if self.clusters is not None:
                data.update({'clusters': self.clusters})
             # TODO addphylodists
            if self.similarity_scores is not None:
                data.update({'similarity_scores': self.similarity_scores})
            if self.KaKs_Scores is not None:
                data.update({'KaKs_Scores': self.KaKs_Scores})
            if self.nans is not None:
                data.update({'nans': self.nans})
                # here use wether or not it is a new protein
            if self.ecai is not None:
                data.update({'ecai': self.ecai})
            if self.cai is not None:
                data.update({'cai': self.cai})
            if self.lenmat is not None:
                data.update({'lengths': self.lenmat.sum(1)})
            if self.GCcount is not None:
                data.update({'gc': self.GCcount})
            source = ColumnDataSource(data=data)
            output_notebook()
            labe = ["show Cluster", "show Doublon", "show Nans", "show KaKs_Scores",
                    "show similarity_scores", "show Length", "show gc", "show ecai", "show cai"]  # 9
            callback = CustomJS(args=dict(source=source), code=str(utils.callback))
            radio_button_group = widgets.RadioButtonGroup(
                labels=labe, callback=callback, active=0)
            hover = HoverTool(tooltips=[("species", "@species"), ("nan_nb", "@nans"),
                                        ("mean_entr", "@meanentropy"), ("length", "@lengths"), ("GCcount", "@gc")])
            p = figure(title="Plot of the homology " + self.homocode,
                       tools=[hover, WheelZoomTool(), PanTool(), SaveTool(), ResetTool()],
                       plot_width=800, plot_height=600)
            p.circle(x='x', y='y', source=source, color='color', size=10)
            if self.centroids is not None:
                if self.centroids.shape[1] > 2:
                    raise UnboundLocalError("you need to have dimensionality reduced to 2D \
                        to have a right plotting of the centroids")
                p.square(self.centroids[:, 0], self.centroids[:, 1], size=12, color="olive", alpha=0.3)
            save(column(radio_button_group, p), "utils/templot/homoplot_interactive_.html")
            show(column(radio_button_group, p))
        else:
            fig = plt.figure(figsize=(40, size))
            if D == 2:
                ax = fig.add_subplot(111)
                ax.scatter(self.reduced[:, 0], self.reduced[:, 1], s=90, c=colors)
            elif D == 3:
                ax = Axes3D(fig)
                ax.scatter(self.reduced[:, 0], self.reduced[:, 1], self.reduced[:, 2], s=140, c=colors)
            elif D == 4:
                ax = fig.add_subplot(111, projection='3d')
                ax.scatter(self.reduced[:, 0], self.reduced[:, 1],
                           self.reduced[:, 2], s=self.reduced[:, 3] * 200, c=colors)
            else:
                raise AttributeError("please choose a D between 2 and 4")
            plt.show()
            plt.savefig("utils/templot/homoplot_matplotlib_.pdf")
        print "homology: " + self.homocode
        print "------------------------------------"
        for key, val in self.metrics.iteritems():
            print key + ': ' + str(val)
        print "------------------------------------"
        mea = self.full.mean(axis=0)
        var = self.full.var(axis=0)
        print "Amino, CUB value mean, variance"
        amino = list(utils.amino)
        for i in range(len(amino)):
            print amino[i] + ": %.2f" % mea[i] + ", %.3f" % var[i]
        if self.KaKs_Scores is not None:
            print "avg KaKs Scores"
            print self.KaKs_Scores.mean()
        if self.similarity_scores is not None:
            print "avg similarity Scores"
            print self.similarity_scores.mean()
        print "------------------------------------"
        if self.meanecai is not None:
            print "mean ecai: " + str(self.meanecai)
        if self.meancai is not None:
            print "mean cai: " + str(self.meancai)
        if self.isrecent is not None:
            if self.isrecent:
                print "fairly recent protein"
        if self.ishighpreserved is not None:
            if self.ishighpreserved:
                print "high preserved protein"
        if self.protein_abundance is not None:
            print "abundance of the protein in a cell: " + str(self.protein_abundance)
        if self.weight is not None:
            print "weight of protein: " + str(self.weight)
        if self.mRNA_abundance is not None:
            print "abundance of mRNA: " + str(self.mRNA_abundance)
        if self.cys_elements:
            print "number of cys_elements: " + str(self.cys_elements)
        if self.is_secreted:
            print "the protein is secreted out of the cell"
        if self.decay_rate is not None:
            print "halflife (in mn): " + str(self.decay_rate)
        if self.tot_volume is not None:
            print "total volume of protein: " + str(self.tot_volume)
        if self.mean_hydrophobicity is not None:
            print "the pseudo hydrophobicity is of: " + str(self.mean_hydrophobicity)
        if self.glucose_cost is not None:
            print "the glucose cost: " + str(self.glucose_cost)
        if self.synthesis_steps is not None:
            print "synthesis steps: " + str(self.synthesis_steps)
        if self.isoelectricpoint is not None:
            print "the isoelectricpoint(Pi) is : " + str(self.isoelectricpoint)

    def clusterize_(self, clustering='gaussian', eps=0.8, homogroupnb=None, assess=True, verbose=True):
        """
        will clusterize the homology using gaussian mixture clustering or DBSCAN

        and order them according
        to the density of each cluster (we are interested in the dense ones)
        and assess the quality using 3 criterion:
        BIC, AIC ,silhouette, cal_hara, phylodistance.

        Args:
            clustering: str flag the clustering algorithm [gaussian, dbscan]
            eps: float, hyperparam of the max size of the nsphere of each cluster
            homogroupnb: int hyperparam of gaussian for the number of clusters
            assess: wether or not to assess the quality of the clustering
            verbose: wether or not to show clustering quality information

        Returns:
            The clusters for each datapoint of the homology as a list[int]
        Raises:
            AttributeError: "Hey, please use gaussian or dbscan"
        """
        self.metrics = {}
        if clustering == 'gaussian' and homogroupnb is not None:
            # http://scikit-learn.org/stable/modules/generated/sklearn.mixture.GaussianMixture.html
            alg = mixture.GaussianMixture(n_components=homogroupnb, n_init=2, init_params='random')
            alg.fit(self.full)
            if assess:
                aic = alg.aic(self.full)
                bic = alg.bic(self.full)
                self.metrics.update({'aic': aic, 'bic': bic})
                if verbose:
                    print "the BIC scores for the GMM is"
                    print aic
                    print "the AIC scores for the GMM is"
                    print bic
                self.metrics.update({'bic': bic, 'aic': aic})
                self.metrics.update
            self.clusters = alg.predict(self.full).tolist()
            cov = alg.covariances_
            self.centroids = alg.means_
            dists = np.zeros((len(self.centroids), 2))
            for i in range(len(self.centroids)):
                dist = np.zeros(len(cov[i]))
                for j in range(len(cov[i])):
                    dist[j] = euclidean(self.centroids[i], cov[i][j])
                dists[i] = [dist.mean(), dist.var()**(0.5)]
            print "the mean,var of the distances of the covariances to the means of each clusters are: "
            print dists
            n_clusters_ = homogroupnb
        elif clustering == 'dbscan':
            # http://scikit-learn.org/stable/modules/generated/sklearn.cluster.DBSCAN.html

            alg = cluster.DBSCAN(eps=eps, min_samples=7,
                                 algorithm='auto', n_jobs=1)
            self.clusters = alg.fit_predict(self.full).tolist()
            n_clusters_ = len(set(self.clusters)) - (1 if -1 in self.clusters else 0)
            if verbose:
                print "Estimated number of clusters using DBscan: " + str(n_clusters_)
        else:
            raise AttributeError("Hey, please use gaussian or dbscan")
        if assess:
            try:
                silhouette = metrics.silhouette_score(self.full, self.clusters).item()
                cal_hara = metrics.calinski_harabaz_score(self.full, self.clusters).item()
            except ValueError:
                silhouette = 0
                cal_hara = 0

            if verbose:
                print 'silhouette_score ' + str(silhouette)
                print 'cal_hara ' + str(cal_hara)
            self.metrics.update({'silhouette': silhouette,
                                 'cal_hara': cal_hara})
            if utils.phylo_distances is not None:
                speciestable = dict(utils.speciestable)
                avg_phylodistance = []
                hastaxons = utils.phylo_distances.index.tolist()
                spe = list(set([speciestable[j] for j in self.names if speciestable[j] in hastaxons]))
                if len(spe) < 2:
                    self.metrics.update({'cluster_phylodistance': [1.]})
                    if verbose:
                        print 'not enough species'
                    return self.clusters
                div = float(utils.phylo_distances[spe].loc[spe].sum().sum()) / (len(spe)**2 - len(spe))
                for i in range(-1, n_clusters_):
                    # Here the first value is for the outliers and the second for the unclusterized
                    # data points
                    ind = np.argwhere(np.array(self.clusters) == i).T[0]
                    species = list(set([speciestable[self.names[j]] for j in ind if speciestable[self.names[j]] in hastaxons]))
                    if len(species) < 2:
                        continue
                    avg_phylodistance.append((float(utils.phylo_distances[species].
                                                    loc[species].sum().sum()) / (len(species)**2 - len(species))) / div)
                self.metrics.update({'cluster_phylodistance': avg_phylodistance})
                if verbose:
                    print 'avg phylo dist of the clusters: ' + str(np.array(self.metrics['cluster_phylodistance']))

        return self.clusters

    # TODO: add a classifier to high or low CUB given species info

    def _dictify(self):
        """
        Used by the saving function. transform the object into a dictionary that can be json serializable

        Args:
            None

        Returns:
            A dict holding every element to be jsonized
        """
        return {"reduced": self.reduced.tolist() if not (self.reduced is None) else None,
                "clusters": self.clusters,
                "full": self.full.tolist() if not (self.full is None) else None,
                "names": self.names,
                "homocode": self.homocode,
                "centroids": self.centroids.tolist() if self.centroids is not None else None,
                "metrics": self.metrics,
                "KaKs_Scores": self.KaKs_Scores.tolist() if self.KaKs_Scores is not None else None,
                "similarity_scores": self.similarity_scores.tolist() if self.similarity_scores is not None else None,
                "proteinids": self.proteinids,
                "geneids": self.geneids,
                "ecai": self.ecai.tolist() if self.ecai is not None else None,
                "meanecai": self.meanecai,
                "cai": self.cai.tolist() if self.cai is not None else None,
                "meancai": self.meancai,
                "nans": self.nans.tolist() if self.nans is not None else None,
                "doub": self.doub.tolist() if self.doub is not None else None,
                "var": self.var.tolist() if self.var is not None else None,
                "mean": self.mean.tolist() if self.mean is not None else None,
                "lenmat": self.lenmat.tolist() if self.lenmat is not None else None,
                "GCcount": self.GCcount.tolist() if self.GCcount is not None else None,
                "reduced_algo": self.reduced_algo,
                "isrecent": self.isrecent,
                "ishighpreserved": self.ishighpreserved,
                "ref": self.ref.tolist() if self.ref is not None else None,
                "refprot": self.refprot,
                "refgene": self.refgene,
                "protein_abundance": self.protein_abundance,
                "weight": self.weight,
                "mRNA_abundance": self.mRNA_abundance,
                "cys_elements": self.cys_elements,
                "is_secreted": self.is_secreted,
                "decay_rate": self.decay_rate,
                "tot_volume": self.tot_volume,
                "othercods": self.othercods,
                "conservation": self.conservation,
                "mean_hydrophobicity": self.mean_hydrophobicity,
                "glucose_cost": self.glucose_cost,
                "synthesis_steps": self.synthesis_steps,
                "isoelectricpoint": self.isoelectricpoint}
