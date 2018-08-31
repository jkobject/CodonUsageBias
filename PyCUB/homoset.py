""""
Created by Jeremie KALFON
Date : 21 FEV 2018
University of Kent, ECE paris
jkobject.com

"""


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import holoviews as hv
from holoviews.operation.datashader import datashade, dynspread, shade
from bokeh.io import save, show
from bokeh.plotting import *
from bokeh.models import *
from bokeh.layouts import column
import utils
import homology as h

from sklearn.cluster import KMeans
from sklearn.cluster import MiniBatchKMeans
from sklearn.metrics.pairwise import cosine_similarity
from kmodes.kmodes import KModes
from sklearn import manifold as man
from sklearn import metrics
from sklearn.decomposition import PCA
from sklearn.preprocessing import normalize
from sklearn.metrics.pairwise import paired_distances as prdist
import tsne
from MulticoreTSNE import MulticoreTSNE as TSNE
from scipy import random
from scipy import sparse

import collections

import pdb


class HomoSet(collections.MutableMapping):
    """HomoSet is the object containing evrey homology as a dictionnary according to thie rhomology code

        from Homoset you can do much computation that requires the set of homologies
        Object where we store an homology group basically where we do our entire
        Computation from.
        Be aware that even if you use str, the keys will be stored as unicode as jsonized dict will
        output unicode in python < 3

        Attributes:
            hashomo_matrix : a np.array[boolean] (species, homologies)
                that store the matrix of gene presence in species
            homo_matrix : np.array[float] (homologies*species(inthehomologies),aminoacids)
                but containing the codon entropy vectors instead
            homodict : dictionnary of homology object containing codon usage per species
                for the gene coresponding to this homology
            homo_namelist : list[str] of all the homology names
            species_namelist: list[str] of all the species names
            clusters: list[int] of clusters for the homoset clustering can be of size of nb of species
                or of nb of homologies
            homogroupnb: int of group for the homology clustering
            red_homomatrix: np.array[float] (homologies*species(inthehomologies),x*y)the reduced 2D version
                of the homology matrix for the full homology matrix
            wasclusterized: bool if the homologies have been clustered or not. usefull for processing requiring clusters
            homocluster: np.array[int] homologies*species of int of cluster number
            homo_matrixnames: np.array[int] corresponding to names in PyCUB.utils.speciestable
            fulleng: np.array[int] number of coding for each amino for each gene for each homology
            datatype: str of the CUB value type [entropy, A_value(entropyLocation), frequency]
            averagehomo_matrix: np.array[float] of the avg CUB values per homologies
            stats: dict of statistics on the clusterings (see get_clusterstats())

    """

    hashomo_matrix = None
    homo_matrix = None
    homo_matrixnames = None
    fulleng = None
    homodict = {}
    homo_namelist = []
    species_namelist = []
    clusters = []
    homogroupnb = 2
    red_homomatrix = None
    wasclusterized = False
    homo_clusters = None
    datatype = ''
    averagehomo_matrix = None
    stats = {}

    def __init__(self, **kwargs):
        """
        will initialize the object with the different values you might have from another project use the data dictionnary to add any type of data

        kwargs:
            a dictionary to any values present in the homoset
        """
        data = kwargs.get("data", None)
        if data is not None:
            utils.CUBD = data.get("CUBD", 18)
            self.homo_matrix = np.asarray(data["homo_matrix"]) if data.get(
                "homo_matrix", None) is not None else None
            self.homo_namelist = [str(i) for i in data.get("homo_namelist", [])]
            self.species_namelist = [str(i) for i in data.get("species_namelist", [])]
            self.homogroupnb = data.get("homogroupnb", 2)
            self.clusters = data.get("clusters", [])
            self.datatype = data.get("datatype", '')
            self.phylo_distances = pd.read_json(data["phylo_distances"], orient='split') if data.get(
                "phylo_distances", None) is not None else None
            self.red_homomatrix = np.asarray(data["red_homomatrix"]) if data.get(
                "red_homomatrix", None) is not None else None
            self.homo_matrixnames = np.asarray(data["homo_matrixnames"]) if data.get(
                "homo_matrixnames", None) is not None else None
            self.best_eps = data.get("best_eps", None)
            a = data.get("speciestable", {})
            if a:
                utils.speciestable = {}
                for key, val in a.iteritems():
                    utils.speciestable.update({int(key): val})
            if data.get("phylo_distances", False):
                utils.phylo_distances = pd.read_json(data.get("phylo_distances"), orient='split')
            utils.meandist = utils.phylo_distances.sum().sum() / (len(
                utils.phylo_distances)**2 - len(utils.phylo_distances)) if utils.phylo_distances is not None else None
            self.hashomo_matrix = np.asarray(data["hashomo_matrix"]) if data.get(
                "hashomo_matrix", None) is not None else None
            self.fulleng = np.asarray(data["fulleng"]) if data.get(
                "fulleng", None) is not None else None
            self.wasclusterized = data.get('wasclusterized', False)
            self.homo_clusters = np.asarray(data["homo_clusters"]) if data.get(
                "homo_clusters", None) is not None else None
            self.averagehomo_matrix = np.asarray(data["averagehomo_matrix"]) if data.get(
                "averagehomo_matrix", None) is not None else None
            self.stats = data.get('stats', {})
            self.homodict = {}
            for key, val in data["homodict"].iteritems():
                self.homodict.update({key: h.homology(data=val)})
            self.cluster_similarity = np.asarray(data['cluster_similarity']) if data.get("cluster_similarity", False) else None
            self.cub_similarity = np.asarray(data['cub_similarity']) if data.get("cub_similarity", False) else None
        else:
            self.homo_matrix = kwargs.get("homo_matrix", None)
            self.homo_namelist = kwargs.get("homo_namelist", [])
            self.species_namelist = kwargs.get("species_namelist", [])
            self.homogroupnb = kwargs.get("homogroupnb", 2)
            self.cluster_similarity = kwargs.get("cluster_similarity", None)
            self.cub_similarity = kwargs.get("cub_similarity", None)
            self.clusters = kwargs.get("clusters", [])
            self.datatype = kwargs.get("datatype", '')
            self.phylo_distances = kwargs.get("phylo_distances", None)
            self.red_homomatrix = kwargs.get("red_homomatrix", None)
            if kwargs.get("speciestable", {}):
                utils.speciestable = kwargs.get("speciestable")
            if kwargs.get("phylo_distances", False):
                utils.phylo_distances = kwargs.get("phylo_distances", None)
            utils.meandist = utils.phylo_distances.sum().sum() / (len(
                utils.phylo_distances)**2 - len(utils.phylo_distances)) if utils.phylo_distances is not None else None
            self.hashomo_matrix = kwargs.get("hashomo_matrix", None)
            self.fulleng = kwargs.get("fulleng", None)
            self.wasclusterized = kwargs.get('wasclusterized', False)
            self.homo_clusters = kwargs.get("homo_clusters", None)
            self.averagehomo_matrix = kwargs.get("averagehomo_matrix", None)
            self.homodict = kwargs.get("homodict", {})
            self.stats = kwargs.get('stats', {})
            self.best_eps = kwargs.get("best_eps", None)
            self.homo_matrixnames = kwargs.get("homo_matrixnames", None)
            if kwargs.get("CUBD", False):
                utils.CUBD = kwargs.get("CUBD", 18)
    # magic methods https://rszalski.github.io/magicmethods/

    def __getitem__(self, key):
        """
        get the homology at the corresponding key

        Args:
            key: str, unicode the key at which to get the homology

        Raises:
            TypeError "the type you should enter is int or unicode or str"
            KeyError
        """
        if type(key) is unicode:
            return self.homodict[key]
        elif type(key) is int:
            return self.homodict[self.homo_namelist[key]]
        elif type(key) is str:
            return self.homodict[unicode(key)]
        else:
            raise TypeError("the type you should enter is int or or unicode or str")

    def __setitem__(self, key, value):
        """
        add an item at the corresponding key

        Args:
            key: str, unicode the key at which to add
            val: PyCUB.homology to add at this key

        Raises:
            TypeError "the type you should enter is int or unicode or str"
            KeyError
        """
        if type(key) is unicode:
            self.homodict[key] = value
        elif type(key) is int:
            self.homodict[self.homo_namelist[key]] = value
        elif type(key) is str:
            self.homodict[unicode(key)] = val
        else:
            raise TypeError("the type you should enter is int or unicode or str")

    def __delitem__(self, key):
        """
        same as in dict()

        Args:
            key: str, unicode the key at which to add

        Raises:
            TypeError "the type you should enter is int or unicode or str"
            KeyError
        """
        if type(key) is unicode:
            del self.homodict[key]
        elif type(key) is int:
            del self.homodict[self.homo_namelist[key]]
        elif type(key) is str:
            del self.homodict[unicode(key)]
        else:
            raise TypeError("the type you should enter is int or unicode or str")

    def __iter__(self):
        """
        same as in dict()
        """
        return iter(self.homodict)

    def iteritems(self):
        """
        same as in dict()
        """
        return self.homodict.iteritems()

    def __len__(self):
        """
        gives you the length of the homoset ( number of homologies)
        """
        return len(self.homodict)

    def update(self, val):
        """
        a function to update the dictionnary that homoset is

        Args:
            val: the dict to append tot this one
        """
        for key, val in val.iteritems():
            self.homodict.update({unicode(key): val})


# CpCpGpApApTpApTpApTpTpCpCpGpApApTpApTpApTpTpCpCpGpApApTpApTpApTpTpCpCpGpApApTpApTpApTpTpTpTpCpCpGpApApTpApTpApTpTp
# GbGbCbTbTbAbTbAbTbAbAbGbGbCbTbTbAbTbAbTbAbAbGbGbCbTbTbAbTbAbTbAbAbGbGbCbTbTbAbTbAbTbAbAbAbGbGbCbTbTbAbTbAbTbAbAbAb

    def plot_all(self, With='tsne', perplexity=60, interactive=False, bins=100, offset=20, iteration=400, redo=False,
                 bypasstsne=False, dotsize=7, inpar=True):
        """
        will plot all the homologies in the full_homo_matrix (and compute it)

        (sometimes around 800 000 datapoints) to look at any kind of relationships
        as there is too much datapoints, the plots are density ones.

        Args:
            With: flag the dim reduction algorithm to use (tsne: need >16gigs of RAM,now use another version
                of tsne for large datasets )(PCA: works well)(lsta/hessian:untested)
            perplexity,iteration: ints of basic tsne hyperparams
            interactive: bool if true should use bokeh else matplotlib

        Returns:
            the desired plot if the size is high and we are interactive

        Raises:
            AttributeError: "not a good red algorithm"
        """
        if self.homo_matrix is None:
            self.loadfullhomo()
        size = len(self.homo_matrix)
        if size < 4000:
            size = 'small'
        elif size < 10000:
            size = 'medium'
        else:
            size = 'big'
        if self.red_homomatrix is None or redo:
            if With == 'tsne':
                if size == 'small' or bypasstsne:
                    if inpar:
                        red = TSNE(n_components=2, perplexity=30.0, verbose=2, n_jobs=-1,
                                   n_iter=iteration).fit_transform(self.homo_matrix)
                        # https://github.com/DmitryUlyanov/Multicore-TSNE}},
                    else:
                        red = man.TSNE(n_components=2, perplexity=30.0, verbose=2,
                                       n_iter=iteration).fit_transform(self.homo_matrix)
                else:
                    red = tsne.tsne(self.homo_matrix, no_dims=2, initial_dims=self.homo_matrix.shape[1],
                                    perplexity=30.0)
            elif With == 'PCA':
                red = PCA(n_components=2).fit_transform(self.homo_matrix)
            elif With == 'ltsa' or With == 'hessian':
                red = man.LocallyLinearEmbedding(self.homo_matrix.shape[0] / 100, 2,
                                                 eigen_solver='auto', method=With).fit_transform(self.homo_matrix)
            else:
                raise AttributeError("not a good red algorithm")
            self.red_homomatrix = red
        dotsize = 2 if size == 'big' else 5 if size == 'medium' else 10
        colormap = list(utils.colormap)
        speciestable = dict(utils.speciestable)
        if size == "small":
            colors = []
            homonames = []
            for i, val in enumerate(self.homodict.values()):
                colors.extend([colormap[i]] * len(val.names))
                homonames.extend([val.homocode] * len(val.names))
        if not interactive:
            fig = plt.figure(figsize=(40, 40))
            ax = fig.add_subplot(111)
            ax.scatter(self.red_homomatrix[:, 0], self.red_homomatrix[:, 1], s=dotsize, color=colors)
            plt.show()
            plt.savefig('utils/templot/full_homoset.pdf')
        else:
            if size == "small":
                data = dict(x=red[:, 0], y=red[:, 1],
                            species=[str(speciestable[n]) for n in self.homo_matrixnames],
                            meanentropy=["%.2f" % self.homo_matrix[i].mean() for i in range(len(self.homo_matrix))],
                            color=colors,
                            homologies=colors,
                            homonames=homonames,
                            lengths=self.fulleng.sum(1))
                source = ColumnDataSource(data=data)
                output_notebook()
                labe = ["show Homologies", "show Length", "show mean"]
                callback = CustomJS(args=dict(source=source), code=str(utils.callback_plotall))
                radio_button_group = widgets.RadioButtonGroup(
                    labels=labe, callback=callback, active=0)
                hover = HoverTool(tooltips=[("species", "@species"), ("mean_entr", "@meanentropy"),
                                            ("length", "@lengths"), ("homology: ", "@homonames")])
                p = figure(title="Full plot of the homoset",
                           tools=[hover, WheelZoomTool(), PanTool(), SaveTool(), ResetTool()],
                           plot_width=800, plot_height=600)
                p.circle(x='x', y='y', source=source, color='color', size=dotsize)
                save(column(radio_button_group, p), "utils/templot/full_homoset.html")
                show(column(radio_button_group, p))
            else:
                hv.notebook_extension('bokeh', 'matplotlib')
                dynspread.max_px = 200
                dynspread.threshold = 0.5
                shade.cmap = "#30a2da"  # to match HV Bokeh default
                points = hv.Points(self.red_homomatrix, label="all " + str(self.red_homomatrix.shape[0]) + " homologies")

                def heatmap(coords, bins=10, offset=20.0, transform=lambda d, m: d, label=None):
                    """
                    Given a set of coordinates, bins them into a 2d histogram grid
                    of the specified size, and optionally transforms the counts
                    and/or compresses them into a visible range starting at a
                    specified offset between 0 and 1.0.
                    """
                    hist, xs, ys = np.histogram2d(coords.T[0], coords.T[1], bins=bins)
                    counts = hist[:, ::-1].T
                    transformed = transform(counts, counts != 0)
                    span = transformed.max() - transformed.min()
                    compressed = np.where(counts != 0, offset + (1.0 - offset) * transformed / span, 0)
                    args = dict(label=label) if label else {}
                    return hv.Image(compressed, bounds=(xs[-1], ys[-1], xs[1], ys[1]), **args)
                print 'please write %output size=200" before calling this function'
                datashader = datashade(points)
                heatmaper = heatmap(self.red_homomatrix, bins=bins, offset=offset)(style=dict(cmap="fire"))
                save(hv.renderer('bokeh').get_plot(datashader).state, 'utils/templot/full_homoset_datashade.html')
                save(hv.renderer('bokeh').get_plot(heatmaper).state, 'utils/templot/full_homoset_heatmap.html')
                return heatmaper + datashader

    def loadfullhomo(self):
        """
        function to concatenate all the homologies in one big array(practicle for certain computations)

        Args:
            None
        """
        homo_matrix = self[self.homo_namelist[0]].full.copy()
        homo_matrixnames = list(self[self.homo_namelist[0]].names)
        fulleng = self[self.homo_namelist[0]].lenmat.copy()
        for x, homo in enumerate(self.homo_namelist[1:]):
            try:
                homo_matrix = np.vstack((homo_matrix, self[homo].full))
                homo_matrixnames.extend(self[homo].names)
                fulleng = np.vstack((fulleng, self[homo].lenmat))
            except ValueError:
                print x, self[homo].full.shape
                pdb.set_trace()
            print '{0}%\r'.format((x * 100) / len(self.homo_namelist)),
        self.homo_matrix = homo_matrix
        self.fulleng = fulleng
        self.homo_matrixnames = np.asarray(homo_matrixnames)
        if np.isinf(self.homo_matrix).any():
            self.homo_matrix[np.isinf(self.homo_matrix)] = 20
            for val in self.values():
                val.full[np.isinf(val.full)] = 20
        print "loaded"

    def loadhashomo(self, withnames=False):
        """
        function to compute the matrix of bool saying wether species X has a gene or more in homology Y

        Args:
            None
        """
        hashomo = np.zeros((len(self.homo_namelist),
                            len(self.species_namelist)), dtype=np.bool)
        for i, homo in enumerate(self.homo_namelist):
            ind = [na for y, na in enumerate(self[homo].names) if not self[homo].doub[y]]
            # get the indeces from each species name
            hashomo[i][ind] = True
        self.hashomo_matrix = hashomo
        print "loaded"

    def size(self):
        """
        the size of the homoset (number of genes)

        Args:
            None

        Returns:
            int the number of genes
        """
        if self.hashomo_matrix is None:
            self.loadhashomo()
        return np.count_nonzero(self.hashomo_matrix)

    def add_random_homology(self):
        """
        a function to populate a homology with random values

        Args:
            None

        Returns:
            the name of the random homology (str)
        """
        full = np.random.rand(len(self.species_namelist), utils.CUBD) * 20
        nans = np.zeros(len(self.species_namelist))
        names = []
        proteinids = []
        for _ in range(len(self.species_namelist)):
            names.append("species_" + str(int(random.rand() * 3000)))
            proteinids.append("proteinids_" + str(int(random.rand() * 4)))
        similarity_scores = np.random.rand(len(self.species_namelist))
        Kaks_Scores = np.random.rand(len(self.species_namelist)) * 1.5
        lenmat = np.random.rand(len(self.species_namelist), 18) * 100
        lenmat = lenmat.astype(int)
        homo = h.homology(full=full, nans=nans, names=names, similarity_scores=similarity_scores,
                          proteinids=proteinids, Kaks_Scores=Kaks_Scores, lenmat=lenmat)
        name = 'random_' + str(int(random.rand() * 1000))
        self.homo_namelist.append(name)
        self.homodict.update({unicode(name): homo})
        return name

    def preprocessing(self, withtaxons=False, withnames=True):
        """
        will compute the full list of names, 

        find doublons, and set the names to ints instead
        of strings. called after loading from ensembl and associate namelist in each homologies to a number
        in utils.speciestable

        Args:
            withtaxons: bool to true calls preprocessing_taxons() else one of the other two
            withnames: bool to true to call preprocessing_names() else calls preprocessing_namelist()

        Returns:
            taxons,species : list, if the names contains an additional list of taxon ids.
                and the corresponding species in the same order only if withtaxons
        """
        print "preprocessing..."
        if self.homodict is not None:
            if withtaxons:
                return self.preprocessing_taxons()
            elif withnames:
                self.preprocessing_names()
            else:
                return self.preprocessing_namelist()

    def preprocessing_taxons(self):
        """
        preprocess the data by computing the names as ints

        creates a speciestable to find the corresponding names
        computes the doublons and compute the species_namelist of this homoset
        will returns the species and the corresponding taxons  extracted from the homology.names

        Args:
            None
        """
        i = 0
        helper = {}
        species_namelist = set([])
        speciestable = {}
        taxons = []
        species = []
        for key, val in self.iteritems():
            doub = np.zeros(len(val.names[0]), dtype=np.bool)
            names = []
            temp = ''
            for j, name in enumerate(val.names[0]):
                if name == temp:
                    doub[j] = True
                else:
                    temp = name
                    old_l = len(species_namelist)
                    species_namelist.add(name)
                    if len(species_namelist) != old_l:
                        species.append(val.names[0][j])
                        taxons.append(val.names[1][j])
                        speciestable.update({i: name})
                        helper.update({name: i})
                        i += 1
                names.append(helper[name])
            val.names = names
            val.doub = doub
        utils.speciestable = dict(speciestable)
        self.species_namelist = species
        return taxons, species

    def preprocessing_names(self):
        """
        same as preprocessing_taxons() but admiting there is no taxon information (Yun's data for example)

        Args:
            None
        """
        species_namelist = set([])
        species = []  # we need to have another species list for the ordering
        # to be kept as it should (a set is ordered to be accessed in log time)
        speciestable = {}
        i = 0
        helper = {}
        for key, val in self.iteritems():
            doub = np.zeros(len(val.names), dtype=np.bool)
            names = []
            temp = ''
            for j, name in enumerate(val.names):
                if name == temp:
                    doub[j] = True
                else:
                    temp = name
                    old_l = len(species_namelist)
                    species_namelist.add(name)
                    if len(species_namelist) != old_l:
                        species.append(val.names[j])
                        speciestable.update({i: name})
                        helper.update({name: i})
                        i += 1
                names.append(helper[name])
            val.names = names
            val.doub = doub
        utils.speciestable = dict(speciestable)
        self.species_namelist = species

    def preprocessing_namelist(self):
        """
        same as preprocessing_names() but without preprocessing the names of each homologies only updating the species_namelist

        Args:
            None
        """
        # as a dict of int is ordered
        species_namelist = set([])
        speciestable = dict(utils.speciestable)
        for key, val in self.iteritems():
            for j, name in enumerate(val.names):
                species_namelist.add(speciestable[name])
        return list(species_namelist)

    def compute_ages(self, preserved=True, minpreserv=0.9, minsimi=0.85):
        """
        will compute whether or not a coding gene is highly preserved and if not how recent it is

        it uses the phylogenetic distances and the similarities amongst homlogies to try to find a
        good proxy

        Args:
            preserved: bool to true if we should find highly preserved genes or not
            minpreserv: float minimal percentage of homologous species that have this homology
            minsimi: float minimal avg similarity between genes to consider them highly preserved

        Raises:
            UnboundLocalError: "you need to have the similarity_scores"
        """
        if self.homodict[self.keys()[-1]].similarity_scores is None:
            raise UnboundLocalError("you need to have the similarity_scores")
        if self.homodict[self.keys()[-1]].isrecent is None:
            phylo_distances = utils.phylo_distances.copy()
            allowed = phylo_distances.index.tolist()
            speciestable = dict(utils.speciestable)
            for _, homo in self.homodict.iteritems():
                # if the homology is one of species that are closely related solely (meaning it is a recent one)
                if preserved:
                    # a homology is highly preserved if it is something that is shared by all
                    # and that has not changed in evolution (high similarity)
                    if len(homo.full) > minpreserv * len(self.species_namelist) and homo.similarity_scores.mean() > minsimi:
                        homo.ishighpreserved = True
                        homo.isrecent = False
                        continue
                    else:
                        homo.ishighpreserved = False
                    homo.isrecent = False
                    species = [speciestable[n] for i, n in enumerate(homo.names) if not homo.doub[i] and speciestable[n] in allowed]
                    if len(species) < 2:
                        # we can't decide ...
                        continue
                    dist = float(phylo_distances[species].loc[species].sum().sum()) / (len(species)**2 - len(species))
                    if dist < utils.meandist * minsimi:
                        # print "found a homology with mean " + str(dist)
                        homo.isrecent = dist / utils.meandist
        else:
            print "it was already loaded"

    def compute_entropyloc(self, using='computejerem'):
        """
        called if need entropy location and used ensembl data. you can always compute entropy location from entropy data.

        Will be much faster than doing it directly when calling
        ensembl's data as it computes the partition function
        only one for each lengths

        Args:
            using: str flags the partition algorithm to use

        Raises:
            UnboundLocalError: "you need to have your CUB values as entropy"
        """
        if self.datatype == 'entropy':
            if self.homo_matrix is None:
                self.loadfullhomo()
            self.homo_matrix = utils.getloc(self.homo_matrix, self.fulleng, using=using)
            self.red_homomatrix = None
            self.homo_clusters = None
            self.datatype = 'entropyloc'
            self.averagehomo_matrix = None
            pos = 0
            for x, homo in enumerate(self.homo_namelist):
                self[homo].full = self.homo_matrix[pos:pos + len(self[homo].full)]
                self[homo].var = self[homo].full.var(0)**(0.5)
                self[homo].mean = self[homo].full.mean(0)
                self[homo].clusters = None
                self[homo].centroids = None
                self[homo].metrics = {}
                self[homo].reduced = None
                self[homo].reduced_algo = None
                pos = pos + len(self[homo].full)
                print '{0}%\r'.format((x * 100) / len(self.homo_namelist)),
        else:
            raise UnboundLocalError("you need to have your CUB values as entropy")

    def remove(self, species):
        """
        remove this list of species from the homologies

        Args:
            species: list[str] the species to remove from all homologies
        """
        for key in self.keys():
            self[key].remove(species)

    def clean_species(self, thresh=0.3):
        """
        will remove a species from all the homologies of this homoset

        Warning, as the all/working homosets share the ref to homologies, deleting some species in working_homoset
        will result in removing some species in some homologies of all_homoset

        Args:
            thresh: float the threshold of avg presence in homologies below which the species are remove
        """
        spe = []
        perc = self.hashomo_matrix.sum(0) / len(self.hashomo_matrix)
        for i in range(len(perc)):
            if perc[i] < thresh:
                spe.append(self.species_namelist[i])
        self.remove(spe)

    def plot_homoperspecies(self):
        """
        will plot the number of homology per spcies in this homology group
        """
        sumed = np.sum(self.hashomo_matrix, axis=1)
        plt.figure(figsize=(40, 10))
        plt.title('number of homologies per species')
        plt.bar(range(len(sumed)), sumed)
        plt.savefig("utils/templot/clustersimilarity.pdf")
        print "you can always look at a particular range of species with 'homo_namelist' "

    def plot_speciesperhomo(self):
        """
        will plot the number of species per homologies in this homology group
        """
        sumed = np.sum(self.hashomo_matrix, axis=0)
        plt.figure(figsize=(40, 10))
        plt.title('number of species per homologies')
        plt.bar(range(len(sumed)), sumed)
        plt.savefig("utils/templot/clustersimilarity.pdf")
        print "you can always look at a particular range of species with 'homo_namelist' "

# CpCpGpApApTpApTpApTpTpCpCpGpApApTpApTpApTpTpCpCpGpApApTpApTpApTpTpCpCpGpApApTpApTpApTpTpTpTpCpCpGpApApTpApTpApTpTp
# GbGbCbTbTbAbTbAbTbAbAbGbGbCbTbTbAbTbAbTbAbAbGbGbCbTbTbAbTbAbTbAbAbGbGbCbTbTbAbTbAbTbAbAbAbGbGbCbTbTbAbTbAbTbAbAbAb

    def cluster_homologies(self, clustering='kmeans', byspecie=False, order=True,
                           plot_ordering=True, homogroupnb=2, findnb=False):
        """
        Compute an homology group :

        from matrix computation using the homo_matrix
        (or from network computation in homologize_from_network)
        Can be computed many times and will updata homoset with the most recent homoset found
        if homoset exists, it will save it.

        Args:
            clustering: str flags to 'kmeans', 'kmodes', 'fast' to use different sk-learn algorithms
            plot: bool flags to true for the function to output ploting
                of the affinity matrix with and without the clusters
            homogroupnb: int nb of groups you want to extract
            byspecie: bool to true if we cluster by species instead of homologies
            order: bool whether or not to order
            plot_ordering: bool to true to plot this ordering
            findnb: bool to true to find the right number of clusters (homogroupnb)

        Returns:
            if findnb, will return the clusters for each homogroupnb up to 9

        Raises:
            AttributeError: "you entered a wrong clustering algorithm"

        """
        if findnb is True:
            homogroupnb = 2
        clusts = []
        while True:
            if clustering == "kmeans":
                # The parallel version of K-Means is broken on OS X
                # when numpy uses the Accelerate Framework.
                # This is expected behavior: Accelerate can be called after a fork
                # but you need to execv the subprocess with the Python binary
                # (which multiprocessing does not do under posix)
                alg = KMeans(n_clusters=homogroupnb, n_jobs=-1)
            elif clustering == "fast":
                alg = MiniBatchKMeans(n_clusters=homogroupnb)

            elif clustering == "kmodes":
                # https://github.com/nicodv/kmodes/blob/master/kmodes/kmodes.py
                alg = KModes(n_clusters=homogroupnb,
                             init='Huang', n_init=2, verbose=1)

            else:
                raise AttributeError("you entered a wrong clustering algorithm")

            self.hashomo_matrix = self.hashomo_matrix.T if byspecie else self.hashomo_matrix
            # BYSPECIES is a special case where one would like to keep all homologies and
            # remove species instead. we would advise the opposite though
            alg.fit(self.hashomo_matrix)
            clust = alg.labels_.astype(int)
            metricA = metrics.silhouette_score(self.hashomo_matrix, clust, metric='euclidean')
            metricB = metrics.calinski_harabaz_score(self.hashomo_matrix, clust)
            if findnb is True:
                # TODO: totest
                clusts[homogroupnb] = clust
                if homogroupnb < 10:
                    homogroupnb += 1
                    print "for " + homogroupnb + " cluster, the metric is " + metricA, metricB
                else:
                    print "you just select the clust by doing PyCUB.homo.clusters =\
                    output['groupnbyouwant']"
                    print "you an run the 'orderfromclust' also to have the plottings and orderings"
                    return clusts
            else:
                break
        print "the quality of the clustering is: [silhouette_score,calinski_harabaz_score]"
        print metricA
        print metricB
        if order:
            self.orderfromclust(homogroupnb, clust, byspecie=byspecie,
                                plot_ordering=plot_ordering)
        else:
            self.hashomo_matrix = self.hashomo_matrix.T if byspecie else self.hashomo_matrix

    def orderfromclust(self, homogroupnb, clust, byspecie=False, plot_ordering=True):
        """
        creates an ordering of every elements 

        (names, homologies according to the found clusters)
        from an ordered cluster of species or of homologies
        self.hashomomatrix should reflect this orientation as well

        Args:
            homogroupnb: int the number of clusters
            clust: np.array(int) the clusters
            byspecie: bool to true to order by species
            plot_ordering: bool to true to plot the ordering

        Returns:
            the ordered hashomomatrix (np.array[bool])
        """
        orderedhas = np.zeros(self.hashomo_matrix.shape, dtype=bool)
        ltemp = [0] * len(self.homo_namelist) if not byspecie else [0] * len(self.species_namelist)
        self.clusters = [0] * len(clust)
        begin = 0
        # reorder all matrices
        for i in range(homogroupnb):
            ind = np.argwhere(clust == i)[:, 0]
            orderedhas[begin:begin + len(ind)] = self.hashomo_matrix[ind]
            # the list as well
            sublist = [self.species_namelist[e] for e in ind] if byspecie else [self.homo_namelist[e] for e in ind]
            ltemp[begin:begin + len(ind)] = sublist
            self.clusters[begin:begin + len(ind)] = clust[ind]
            begin += len(ind)
        if byspecie:
            self.species_namelist = ltemp
        else:
            self.homo_namelist = ltemp
        if plot_ordering:
            self._plot_clust(self.hashomo_matrix, orderedhas)

        if byspecie:
            self.hashomo_matrix = orderedhas.T
        else:
            self.hashomo_matrix = orderedhas
        self.homogroupnb = homogroupnb
        return orderedhas

    def get_clusterstats(self, sort=True, interactive=True, redo=False):
        """
        will find the number of cluster i per homologies and per species

        plot for each species, how much its genes are outliers, how much are belonging
        to a secondary cluster and how much are belonging to the principal cluster.
        --> create a long stacked bar plot with these values

        Args:
            sort: bool to true to sort everyhting according to the statistics differences
            interactive: bool to true to have an interactive barplot
            redo: bool to true not to reused cached data

        Raises:
            UnboundLocalError: "you need to find clusters first"
        """
        homoclusters = np.zeros((len(self.homodict), 3))
        speciestable = dict(utils.speciestable)
        specluster = np.zeros((len(self.species_namelist), 3))
        if self.wasclusterized:
            if not self.stats or redo:
                for j, homo in enumerate(self.homo_namelist):
                    x = -1
                    val = self[homo]
                    hprop = np.zeros(3)
                    for i in range(len(val.clusters)):
                        if not val.doub[i]:
                            x += 1
                        if val.clusters[i] == 0:  # primary clusters
                            hprop[0] += 1
                            specluster[val.names[x]][0] += 1
                        elif val.clusters[i] == -1:  # outliers clusters
                            hprop[1] += 1
                            specluster[val.names[x]][1] += 1
                        else:   # secondary clusters
                            hprop[2] += 1
                            specluster[val.names[x]][2] += 1

                    homoclusters[j] = np.divide(hprop, 1 + i) if i != 0 else np.array([1, 0, 0])
                specluster = np.divide(specluster.T, specluster.sum(1)).T
                specluster = specluster[~np.isnan(specluster).any(axis=1)]
                self.stats = {'homologies': homoclusters.tolist(), 'species': specluster.tolist()}
                if sort:
                    ind = sorted(range(len(self.stats['homologies'])), key=lambda k: self.stats['homologies'][k])
                    self.stats['homologies'].sort()
                    if self.hashomo_matrix is not None:
                        self.hashomo_matrix[:] = self.hashomo_matrix[ind]
                    if self.homo_namelist is not None:
                        self.homo_namelist[:] = [self.homo_namelist[i] for i in ind]
                    if self.clusters:
                        self.clusters[:] = [self.clusters[i] for i in ind]
                    if self.homo_matrix is not None:
                        self.homo_matrix = self.homo_matrix[ind]
                    if self.homo_matrixnames is not None:
                        self.homo_matrixnames = self.homo_matrixnames[ind]
                    if self.fulleng is not None:
                        self.fulleng = self.fulleng[ind]
                    if self.homo_clusters is not None:
                        self.homo_clusters = self.homo_clusters[ind]
                    ind = sorted(range(len(self.stats["species"])), key=lambda k: self.stats["species"][k])
                    self.species_namelist = speciestable.values()
                    self.species_namelist = [self.species_namelist[i] for i in ind]
                    for val in speciestable.values():
                        if val not in self.species_namelist:
                            self.species_namelist.append(val)
                    self.stats["species"].sort()
                    print "beware, species_namelist does not represent speciesdict/loadhashomo/loadfullhomo/etc"
                    print "just do homoset.species_namelist = cub.speciesdict().values() to retrieve the correspondance"
            else:
                print "using cached data to save computation"
            return self._barplot()
        else:
            raise UnboundLocalError("you need to find clusters first")

    def compare_clusters(self, cubdistance_matrix=True, plot=True, interactive=True, size=40):
        """
        for each clusters in homologies, will compare them with a similarity matrix and a distance matrix

        compare amongst the working homoset homologies, the clusters together,
        by what species they contains by creating a new vector of species presence
        in each cluster and plotting the similarity matrix of each of those vectors.
        --> create a compare function in homoset of
        homologies clusters similarity matrix and ordering.
        basically the distance should be nan if it has not the species,
        -1 if outlier to other and 1 if one cluster to another and zeros if the same to the same

        Args:
            cubdistance_matrix: bool to true if want to compute the matrix of the averageCUB value distances
                summed for each cluster amongst the homologies
            plot: bool to true to plot
            size: int the size of the plot
        """
        j = 0
        simimatrix = np.zeros((len(self.homodict), len(self.homodict)))
        for _, homo in self.iteritems():
            a = set(homo.names)
            se = np.zeros(len(self.species_namelist), dtype=int) - 2
            se[np.array(homo.names)[np.invert(homo.doub)]] = np.array(homo.clusters)[np.invert(homo.doub)]
            i = 0
            for _, homocmp in self.iteritems():
                if i < j:
                    simimatrix[j, i] = simimatrix[i, j]
                    i += 1
                    continue
                if i == j:
                    simimatrix[i, i] = 1.
                    i += 1
                    continue
                similar = len(a.intersection(homocmp.names))
                if not similar:
                    i += 1
                    continue
                comp = np.zeros(len(self.species_namelist), dtype=int) - 3
                comp[np.array(homocmp.names)[np.invert(homocmp.doub)]] = np.array(homocmp.clusters)[np.invert(homocmp.doub)]
                # for each clusters find how much their clusters share same species,
                # normalized by the amount of species they have in common
                # as it will count as similar species they do
                simimatrix[j, i] = float((se == comp).sum()) / similar
                i += 1
            j += 1
        self.cluster_similarity = simimatrix.copy()
        if plot:
            simiclust = self.plot_simiclust(interactive=interactive, size=size)
        if cubdistance_matrix:
            j = 0
            simimatrix = np.zeros((len(self.homodict), len(self.homodict)))
            for _, homo in self.iteritems():
                a = set(homo.names)
                se = np.zeros((len(self.species_namelist), utils.CUBD))
                se[np.array(homo.names)[np.invert(homo.doub)]] = homo.full[np.invert(homo.doub)]
                i = 0
                for _, homocmp in self.iteritems():
                    if i < j:
                        simimatrix[j, i] = simimatrix[i, j]
                        i += 1
                        continue
                    if i == j:
                        simimatrix[i, i] = 0.
                        i += 1
                        continue
                    similar = len(a.intersection(homocmp.names))
                    if not similar:
                        simimatrix[j, i] = np.NaN
                        i += 1
                        continue
                    comp = se.copy()
                    # we need to have the values of the homo we compare to that are not in homo we compare
                    # from, to be the same in order for the distance to be zero at these values
                    comp[np.array(homocmp.names)[
                        np.logical_and(np.invert(homocmp.doub),
                                       se.any(1)[homocmp.names])]] = homocmp.full[
                        np.logical_and(np.invert(homocmp.doub),
                                       se.any(1)[homocmp.names])]
                    simimatrix[j, i] = prdist(se, comp).sum() / similar
                    i += 1
                j += 1
            simimatrix = simimatrix / simimatrix.max()
            self.cub_similarity = 1 - simimatrix
            self.cub_similarity = np.nan_to_num(self.cub_similarity)
            if plot:
                distcub = self.plot_distcub(interactive=interactive, size=size)
                if interactive:
                    show(column(distcub, simiclust))
        else:
            if interactive:
                show(simiclust)

    def find_clusters(self, clustering='dbscan', homogroupnb=None,
                      assess=True, eps=0.8, best_eps=True, trainingset=30,
                      epoch=20, ranges=[0.2, 0.9], size=10, redo=False):
        """
        Finds, for each homologies in the working homoset, groups that are part of compact clusters

        it will be using gaussian mixture clustering or DBSCAN and order them according
        to the density of each cluster (we are interested in the densest ones) and assess
        the quality using 3 criterion:BIC, number of outliers,
        also: - find if we are close to ancestry tree,
        here we need to represent a comparison of the closeness
        in a phylogenetic tree to a cluster of species
        --> given a grouping of phylogenetic tree, what cluster is the most similar to it

        Args:
            clustering: method (DBSCAN, gaussian mixture)
            homogroupnb: the number of groups can be a number or else will look for the better number of
                cluster according to assessments.
            assess: plot or not the assessments
            eps: the hyperparams
            best_eps: bool to true wether or not to do an hyperparam greedy search
            trainingset: int if best_eps, the size of the training ser
            epoch: int the number of increasing trials
            ranges: tuple[float] the two min and max values to use
            size: int the x size of the plot
        """
        if best_eps:
            eps = self.findbest_eps(trainingset, clustering, epoch=epoch, ranges=ranges, size=size, redo=redo)
        for val in self.homo_namelist:
            print val
            self[val].clusterize_(clustering=clustering, eps=eps,
                                  homogroupnb=homogroupnb, assess=assess)
            print "-----------------------------"
        self.wasclusterized = True
        if best_eps:
            score = 0
            for _, homo in self.iteritems():
                if len(homo.metrics["cluster_phylodistance"]) == 1:
                    score += 1
                score += np.array(homo.metrics["cluster_phylodistance"][1:]).mean() - (homo.metrics["cluster_phylodistance"][0] - 1)
            print "-----------------TOT------------------"
            print score

    def findbest_eps(self, trainingset=400, clustering="dbscan",
                     epoch=20, ranges=[0.2, 0.9], size=10, redo=False):
        """
        will find the best eps hyperparameter (the one that minimizes the evolutionary distance within its clusters)

        Args:
            trainingset: int the number of homologies in your training set (should be 20% of the total)
            clustering: str flags the clustering algorithm from which to find the best hyperparam
            epoch: int the number of increasing trials
            ranges: tuple[float] the two min and max values to use
            size: int the x size of the plot
            redo: wether or not to recompute everyhting if it has already been computed once

        Returns;
            the best value (often a float)

        Raises:
            LookupError: "you need to compute the phylogenetic distances first to have some form of labels"
        """
        if self.best_eps is None or redo:
            scores = []
            vals = []
            best_score = 10000000
            best_eps = 0
            val = ranges[0]
            if utils.phylo_distances is None:
                raise LookupError("you need to compute the phylogenetic distances first to have some form of labels")
            if clustering == "gaussian":
                epoch = 8
                val = 1
            # random draw
            homolist = [self.homo_namelist[int(random.rand() * len(self.homo_namelist))] for _ in range(trainingset)]
            for epochstep in range(epoch):
                score = 0
                print '\repoch: ' + str(epochstep),
                for ste in range(trainingset):
                    homo = self[homolist[ste]]
                    homo.clusterize_(clustering=clustering,
                                     eps=val, homogroupnb=val, assess=True, verbose=False)
                    sc = 0
                    for i, _ in enumerate(homo.metrics["cluster_phylodistance"]):
                        # weighted mean of the phylo distances inter clusters by the proportion of datapoint
                        # in each of them
                        if len(homo.metrics["cluster_phylodistance"]) > 1:
                            weight = (float(len(
                                np.argwhere(np.array(homo.clusters) == i - 1).T[0])) / len(homo.clusters))
                        else:
                            weight = 1
                        sc += homo.metrics["cluster_phylodistance"][i] * weight
                    score += sc
                # get the best eps of all
                # compute its similarity on the full set of homologies
                score = score / ste
                scores.append(score)
                vals.append(val)
                if score < best_score:
                    best_score = score
                    best_eps = val
                val += float((ranges[1] - ranges[0])) / (epoch - 1) if clustering == "dbscan" else 1
            self.best_eps = best_eps
            plt.figure(figsize=(size, 20))
            plt.title('scores for best clustering hyperparameter')
            plt.plot(vals, scores)
            plt.savefig("utils/templot/scores_" + clustering + "_" + str(epoch) + "epochs.pdf")
            plt.show()
        print "the best value is : " + str(self.best_eps)
        return self.best_eps

# CpCpGpApApTpApTpApTpTpCpCpGpApApTpApTpApTpTpCpCpGpApApTpApTpApTpTpCpCpGpApApTpApTpApTpTpTpTpCpCpGpApApTpApTpApTpTp
# GbGbCbTbTbAbTbAbTbAbAbGbGbCbTbTbAbTbAbTbAbAbGbGbCbTbTbAbTbAbTbAbAbGbGbCbTbTbAbTbAbTbAbAbAbGbGbCbTbTbAbTbAbTbAbAbAb

    def _barplot(self, interactive=True):
        """
        called by statistics function to plot a barplot of the proportion of different cluster per homologies and per species

        Args:
            interactive: bool to true if you want to use the bokeh interactive version
        Returns:
            the barplot holoviews object (will directly render if in a notebook)
        """
        hv.notebook_extension('bokeh') if interactive else hv.notebook_extension('matplotlib')
        print 'please write %output size=100" before calling this function'
        label = ["outliers", "primary", "secondary"]
        dims = dict(kdims='homologies', vdims='props')
        plot = hv.Overlay([hv.Area(np.array(self.stats["homologies"]).T[i], label=label[i], **dims) for i in range(3)])
        overlay = hv.Area.stack(plot).relabel("cluster info: primary_clust +outlier +secondary_clust for homologies")
        save(hv.renderer('bokeh').get_plot(overlay).state, 'utils/templot/barplot_homoset.html')
        dims = dict(kdims='species', vdims='props')
        plotspe = hv.Overlay([hv.Area(np.array(self.stats["species"]).T[i], label=label[i], **dims) for i in range(3)])
        overlayspe = hv.Area.stack(plotspe).relabel("cluster info: primary_clust +outlier +secondary_clust for species")
        save(hv.renderer('bokeh').get_plot(overlayspe).state, 'utils/templot/barplot_species.html')
        return overlay + plot + overlayspe + plotspe

    def _dictify(self, savehomos=False):
        """
        Used by the saving function. transform the object into a dictionary that can be json serializable

        Args:
            savehomos: bool to true if you consider this homology as containing all the information
            (the other ones only have references to a subset of the data of this one)

        Returns:
            the dictionnary of all the data in this Object in the correct format to be jsonized
        """
        dictihomo = {}
        if savehomos:
            for key, val in self.iteritems():
                dictihomo.update({key: val._dictify()})

        return {"hashomo_matrix": self.hashomo_matrix.tolist() if self.hashomo_matrix is not None else None,
                "homo_matrix": self.homo_matrix.tolist() if self.homo_matrix is not None else None,
                "clusters": self.clusters,
                "best_eps": self.best_eps,
                "homo_namelist": self.homo_namelist,
                "cluster_similarity": self.cluster_similarity.tolist() if self.cluster_similarity is not None else None,
                "cub_similarity": self.cub_similarity.tolist() if self.cub_similarity is not None else None,
                "homodict": dictihomo,
                "species_namelist": self.species_namelist,
                "speciestable": dict(utils.speciestable) if savehomos else None,
                "homogroupnb": self.homogroupnb,
                "homo_matrixnames": self.homo_matrixnames.tolist() if self.homo_matrixnames is not None else None,
                "stats": self.stats,
                "datatype": self.datatype,
                "phylo_distances": utils.phylo_distances.to_json(orient='split') if
                utils.phylo_distances is not None and savehomos else None,
                "CUBD": utils.CUBD,
                "red_homomatrix": self.red_homomatrix.tolist() if self.red_homomatrix is not None else None,
                "fulleng": self.fulleng.tolist() if self.fulleng is not None else None,
                "wasclusterized": self.wasclusterized,
                "homo_clusters": self.homo_clusters.tolist() if self.homo_clusters is not None else None,
                "averagehomo_matrix": self.averagehomo_matrix.tolist() if self.averagehomo_matrix is not None else None}

    def plot_hashomo(self, invert=False, size=40, interactive=False, rangeto=None):
        """
        plot the has homo matrix

        the interactive version allows you to see each particular datapoint with much more precision

        Args:
            interactive: bool to use the bokeh interactive version
            size: int size of the matrix
            invert: bool flag to true to invert the plot

        """
        if interactive:
            print "no invert here"
            colormap = list(utils.colormap)
            xname = []
            yname = []
            if rangeto is None:
                rangeto = range(len(self.homo_namelist))
            color = []
            homoclust = len(self.clusters) == len(self.homo_namelist)
            for i in range(len(self.species_namelist)):
                for j in rangeto:
                    xname.append(self.homo_namelist[j])
                    yname.append(self.species_namelist[i])
                    if homoclust:
                        color.append(colormap[self.clusters[j]] if self.hashomo_matrix[j, i] else '#999999')
                    else:
                        color.append(colormap[0] if self.hashomo_matrix[j, i] else '#999999')
            data = dict(
                xname=xname,
                yname=yname,
                colors=color,
            )
            output_notebook()
            hover = HoverTool(tooltips=[('names: ', '@yname, @xname')])
            p = figure(title="presence of homologies between species",
                       x_range=list(reversed([self.homo_namelist[y] for y in rangeto])),
                       y_range=self.species_namelist,
                       x_axis_location="above", tools=[hover, WheelZoomTool(), PanTool(), SaveTool(), ResetTool()])

            p.plot_width = 2000
            p.plot_height = 2000
            p.grid.grid_line_color = None
            p.axis.axis_line_color = None
            p.axis.major_tick_line_color = None
            p.axis.major_label_text_font_size = "5pt"
            p.axis.major_label_standoff = 0
            p.xaxis.major_label_orientation = np.pi / 3
            p.rect('xname', 'yname', 0.9, 0.9, source=data,
                   color='colors', line_color=None,
                   hover_line_color='black', hover_color='colors')
            save(p, 'utils/templot/hashomomatrix.html')
            show(p)  # show the plot
        else:
            plt.figure(figsize=(size, 200))
            plt.title('the regular matrix')
            plt.imshow(self.hashomo_matrix.T if invert else self.hashomo_matrix)
            plt.savefig("utils/templot/hashomomatrix.pdf")
            plt.show()

    def plot_simiclust(self, interactive=True, size=40):
        """
        plot the similarity matrix of each homologies from its clusters

        the interactive version allows you to see each particular datapoint with much more precision

        Args:
            interactive: bool to use the bokeh interactive version
            size: size of the matrix
        """
        if self.cluster_similarity is not None:
            if interactive:
                colormap = list(utils.colormap)
                xname = []
                yname = []
                xname.extend(self.homo_namelist * len(self.homo_namelist))
                for i in self.homo_namelist:
                    yname.extend([i] * len(self.homo_namelist))
                data = dict(
                    xname=xname,
                    yname=yname,
                    colors=[colormap[1]] * (len(self.homo_namelist)**2),
                    alphas=self.cluster_similarity.flatten(),
                )
                output_notebook()
                hover = HoverTool(tooltips=[('names: ', '@yname, @xname'), ('similarity: ', '@alphas')])
                p = figure(title="similarity matrix of each homologies from its clusters",
                           x_range=list(reversed(self.homo_namelist)), y_range=self.homo_namelist,
                           x_axis_location="above", tools=[hover, WheelZoomTool(), PanTool(), SaveTool(), ResetTool()])
                p.plot_width = 800
                p.plot_height = 800
                p.grid.grid_line_color = None
                p.axis.axis_line_color = None
                p.axis.major_tick_line_color = None
                p.axis.major_label_text_font_size = "5pt"
                p.axis.major_label_standoff = 0
                p.xaxis.major_label_orientation = np.pi / 3
                p.rect('xname', 'yname', 0.9, 0.9, source=data,
                       color='colors', alpha='alphas', line_color=None,
                       hover_line_color='black', hover_color='colors')
                save(p, 'utils/templot/simicluster.html')
                return p  # show the plot
            else:
                plt.figure(figsize=(size, 200))
                plt.title('the homologies similarity of their clusters')
                plt.imshow(self.cluster_similarity)
                plt.savefig("utils/templot/simiclust.pdf")
                plt.show()

    def plot_distcub(self, interactive=False, size=40):
        """
        plot the distance matrix of each homologies from the average of their CUB distances

        the interactive version allows you to see each particular datapoint with much more precision

        Args:
            interactive: bool to use the bokeh interactive version
            size: size of the matrix

        """
        if self.cub_similarity is not None:
            if interactive:
                colormap = list(utils.colormap)
                xname = []
                yname = []
                xname.extend(self.homo_namelist * len(self.homo_namelist))
                for i in self.homo_namelist:
                    yname.extend([i] * len(self.homo_namelist))
                data = dict(
                    xname=xname,
                    yname=yname,
                    colors=[colormap[2]] * (len(self.homo_namelist)**2),
                    alphas=self.cub_similarity.flatten(),
                )

                output_notebook()
                hover = HoverTool(tooltips=[('names: ', '@yname, @xname'), ('distcub: ', '@alphas')])
                p = figure(title="the distance matrix of each homologies from the average of their CUB distances",
                           x_range=list(reversed(self.homo_namelist)), y_range=self.homo_namelist,
                           x_axis_location="above", tools=[hover, WheelZoomTool(), PanTool(), SaveTool(), ResetTool()])
                p.plot_width = 800
                p.plot_height = 800
                p.grid.grid_line_color = None
                p.axis.axis_line_color = None
                p.axis.major_tick_line_color = None
                p.axis.major_label_text_font_size = "5pt"
                p.axis.major_label_standoff = 0
                p.xaxis.major_label_orientation = np.pi / 3
                p.rect('xname', 'yname', 0.9, 0.9, source=data,
                       color='colors', alpha='alphas', line_color=None,
                       hover_line_color='black', hover_color='colors')
                save(p, 'utils/templot/distcub.html')
                return p  # show the plot
            else:
                plt.figure(figsize=(size, 200))
                plt.title('the homologies similarities to the avg dist of their CUB')
                plt.imshow(self.cub_similarity)
                plt.savefig("utils/templot/distcub.pdf")
                plt.show()

    def _plot_clust(self, mat, orderedmat):
        """
        will plot the correlation matrix of the has_homomatrix before and after ordering allowing one to show its effect

        Args:
            mat: np.array[bool] the current homology matrix
            orderedmat: np.array[bool] the ordered homology matrix
        """

        # the regular matrix
        plt.figure(figsize=(40, 40))
        plt.title('the regular matrix')
        plt.imshow(mat.T)
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
        plt.imshow(orderedmat.T)
        plt.show()
        plt.savefig("utils/templot/similarity_matrix.pdf")
        # affinity of the ordered matrix
        mat_sparseo = sparse.csr_matrix(orderedmat)
        similaritieso = cosine_similarity(mat_sparseo)
        plt.figure(figsize=(40, 40))
        plt.title('the affinity of the ordered matrix')
        plt.imshow(similaritieso)
        plt.show()
