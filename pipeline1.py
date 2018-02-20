import numpy as np
import scipy as sp
from sklearn import preprocessing as prep
from sklearn import manifold as man
from sklearn import cluster
import pandas as pd
import networkx as netx
import glob
import os
import cpickle as pickle  # usefull to save our object

import matplotlib.pyplot as plt
from bokeh.plotting import *
from bokeh.models import *

CodonsDict = {
    'TTT': 0, 'TTC': 0, 'TTA': 0, 'TTG': 0, 'CTT': 0,
    'CTC': 0, 'CTA': 0, 'CTG': 0, 'ATT': 0, 'ATC': 0,
    'ATA': 0, 'ATG': 0, 'GTT': 0, 'GTC': 0, 'GTA': 0,
    'GTG': 0, 'TAT': 0, 'TAC': 0, 'TAA': 0, 'TAG': 0,
    'CAT': 0, 'CAC': 0, 'CAA': 0, 'CAG': 0, 'AAT': 0,
    'AAC': 0, 'AAA': 0, 'AAG': 0, 'GAT': 0, 'GAC': 0,
    'GAA': 0, 'GAG': 0, 'TCT': 0, 'TCC': 0, 'TCA': 0,
    'TCG': 0, 'CCT': 0, 'CCC': 0, 'CCA': 0, 'CCG': 0,
    'ACT': 0, 'ACC': 0, 'ACA': 0, 'ACG': 0, 'GCT': 0,
    'GCC': 0, 'GCA': 0, 'GCG': 0, 'TGT': 0, 'TGC': 0,
    'TGA': 0, 'TGG': 0, 'CGT': 0, 'CGC': 0, 'CGA': 0,
    'CGG': 0, 'AGT': 0, 'AGC': 0, 'AGA': 0, 'AGG': 0,
    'GGT': 0, 'GGC': 0, 'GGA': 0, 'GGG': 0}


# this dictionary shows which codons encode the same AA
SynonymousCodons = {
    'CYS': ['TGT', 'TGC'],
    'ASP': ['GAT', 'GAC'],
    'SER': ['TCT', 'TCG', 'TCA', 'TCC', 'AGC', 'AGT'],
    'GLN': ['CAA', 'CAG'],
    #'MET': ['ATG'],
    'ASN': ['AAC', 'AAT'],
    'PRO': ['CCT', 'CCG', 'CCA', 'CCC'],
    'LYS': ['AAG', 'AAA'],
    'STOP': ['TAG', 'TGA', 'TAA'],
    'THR': ['ACC', 'ACA', 'ACG', 'ACT'],
    'PHE': ['TTT', 'TTC'],
    'ALA': ['GCA', 'GCC', 'GCG', 'GCT'],
    'GLY': ['GGT', 'GGG', 'GGA', 'GGC'],
    'ILE': ['ATC', 'ATA', 'ATT'],
    'LEU': ['TTA', 'TTG', 'CTC', 'CTT', 'CTG', 'CTA'],
    'HIS': ['CAT', 'CAC'],
    'ARG': ['CGA', 'CGC', 'CGG', 'CGT', 'AGG', 'AGA'],
    #'TRP': ['TGG'],
    'VAL': ['GTA', 'GTC', 'GTG', 'GTT'],
    'GLU': ['GAG', 'GAA'],
    'TYR': ['TAT', 'TAC']}


class Genes(object):
    """docstring for Genes

        the only restriction on the file we want is any thing readable by pandas (txt, csv, xls,...)
        with the row 'species' and 'entropylocation'. the file should look like : 
            GENE-AMINOACID.*
            the - can be anything else...


    """

    def download_data(self, name='first500'):
    """download a file from the file list with the url of its location
 
    using urllib, you can add you own name and location in this global parameter
 
    Parameters:
    -----------
 
    name: str
        the path of the file correspondong to a file in the filelist 
        (''Sue_2x_3000_40_-46.tif' or 'demoMovieJ.tif')
    Raise:
    ---------
        WrongFolder Exception
    """

    #\bug
    #\warning

    file_dict = {'first500': 'https://www.dropbox.com/s/9j48pg1eixpnub4/homology601t1000.zip?dl=1'}
    path_data = './' + name
        if not os.path.exists(path_movie):
            url = file_dict[name]
            print("downloading " + name + "with urllib")
            f = urlopen(url)
            data = f.read()
            with open(path_movie, "wb") as code:
                code.write(data)
        else print("File already downloaded")

    def __init__(self, folder='first500', doAll=False, genelist=["YAL019W"], aminonb=18,
                 minspecies=1, separation="homology"):
        """
        initialize ou codonClass with 


        Returns:
        ---------
        genevect : dict{genename : gene}

        specieslist : list[list[speciesname]]

        folder : nameofthefolder

        genelist = list[genename]

        score = dict{species : #gene}

        gene = pandas.dataframe(row = speciesname, column = aminoname, data = entropyvalue)
        """
        super(Genes, self).__init__()
        nameA = "empty"
        nameB = "empty"
        # we create all our datastructures
        self.genevect = {}
        self.specieslist = []
        self.folder = folder
        self.doAll = doAll
        self.score = {}
        self.isminspecified = minspecies > 1  # we are going to keep that for further usage when dev19.
        self.saved = False

        if not os.path.isdir(folder):
            download_data()

        if doAll:
            # we verbose a bit in this class
            print "we are doing all the " + len(os.listdir(self.folder)) + " files"
            for f in os.listdir(self.folder):
                nameA = f.split(separation)[0]
                if(nameA != nameB):
                    nameB = nameA
                    self.genelist.append(nameB)
        else:
            self.genelist = genelist
        for gene in self.genelist:
            try:
                genDF, species = self.readcods_gene(folder, separation, gene, aminonb, separation)
                if len(species) > minspecies:
                    self.genevect.update({gene: genDF})
                    self.specieslist.append(species)
                else print gene + "contains less species than minspecies"
            except OSError:
                print "you do not have the files here"
            except ValueError:
                print gene + " has non matching components for a same gene.."

    def readcods_gene(self, separation, folder="kalfonDTA", gene="YAL019W", aminonb=18):
        """
        read the things and gives you back a nice dict

        :param folder: the folder you wanan look onto
        :param type: the type of gene you are looking for
        :return two dictionaries : one being a set of dictionaries looking like the file we have
                                    the other a dictionnary of homologous genes and their
                                    18 components values for each species
                a list of species



        """
        gendict = {}
        # We want to get the basic information from the files
        first_file = glob.glob(folder + "/" + gene + separation + "His.*")[0]
        print gene
        meta = pd.read_csv(first_file).dropna()
        rows = meta.shape[0]
        gentab = np.zeros((aminonb, rows))
        amino = []
        i = 0
        for file in glob.glob(folder + "/" + gene + "*.*"):
            if file[-7:-4] != 'ror':  # TODO: write if it belongs to a list of amino
                                        # acid (given by the user)
                amino.append(file[-7:-4])
                struct = pd.read_csv(file).dropna().reset_index()
                gentab, species = check_nans(struct, species, gentab, i)
                i += 1
        genDF = pd.DataFrame(data=gentab, index=species, column=amino)
        return genDF, species

    def check_nans(struct, species, gentab, i):
        """
        we check if there is any nans in the structure given by reading a csv object 
        of the format we want
        """

        newspecies = struct['species'].values.tolist()
        nans = struct['entropylocation'].values.tolist()
        if newspecies != species:  # we might have to remove some species
            common = newspecies & species
            # those two are now just a mask of what we need to remove
            species -= common
            newspecies -= common
            # we remove the species
            nans.remove(nans[newspecies.index()])
            for gene in gentab:
                gene.remove(gene[species.index()])
        gentab[:, i] = nans
        return gentab, common

    def preprocess(self, min_gene=4):
        """
        preprocess all the data by getting it
        normalized and centerized and 
        remove the species which have less than the min_sim_gene genes
        /!\ all the gentabs in genereg should now be the same size
        add a number of gene to each species in the score datastructure
        """
        totspecies = []
        score = {}
        same = []
        notsame = []

        if min_gene > len(self.genelist):
            print "your minspecie is higher than your genelist size... setting minspecies to 1"
            min_gene = 1

        if min_gene > 1:
            if not bool(self.score):  # if it exists
                print "we are now computin the score of each species...."
                print "saccharomyces_cerevisiae should be the top one"
                i = 0

                for species in self.specieslist:  # we go throught the list of lists
                    i += 1
                    print i if i % 30 == 0  # just verbosing
                    same = species & totspecies  # we compare them
                    notsame = species - same
                    for spe in same:
                        score[spe] += 1
                    for spe in notsame:  # for all the differents we update the score list
                        score.update({spe: 1})
                    totspecies = totspecies | species
            for key, value in score.iteritems():
                if value < min_gene:
                    for gene in

        scaler = prep.StandardScaler()
        for gene in self.genevect.values():
            scaler.partial_fit(gene)
        for gene in self.genevect.values():
            gene = scaler.transform(gene)

    def preprocess_gene(self, gene):
        """
        preprocess one gene by getting
        normalized and centerized it

        :param gene:  a matrix of gene codon usage per species
        """
        stdscal = prep.StandardScaler().fit(gene)
        self.scaled = stdscal.transform(gene)

    def restoretofile(filename):
        """
        sometimes in the course of the pipeline, modification are made to our dataframes of 
        genevect or genetab.
        we should be able to save them and restore our Gene object to any moment we want
        """

    def reduce_dim(self, gene, n=2, perplexity=40):
        """
        reduce the dimensionality of your gene dataset to a defined dimension 
        using the t-SNE algorithm

        :param gene:  a matrix of gene codon usage per species
                  n:  the desired dimension
        perplexity :  an optional value when you know about tsne 

        :return tsned: the reduced dataset
        """
        tsned = man.TSNE(n_components=n, perplexity=perplexity).fit_transform(self.scaled)
        return tsned

    def clusterize(self, gene, n_clusters=4):
        """
        groups the gene dataset by using the kmeans algorithm (according to a value of closeness)

        :param gene:  a matrix of gene codon usage per species
        n_clusters : the number of groups you want to create

        :return centroids: points in HighDimension corresonding to the center of each cluster
                labels   : numerical values for differentiating each values accroding to its 
                           cluster
        """
        kmean = cluster.MiniBatchKMeans(n_clusters=n_clusters)
        kmean.fit(gene)
        labels = kmean.labels_
        centroids = kmean.cluster_centers_
        return centroids, labels

    def plot_gene(self, tsnedgene, species, getimage=False,
                  labels=False, centroids=False):
        """
        use bokeh to create nice interactive plots (either on jupyter notebook or as
         html-javascripts pages)
        in these plots of your gene dataset you can have a look at your previously 
        dimensionality reduced
        gene dataset and hover over points to have a look at the species or display 
        colors and centroids according
        to the clusterize function's output

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

    def getAllWithMoreThan(self, num=4, type="Genes"):
        """Get all the Genes/Species containing more than X Species/Genes
           /!\ it will remove the surplux of species/Genes in list containing more 
           than the specified amount ( for comparison purposes you don't want to have 
           different amounts)
           It outputs a reduced list similar to GeneVect

            :param 
            :return: 

        """

    def getGeneSimilarity(self):
        """Output a similarity matrix displaying how close each species are 
        to one another according to their number of Shared genes 
        (should be similar to the phylogenetic tree)
        Also Outputs the most similar species in a list (without containing the ones
        that share no gene)

            :param 
            :return: 

        """

    def sameSame(self):
        """Output a really reduced Genevect-like list containing only
        a group of species with exactly same genes. 

            :param 
            :return: 

        """

    def batchPlot
