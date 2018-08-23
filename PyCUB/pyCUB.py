"""
Created by Jeremie KALFON
Date : 21 FEV 2018
University of Kent, ECE paris
jkobject.com

PyCUB is a Project which goal is to understand the particular dynamics of the codon usage
bias.

you can find more about pycub here ../README.md

if you find you need to add anything please contact me @jkobject or directly make pull requests.
"""
import os
import json
import zipfile
import shutil
from ftplib import FTP
import gzip
import copy
import requests
from sklearn.preprocessing import normalize
try:
    from urllib2 import urlopen as urlopen
except:
    from urllib.request import urlopen as urlopen

from sklearn.neural_network import MLPRegressor
from joblib import Parallel, delayed
import multiprocessing
from functools32 import lru_cache

from rpy2.robjects.packages import importr
from ete2 import NCBITaxa
from rpy2 import robjects
import rpy2.robjects.packages as rpackages
from rpy2.robjects.vectors import StrVector


import pandas as pd
import numpy as np

from scipy.spatial.distance import euclidean
from scipy.sparse.csgraph import dijkstra
from scipy.stats import spearmanr
from sklearn import manifold as man
from sklearn.decomposition import PCA
from sklearn.linear_model import MultiTaskLassoCV, LassoCV
from sklearn import cluster

import espece as spe
import homoset as hset
import utils
import homology as h

from bokeh.plotting import *
from bokeh.models import *
from bokeh.io import save, show
from bokeh.layouts import column
import matplotlib.pyplot as plt

from Bio import SeqIO
# Lasso, cross Validation, Principal Component Analysis, T-SNE, spearman's rho,
# djikstra alg with fbonacci's boosting, multinomial distribution,
# multivariate normal approximation to multinomial, multiprocessing,
# parallel computing, entropy, AUC_entropy, pseudo phylodistance
# by cophenetic matrix of a dendogram, cosine similarity, hyperparams grid search,
# KMeans, MiniBatchKMeans, KModes, silhouette_score, calinski_harabaz_score, akaike information criterion,
# bayesian information criterion  binary search, recurssive function,
# dynamic programming, endres distance, kulback leiber divergence, gaussian mixture clustering
# neural networks, Perceptron, DBscan
# python, js, REST, ftp, json, doxygen,gzip,
# CAI, tRNA, CUB, 3D conformation of DNA, CUF, fungi, animals, plants, GCcount, KaKs_Scores, hydrophob, synthcost, isoepoint
# HiC data, fasta, fast,

# biology, genomics, genetics, population genomics, comparative genomics, computational biology,
# bioinformatics, machine learning, statistics, statistical learning, informatics, computer science
# knowledge discovery, datascience, big data, cloud computing, scientific computing

import pdb


class PyCUB(object):
    """PyCUB is the main object of the project that allows the user to access most of the functions

        When using it, please follow the documentation and examples on notebooks thought you can
        still use it as you please and use some of the nice tricks provided here and in python

        Args:
            species: dictionary of Espece objects from the name of the species.
                (see espece.py)
            working_homoset: PyCUB.homoset object that stores a subset of the homologies
                you want to work on
            all_homoset PyCUB.homoset that stores the all the homologies
            session: str the session name you want to use (will appear in the savings for example
            _is_saved : bool trivial system only boolean
            links: dict of all the links readily available in PyCUB.
                for the project of Jeremie KALFON please use whatever datasets you may find usefull
                (you can also download from Ensembl)

            coeffgenes: np.array regressing values for each attributes
            scoregenes: the score of the regressor
            scorespecies: the score of the regressor
            coeffspecies: np.array regressing values for each attributes

    """

    links = {'yun': {
        'homology1t500.zip': 'https://www.dropbox.com/s/fmh0ljf02twn4vw/homology1t500.zip?dl=1',
        'homology501t1000.zip': 'https://www.dropbox.com/s/ld4ar5pnh0f1w1w/homology501t1000.zip?dl=1',
        'homology1001t2000.zip': 'https://www.dropbox.com/s/he39xu9c0n2jw8n/homology1001t2000.zip?dl=1',
        'homology2001t2500.zip': 'https://www.dropbox.com/s/8w73jbs3r0ugqb8/homology2001t2500.zip?dl=1',
        'homology2501t3000.zip': 'https://www.dropbox.com/s/86d23iaetw3hmzy/homology2501t3000.zip?dl=1',
        'homology3001t3500.zip': 'https://www.dropbox.com/s/mr1tefylq11l3ee/homology3001t3500.zip?dl=1',
        'first50.zip': 'https://www.dropbox.com/s/m3vob12ztfqs8gh/first50.zip?dl=1'},
        'mymeta': {
        'Amino Acid Properties README.txt': 'https://www.dropbox.com/s/3tb2j69l0acirt0/\
        Amino%20Acid%20Properties%20README.txt?dl=1',
        'Amino Acid Properties.csv':
            'https://www.dropbox.com/s/g157emzyid2qi83/Amino%20Acid%20Properties.csv?dl=1',
        'cerevisae_prot_abundance.csv':
            'https://www.dropbox.com/s/t77016m5fqzb2fc/cerevisae_prot_abundance.csv?dl=1',
        'names_with_links.csv':
            'https://www.dropbox.com/s/voj26r0onvvqvx2/names_with_links.csv?dl=1',
        'order_name461.csv':
            'https://www.dropbox.com/s/0708046ld1pcju4/order_name461.csv?dl=1',
        'Yun_Species_Context':
            'https://www.dropbox.com/s/rdse1rco04hmuwf/Yun_Species_Context.csv?dl=1',
        'homolist.json':
            'https://www.dropbox.com/s/5a3h8hps9eozd8g/homolist.json?dl=1'
    },
        'meta': {
        'fungi':
        'ftp://ftp.ensemblgenomes.org/pub/release-39/fungi/species_metadata_EnsemblFungi.json',
        'bacteria':
        'ftp://ftp.ensemblgenomes.org/pub/release-39/bacteria/species_metadata_EnsemblBacteria.json',
        'plants':
        'ftp://ftp.ensemblgenomes.org/pub/release-39/plants/species_metadata_EnsemblPlants.json'
    }
    }
    species = {}
    working_homoset = None
    all_homoset = None
    _is_saved = False
    session = None
    coeffgenes = None
    scoregenes = None
    scorespecies = None
    coeffspecies = None

    def __init__(self, species={}, _is_saved=False,
                 _is_loaded=False, working_homoset=None, all_homoset=None, session='session1'):
        """
        will initialize the object with the different values you might have from another project

        Args:
            species: dictionary of Espece objects from the name of the species.
                (see espece.py)
            working_homoset : PyCUB.homoset object that stores a subset of the homologies
                you want to work on
            all_homoset PyCUB.homoset that stores the all the homologies
            session : str the session name you want to use (will appear in the savings for example
            _is_saved : bool trivial system only boolean

        """
        self.species = species
        self.working_homoset = working_homoset
        self.all_homoset = all_homoset
        self._is_saved = _is_saved
        self._is_loaded = _is_loaded
        self.session = session
        self.homolist = None
        print "working on session: " + self.session
        if os.path.isdir('utils/save/' + session):
            print 'you already have a session here (just a warning)'
        else:
            os.mkdir('utils/save/' + session)

    # create a function to find all homologies from a species

# CpCpGpApApTpApTpApTpTpCpCpGpApApTpApTpApTpTpCpCpGpApApTpApTpApTpTpCpCpGpApApTpApTpApTpTpTpTpCpCpGpApApTpApTpApTpTp
# GbGbCbTbTbAbTbAbTbAbAbGbGbCbTbTbAbTbAbTbAbAbGbGbCbTbTbAbTbAbTbAbAbGbGbCbTbTbAbTbAbTbAbAbAbGbGbCbTbTbAbTbAbTbAbAbAb

    def getHomologylist(self, species='saccharomyces_cerevisiae', kingdom='fungi'):
        """
        A function to retrieve the homologies directly from a given species 

        (it is better to use
        one of the key species for the different kingdoms (sacharomyces, HS, Arabidopsis..))

        Args:
            specie: str the name of the specie to get the homology from
            kingdom: str the kingdom where we can find this specie
        """
        location = 'ftp.ensemblgenomes.org' if kingdom != 'vertebrate' else 'ftp.ensembl.org'
        release = 'release-40/' if kingdom != 'vertebrate' else 'release-93'
        ftp = FTP(location)
        ftp.login()
        if kingdom == 'vertebrate':
            kingdom = ''
        ftp.cwd('pub/' + release + kingdom + '/fasta/')
        data = []
        name = []
        ftp.retrlines('NLST', data.append)
        for d in data:
            if d == species:
                ftp.cwd(d)
                link = []
                ftp.cwd('cds')
                ftp.retrlines('NLST', link.append)
                with open("utils/data/temp.fa.gz", "wb") as file:
                    for i in link:
                        if i[-9:] == "all.fa.gz":
                            ftp.retrbinary("RETR " + i, file.write)
                with gzip.open("utils/data/temp.fa.gz", "rt") as handle:
                    for record in SeqIO.parse(handle, "fasta"):
                        name.append(record.name)
        self.homolist = name

    def get_data(self, From='yun', homonames=None, kingdom='fungi', sequence='cdna',
                 additional='type=orthologues', saveonfiles=False, normalized=True, setnans=False,
                 by="entropy", using="normal", tRNA=True, getCAI=True, first=20, inpar=True):
        """
        Download the data from somewhere on the web (Ensembl, Yun(with links))

        you can provide a lot of different values to scrape Ensembl's datasets
        it will compute from ensembl to retrieve a similar dataset as what yun's
        data is.

        Args:
            From: str flag 'yun' or 'ensembl':
            homonames: list[str] what particular homologies you want to scrap if 'all' and you have used the
                getHomologylist() function, will get the homologies from there
            kingdom: str same for kingdoms
            sequence: str the type of sequences you want to use
            additional: str additional information about the scrapping
            saveonfiles: bool save the unprocessed data before populating working homoset
            normalized: bool if you want the values to be normalized by the length of the codons
                    (lengths are always saved)
            setnans: bool if you want to save the nans as metadata
            by: str flag 'entropy', 'entropyLocation' (entropy location), 'frequency'
            using: str flag 'random' 'normal' 'permutation' 'full'
            inpar: bool or int for parallel computing and number of core
            tRNA: bool whether or not to compute tRNA data
            getCAI: bool flag to true to retrieve the CAI as well
            first: int the first most expressed genes to compute the CAI ref statistics

        Raises:
            AttributeError: "you can't compute codon frequency with Yun's data...", 'not the right From'
            UnboundLocalError: "you have to load the homologies first"
            AttributeError: 'not the right From'

        http://rest.ensemblgenomes.org/
        """
        if by == 'frequency':
            print "you will have a larger dimensional matrix (59D)"
        if type(inpar) is int:
            num_cores = inpar
        else:
            num_cores = -1 if inpar else 1
        if From == 'yun':
            if by is 'frequency':
                raise AttributeError("you can't compute codon frequency with Yun's data...")
            Parallel(n_jobs=8)(delayed(utils.getyun)(key, val) for key, val in
                               self.links['yun'].iteritems())
            self.load(All=False if homonames is not None else True, filename=homonames,
                      From=From, by=by, inpar=inpar)
        elif From == 'ensembl':
            if homonames == 'all' or homonames is None:
                if self.homolist is None and kingdom == 'fungi':
                    with open('utils/meta/homolist.json', "r") as f:
                        self.homolist = json.loads(f.read())
                else:
                    if self.homolist is None:
                        raise UnboundLocalError("you have to load the homologies first")
                    print "using the loaded homolist from ensembl"
            else:
                self.homolist = homonames

            self.all_homoset = hset.HomoSet(datatype=by)
            print "doing all " + str(len(self.homolist)) + " homologies"
            print ' '
            homonamelist = []
            getCAI = self.createRefCAI(first=first, kingdom=kingdom) if getCAI else None
            if bool(inpar):
                values = Parallel(n_jobs=num_cores)(delayed(utils.loadfromensembl)(
                    name, kingdom, sequence,
                    additional, saveonfiles,
                    normalized, setnans, i, by, using, getCAI) for i, name in enumerate(self.homolist))
                for i, val in enumerate(values):
                    if val is not None:
                        homonamelist.append(self.homolist[i])
                        self.all_homoset.update({self.homolist[i]: val})

            else:
                for i, name in enumerate(self.homolist):
                    homo = utils.loadfromensembl(name, kingdom, sequence,
                                                 additional, saveonfiles,
                                                 normalized, setnans, i, by, using, getCAI)
                    if homo is not None:
                        homonamelist.append(name)
                        self.all_homoset.update({name: homo})
            self.all_homoset.datatype = by
            self.all_homoset.homo_namelist = homonamelist
            # TODO: test full pipeline with frequency/entropy/entropylocation
            taxons, species = self.all_homoset.preprocessing(withtaxons=True)
            if tRNA:
                print "computing tRNA copy numbers"
            for i, spece in enumerate(species):
                espece_val = spe.Espece(name=spece, taxonid=taxons[i])
                if tRNA:
                    espece_val.get_tRNAcopy(by=by, setnans=setnans)
                self.species.update({spece: espece_val})
            self.all_homoset.loadhashomo()
        else:
            raise AttributeError('not the right From')

    def get_metadata_Ensembl(self, kingdoms):
        """
        download it and put it where it belongs in the Espece object

        parse the server https://fungi.ensembl.org/info/website/ftp/index.html
        will also get the metadata from the kingdoms that you are analysing

        Args:
            kingdoms: str flag the type of kingdoms you wanna have 'fungi' 'bacteria' 'plants' 'animals'
        """
        if not os.path.exists('utils/meta'):
            os.mkdir('utils/meta')
        url = self.links['meta'][kingdoms]
        print "downloading " + kingdoms + " with urllib"
        if not os.path.exists('utils/meta/' + kingdoms + '.json'):
            f = urlopen(url)
            data = f.read()
            with open('utils/meta/' + kingdoms + '.json', "wb") as code:
                code.write(data)
            print "downloaded"

    def get_mymetadata(self, From='jerem', inpar=True):
        """
        Go ahead and design your own metadata retrieval here. 

        obviously you woud need to change some other functions.
        for me it is mean protein abundances in cerevisiae cells.

        Args:
            From: str flag designer of the function to load metadatas
            inpar: bool for parallel processing
        """
        if not os.path.exists('utils/meta'):
            os.mkdir('utils/meta')
        if From == 'jerem':
            num_cores = -1 if inpar else 1
            Parallel(n_jobs=num_cores)(delayed(utils.mymeta)(key, val) for key, val in
                                       self.links['mymeta'].iteritems())

        if From == 'tobias':
            self.import_metadataTobias()

    def import_metadataTobias(self):
        """
        will import the metadata obtained from tobias for the fungi species affiliated to cerevisiae to each species for further diagnostics.

        Populates metadata[num_genes, plant_pathogen, animal_pathogen, genome_size, plant_symbiotic, brown_rot, white_rot]
        for each species
        and weight, mRNA_abundance, is_secreted, protein_abundance, cys_elements, decay_rate for each homology

        Args:
            None
        """
        # species metadata
        data = pd.read_csv("utils/meta/Yun_Species_Context.csv")
        for i, species in enumerate(data["Genome"]):
            if species in self.species:
                if self.species[species].metadata is None:
                    self.species[species].metadata = {
                        "isplant_pathogen": False,
                        "isanimal_pathogen": False,
                        "isplant_symbiotic": False,  # endophyte or mycorrhizal
                        "isbrown_rot": False,
                        "iswhite_rot": False
                    }
                self.species[species].num_genes = int(data["No_Genes"][i])
                self.species[species].metadata["isplant_pathogen"] = bool(data["plant_pathogen"][i])
                self.species[species].metadata["isanimal_pathogen"] = bool(data["animal_pathogen"][i])
                self.species[species].genome_size = int(data["Genome_Size"][i])
                self.species[species].metadata["isplant_symbiotic"] = bool(data["mycorrhizal"][i] or data["endophyte"][i])
                self.species[species].metadata["isbrown_rot"] = bool(data["brown_rot"][i])
                self.species[species].metadata["iswhite_rot"] = bool(data["white_rot"][i])
        # protein metadata
        data = pd.read_csv("utils/meta/protdata/tob_currated.csv")
        for i, homo in enumerate(data["ORF"]):
            if unicode(homo) in self.all_homoset.keys():
                self.all_homoset[homo].weight = data["Molecular Weight (Da)"][i]
                self.all_homoset[homo].protein_abundance = data["Protein Abundance (molecules per cell)"][i]
                self.all_homoset[homo].mRNA_abundance = float(data["mRNA Abundance (molecules per cell)"][i].replace(',', '.'))\
                    if type(data["mRNA Abundance (molecules per cell)"][i]) is str else data["mRNA Abundance (molecules per cell)"][i]
                self.all_homoset[homo].decay_rate = float(data["Protein decay rate (min-1)"][i].replace(',', '.'))\
                    if type(data["Protein decay rate (min-1)"][i]) is str else data["Protein decay rate (min-1)"][i]
        data = pd.read_csv("utils/meta/protdata/PDI_substrates.csv")
        for i, homo in enumerate(data["ORF"]):
            if unicode(homo) in self.all_homoset.keys():
                self.all_homoset[homo].is_secreted = True
                self.all_homoset[homo].protein_abundance = data["proteins_per_cell"][i]
                self.all_homoset[homo].cys_elements = data["Cys"][i]
                self.all_homoset[homo].decay_rate = data["degradation"][i]

# LOADINGS AND SAVINGS

    def load(self, session=None, All=False, filename='first500', From=None, by='entropy', tRNA=True, inpar=True):
        """
        Get the data that is already present on a filename

        Either load from Yun's datasets or from an already saved session.
        Is being called by get_data. But you can call it to just use one of Yun's files
        as well

        Args:
            From: str if this flag is set to 'yun' it means that the filename is made of Yundata
                in which case we will create directly the homology map in the same time as the rest
                of the PyCUB object.
            All: bool set to true if load everything from Yun
            by: str same flag as get_data (for Yun's files here).
            filename: str the particular filename when not loading them all
            session: str if a session name is provided, then will load a zip file from
                this session's folder

        Returns:
            May return additionals if loading from a session where one decided to save more than the two All/working
            homologies. to be handled separately
        """
        separation = "homology"
        folder = "utils/data/" if session is None else "utils/save/" + session + '/'
        if not self._is_loaded:
            if From == 'yun' and session is None:
                if All:
                    # if we want them all we will take the first and then use loadmore
                    # function
                    filename = 'homology1t500'
                # then we process it according to how Yun Displays its data, trying to fill in
                # as much as we can
                homo_namelist = []
                nameB = ""
                self.all_homoset = hset.HomoSet(datatype=by)
                if not os.path.exists(folder + filename):
                    # we the file has not been zipped already
                    zip_ref = zipfile.ZipFile(folder + filename + '.zip', 'r')
                    zip_ref.extractall(folder + filename)
                    zip_ref.close()
                    note = False
                else:
                    note = True
                file = folder + filename
                if len(os.listdir(file)) == 2:
                    # the case with macos and special zipping... lame
                    filename = file + '/' + filename
                else:
                    filename = file

                # getting all the homology names
                print "Reviewing all the " + str(len(os.listdir(filename))) + " files"
                for f in sorted(os.listdir(filename)):
                    if f.endswith(".txt"):
                        nameA = f.split(separation)[0]
                        if(nameA != nameB):
                            nameB = nameA
                            homo_namelist.append(nameB)
                self.all_homoset.homo_namelist = homo_namelist

                # getting all the species names and instanciating the species object
                df, dflink = utils.retrievenames()
                i = 0
                df = df.sort_values(by='name')
                for _, row in df.iterrows():
                    espece_val = spe.Espece(name=row['name'],
                                            link=dflink.loc[dflink['name'] == row['name'], 'b'].tolist()[0])
                    if tRNA:
                        espece_val.get_tRNAcopy(by=by)
                    self.species.update({
                        row['name']: espece_val})
                    utils.speciestable.update({i: row['name']})
                    i += 1
                self.all_homoset.species_namelist = df['name'].tolist()

                # getting the homologies now
                dou = 0
                if by == "entropy":
                    by = "entropyValue"
                if inpar:
                    values = Parallel(n_jobs=-1)(delayed(utils.homoyun)(
                        separation, filename,
                        homology, by=by) for homology in homo_namelist)
                    for i, val in enumerate(values):
                        self.all_homoset.update({homo: h.homology(
                            full=val[0], names=val[1].tolist(),
                            nans=val[2],
                            homocode=homo,
                            lenmat=val[3], doub=val[4])})
                        dou += np.count_nonzero(val[4])
                else:
                    for homo in homo_namelist:
                        val = utils.homoyun(separation, filename, homo, by=by)
                        self.all_homoset.update({homo: h.homology(
                            full=val[0], names=val[1].tolist(),
                            nans=val[2],
                            homocode=homo,
                            lenmat=val[3], doub=val[4])})
                        dou += np.count_nonzero(val[4])
                # create the hashomomatrix
                self.all_homoset.preprocessing(withnames=True)
                self.all_homoset.loadhashomo()
                self.all_homoset.datatype = by
                print "you had " + str(dou) + " same species homologies"
                print "reviewed " + str(len(homo_namelist)) + " homologies "

                # if we haven't change the working with processing
                self._is_saved = False
                if not note:
                    shutil.rmtree(file)
                self._is_loaded = True
                if All:
                    self.loadmore('homology4501t5000', by=by)
                    self.loadmore('homology3501t4000', by=by)
                    self.loadmore('homology2501t3000', by=by)
                    self.loadmore('homology601t1000', by=by)
                    self.loadmore('homology4001t4500', by=by)
                    self.loadmore('homology3001t3500', by=by)
                    self.loadmore('homology2001t2500', by=by)
                    self.loadmore('homology1001t2000', by=by)
                print "you now have " + str(np.count_nonzero(self.all_homoset.hashomo_matrix)) +\
                    " genes in total"
            elif From is None:
                if not os.path.isfile(folder + filename + '.json'):
                    print "unzipping " + folder + filename + '.json.gz'
                    os.system("gzip -d " + folder + filename + '.json.gz')
                with open(folder + filename + ".json", "r") as f:
                    print "loading from " + filename
                    additionals = self._undictify(json.loads(f.read()))
                print "it worked !"
                os.system("gzip " + folder + filename + '.json')
                print "you now have " + str(np.count_nonzero(self.all_homoset.hashomo_matrix)) +\
                    " genes in total"
                return additionals

        else:
            print "hey, it looks like this object has already loaded some things"
            print "please use loadmore or use another object"
            print "you can delete this one with 'del' "

    def loadmore(self, filename='first500', by='entropyLocation'):
        """
        Get the data that is already present on a filename when you already have data

        is usefull to load more of Yun's datasets.
        is called when load is set to All


        Args:
            From: if this flag is set to 'yun' it means that the filename is made of Yundata
                in which case we will create directly the homology map in the same time as the rest
                of the PyCUB object. Here it is the only available option.
            filename: str the filename to additionaly load
            by: flag same as before

        Raises:
            UnboundLocalError: "You should try load first, this object is empty"

        """
        folder = "utils/data/"
        separation = "homology"
        if self._is_loaded:
            # then we process it according to how Yun Displays its data, trying to fill in
            # as much as we can
            homo_namelist = []
            nameB = ""
            if not os.path.exists(folder + filename):
                print "unzipping " + filename
                zip_ref = zipfile.ZipFile(folder + filename + '.zip', 'r')
                zip_ref.extractall(folder + filename)
                zip_ref.close()
                note = False
            else:
                note = True
            file = folder + filename
            if len(os.listdir(file)) < 3:
                file = file + '/' + filename
            print "Reviewing all the " + str(len(os.listdir(file))) + " files"
            for f in sorted(os.listdir(file)):
                if f.endswith(".txt"):
                    nameA = f.split(separation)[0]
                    if(nameA != nameB):
                        nameB = nameA
                        homo_namelist.append(nameB)

            # comparing two lists
            notdup = [item for item in homo_namelist if not (item in self.all_homoset.homo_namelist)]
            dup = len(homo_namelist) - len(notdup)
            if dup != 0:
                print "there is " + str(dup) + " duplicate from previous loads.. not cool"
                homo_namelist = notdup
            # update homonamelist
            self.all_homoset.homo_namelist.extend(homo_namelist)
            dou = 0
            print "start the iteration process, I hope you haven't clusterized \
            your data yet..else it won't work (for now)"
            num_cores = multiprocessing.cpu_count()
            values = Parallel(n_jobs=num_cores)(delayed(utils.homoyun)(
                i, homology, separation, file, self.all_homoset.species_namelist, by=by) for i,
                homology in enumerate(homo_namelist))
            for val in values:
                self.all_homoset.update(val[0])
                self.all_homoset.hashomo_matrix = np.vstack((self.all_homoset.hashomo_matrix, val[1]))\
                    if self.all_homoset.hashomo_matrix is not None else val[1]
                dou += val[2]
            print "you had " + str(dou) + " same species homologies (it can't be processed! for now)"
            print "reviewed " + str(len(homo_namelist)) + " homologies "
            self._is_saved = False
            if not note:
                shutil.rmtree(folder + filename)
            self._is_loaded = True
            print "you now have " + str(np.count_nonzero(self.all_homoset.hashomo_matrix)) +\
                " genes in total"
        else:
            raise UnboundLocalError("You should try load first, this object is empty")

    def save(self, name, save_workspace=True, save_homo=True, add_homosets={}, cmdlinetozip="gzip"):
        """
        call to save your work. you should call save on specific data structure if this is what you want to save.

        Will call other object's save, will transform all the variable into dict and save the dicts
        as json files. will save the df also as json files. PyCUB and homoset have
        their own json file.
        adding some params because else the object may be too big

        Args:
            name: str the name of the particular save on this session
            save_workspace: bool to fale not to save working_homoset
            save_homo: bool to false not to save all_homoset
            cmdlinetozip: str you need to tell the platform how to zip on your system uses gzip by default
                but it needs to be installed

        """
        filename = "utils/save/" + self.session + '/' + name + ".json"
        print "writing in " + name
        dictify = self._dictify(save_workspace, save_homo, add_homosets)
        data = json.dumps(dictify, indent=4, separators=(',', ': '))
        dirname = os.path.dirname(filename)
        if not os.path.exists(dirname):
            os.makedirs(dirname)
        with open(filename, 'w') as f:
            f.write(data)
            print "it worked !"
        # now we zip to save 90% memory space
        if cmdlinetozip == 'mac':
            os.system("ditto -c -k --sequesterRsrc " + filename + ' ' + filename + '.zip')
            os.remove(filename)
        if cmdlinetozip == 'gzip':
            os.system("gzip " + filename)
        self._is_saved = True


# PREPROCESSINGS
# CpCpGpApApTpApTpApTpTpCpCpGpApApTpApTpApTpTpCpCpGpApApTpApTpApTpTpCpCpGpApApTpApTpApTpTpTpTpCpCpGpApApTpApTpApTpTp
# GbGbCbTbTbAbTbAbTbAbAbGbGbCbTbTbAbTbAbTbAbAbGbGbCbTbTbAbTbAbTbAbAbGbGbCbTbTbAbTbAbTbAbAbAbGbGbCbTbTbAbTbAbTbAbAbAb

    def get_working_homoset(self, clusternb=None, species=None, homologies=None, cleanhomo=None):
        """
        create a subset of all_homoset on which you would like to do further computation

        To use once you have clustered homology groups, else takes everything.
        Can also be used just to get a subset of the all homosets.


        Args:
            clusternb: int set the cluster of the group you want to get need to be between 1 and homogroupnb
            homologies: list[str] the subset as a list you want to get from all_homoset
                (can be additional to a clusternb)
            species: list[str] the subset as a list you want to get from all_homoset
                (can be additional to a clusternb)
            cleanhomo: float if the homology is only shared by less than this amount amongst the species present
                in this homoset, removes them.

        Returns:
            a HomoSet object (see homoset.py)

        Raises:
            UnboundLocalError: "you have not clusterized you 'all_homoset'. you want to just use 'order_from_matrix' on it."

        """
        leng = len(self.all_homoset.clusters)
        if clusternb is not None and leng > 0:
            # if clustering has been done
            if leng == len(self.all_homoset.species_namelist):
                # version by species
                homo = hset.HomoSet()
                ind = np.argwhere(np.asarray(self.all_homoset.clusters) == clusternb - 1)[:, 0]
                species_name = [self.all_homoset.species_namelist[i] for i in ind]
                homo.homodict = dict(self.all_homoset.homodict)
                homo.homodict.remove(species_name)
                homo.homo_namelist = self.all_homoset.homo_namelist
                homo.species_namelist = species_name

            elif leng == len(self.all_homoset.homodict):
                # version by homologies
                homo = hset.HomoSet()
                ind = np.argwhere(np.asarray(self.all_homoset.clusters) == clusternb - 1)[:, 0]
                if cleanhomo is not None:
                    perc = self.all_homoset.hashomo_matrix.sum(1).astype(float) / self.all_homoset.hashomo_matrix.shape[1]
                    homo_name = [self.all_homoset.homo_namelist[i] for i in ind if perc[i] > cleanhomo]
                else:
                    homo_name = [self.all_homoset.homo_namelist[i] for i in ind]
                for x in homo_name:
                    homo.update({x: self.all_homoset.homodict[x]})
                homo.homo_namelist = homo_name
                homo.species_namelist = self.all_homoset.species_namelist

        else:
            raise UnboundLocalError("you have not clusterized you 'all_homoset'. you want to just use" +
                                    "'order_from_matrix' on it.")
            return False
        homo.datatype = self.all_homoset.datatype
        if homologies is not None:
            homo.homodict = {k: homo.homodict[k] for k in homologies}
        if species is not None:
            other = [item for item in homo.species_namelist if item not in species]
            homo.remove(sepcies=other)
        homo.loadhashomo()
        # big mistake to call this as it removes also species from the all homology as they both share the same objects
        # if cleanspecies is not None:
        #    homo.clean_species(thresh=cleanspecies)
        #    homo.loadhashomo()
        homo.loadfullhomo()
        homo.datatype = self.all_homoset.datatype
        self.working_homoset = homo
        return homo

    def get_subset(self, homoset, withcopy=False, clusternb=None, species=None, homologies=None):
        """
        either changes or returns a subset of the provided homoset

        To use once if you want to further refine a set of homologies


        Args:
            homoset: PyCUB.homoset to get a subset from
            withcopy: bool to true if we don't want to change the homoset object but create a copy from it
            clusternb: int set the cluster of the group you want to get need to be between 1 and homogroupnb
            homologies: list[str] the subset as a list or a tuple of int
            species: list[the subset as a list, or a list of int

        Returns:
            a HomoSet object (see homoset.py)
        """
        if withcopy:
            homo = copy.deepcopy(homoset)
        else:
            homo = homoset
        if homologies is not None:
            homo.homodict = {k: homoset[k] for k in homologies} if type(homologies[0]) is str else \
                {homoset.homo_namelist[k]: homoset[k] for k in homologies}
            homo.homo_namelist = homo.homodict.keys()
        homo.hashomo_matrix = None
        homo.homo_matrix = None
        homo.homo_matrixnames = None
        homo.fulleng = None
        homo.red_homomatrix = None
        homo.homo_clusters = None
        homo.averagehomo_matrix = None
        homo.stats = {}
        if species is not None:
            other = [item for item in homoset.species_namelist if not (item in species)]
            homo.remove(sepcies=other)
        homo.loadhashomo()
        homo.loadfullhomo()
        return homo

    def get_full_genomes(self, kingdom='fungi', seq='cds', avg=True, by="entropy", normalized=False):
        """
        go trought all full genome fasta files in the ftp server of ensemblgenomes and

        download then parse them to get the full entropy of the genome.
        usefull for futher comparison steps.
        will populate the fullentropy, fullvarentropy, fullGCcount,
        varGCcount of each species where the full sequence is known

        Args:
            kingdom: str flags the relevant kingdom of you current session [fungi,plants,bacteria, animals]
            seq: str flags the type of sequence you consider the full genome is (coding or non coding or full) [cds, all, cda]
            avg: bool to true if we average over each gene or get the full dna in one go.
            by: str flags what type of computation should be done [entropy,frequency]
            normalized: should we normalize the entorpy by length

        Raises:
            MemoryError: "this sequence is too long to be computed (> 1 billion bp)"

        """
        def _compute_full_entropy(handle, by='entropy', avg=True, normalized=False, setnans=False):
            """
            called by get full genomes, either calls utils.computeYun or process the full coding sequence in one go
            Private method. see 'get_full_genomes()'
            """
            GCcount = []
            val = []
            if avg:
                for record in SeqIO.parse(handle, "fasta"):
                    codseq = [record.seq._data[i:i + 3] for i in range(0, len(record.seq._data), 3)]
                    valH, _, _ = utils.computeyun(codseq, setnans=setnans, normalized=normalized, by=by)
                    val.append(valH)
                    GCcount.append(float(record.seq._data.count('G') + record.seq._data.count('C')) / len(record.seq._data))
                return np.array(val), np.array(GCcount)
            else:
                print "not working, overflow ..."
                """
                pdb.set_trace()
                c = []
                amino = list(utils.amino)
                codons = dict(utils.codons)
                GCcount = 0.
                for x, record in enumerate(SeqIO.parse(handle, "fasta")):
                    if len(c) > 1000000000:
                        raise MemoryError("this sequence is too long to be computed (> 1 billion bp)")
                    GCcount += (record.seq._data.count('G') + record.seq._data.count('C'))
                    for i in range(0, len(record.seq._data), 3):
                        c.append(record.seq._data[i:i + 3])
                valH = np.zeros(len(amino)) if by != 'frequency' else np.zeros(59)  # the number of codons usefull
                utils.CUBD = len(amino) if by != 'frequency' else 59
                pos = 0
                GCcount = float(GCcount) / len(c)
                for k, amin in enumerate(amino):
                    subcodons = codons[amin]
                    nbcod = len(subcodons)  # replace Cleng
                    count = np.zeros(nbcod)
                    X = np.zeros(nbcod)
                    mn = np.ones(nbcod) / nbcod
                    for j, val in enumerate(c):
                        for i, cod in enumerate(codons[amin]):
                            if val == cod:
                                count[i] += 1
                                c.pop(j)
                                break
                    lengsubseq = count.sum()  # replace subSlength
                    if by == 'frequency':
                        E = count / lengsubseq
                        valH[pos:pos + nbcod] = E
                        pos += nbcod
                    elif by == "entropy":
                        Yg = multinomial.pmf(x=count, n=lengsubseq, p=mn)
                        # efor part
                        div, i = divmod(lengsubseq, nbcod)
                        X[:int(i)] = np.ceil(div) + 1
                        X[int(i):] = np.floor(div)
                        Eg = multinomial.pmf(x=X, n=lengsubseq, p=mn)
                        # end here
                        valH[k] = -np.log(Yg / Eg) / lengsubseq if normalized else -np.log(Yg / Eg)
                print "missed codons: "+str(len(c))
                return valH, GCcount
                """
        location = 'ftp.ensemblgenomes.org' if kingdom != 'vertebrate' else 'ftp.ensembl.org'
        release = 'release-40/' if kingdom != 'vertebrate' else 'release-93'
        ftp = FTP(location)
        ftp.login()
        if kingdom == 'vertebrate':
            kingdom = ''
        ftp.cwd('pub/' + release + kingdom + '/fasta/')
        data = []
        ftp.retrlines('NLST', data.append)
        species_namelist = self.species.keys()
        for d in data:
            ftp.cwd(d)
            if d[-10:] == 'collection':
                subdata = []
                ftp.retrlines('NLST', subdata.append)
                for sub in subdata:
                    if sub in species_namelist:
                        link = []
                        ftp.cwd(sub + '/' + seq)
                        ftp.retrlines('NLST', link.append)
                        with open("utils/data/temp.fa.gz", "wb") as file:
                            for i in link:
                                if i[-9:] == "all.fa.gz":
                                    ftp.retrbinary("RETR " + i, file.write)
                        with gzip.open("utils/data/temp.fa.gz", "rt") as handle:
                            val, gccount = _compute_full_entropy(handle, by, avg)
                        self.species[sub].fullentropy = val.mean(1) if avg else val
                        self.species[sub].fullvarentropy = (val.var(1)**(0.5)).mean() if avg else None
                        self.species[sub].fullGCcount = gccount.mean() if avg else gccount
                        self.species[sub].varGCcount = gccount.var()**(0.5) if avg else None
                        ftp.cwd('../..')
                        os.remove("utils/data/temp.fa.gz")

            else:
                if d in species_namelist:
                    link = []
                    ftp.cwd(seq)
                    ftp.retrlines('NLST', link.append)
                    with open("utils/data/temp.fa.gz", "wb") as file:
                        for i in link:
                            if i[-9:] == "all.fa.gz":
                                ftp.retrbinary("RETR " + i, file.write)
                                break
                    with gzip.open("utils/data/temp.fa.gz", "rt") as handle:
                        val, gccount = _compute_full_entropy(handle, by, avg)
                        pdb.set_trace()
                    self.species[d].fullentropy = val.mean(1) if avg else val
                    self.species[d].fullvarentropy = (val.var(1)**(0.5)).mean() if avg else None
                    self.species[d].fullGCcount = gccount.mean() if avg else gccount
                    self.species[d].varGCcount = gccount.var()**(0.5) if avg else None
                    ftp.cwd('..')
                    os.remove("utils/data/temp.fa.gz")
            ftp.cwd('..')

    def get_taxons(self):
        """
        find the taxons of each referenced species (see PyCUB.Espece.gettaxons())

        Args:
            None
        """
        for key, val in self.species.iteritems():
            try:
                val.gettaxons()
            except:
                print key + " has no referenced taxons"
        print "got taxons"

    def get_evolutionary_distance(self, display_tree=False, size=40):
        """
        uses metadata of the ancestry tree and computes a theoretical evolutionary distance matrix between each species

        can optionaly take any hierarchical evolutionary file between a group of species
        will populate utils.phylo_distances with a pandas.df of the phylodistance
        and meandist with the average distance amongst species in the df, species are referenced
        by their taxon ids. you have to have taxons in your species. will also plot the distance matrix

        Args:
            display_tree: bool to true to print the phylogenetic tree as a txt
                (may be quite big)
            size: int the x size of the plot

        Raises:
            EnvironmentError: "you need to have R installed to compute the distance"
        """

        ncbi = NCBITaxa()
        taxons = []
        for key, val in self.species.iteritems():
            if val.taxonid is not None and val.taxonid != '':
                taxons.append(val.taxonid)
        tree = ncbi.get_topology(taxons)  # taxons
        # finding what this tree looks like
        if display_tree:
            print tree.get_ascii(attributes=["sci_name", "rank"])

        with open('utils/meta/metaphylo/temp_tree.phy', 'w') as f:  # maybe will be newick format...
            f.write(tree.write())
        # """
        try:
            # https://stackoverflow.com/questions/19894365/running-r-script-from-python
            base = importr('base')
            utiles = importr('utils')
        except:
            print EnvironmentError("you need to have R installed to compute the distance")
            return
        if not rpackages.isinstalled('treeio'):
            robjects.r('''
                source("https://bioconductor.org/biocLite.R")
                biocLite("treeio")
                ''')
        robjects.r('''
            treeText <- readLines("utils/meta/metaphylo/temp_tree.phy")
            treeText <- paste0(treeText, collapse="")
            library(treeio)
            tree <- read.tree(text = treeText) ## load tree
            distMat <- cophenetic(tree)
            write.table(distMat,"utils/meta/metaphylo/phylodistMat_temp.csv")
        ''')
        # """
        df = pd.read_csv("utils/meta/metaphylo/phylodistMat_temp.csv", delim_whitespace=True)
        dcol = {}
        dind = {}
        for name, species in self.species.iteritems():
            if species.taxonid:
                dind.update({int(species.taxonid): name})
                dcol.update({unicode(species.taxonid): name})
        df = df.rename(index=dind, columns=dcol)
        utils.phylo_distances = df
        utils.meandist = df.sum().sum() / (len(df)**2 - len(df))
        self.plot_distances(size=size)

    def createRefCAI(self, speciestocompare='saccharomyces_cerevisiae', kingdom='fungi', first=20):
        """
        do a compute CAI

        where we get Tobias' data to find highly expressed genes and
        use them to compute codon frequency for the reference set and use it to compute
        the CAI and mean CAI for each ho
        mology.

        Args:
            speciestocompare: str the name of the species to retrieve the genes from
            kingdom: str the kingdom where we can find it
            first: the number of highly expressed genes to retrieve

        """
        if kingdom != 'fungi':
            print "if kingdom is not fungi, need to provide another file"
            return
        data = pd.read_csv("utils/meta/protdata/tob_currated.csv")
        homonames = data["ORF"].values
        expres = data["Protein Abundance (molecules per cell)"].values
        ind = expres.argsort()
        highlyexpressed = [homonames[i] for i in ind[:first]]

        location = 'ftp.ensemblgenomes.org' if kingdom != 'vertebrate' else 'ftp.ensembl.org'
        release = 'release-40/' if kingdom != 'vertebrate' else 'release-93'
        ftp = FTP(location)
        ftp.login()
        if kingdom == 'vertebrate':
            kingdom = ''
        ftp.cwd('pub/' + release + kingdom + '/fasta/')
        data = []
        ftp.retrlines('NLST', data.append)
        for d in data:
            if d == speciestocompare:
                ftp.cwd(d)
                link = []
                ftp.cwd('cds')
                ftp.retrlines('NLST', link.append)
                with open("utils/data/temp.fa.gz", "wb") as file:
                    for i in link:
                        if i[-9:] == "all.fa.gz":
                            ftp.retrbinary("RETR " + i, file.write)
                codseq = []
                with gzip.open("utils/data/temp.fa.gz", "rt") as handle:
                    for record in SeqIO.parse(handle, "fasta"):
                        if record.id in highlyexpressed:
                            codseq.extend([record.seq._data[i:i + 3] for i in range(0, len(record.seq._data), 3)])
                os.remove("utils/data/temp.fa.gz")
        return utils.reference_index(codseq, forcai=True)

    def speciestable(self):
        """
        a copy of the utils.speciestable

        Args:
            None

        Returns:
            a copy of the utils.speciestable (dict[int,str] of species to their PyCUB coded value
        """
        return dict(utils.speciestable) if utils.speciestable is not None else None

    def phylo_distances(self):
        """
        a copy of the phylodistances dataframe see (get_evolutionary_distance())

        Args:
            None

        Returns:
            a copy of the phylodistances dataframe see (get_evolutionary_distance())
        """
        return utils.phylo_distances.copy() if utils.phylo_distances is not None else None

    def compute_averages(self, homoset):
        """
        compute the average entropy

        Will add species related averages gotten from this homoset in the species
        container, and in the homoset everytime you compute averages from a set, they will get erased !

        Args:
            homoset: PyCUB.homoset from which to compute the averages

        """
        if homoset.homo_matrix is None:
            homoset.loadfullhomo()
        if homoset.hashomo_matrix is None:
            homoset.loadhashomo()
        # we get all CUB values pertaining to one specific species from
        # the homoset with the full homo matrix
        _, counts = np.unique(homoset.homo_matrixnames, return_counts=True)
        ind = homoset.homo_matrixnames.argsort()
        GCmat = np.zeros((len(homoset.homo_namelist), len(homoset.species_namelist)))
        for i, val in enumerate(homoset.homo_namelist):
            GCmat[i, homoset[val].names] = homoset[val].GCcount
        GCmat = GCmat.sum(0) / np.count_nonzero(GCmat, 0)
        for i, spece in enumerate(homoset.species_namelist):
            self.species[spece].meanGChomo = GCmat[i]
        pos = 0
        speciestable = dict(utils.speciestable)
        for i, un in enumerate(counts):
            aslicetoavg = homoset.homo_matrix[ind[pos:pos + un]]
            bslicetoavg = homoset.fulleng[ind[pos:pos + un]]

            self.species[speciestable[i]].average_entropy = aslicetoavg.mean(axis=0)
            self.species[speciestable[i]].average_size = bslicetoavg.sum(0).mean()
            # variances are mean variances over all values
            self.species[speciestable[i]].var_entropy = (aslicetoavg.var(axis=0)**(0.5)).mean()
            pos += un
            self.species[speciestable[i]].tot_homologies = un
        print "homology averages : " + str(homoset.homo_matrix.mean(axis=0))
        homoset.averagehomo_matrix = np.array([homoset[homo].mean for homo in homoset.homo_namelist])
        for _, val in homoset.iteritems():
            val.compute_averages()

# CpCpGpApApTpApTpApTpTpCpCpGpApApTpApTpApTpTpCpCpGpApApTpApTpApTpTpCpCpGpApApTpApTpApTpTpTpTpCpCpGpApApTpApTpApTpTp
# GbGbCbTbTbAbTbAbTbAbAbGbGbCbTbTbAbTbAbTbAbAbGbGbCbTbTbAbTbAbTbAbAbGbGbCbTbTbAbTbAbTbAbAbAbGbGbCbTbTbAbTbAbTbAbAbAb

    def compare_species(self, showvar=True, reducer='tsne', perplexity=40, eps=0.3, size=10):
        """
        compare the species according to their mean CUB,

        plot the mean CUB
        to their full CUB, to their tRNA copy numbers, to the euclidean distance of their CUB
        to the one of their phylogenetic matrix.
        in this plot, the mean entropy value is plotted as a regular homology plot but each dot is a species
        thus we can compare them together, moreover, the size of the dots informs oneself of the variance
        in entropy per species. the color intensity informs on how much this is close to what is given by the
        tRNA values. (additional information such as name of the species, number of tRNA values,metadata and the point
        from its tRNA value is plotted when hovering the dot)
        then we can also compare the mean of homologies or else, the full entropy of the cdna sequence per species
        is also computed the euclidean distance amongst the species for the full entropy to see if a difference can be linked
        to some evolutionary for one codon

        Args:
            size: the average size of the datapoints in the pointcloud representation of this dataset
            showvar: bool to true, show the mean variance in CUB values accros this homology as a variation in dot sizes
            eps: float the hyperparamter of the clustering algorithm applied to this dataset
            homoset: PyCUB.homoset the homoset to use
            reducer: str the reducer to use 'tsne' or 'PCA'
            perplexity: int the perplexity hyperparam for tSNE


        Raises:
            UnboundLocalError: "you need to compute the averages of the all_homoset. use PyCUB.compute_averages(homoset)"
            UnboundLocalError: "you have to compute tRNA values, use PyCUB.get_tRNAcopy()"
            AttributeError: "try avg or full"

        """
        e_subspecies = np.zeros((len(self.species), utils.CUBD))
        vare_subspecies = np.zeros((len(self.species), 1))
        tRNAentropydist = np.zeros(len(self.species))
        tRNA_number = np.zeros((len(self.species), 1), dtype=int)
        genome_size = np.zeros((len(self.species), 1), dtype=int)
        efulldiff = np.zeros((len(self.species), 1))
        fullvarentropy = np.zeros((len(self.species), 1))
        varGCcount = np.zeros((len(self.species), 1))
        gcfulldiff = np.zeros((len(self.species), 1))
        num_genes = np.zeros((len(self.species), 1), dtype=int)
        meanGChomo = np.zeros((len(self.species), 1))
        suff = 0
        phylo_distances = self.phylo_distances()
        for i, (name, specie) in enumerate(self.species.iteritems()):
            if specie.average_entropy is not None:
                e_subspecies[i] = specie.average_entropy
                vare_subspecies[i] = specie.var_entropy
                meanGChomo[i] = specie.meanGChomo
            else:
                raise UnboundLocalError("you have to compute averages, use PyCUB.compute_averages(homoset)")
            if specie.fullentropy is not None:
                efulldiff[i] = euclidean(specie.average_entropy, specie.fullentropy)
                gcfulldiff[i] = euclidean(specie.meanGChomo, specie.fullGCcount)

            if specie.fullvarentropy is not None:
                varGCcount[i] = specie.varGCcount
                fullvarentropy[i] = specie.fullvarentropy
            if specie.copynumbers is not None:
                if specie.tRNAentropy is not None:
                    if to == 'full' and specie.fullentropy is not None:
                        tRNAentropydist[i] = euclidean(specie.tRNAentropy, specie.fullentropy)
                    elif to == 'avg' or specie.fullentropy is None:
                        tRNAentropydist[i] = euclidean(specie.tRNAentropy, specie.average_entropy)
                    else:
                        raise AttributeError("try avg or full")
                    suff += 1
                    # if is zero, will be black colored
                if specie.copynumbers.get("datapoints", False):
                    tRNA_number[i] = specie.copynumbers["datapoints"]
                    # if is zero, will be black colored
            else:
                raise UnboundLocalError("you have to compute tRNA values, use PyCUB.get_tRNAcopy()")
            if specie.genome_size:
                genome_size[i] = specie.genome_size
            if specie.num_genes:
                num_genes[i] = specie.num_genes
            # if is zero, will be black colored
        print "we have " + str(suff) + " species with sufficient statistics in their tRNA values"
        if reducer == 'tsne':
            red = man.TSNE(n_components=2, perplexity=perplexity).fit_transform(e_subspecies)
        elif reducer == 'pca':
            red = PCA(n_components=2).fit_transform(e_subspecies)
        alg = cluster.DBSCAN(eps=eps, min_samples=6, algorithm='auto', n_jobs=-1)
        clusters = alg.fit_predict(e_subspecies).tolist()
        colormap = list(utils.colormap)
        colors = [colormap[0] if not dist else
                  utils.rgb2hex((26, 188, 10 + np.floor(246 * dist))) for dist in tRNAentropydist]
        data = dict(x=red[:, 0], y=red[:, 1],
                    species=self.species.keys(),
                    meanentropy=["%.2f" % i.mean() for i in e_subspecies],
                    color=colors,
                    vare_subspecies=vare_subspecies,
                    tRNA_number=tRNA_number,
                    recent=colors,
                    tRNAentropydist=tRNAentropydist,
                    genome_size=genome_size,
                    num_genes=num_genes,
                    efulldiff=efulldiff,
                    gcfulldiff=gcfulldiff,
                    gccount=meanGChomo,
                    varGCcount=varGCcount,
                    fullvarentropy=fullvarentropy,
                    clusters=clusters,
                    size=[(size / 2) + (size * 25 * (val.var_entropy))
                          for _, val in self.species.iteritems()] if showvar else size)
        # compare it to the phylogenetic distance matrix and the sizes
        if phylo_distances is not None:
            names = list(phylo_distances.index)
            distances = np.zeros(len(self.species))
            for i, val in enumerate(self.species.keys()):
                if val in names:
                    distances[i] = phylo_distances[val].values.mean()
            data.update({'distances': distances})
        labe = ["show clusters", "show num_genes", "show genome_size",
                "show full/mean CUB diff", "show GCcount", "show full/mean GC diff",
                "show GC variance", "show distance to tRNA UB", "show tRNA number", "show avg phylodistance",
                "show full phylo distance"]  # 11
        for i, (key, val) in enumerate(self.species[self.species.keys()[0]].metadata.iteritems()):
            labe.append("show if " + key)
            data.update({str(i + 11): [espe.metadata[key] if espe.metadata is not None else False for _, espe in self.species.iteritems()]})

        """
        to implement full phylo distance, we need to add onHover, and put in data a dict which associate
        for each species name, a species ordered list of normalized distance values from the phylodistance
        dataframe
        """

        source = ColumnDataSource(data=data)
        output_notebook()

        callback = CustomJS(args=dict(source=source), code=utils.callback_allgenes)
        radio_button_group = widgets.RadioButtonGroup(
            labels=labe, callback=callback, active=7)
        hover = HoverTool(tooltips=[("species: ", "@species"), ("mean entr: ", "@meanentropy"),
                                    ("phylodistances: ", "@distances"), ("tRNA CN: ", "@tRNA_number"),
                                    ("dist2tRNA Bias Value: ", "@tRNAentropydist"), ("GCbias: ", "@gccount")])
        p = figure(title="exploration of every homologies",
                   tools=[hover, WheelZoomTool(), PanTool(), SaveTool(), ResetTool()],
                   plot_width=800, plot_height=600)
        p.circle(x='x', y='y', source=source, color='color', size='size')
        save(column(radio_button_group, p), "utils/templot/homology_compare.html")
        show(column(radio_button_group, p))

    def compute_ages(self, homoset, preserved=True, minpreserv=0.9, minsimi=0.85):
        homoset.compute_ages(preserved=preserved, minpreserv=minpreserv, minsimi=minsimi)

    def regress_on_species(self, without=[""], full=True, onlyhomo=False, perctrain=0.8, algo="lasso",
                           eps=0.001, n_alphas=100):
        """
        Will fit a regression curve on the CUB values of the different species according to the metadatas available for each of them.

        It will try to see if there is enough information in the metadata to retrieve CUB values. and if there is,
        how much for each metadata (if we constraint the number of regressors) is it better for mean homology CUB
        or full genome CUB ?
        or raw frequency, should we remove some data?

        Args:
            without: list[str] of flags [similarity_scores, KaKs_Scores, nans, lenmat, GCcount, weight,
                protein_abundance, mRNA_abundance, decay_rate, cys_elements, tot_volume, mean_hydrophobicity,
                glucose_cost, synthesis_steps, is_recent, meanecai]
            full: bool flags to true to use full CUB values or meanCUB values,  as regressee
            homoset: PyCUB.homoset the homoset to use
            perctrain: the percentage of training set to total set ( the rest is used as test set)
            algo: str flag to lasso or nn to use either Lasso with Cross Validation, or a 2 layer  neural net
            eps: the eps value for the Lasso
            n_alphas: the number of alphas for the lasso

        Returns:
            scoregenes: float, the score of the regression performed
            coeffgenes: the coefficient applied to each category (for each CUB value if using full)
            attrlist: the corresponding list[str] of attribute used

        Raises:
            UnboundLocalError: "wrong params"
        """
        params = []
        phylo_distances = self.phylo_distances()
        espece = self.species.values()[0]
        attrlist = ["average_size", "num_genes", "genome_size", "average_size", "fullGCcount",
                    "varGCcount", "tot_homologies"]
        for attr in attrlist:
            if getattr(espece, attr) is not None and attr not in without:
                arr = np.array([getattr(spece, attr) for spece in self.species.values()]).astype(float)
                arr = arr / arr.max()
                params.append(arr)

        if espece.copynumbers is not None and "copynumbers" not in without:
            attrlist.append("copynumbers")
            arr = np.array([spece.copynumbers.get('datapoints', 0) for spece in self.species.values()]).astype(float)
            arr = arr / arr.max()
            params.append(arr)
        for key, val in espece.metadata.iteritems():
            if str(key) not in without:
                attrlist.append(key)
                params.append(np.array([spece.metadata[key] if spece.metadata is not None else False
                                        for spece in self.species.values()]).astype(int))
        if phylo_distances is not None and "phylo_distances" not in without:
            attrlist.append("phylo_distances")
            names = list(phylo_distances.index)
            arr = np.array([phylo_distances[val].values.mean() if val in names else 0 for val in self.species.keys()])
            arr = arr / arr.max()
            params.append(arr)
        dataset = np.zeros((len(self.species), utils.CUBD))
        for i, (_, v) in enumerate(self.species.iteritems()):
            dataset[i] = v.average_entropy if onlyhomo else v.fullentropy
        if not full:
            dataset = dataset.mean(1)

        if algo == "lasso":
            # http://scikit-learn.org/stable/modules/generated/sklearn.linear_model.LassoCV.html
            model = MultiTaskLassoCV(eps=eps, n_alphas=n_alphas,
                                     alphas=None, fit_intercept=True, normalize=False,
                                     max_iter=1000, tol=0.0001, copy_X=False, cv=None,
                                     verbose=False, n_jobs=1, random_state=None, selection='cyclic') \
                if full else LassoCV(eps=eps, n_alphas=n_alphas,
                                     alphas=None, fit_intercept=True, normalize=False, precompute='auto',
                                     max_iter=1000, tol=0.0001, copy_X=False, cv=None, verbose=False, n_jobs=-1,
                                     positive=False, random_state=None, selection='cyclic')
        elif algo == "nn" and not full:
            # http://scikit-learn.org/stable/modules/generated/sklearn.linear_model.Perceptron.html
            model = MLPRegressor(hidden_layer_sizes=(len(attrlist), len(attrlist)), activation='relu', solver='adam', alpha=0.0001,
                                 batch_size='auto', learning_rate='constant', learning_rate_init=0.001,
                                 power_t=0.5, max_iter=200, shuffle=True, random_state=None, tol=0.0001,
                                 verbose=1, warm_start=False, momentum=0.9, nesterovs_momentum=True,
                                 early_stopping=False, validation_fraction=0.1, beta_1=0.9, beta_2=0.999, epsilon=1e-08)
        else:
            raise UnboundLocalError("wrong params")
        pdb.set_trace()
        params = np.vstack(params).T
        model.fit(params[:int(len(self.species) * perctrain)], dataset[:int(len(self.species) * perctrain)])
        self.scorespecies = model.score(params[int(len(self.species) * perctrain):],
                                        dataset[int(len(self.species) * perctrain):], sample_weight=None)
        self.coeffspecies = model.coef_.tolist() if algo == "lasso" else model.coefs_
        print "the R^2 score is of: " + str(self.scorespecies)
        print "-------------------------------"
        if model == "lasso":
            for i, val in enumerate(attrlist):
                print val + ": " + str(self.coeffspecies[i])
        return self.scorespecies, self.coeffspecies, attrlist

    # TODO: create a model that can say if the species has a high  CUB or low, given data on the species and on the protein

    def compare_homologies(self, homoset, homosapiens=False, mindistance=10, preserved=True, size=10,
                           minpreserv=0.9, minsimi=0.9, showvar=True, eps=0.28, reducer='tsne', perplexity=40):
        """
        finds for species with a common ancester separated by a pseudo phylogenetic distance X,

        genes/functions that are novel to only a subset.
        plot the two with their differences the differences between the two
        also considers an homology as highly preserved if it is shared amongst most of the species
        and if the average similarity score is high amongst this homology.
        also shows if there is a relationship between the number of amino acids the sequence does not
        encode for and the codon usage bias
        We could have used the sequence dating of ensembl but it only works for homo sapiens for now
        maybe use it for homosapiens later

        Args:
            homosapiens: bool to true if we should use homosapiens dataset on gene dates
            mindistance: int the minimal phylogenetic distance between in average in this homology to consider it
                highly conserved
            preserved: bool to true if we should find highly preserved genes or not
            size: the average size of the datapoints in the pointcloud representation of this dataset
            minpreserv: float minimal percentage of homologous species that have this homology
            minsimi: float minimal avg similarity between genes to consider them highly preserved
            showvar: bool to true, show the mean variance in CUB values accros this homology as a variation in dot sizes
            eps: float the hyperparamter of the clustering algorithm applied to this dataset
            homoset: PyCUB.homoset the homoset to use
            reducer: str the reducer to use 'tsne' or 'PCA'
            perplexity: int the perplexity hyperparam for tSNE


        Raises:
            UnboundLocalError: "you need to compute the averages of the all_homoset. use PyCUB.compute_averages(homoset)"
        """
        if not homosapiens:
            if homoset[-1].isrecent is None:
                homoset.compute_ages(preserved=preserved, minpreserv=minpreserv, minsimi=minsimi)
        else:
            pass
            # TODO: code the version for homo sapiens where we know exactly this distance with more data
            # and better inference metrics

        # display the differences between recent homologies and older ones
        if homoset.averagehomo_matrix is None:
            raise UnboundLocalError("you need to compute the averages of homoset. use PyCUB.compute_averages(homoset)")
        if reducer == 'tsne':
            red = man.TSNE(n_components=2, perplexity=perplexity).fit_transform(homoset.averagehomo_matrix)
        elif reducer == 'pca':
            red = PCA(n_components=2).fit_transform(homoset.averagehomo_matrix)
        else:
            raise AttributeError("wrong algorithm")
        alg = cluster.DBSCAN(eps=eps, min_samples=7, algorithm='auto', n_jobs=-1)
        clusters = alg.fit_predict(homoset.averagehomo_matrix).tolist()
        n_clusters_ = len(set(clusters))
        if n_clusters_ > 10:
            print "ooups you have more than 10 clusters"
        colormap = list(utils.colormap)
        colors = [colormap[0] if not homoset[homo].isrecent else
                  utils.rgb2hex((26, 188, np.floor(156 * homoset[homo].isrecent))) for homo in homoset.homo_namelist]
        data = dict(x=red[:, 0], y=red[:, 1],
                    homologies=homoset.homo_namelist,
                    meanentropy=["%.2f" % homoset.averagehomo_matrix[i].mean()
                                 for i in range(len(homoset.averagehomo_matrix))],
                    color=colors,
                    recent=colors,
                    clusters=clusters,
                    size=[(size / 2) + (size * 25 * self.all_homoset[homo].var.mean()) if showvar else size for homo in self.all_homoset.homo_namelist])
        # add average of similar protein name
        values = ["similarity_scores", "KaKs_Scores", "nans", "lenmat", "GCcount", "weight",
                  "protein_abundance", "mRNA_abundance", "decay_rate", "is_secreted", "cys_elements",
                  "tot_volume", "mean_hydrophobicity", "glucose_cost", "synthesis_steps", "isoelectricpoint", "meanecai", "meancai", "conservation"]

        labe = ["show Recent/preserved", "showclusters", "show avg similarity_scores", "show avg KaKs_Scores", "show Nans avg",
                "show avg Length", "show avg GCcount", "Show weight", "Show prot abundance", "Show mRNA abundance",
                "Show half life", "Show secreted", "Show num of cys", "Show volume", "Show hydrophobicity", "show cost (glucose)",
                "Show synthesis cost", "Show Pi", "Show ECAI", "Show CAI", "show amino Conservation"]  # 21
        templabe = labe[:2]
        i = 2
        for val in values[:5]:
            if getattr(homoset[0], val) is not None:
                data.update({val: [getattr(homoset[homo], val).mean() for homo in homoset.homo_namelist]})
                templabe.append(labe[i])
            else:
                templabe.append(" ")
            i += 1
        for val in values[5:]:
            if getattr(homoset[0], val) is not None:
                data.update({val: [getattr(homoset[homo], val) for homo in homoset.homo_namelist]})
                templabe.append(labe[i])
            else:
                templabe.append(" ")
            i += 1

        source = ColumnDataSource(data=data)
        output_notebook()
        callback = CustomJS(args=dict(source=source), code=utils.callback_allhomo)
        radio_button_group = widgets.RadioButtonGroup(
            labels=templabe, callback=callback, active=0)
        hover = HoverTool(tooltips=[("homologies: ", "@homologies"), ("avg nans: ", "@nans"), ("similarity scores: ", "@similarity_scores"),
                                    ("mRNA abundance: ", "@mRNA_abundance"), ("mean ecai: ", "@meanecai"), ("amino conservation: ", "@conservation"),
                                    ("mean_entr: ", "@meanentropy"), ("length: ", "@lengths"), ("GCcount: ", "@gc")])
        p = figure(title="exploration of every homologies",
                   tools=[hover, WheelZoomTool(), PanTool(), SaveTool(), ResetTool()],
                   plot_width=800, plot_height=600)
        p.circle(x='x', y='y', source=source, color='color',
                 size='size')
        save(column(radio_button_group, p), "utils/templot/homology_compare.html")
        show(column(radio_button_group, p))

    def regress_on_genes(self, homoset, full=True, without=['meanecai', 'meancai'], perctrain=0.8, algo="lasso", eps=0.001, n_alphas=100):
        """
        Will fit a regression curve on the CUB values of the different homologies according to the metadatas available for each of them.

        It will try to see if there is enough information in the metadata to retrieve CUB values. and if there is,
        how much for each metadata (if we constraint the number of regressors) is it better for entropy values, mean entropy
        or ECAI values
        or raw frequency, should we remove some data

        Args:
            without: list[str] of flags [similarity_scores, KaKs_Scores, nans, lenmat, GCcount, weight,
                protein_abundance, mRNA_abundance, decay_rate, cys_elements, tot_volume, mean_hydrophobicity,
                glucose_cost, synthesis_steps, is_recent, meanecai]
            full: bool flags to true to use full CUB values or meanCUB values,  as regressee
            homoset: PyCUB.homoset the homoset to use
            perctrain: the percentage of training set to total set ( the rest is used as test set)
            algo: str flag to lasso or nn to use either Lasso with Cross Validation, or a 2 layer  neural net
            eps: the eps value for the Lasso
            n_alphas: the number of alphas for the lasso

        Returns:
            scoregenes: float, the score of the regression performed
            coeffgenes: the coefficient applied to each category (for each CUB value if using full)
            attrlist: the corresponding list[str] of attribute used

        Raises:
            UnboundLocalError: "wrong params"
        """
        params = []
        values = ["similarity_scores", "KaKs_Scores", "nans", "lenmat", "GCcount", "weight",
                  "protein_abundance", "mRNA_abundance", "decay_rate", "is_secreted", "cys_elements",
                  "tot_volume", "mean_hydrophobicity", "glucose_cost", "synthesis_steps",
                  "isoelectricpoint", "meanecai", "meancai", "conservation"]
        attrlist = []
        for val in values[:5]:
            if val not in without:
                if getattr(homoset[0], val) is not None:
                    arr = np.nan_to_num(np.array([getattr(homoset[homo], val).mean() for homo in homoset.homo_namelist])).astype(float)
                    arr = arr / arr.max()
                    params.append(arr)
                    attrlist.append(val)
        pdb.set_trace()
        for val in values[5:]:
            if val not in without:
                if getattr(homoset[0], val) is not None:
                    arr = np.nan_to_num(np.array([getattr(homoset[homo], val) for homo in homoset.homo_namelist])).astype(float)
                    arr = arr / arr.max()
                    params.append(arr)
                    attrlist.append(val)
        dataset = homoset.averagehomo_matrix if full else homoset.averagehomo_matrix.mean(1)
        if algo == "lasso":
            # http://scikit-learn.org/stable/modules/generated/sklearn.linear_model.LassoCV.html
            model = MultiTaskLassoCV(eps=eps, n_alphas=n_alphas,
                                     alphas=None, fit_intercept=True, normalize=False,
                                     max_iter=1000, tol=0.0001, copy_X=False, cv=None,
                                     verbose=False, n_jobs=1, random_state=None, selection='cyclic')\
                if full else LassoCV(eps=eps, n_alphas=n_alphas,
                                     alphas=None, fit_intercept=True, normalize=False, precompute='auto',
                                     max_iter=1000, tol=0.0001, copy_X=False, cv=None, verbose=False, n_jobs=-1,
                                     positive=False, random_state=None, selection='cyclic')
        elif algo == "nn" and not full:
            # http://scikit-learn.org/stable/modules/generated/sklearn.linear_model.Perceptron.html
            model = MLPRegressor(hidden_layer_sizes=(len(attrlist), len(attrlist)), activation='relu', solver='adam', alpha=0.0001,
                                 batch_size='auto', learning_rate='constant', learning_rate_init=0.001,
                                 power_t=0.5, max_iter=200, shuffle=True, random_state=None, tol=0.0001,
                                 verbose=1, warm_start=False, momentum=0.9, nesterovs_momentum=True,
                                 early_stopping=False, validation_fraction=0.1, beta_1=0.9, beta_2=0.999, epsilon=1e-08)
        else:
            raise UnboundLocalError("wrong params")
        params = np.vstack(params).T
        model.fit(params[:int(len(homoset.homo_namelist) * perctrain)], dataset[:int(len(homoset.homo_namelist) * perctrain)])
        self.scoregenes = model.score(params[int(len(homoset.homo_namelist) * perctrain):],
                                      dataset[int(len(homoset.homo_namelist) * perctrain):], sample_weight=None)
        self.coeffgenes = model.coef_.tolist() if algo == "lasso" else model.coefs_
        print "the R^2 score is of: " + str(self.scoregenes)
        print "-------------------------------"
        if model == "lasso":
            for i, val in enumerate(attrlist):
                print val + ": " + str(self.coeffgenes[i])
        return self.scoregenes, self.coeffgenes, attrlist

    def getRelation2G3DD(self, species_name='saccharomyces_cerevisiae', kingdom='fungi',
                         intrachromosome="utils/meta/3Dmodel/interactions_HindIII_fdr0.01_intra_cerevisiae.csv",
                         interchromose=["utils/meta/3Dmodel/cerevisiae_inter1.csv",
                                        "utils/meta/3Dmodel/cerevisiae_inter2.csv",
                                        "utils/meta/3Dmodel/cerevisiae_inter3.csv",
                                        "utils/meta/3Dmodel/cerevisiae_inter4.csv",
                                        "utils/meta/3Dmodel/cerevisiae_inter5.csv"], bins=2000, seq='cds', use='diament2', compute="jerem"):
        """
        https://www.nature.com/articles/ncomms6876

        retrieve the data for the species sacharomyces cerevisiae and Schizosaccharomyces pombe
        and find if similarity distances of CUB using entropy between genes of this species is predictive
        of closeness of genes in the nucleus.

        Used to confirm a work on nature and see if we can have some similar results by only looking at the
        CUB

        Args:
            species_name: str the name of the species to look for
            kingdom: str the kingdom in which to find the species
            intrachromosome: str the location of the csv interaction data for intrachromosome respecting the format
                of the default file
            interchromose: str the location of the csv interaction data for interchromose respecting the format
                of the default file
            bins: int, the number of bin to use (a power of 2)
            seq: the type of sequence to compare to. (to compute the CUB from)
        """
        # get gene distance matrix from entropy value distance or Andres Schindelin metrics
        # compare to see how much the distance between one can explain the distance between another by
        # regression
        # retrieve the data.
        intra = pd.read_csv(intrachromosome, delim_whitespace=True).drop(columns=["qvalue", "freq"])
        inter = pd.concat([pd.read_csv(interchro) for interchro in interchromose]).drop(columns=[
            "P value", "Q value", "sequence frequency"])
        # getting all the genes
        torname = {"HindIII fragment": "locus1",
                   "HindIII fragment.1": "locus2",
                   "chromosome": "chr1",
                   "chromosome.1": "chr2"}
        inter = inter.rename(torname, axis="columns")
        df = pd.concat([intra, inter])
        df = df.sort_values(by=["chr1", "locus1"])
        df = df.reset_index()
        chrom2int = {
            'I': 1,
            'II': 2,
            'III': 3,
            'IV': 4,
            'V': 5,
            'VI': 6,
            'VII': 7,
            'VIII': 8,
            'IX': 9,
            'X': 10,
            'XI': 11,
            'XII': 12,
            'XIII': 13,
            'XIV': 14,
            'XV': 15,
            'XVI': 16,
            'XVII': 17,
            'XVIII': 18,
            'XIX': 19,
            'XX': 20,
            'XXI': 21
        }

        location = 'ftp.ensemblgenomes.org' if kingdom != 'vertebrate' else 'ftp.ensembl.org'
        release = 'release-40/' if kingdom != 'vertebrate' else 'release-93'
        ftp = FTP(location)
        ftp.login()
        if kingdom == 'vertebrate':
            kingdom = ''
        ftp.cwd('pub/' + release + kingdom + '/fasta/')
        data = []
        ftp.retrlines('NLST', data.append)
        pdb.set_trace()
        for d in data:
            if d in species_name:
                ftp.cwd(d)
                link = []
                ftp.cwd(seq)
                ftp.retrlines('NLST', link.append)
                with open("utils/data/temp.fa.gz", "wb") as file:
                    for i in link:
                        if i[-9:] == "all.fa.gz":
                            ftp.retrbinary("RETR " + i, file.write)
                vals = []
                cufs = []
                cubs = []
                names = []
                positions = []
                nb = 0
                codons = list(utils.codamino.keys())
                with gzip.open("utils/data/temp.fa.gz", "rt") as handle:
                    for record in SeqIO.parse(handle, "fasta"):
                        uncounted = len(record.seq._data) - (record.seq._data.count('A') + record.seq._data.count('T') +
                                                             record.seq._data.count('C') + record.seq._data.count('G'))
                        if uncounted:
                            print "ref uncounted = " + str(uncounted)
                            codseq = record.seq._data.replace("Y", "T").replace("R", "G").replace("K", "G")\
                                .replace("M", "A").replace('S', 'C').replace("W", "A").replace("B", "C").replace("D", "T")\
                                .replace("H", "T").replace("V", "G").replace("N", "C")
                        else:
                            codseq = record.seq._data
                        codseq = [codseq[i:i + 3] for i in range(0, len(codseq), 3)]
                        leng = len(codseq)
                        CuF = np.zeros(64)
                        for i, val in enumerate(codons):
                            CuF[i] = codseq.count(val)
                        nb += 1
                        valH, CuB, _, _ = utils.computeyun(codseq, setnans=False, normalized=False,
                                                           by="entropy" + "frequency")
                        server = "http://rest.ensemblgenomes.org" if kingdom != 'vertebrate' else "http://rest.ensembl.org"
                        names.append(record.id)
                        vals.append(valH)
                        cufs.append((CuF / leng).tolist())
                        cubs.append(CuB)
                        ext = "/lookup/id/" + record.id + "?expand=1"
                        r = requests.get(server + ext, headers={"Content-Type": "application/json"})
                        if not r.ok:
                            r.raise_for_status()
                            sys.exit()
                        print '{0} genes\r'.format(nb),
                        decoded = r.json()
                        chrom = chrom2int.get(str(decoded["seq_region_name"]), False)
                        if not chrom:
                            names.pop()
                            vals.pop()
                            cufs.pop()
                            cubs.pop()
                            continue
                        # we associate valH to a position retrieve the positions
                        positions.append([chrom, (decoded["start"] + decoded["end"]) / 2])
                ind = sorted(range(len(positions)), key=lambda k: positions[k])
                vals = np.array(vals)[ind]
                cufs = np.array(cufs)[ind]
                cubs = np.array(cubs)[ind]
                positions = np.array(positions)[ind]
                names = [names[i] for i in ind]

                tempchrom = 0
                tempind = 0
                tempdf = df.loc[df['chr1'] == 1]
                if compute == 'jerem':
                    gene2pos = {}
                    pos2gene = {}
                    dists = np.zeros(len(positions))
                    tempdist = 1000000
                    for n, val in enumerate(positions):
                        if val[0] >= tempchrom + 1:
                            tempdf = df.loc[df['chr1'] == val[0]]
                            tempchrom = val[0]
                            tempind = tempdf.index[0]
                            maxind = tempdf.index[-1]
                        # Here we could use a modified binar search instead
                        while abs(tempdf['locus1'][tempind] - val[1]) <= tempdist:
                            tempdist = abs(tempdf['locus1'][tempind] - val[1])
                            tempind += 1
                            if tempind >= maxind:
                                break
                        dists[n] = tempdist
                        tempind -= 1
                        # we found a position
                        tempdist = 10000000
                        gene2pos.update({n: [tempchrom, tempdf['locus1'][tempind]]})
                        if pos2gene.get((tempchrom, tempdf['locus1'][tempind]), False):
                            pos2gene[(tempchrom, tempdf['locus1'][tempind])].append(n)
                        else:
                            pos2gene.update({(tempchrom, tempdf['locus1'][tempind]): [n]})
                    # for each gene positions we look at the closest point in the contact map (list of positison)
                    # we create a mapping dict for that.
                    tempchrom = 0
                    print "average distance is" + str(dists.mean())
                    missedrelation = 0
                    tempdf = df.loc[df['chr1'] == 1]
                    dist3D = np.zeros((len(positions), len(positions)), dtype=int)
                    dist3D += 1000000
                    np.fill_diagonal(dist3D, 0)
                    # Doing the first one for efficiency
                    tempdf = df.loc[df['chr1'] == 1]
                    relatedto = tempdf.loc[tempdf['locus1'] == gene2pos[0][1]]
                    chro = list(relatedto["chr2"])
                    loc = list(relatedto["locus2"])
                    for i in range(len(chro)):
                        pos = pos2gene.get((chro[i], int(loc[i])), False)
                        if pos:
                            for p in pos:
                                dist3D[p, 0] = 1
                                dist3D[0, p] = 1
                        else:
                            missedrelation += 1
                    n = 1
                    for val in positions[1:]:
                        dist3D[n, n - 1] = 1
                        dist3D[n - 1, n] = 1
                        if val[0] >= tempchrom + 1:
                            tempdf = df.loc[df['chr1'] == val[0]]
                        relatedto = tempdf.loc[tempdf['locus1'] == gene2pos[n][1]]
                        chro = list(relatedto["chr2"])
                        loc = list(relatedto["locus2"])
                        for i in range(len(chro)):
                            pos = pos2gene.get((chro[i], int(loc[i])), False)
                            if pos:
                                for p in pos:
                                    dist3D[p, n] = 1
                                    dist3D[n, p] = 1
                            else:
                                missedrelation += 1
                        n += 1
                    print "got " + str(missedrelation) + " missed relation"
                else:
                    gene2pos = {}
                    pos2gene = {}
                    dists = np.zeros(len(positions))
                    tempdist = 1000000
                    for n, val in enumerate(positions):
                        if val[0] >= tempchrom + 1:
                            tempdf = df.loc[df['chr1'] == val[0]]
                            tempchrom = val[0]
                            tempind = tempdf.index[0]
                            maxind = tempdf.index[-1]
                        # Here we could use a modified binar search instead
                        while abs(tempdf['locus1'][tempind] - val[1]) <= tempdist:
                            tempdist = abs(tempdf['locus1'][tempind] - val[1])
                            tempind += 1
                            if tempind >= maxind:
                                break
                        dists[n] = tempdist
                        tempind -= 1
                        # we found a position
                        tempdist = 10000000
                        gene2pos.update({n: [tempchrom, tempdf['locus1'][tempind]]})
                        if pos2gene.get((tempchrom, tempdf['locus1'][tempind]), False):
                            pos2gene[(tempchrom, tempdf['locus1'][tempind])].append(n)
                        else:
                            pos2gene.update({(tempchrom, tempdf['locus1'][tempind]): [n]})
                    # for each gene positions we look at the closest point in the contact map (list of positison)
                    # we create a mapping dict for that.
                    tempchrom = 0
                    print "average distance is" + str(dists.mean())
                    missedrelation = 0
                    tempdf = df.loc[df['chr1'] == 1]
                    dist3D = np.zeros((len(positions), len(positions)), dtype=int)
                    dist3D += 1000000
                    np.fill_diagonal(dist3D, 0)
                    # Doing the first one for efficiency
                    tempdf = df.loc[df['chr1'] == 1]
                    relatedto = tempdf.loc[tempdf['locus1'] == gene2pos[0][1]]
                    chro = list(relatedto["chr2"])
                    loc = list(relatedto["locus2"])
                    for i in range(len(chro)):
                        pos = pos2gene.get((chro[i], int(loc[i])), False)
                        if pos:
                            for p in pos:
                                dist3D[p, 0] = 1
                                dist3D[0, p] = 1
                        else:
                            missedrelation += 1
                    n = 1
                    for val in positions[1:]:
                        dist3D[n, n - 1] = 1
                        dist3D[n - 1, n] = 1
                        if val[0] >= tempchrom + 1:
                            tempdf = df.loc[df['chr1'] == val[0]]
                        relatedto = tempdf.loc[tempdf['locus1'] == gene2pos[n][1]]
                        chro = list(relatedto["chr2"])
                        loc = list(relatedto["locus2"])
                        for i in range(len(chro)):
                            pos = pos2gene.get((chro[i], int(loc[i])), False)
                            if pos:
                                for p in pos:
                                    dist3D[p, n] = 1
                                    dist3D[n, p] = 1
                            else:
                                missedrelation += 1
                        n += 1
                    print "got " + str(missedrelation) + " missed relation"
                shortestpath = dijkstra(dist3D, directed=False)  # to set shortest as 1

                # for each genes, we look if there is a contact gene in the contact map with the mapping
                # if the genes are the ones next to each other, we mark them as close to each other
                # true also if on different chrom
                # then we compute 1hop distances on this matrix
                #
                # if tempchrom + 1 == val['chr1']
                # val['chr1']
                # for all values
                # we take n subsets for which we compute the average CUB,
                # then for each we compute a CUB distance value as a big distance
                # matrix of size n x n n should be 32 000
                # we create another distance matrix using an Andres distance metrics on the CUF values
                if use != "jerem1":
                    distcuf = np.zeros((len(positions), len(positions)), dtype=float)
                    distcub = np.zeros((len(positions), len(positions)), dtype=float)
                    distent = np.zeros((len(positions), len(positions)), dtype=float)
                    vals = np.ma.masked_equal(vals, 0)
                    cufs = np.ma.masked_equal(cufs, 0)
                    cubs = np.ma.masked_equal(cubs, 0)
                    for val in vals:
                        i = 0
                        for comp in vals:
                            if i < j:
                                distent[j, i] = distent[i, j]
                                distcub[j, i] = distcub[i, j]
                                distcuf[j, i] = distcuf[i, j]
                            elif i > j:
                                distent[j, i] = euclidean(val, comp)
                                distcub[j, i] = utils.endresdistance(cubs[j], cubs[i])
                                distcuf[j, i] = utils.endresdistance(cufs[j], cufs[i])
                            i += 1
                        j += 1
                        print '\rdistcomputation ' + str(j) + 'over ' + str(len(vals)),

                    if use == "diament1":
                        div, i = divmod(len(names), bins)
                        X = np.zeros(bins)
                        cuf = np.zeros((bins, bins))
                        cub = np.zeros((bins, bins))
                        ent = np.zeros((bins, bins))
                        dist3D = np.zeros((bins, bins), dtype=float)
                        X[:int(i)] = np.ceil(div) + 1
                        X[int(i):] = np.floor(div)
                        start = 0
                        for num, val in enumerate(X):
                            star = 0
                            for nu, va in enumerate(X):
                                dist3D[num, nu] = shortestpath[start:start + int(val), star:star + int(va)].mean()
                                cuf[num, nu] = distcuf[start:start + int(val), star:star + int(va)].mean()
                                cub[num, nu] = distcub[start:start + int(val), star:star + int(va)].mean()
                                ent[num, nu] = distent[start:start + int(val), star:star + int(va)].mean()
                                star += int(va)
                            start += int(val)
                            print '\rbinninng ' + str(num),
                        del distcuf
                        del distcub
                        del distent
                        del shortestpath
                        self.rho_ent, self.pent = spearmanr(ent, dist3D, axis=None)
                        self.rho_cub, self.pcub = spearmanr(cub, dist3D, axis=None)
                        self.rho_cuf, self.pcuf = spearmanr(cuf, dist3D, axis=None)
                    elif use == "diament2":
                        div, i = divmod(len(names)**2, bins)
                        X = np.zeros(bins)
                        cuf = np.zeros(bins)
                        cub = np.zeros(bins)
                        ent = np.zeros(bins)
                        dist3D = np.zeros(bins)
                        X[:int(i)] = np.ceil(div) + 1
                        X[int(i):] = np.floor(div)

                        distent = np.ravel(distent)
                        shortestpath = np.ravel(shortestpath)
                        distcub = np.ravel(distcub)
                        distcuf = np.ravel(distcuf)
                        # for CUF
                        ind = np.argsort(distent)
                        distent = distent[ind]
                        sortshortestpath = shortestpath[ind]
                        start = 0
                        for i, val in enumerate(X):
                            dist3D[i] = sortshortestpath[start:start + int(val)].mean()
                            ent[i] = distent[start:start + int(val)].mean()
                            start += int(val)
                        self.rho_ent, self.pent = spearmanr(ent, dist3D)
                        # for CUB
                        ind = np.argsort(distcuf)
                        distcuf = distcuf[ind]
                        sortshortestpath = shortestpath[ind]
                        start = 0
                        for i, val in enumerate(X):
                            dist3D[i] = sortshortestpath[start:start + int(val)].mean()
                            cuf[i] = distcuf[start:start + int(val)].mean()
                            start += int(val)
                        self.rho_cuf, self.pcuf = spearmanr(cuf, dist3D)
                        # for CUF
                        ind = np.argsort(distcub)
                        distcub = distcub[ind]
                        sortshortestpath = shortestpath[ind]
                        start = 0
                        for i, val in enumerate(X):
                            dist3D[i] = sortshortestpath[start:start + int(val)].mean()
                            cub[i] = distcub[start:start + int(val)].mean()
                            start += int(val)
                        self.rho_cub, self.pcub = spearmanr(cub, dist3D)

                else:
                    # computation jerem 1st
                    div, i = divmod(len(names), bins)
                    X = np.zeros(bins)
                    cuf = np.zeros((bins, 64))
                    cub = np.zeros((bins, 59))
                    ent = np.zeros((bins, 18))
                    dist3D = np.zeros((bins, bins), dtype=float)
                    X[:int(i)] = np.ceil(div) + 1
                    X[int(i):] = np.floor(div)
                    start = 0
                    for num, val in enumerate(X):
                        ent[num] = vals[start:start + int(val)].mean(0)
                        cuf[num] = cufs[start:start + int(val)].mean(0)
                        cub[num] = cubs[start:start + int(val)].mean(0)
                        star = 0
                        for nu, va in enumerate(X):
                            dist3D[num, nu] = shortestpath[start:start + int(val), star:star + int(va)].mean()
                            star += int(va)
                        start += int(val)
                        print '\rbinninng ' + str(num),
                    del vals
                    del cufs
                    del shortestpath
                    distcuf = np.zeros((bins, bins), dtype=float)
                    distcub = np.zeros((bins, bins), dtype=float)
                    distent = np.zeros((bins, bins), dtype=float)
                    ent = np.ma.masked_equal(ent, 0)
                    cuf = np.ma.masked_equal(cuf, 0)
                    cub = np.ma.masked_equal(cub, 0)
                    j = 0
                    pdb.set_trace()
                    for val in ent:
                        i = 0
                        for comp in ent:
                            if i < j:
                                distent[j, i] = distent[i, j]
                                distcub[j, i] = distcub[i, j]
                                distcuf[j, i] = distcuf[i, j]
                            elif i > j:
                                distent[j, i] = euclidean(val, comp)
                                distcub[j, i] = utils.endresdistance(cub[j], cub[i])
                                distcuf[j, i] = utils.endresdistance(cuf[j], cuf[i])

                            i += 1
                        j += 1
                        print '\rdistcomputation ' + str(j),

                    self.rho_ent, self.pent = spearmanr(distent, dist3D, axis=None)
                    self.rho_cub, self.pcub = spearmanr(distcub, dist3D, axis=None)
                    self.rho_cuf, self.pcuf = spearmanr(distcuf, dist3D, axis=None)
                return self.rho_ent, self.pent, self.rho_cuf, self.pcuf, self.rho_cub, self.pcub

        def search(value, arr):
            """
            modified binary search for closest value search in a sorted array

            David Soroko @https://stackoverflow.com/questions/30245166

            Args:
                value: a value to look for 
                arr: array like, an array to skim 

            Returns:
                the closest value present in the array
            """
            if value < arr[0]:
                return arr[0]
            if value > arr[-1]:
                return arr[-1]
            lo = 0
            hi = len(arr) - 1
            while lo <= hi:
                mid = (hi + lo) / 2
                if value < arr[mid]:
                    hi = mid - 1
                elif value > arr[mid]:
                    lo = mid + 1
                else:
                    return arr[mid]
            return arr[lo] if arr[lo] - value < value - arr[hi] else arr[hi]

        # we then do a correlation to the two matrices and see if they are highly correlated
        # we will use the spearman's rho for each bins

    def plot_distances(self, size=40):
        """
        plot the phylogenetic distance matrix

        Args:
            size: int the x size of the plot

        Raises:
            UnboundLocalError: "compute the phylo distance matrix first (look at the doc)"
        """
        if utils.phylo_distances is not None:
            plt.figure(figsize=(size, 200))
            plt.title('evolutionary distances')
            plt.imshow(1. / (1 + np.array(utils.phylo_distances)))
            plt.savefig("utils/templot/evolutionarydistances.pdf")
            plt.show()
        else:
            raise UnboundLocalError("compute the phylo distance matrix first (look at the doc)")
# SPECIAL FUNCTION

    def _dictify(self, save_workspace, save_homo, add_homosets):
        """
        Used by the saving function. transform the workspace object into a dictionary that can be json serializable

        adding some params because else the object may be too big

        Args:
            save_workspace: bool to save working_homoset
            save_homo: bool to save all_homoset
            ass_homosets: PyCUB.homoset instances to add to this dict

        Return:
            A dict holding every element to be jsonized
        """
        dictispecies = {}
        for key, val in self.species.iteritems():
            dictispecies.update({key: val._dictify()})
        di = {"species": dictispecies,
              "all_homoset": self.all_homoset._dictify(savehomos=True) if save_homo and
              (self.all_homoset is not None) else None,
              "working_homoset": self.working_homoset._dictify() if save_workspace and
              (self.working_homoset is not None) else None
              }
        for key, val in add_homosets:
            di.update({key: val})
        return di

    def _undictify(self, data):
        """
        same function but to retransform everything

        Here we don't use other classes undictify functions but we just recreate them by passing it
        to their init methods which is clearer.

        Args:
            data: dict to undictify into the workspace object

        Returns:
            Other PyCUB.homosets that would have been saved as add_homosets
        """
        self.species = {}
        self._is_saved = False
        ret = {}
        if data["all_homoset"] is not None:
            self.all_homoset = hset.HomoSet(data=data["all_homoset"])
        for key, val in data.iteritems():
            if key == 'species':
                for ke, va in val.iteritems():
                    self.species.update({ke: spe.Espece(data=va)})
            elif key == "working_homoset":
                if val is not None:
                    self.working_homoset = hset.HomoSet(data=val)
                    for v in self.working_homoset.homo_namelist:  # as we don't save the homologies twice here
                        self.working_homoset.update({v: self.all_homoset[v]})
            elif key == "all_homoset":
                pass
            else:
                ret.update({key: val})
        return ret

    def loadspeciestable(self):
        """
        short function to retrieve the speciestable from Disk

        Args:
            None
        Raises:
            IOError: "no speciestable file"
        """
        filename = "utils/meta/savings/speciestable.json"
        if not os.path.isfile(filename):
            raise IOError("no speciestable file")
        with open(filename, "r") as f:
            speciestable = json.loads(f.read())
            print "it worked !"
        utils.speciestable = {}
        for key, val in speciestable.iteritems():
            utils.speciestable.update({int(key): str(val)})

    def savespeciestable(self):
        """
        short function to put the speciestable on Disk

        This is done since there may be some memory leakage, probably due to some autoreloading behavior
        of the global data stored on utils.

        Args:
            None
        """
        filename = "utils/meta/savings/speciestable.json"
        data = json.dumps(dict(utils.speciestable), indent=4, separators=(',', ': '))
        dirname = os.path.dirname(filename)
        if not os.path.exists(dirname):
            os.makedirs(dirname)
        with open(filename, 'w') as f:
            f.write(data)
            print "it worked !"
# TODO: [in the end] check that everything is imported, create a req file, check everything is saved,
# retry imports on a new machine, export doc in html and latex, redo the readme (the one for the data
# as well),
# create the file for loading in pip, put it in pip. Create a small medium article. create an html index page
# referencing the project, the documentation, the html plots, the dissertation, the article.


"""
    def offload(self, var='species'):
        ""
        offload when your file is too big (easier than using HDF or something).. for now
        save and

        ""
        print "offload when your file is too big (easier than using HDF or something).. for now"

        if var == 'species'
            dictispecies = {}
            for key, val in self.species.iteritems():
                dictispecies.update({key: val._dictify()})
            filename = "utils/save/" + self.session + '/' + species + ".json"
            print "writing in " + name
            dictify = self._dictify(save_workspace, save_homo)
            data = json.dumps(dictify, indent=4, separators=(',', ': '))
            dirname = os.path.dirname(filename)
            if not os.path.exists(dirname):
                os.makedirs(dirname)
            with open(filename, 'w') as f:
                f.write(data)
                print "it worked !"

            self.species = 'off'
            self.all_homoset = 'off'
        else:
            print "offload only when you have saved your object and got a working homoset"
"""
