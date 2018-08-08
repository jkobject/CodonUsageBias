""""
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
try:
    from urllib2 import urlopen as urlopen
except:
    from urllib.request import urlopen as urlopen

from joblib import Parallel, delayed
import multiprocessing

from rpy2.robjects.packages import importr
from ete2 import NCBITaxa
from rpy2 import robjects
import rpy2.robjects.packages as rpackages
from rpy2.robjects.vectors import StrVector


import pandas as pd
import numpy as np
from scipy.stats import multinomial
from scipy.spatial.distance import euclidean
import espece as spe
import homoset as hset
import utils
import homology as h

from bokeh.plotting import *
from bokeh.models import *
from bokeh.io import save, show
from bokeh.layouts import column

import matplotlib.pyplot as plt
from sklearn import manifold as man
from sklearn.decomposition import PCA
from sklearn.linear_model import Lasso
from Bio import SeqIO


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
        print "working on session: " + self.session
        if os.path.isdir('utils/save/' + session):
            print 'you already have a session here (just a warning)'
        else:
            os.mkdir('utils/save/' + session)

    # create a function to find all homologies from a species

# CpCpGpApApTpApTpApTpTpCpCpGpApApTpApTpApTpTpCpCpGpApApTpApTpApTpTpCpCpGpApApTpApTpApTpTpTpTpCpCpGpApApTpApTpApTpTp
# GbGbCbTbTbAbTbAbTbAbAbGbGbCbTbTbAbTbAbTbAbAbGbGbCbTbTbAbTbAbTbAbAbGbGbCbTbTbAbTbAbTbAbAbAbGbGbCbTbTbAbTbAbTbAbAbAb

    def get_data(self, From='yun', homonames=None, kingdom='compara=fungi', sequence='cdna',
                 additional='type=orthologues', saveonfiles=False, normalized=True, setnans=False,
                 by="entropy", using="normal", tRNA=True, inpar=True):
        """
        Download the data from somewhere on the web (Ensembl, Yun(with links))

        you can provide a lot of different values to scrape Ensembl's datasets
        it will compute from ensembl to retrieve a similar dataset as what yun's
        data is.

        Args:
            From: str flag 'yun' or 'ensembl':
            homonames: list[str] what particular homologies you want to scrap
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

        Raises:
            AttributeError: "you can't compute codon frequency with Yun's data...", 'not the right From'


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
                with open('utils/meta/homolist.json', "r") as f:
                    homonames = json.loads(f.read())
            self.all_homoset = hset.HomoSet(datatype=by)
            print "doing all " + str(len(homonames)) + " homologies"
            homonamelist = []
            if bool(inpar):
                values = Parallel(n_jobs=num_cores)(delayed(utils.loadfromensembl)(
                    name, kingdom, sequence,
                    additional, saveonfiles,
                    normalized, setnans, i, by, using) for i, name in enumerate(homonames))
                for i, val in enumerate(values):
                    if val is not None:
                        homonamelist.append(homonames[i])
                        self.all_homoset.update({homonames[i]: val})

            else:
                for i, name in enumerate(homonames):
                    homo = utils.loadfromensembl(name, kingdom, sequence,
                                                 additional, saveonfiles,
                                                 normalized, setnans, i, by, using)
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
        will import the metadata obtained from tobias for the fungi species affiliated to
        cerevisiae to each species for further diagnostics.

        Populates metadata[num_genes, plant_pathogen, animal_pathogen, genome_size, plant_symbiotic, brown_rot, white_rot]
        for each species
        and weight, mRNA_abundance, is_secreted, protein_abundance, cys_elements, decay_rate for each homology
        """
        # species metadata
        data = pd.read_csv("utils/meta/Yun_Species_Context.csv")
        for i, species in enumerate(data["Genome"]):
            if species in self.species:
                if self.species[species].metadata is None:
                    self.species[species].metadata = {
                            "isplant_pathogen": None,
                            "isanimal_pathogen": None,
                            "isplant_symbiotic": None,  # endophyte or mycorrhizal
                            "isbrown_rot": None,
                            "iswhite_rot": None
                        }
                self.species[species].num_genes = data["No_Genes"][i]
                self.species[species].metadata["plant_pathogen"] = data["plant_pathogen"][i]
                self.species[species].metadata["animal_pathogen"] = data["animal_pathogen"][i]
                self.species[species].genome_size = data["Genome_Size"][i]
                self.species[species].metadata["plant_symbiotic"] = data["mycorrhizal"][i] or data["endophyte"][i]
                self.species[species].metadata["brown_rot"] = data["brown_rot"][i]
                self.species[species].metadata["white_rot"] = data["white_rot"][i]
        # protein metadata
        data = pd.read_csv("utils/meta/protdata/tob_currated.csv")
        for i, homo in enumerate(data["ORF"]):
            if unicode(homo) in self.all_homoset.keys():
                self.all_homoset[homo].weight = data["Molecular Weight (Da)"][i]
                self.all_homoset[homo].protein_abundance = data["Protein Abundance (molecules per cell)"][i]
                self.all_homoset[homo].mRNA_abundance = data["mRNA Abundance (molecules per cell)"][i]
                self.all_homoset[homo].decay_rate = data["Protein decay rate (min-1)"][i]
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
        call to save your work. you should call save on specific data structure if this is what you
        want to save.
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
            homo.homodict = {k: homoset.homodict[k] for k in homologies}
            homo.homo_namelist = homologies
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

        will populate the fullentropy, fullvarentropy, fullGCcount, varGCcount of each species where the full sequence is known

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
            uncounted = 0
            if avg:
                for record in SeqIO.parse(handle, "fasta"):
                    valH, _, _, uncount = utils.computeyun(record.seq._data, setnans=setnans, normalized=normalized, by=by)
                    val.append(valH)
                    uncounted += uncount
                    GCcount.append(float(record.seq._data.count('G') + record.seq._data.count('C')) / len(record.seq._data))
                print "missed codons: "+str(uncounted)
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

        ftp = FTP('ftp.ensemblgenomes.org')
        ftp.login()
        ftp.cwd('pub/release-40/' + kingdom + '/fasta/')
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
                        self.species[sub].fullvarentropy = val.var(1).mean() if avg else None
                        self.species[sub].fullGCcount = gccount.mean() if avg else gccount
                        self.species[sub].varGCcount = gccount.var() if avg else None
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
                    self.species[d].fullvarentropy = val.var(1).mean() if avg else None
                    self.species[d].fullGCcount = gccount.mean() if avg else gccount
                    self.species[d].varGCcount = gccount.var() if avg else None
                    ftp.cwd('..')
                    os.remove("utils/data/temp.fa.gz")
            ftp.cwd('..')

    def get_taxons(self):
        """
        find the taxons of each referenced species (see PyCUB.Espece.gettaxons())
        """
        for key, val in self.species.iteritems():
            try:
                val.gettaxons()
            except:
                print key + " has no referenced taxons"
        print "got taxons"

    def get_evolutionary_distance(self, display_tree=False, size=40):
        """
        uses metadata of the ancestry tree and computes a theoretical evolutionary
        distance matrix between each species

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
            if val.taxonid:
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
            utiles.install_packages(StrVector('treeio'))
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

    def speciestable(self):
        """
        a copy of the utils.speciestable

        Returns:
            a copy of the utils.speciestable (dict[int,str] of species to their PyCUB coded value
        """
        return dict(utils.speciestable) if utils.speciestable is not None else None

    def phylo_distances(self):
        """
        a copy of the phylodistances dataframe see (get_evolutionary_distance())

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
        pos = 0
        speciestable = dict(utils.speciestable)
        for i, un in enumerate(counts):
            aslicetoavg = homoset.homo_matrix[ind[pos:pos + un]]
            bslicetoavg = homoset.fulleng[ind[pos:pos + un]]
            self.species[speciestable[i]].average_entropy = aslicetoavg.mean(axis=0)
            self.species[speciestable[i]].average_size = bslicetoavg.sum(0).mean()
            # variances are mean variances over all values
            self.species[speciestable[i]].var_entropy = aslicetoavg.var(axis=0).mean()
            pos += un
            self.species[speciestable[i]].tot_homologies = un
        print "average all species : " + str(homoset.homo_matrix.mean(axis=0))
        homoset.averagehomo_matrix = np.array([homoset[homo].mean for homo in homoset.homo_namelist])
        for _, val in homoset.iteritems():
            val.compute_averages()
# CpCpGpApApTpApTpApTpTpCpCpGpApApTpApTpApTpTpCpCpGpApApTpApTpApTpTpCpCpGpApApTpApTpApTpTpTpTpCpCpGpApApTpApTpApTpTp
# GbGbCbTbTbAbTbAbTbAbAbGbGbCbTbTbAbTbAbTbAbAbGbGbCbTbTbAbTbAbTbAbAbGbGbCbTbTbAbTbAbTbAbAbAbGbGbCbTbTbAbTbAbTbAbAbAb

    def compare_species(self, showvar=True, to='full'):
        """
        compare the species according to their mean CUB, plot the mean CUB
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
            showvar: bool to true if you want the plot to show the CUB value mean variance as sizes of the dots.
            to: str, flag to 'full' or 'avg' to know to what the tRNA entropy should be compared against (either
                the full entropy values from the get_full_genomes() function or the averaged one from the set of
                homologies we have)

        Raises:
             UnboundLocalError: "you have to compute averages, use PyCUB.compute_averages(homoset)",
                "you have to compute tRNA values, use PyCUB.get_tRNAcopy()"
            AttributeError: "try avg or full"

        """
        # TODO: totest
        pdb.set_trace()
        e_subspecies = np.zeros((len(self.species), utils.CUBD))
        s_subspecies = np.zeros((len(self.species), 1))
        vare_subspecies = np.zeros((len(self.species), 1))
        tRNAentropydist = np.zeros((len(self.species), utils.CUBD))
        tRNA_number = np.zeros((len(self.species), 1))
        genome_size = np.zeros((len(self.species), 1))
        num_genes = np.zeros((len(self.species), 1))
        expectedncgenome = np.zeros((len(self.species), 1))
        suff = 0
        phylo_distances = self.phylo_distances()
        for i, (name, specie) in self.species.iteritems():
            if specie.average_entropy is not None:
                e_subspecies[i] = specie.average_entropy
                s_subspecies[i] = specie.average_size
                vare_subspecies[i] = specie.var_entropy
            else:
                raise UnboundLocalError("you have to compute averages, use PyCUB.compute_averages(homoset)")

            if specie.copynumbers is not None:
                if specie.tRNAentropy is not None:
                    if to == 'full':
                        tRNAentropydist[i] = euclidean(specie.tRNAentropy, specie.fullentropy)
                    elif to == 'avg':
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
                    expectedncgenome[i] = specie.genome_size - (specie.num_genes * specie.average_size)
                # if is zero, will be black colored
            if specie.num_genes:
                num_genes[i] = specie.num_genes
            # if is zero, will be black colored
        print "we have " + str(suff) + " species with sufficient statistics in their tRNA values"
        if alg == 'tsne':
            red = man.TSNE(n_components=n, perplexity=perplexity).fit_transform(e_subspecies)
        elif alg == 'pca':
            red = PCA(n_components=n).fit_transform(e_subspecies)
        alg = cluster.DBSCAN(eps=eps, min_samples=7, algorithm='auto', n_jobs=-1)
        clusters = alg.fit_predict(e_subspecies).tolist()
        colormap = list(utils.colormap)
        colors = [colormap[0] if not dist else
                  utils.rgb2hex([26, 188, 50 + np.floor(106 * dist)]) for dist in tRNAentropydist]
        data = dict(x=red[:, 0], y=red[:, 1],
                    species=self.species.keys(),
                    meanentropy=["%.2f" % i.mean() for i in e_subspecies],
                    color=colors,
                    recent=colors,
                    s_subspecies=s_subspecies,
                    vare_subspecies=vare_subspecies,
                    tRNAentropydist=tRNAentropydist,
                    tRNA_number=tRNA_number,
                    genome_size=genome_size,
                    expectedncgenome=expectedncgenome,
                    num_genes=num_genes,
                    clustes=clusters)

        if self.species[self.species.keys()[0]].fullentropy is not None:
            # we have computed the full entropy for the species.
            data.update({'efull': [specie.fullentropy for _, specie in self.species.iteritems()]})
            data.update({'fullGCcount': [specie.fullGCcount for _, specie in self.species.iteritems()]})
        else:
            print "you could compute the full entropy of the species to \
            see the differences with the regular entropy"
        if self.species[self.species.keys()[0]].fullvarentropy is not None:
            data.update({'varefull': [specie.fullvarentropy for _, specie in self.species.iteritems()]})
            data.update({'varGCcount': [specie.varGCcount for _, specie in self.species.iteritems()]})
        # compare it to the phylogenetic distance matrix and the sizes
        if phylo_distances is not None:
            names = list(phylo_distances.index)
            distances = np.zeros(len(self.species))
            for i, val in enumerate(self.species.keys()):
                if val in names:
                    distances[i] = phylo_distances[val].values.mean()
            data.update({'distances': distances})

        source = ColumnDataSource(data=data)
        output_notebook()
        labe = ["show clusters", "show num_genes", "show expectedncgenome", "show genome_size",
                "show full genome CUB", "show homology only CUB", "show avg GCcount", "show variance of full CUB",
                "show variance of full GC count", "show distance to tRNA UB", "show tRNA number", "show avg phylodistance",
                "show full phylo distance"]
        callback = CustomJS(args=dict(source=source), code=utils.callback_all)
        radio_button_group = widgets.RadioButtonGroup(
            labels=labe, callback=callback, active=9)
        hover = HoverTool(tooltips=[("species: ", "@species"), ("mean entr: ", "@meanentropy"),
                                    ("expectedncgenome: ", "@expectedncgenome"), ("tRNA CN: ", "@tRNA_number"),
                                    ("dist2tRNA Bias Value: ", "@tRNAentropydist"), ("fullGCbias: ", "@fullGCcount")])
        p = figure(title="exploration of every homologies",
                   tools=[hover, WheelZoomTool(), PanTool(), SaveTool(), ResetTool()],
                   plot_width=800, plot_height=600)
        p.circle(x='x', y='y', source=source, color='color',
                 size=[(size / 2) + (size * 2 * (self.all_homoset[homo].var))] if showvar else size)
        save(column(radio_button_group, p), "utils/templot/homology_compare.html")
        show(column(radio_button_group, p))

    def regress_on_species(self, without = [""], full=True,onlyhomo=False, perctrain=0.8):
        """
        """
        # TODO: add regression (what is most explanatory?, How precise (test set), which one can be removed?)
        pdb.set_trace()
        params = []
        phylo_distances = self.phylo_distances()
        spe = self.species.values()[0]        
        if spe.genome_size and "genome_size" not in without:
            params.append([spe.genome_size for spe in self.species.values()])
            if spe.num_genes and "num_genes" not in without:
                params.append([specie.genome_size - (specie.num_genes * specie.average_size) for specie in self.species.values()])
        if spe.isplant_pathogen is not None and "isplant_pathogen" not in without:
            params.append([spe.isplant_pathogen for spe in self.species.values()])
        if spe.isanimal_pathogen is not None and "isanimal_pathogen" not in without:
            params.append([spe.isanimal_pathogen for spe in self.species.values()])
        if spe.isplant_symbiotic is not None and "isplant_symbiotic" not in without:
            params.append([spe.isplant_symbiotic for spe in self.species.values()])
        if spe.isbrown_rot is not None and "isbrown_rot" not in without:
            params.append([spe.isbrown_rot for spe in self.species.values()])
        if spe.iswhite_rot is not None and "iswhite_rot" not in without:
            params.append([spe.iswhite_rot for spe in self.species.values()])
        if spe.copynumbers is not None and "copynumbers" not in without:
            params.append([spe.copynumbers for spe in self.species.values()])
        if spe.average_size is not None and "average_size" not in without:
            params.append([spe.average_size for spe in self.species.values()])
        if spe.fullGCcount is not None and "fullGCcount" not in without:
            params.append([spe.fullGCcount for spe in self.species.values()])
        if spe.varGCcount is not None and "varGCcount" not in without:
            params.append([spe.varGCcount for spe in self.species.values()])
        if spe.tot_homologies is not None and "tot_homologies" not in without:
            params.append([spe.tot_homologies for spe in self.species.values()])
        if phylo_distances is not None and "phylo_distances" not in without:
            names = list(phylo_distances.index)
            params.append([phylo_distances[val].values.mean() for val in self.species.keys() if val in names])
        dataset = np.zeros((len(self.sepcies),utils.CUBD))
        for _,v in self.species.iteritems():
            dataset[i] = v.average_entropy if onlyhomo else v.fullentropy
        if not full:
            dataset = dataset.mean(1)
        if algo == "lasso":
            # http://scikit-learn.org/stable/modules/generated/sklearn.linear_model.LassoCV.html#sklearn.linear_model.LassoCV
            lm = LassoCV(eps=eps, n_alphas=n_alphas, fit_intercept=True, normalize=True,
                         precompute=False, max_iter=1000, verbose=2, tol=0.0001, warm_start=False,
                         positive=False, random_state=None, selection='cyclic')
        elif algo == "nn":
            # http://scikit-learn.org/stable/modules/generated/sklearn.linear_model.Perceptron.html#sklearn.linear_model.Perceptron
            lm = Perceptron(penalty=None, alpha=0.0001, fit_intercept=True, max_iter=None, tol=None,
                            shuffle=True, verbose=2, eta0=1.0, n_jobs=-1, random_state=0, class_weight=None,
                            warm_start=False)
        params = np.asarray(params).T
        lm.fit(params[:len(self.species)*perctrain], dataset[:len(self.species)*perctrain])
        score = lm.score(params[len(self.species)*perctrain:], dataset[len(self.species)*perctrain:], sample_weight=None)
        coeff = lm.get_params()
        coef_ = lm.coeff_
        intercept = lm.intercept_
        print "the R^2 score if of: "+str(score)
        return coeff, coef_, intercept, score

    #TODO: create a model that can say if the species has a high  CUB or low, given data on the species and on the protein

    def compare_homologies(self, homoset, homosapiens=False, mindistance=10, preserved=True, size=10,
                           minpreserv=0.9, minsimi=0.9, showvar=True, eps=0.48):
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
            eps: the hyperparamter of the clustering algorithm applied to this dataset


        Raises:
            UnboundLocalError: "you need to compute the averages of the all_homoset. use PyCUB.compute_averages(homoset)"


        """
        # TODO: totest
        pdb.set_trace()
        if not homosapiens:
            if homoset[0].isrecent is None:
                homoset.compute_ages(preserved=preserved, minpreserv=minpreserv, minsimi=minsimi)
        else:
            pass
            # TODO: code the version for homo sapiens where we know exactly this distance with more data
            # and better inference metrics

        # display the differences between recent homologies and older ones
        if homoset.averagehomo_matrix is None:
            raise UnboundLocalError("you need to compute the averages of homoset. use PyCUB.compute_averages(homoset)")
        if alg == 'tsne':
            red = man.TSNE(n_components=n, perplexity=perplexity).fit_transform(homoset.averagehomo_matrix)
        elif alg == 'pca':
            red = PCA(n_components=n).fit_transform(homoset.averagehomo_matrix)
        else:
            raise AttributeError("wrong algorithm")
        alg = cluster.DBSCAN(eps=eps, min_samples=7, algorithm='auto', n_jobs=-1)
        clusters = alg.fit_predict(homoset.averagehomo_matrix).tolist()
        n_clusters_ = len(set(clusters))
        if n_clusters_ > 10:
            print "ooups you have more than 10 clusters"
        colormap = list(utils.colormap)
        colors = [colormap[0] if not homoset[homo].isrecent else
                  utils.rgb2hex([26, 188, np.floor(156 * homoset[homo].isrecent)]) for homo in homoset.homo_namelist]
        data = dict(x=red[:, 0], y=red[:, 1],
                    homologies=homoset.homo_namelist,
                    meanentropy=["%.2f" % homoset.averagehomo_matrix[i].mean() for i in range(len(homoset.averagehomo_matrix))],
                    color=colors,
                    recent=colors,
                    clusters=clusters)
        # add average of similar protein name
        if homoset[0].similarity_scores is not None:
            data.update({'similarity_scores': [homoset[homo].similarity_scores.mean() for homo in homoset.homo_namelist]})
        if homoset[0].KaKs_Scores is not None:
            data.update({'KaKs_Scores': [homoset[homo].KaKs_Scores.mean() for homo in homoset.homo_namelist]})
        if homoset[0].nans is not None:
            data.update({'nans': [homoset[homo].nans.sum(1).mean() for homo in homoset.homo_namelist]})
        if homoset[0].lenmat is not None:
            data.update({'lengths': [homoset[homo].lenmat.sum(1).mean() for homo in homoset.homo_namelist]})
        if homoset[0].GCcount is not None:
            data.update({'gc': [homoset[homo].GCcount.mean() for homo in homoset.homo_namelist]})
        if homoset[0].weight is not None:
            data.update({'weight': [homoset[homo].weight for homo in homoset.homo_namelist]})
        if homoset[0].protein_abundance is not None:
            data.update({'protein_abundance': [homoset[homo].protein_abundance for homo in homoset.homo_namelist]})
        if homoset[0].mRNA_abundance is not None:
            data.update({'mRNA_abundance': [homoset[homo].mRNA_abundance for homo in homoset.homo_namelist]})
        if homoset[0].decay_rate is not None:
            data.update({'decay_rate': [homoset[homo].decay_rate for homo in homoset.homo_namelist]})
        if homoset[0].is_secreted is not None:
            data.update({'is_secreted': [homoset[homo].is_secreted for homo in homoset.homo_namelist]})
        if homoset[0].cys_elements is not None:
            data.update({'cys_elements': [homoset[homo].cys_elements for homo in homoset.homo_namelist]})
        if homoset[0].tot_volume is not None:
            data.update({'tot_volume': [homoset[homo].tot_volume for homo in homoset.homo_namelist]})
        if homoset[0].mean_hydrophobicity is not None:
            data.update({'mean_hydrophobicity': [homoset[homo].mean_hydrophobicity for homo in homoset.homo_namelist]})
        if homoset[0].glucose_cost is not None:
            data.update({'glucose_cost': [homoset[homo].glucose_cost for homo in homoset.homo_namelist]})
        if homoset[0].synthesis_steps is not None:
            data.update({'synthesis_steps': [homoset[homo].synthesis_steps for homo in homoset.homo_namelist]})
        if homoset[0].meancai is not None:
            data.update({'meancai': [homoset[homo].meancai for homo in homoset.homo_namelist]})
        source = ColumnDataSource(data=data)
        output_notebook()
        labe = ["show Recent/preserved", "show Nans avg", "show avg KaKs_Scores",
                "show avg similarity_scores", "show avg Length", "show avg GCcount", "showclusters",
                "Show weight", "Show prot abundance", "Show mRNA abundance", "Show half life", "Show secreted",
                "Show num of cys", "Show volume", "Show hydrophobicity", "show cost (glucose)",
                "Show synthesis cost", "Show conserved", "Show Pi", "Show CAI"]  # 20
        callback = CustomJS(args=dict(source=source), code=utils.callback_all)
        radio_button_group = widgets.RadioButtonGroup(
            labels=labe, callback=callback, active=0)
        hover = HoverTool(tooltips=[("homologies", "@homologies"), ("avg nans", "@nans"),
                                    ("mean_entr", "@meanentropy"), ("length", "@lengths"), ("GCcount", "@gc")])
        p = figure(title="exploration of every homologies",
                   tools=[hover, WheelZoomTool(), PanTool(), SaveTool(), ResetTool()],
                   plot_width=800, plot_height=600)
        p.circle(x='x', y='y', source=source, color='color',
                 size=[(size / 2) + (3 * self.all_homoset[homo].var)] if showvar else size)
        save(column(radio_button_group, p), "utils/templot/homology_compare.html")
        show(column(radio_button_group, p))

    def regress_on_genes(self, compare="full", without=['meancai'], algo="lasso", eps=0.001, n_alphas=100):
        """
        Will fit a regression curve on the CUB values of the different homologies according to the metadatas
        available for each of them.

        It will try to see if there is enough information in the metadata to retrieve CUB values. and if there is,
        how much for each metadata (if we constraint the number of regressors) is it better for entropy values, mean entropy
        or CAI values
        or raw frequency, should we remove some data

        Args:s
            without: list[str] of flags [similarity_scores, KaKs_Scores, nans, lenmat, GCcount, weight,
                protein_abundance, mRNA_abundance, decay_rate, cys_elements, tot_volume, mean_hydrophobicity,
                glucose_cost, synthesis_steps, is_recent, meancai]
            full: str flags to ["full", "mean", "cai"] full CUB values, meanCUB values, mean CAI as regressee

        """
        # TODO: totest
        pdb.set_trace()
        params = []
        if self.all_homoset[0].similarity_scores is not None and "similarity_scores" not in without:
            params.append([self.all_homoset[homo].similarity_scores.mean() for homo in self.all_homoset.homo_namelist])
        if self.all_homoset[0].KaKs_Scores is not None and "KaKs_Scores" not in without:
            params.append([self.all_homoset[homo].KaKs_Scores.mean() for homo in self.all_homoset.homo_namelist])
        if self.all_homoset[0].nans is not None and "nans" not in without:
            params.append([self.all_homoset[homo].nans.sum(1).mean() for homo in self.all_homoset.homo_namelist])
        if self.all_homoset[0].lenmat is not None and "lenmat" not in without:
            params.append([self.all_homoset[homo].lenmat.sum(1).mean() for homo in self.all_homoset.homo_namelist])
        if self.all_homoset[0].GCcount is not None and "GCcount" not in without:
            params.append([self.all_homoset[homo].GCcount.mean() for homo in self.all_homoset.homo_namelist])
        if self.all_homoset[0].weight is not None and "weight" not in without:
            params.append([self.all_homoset[homo].weight for homo in self.all_homoset.homo_namelist])
        if self.all_homoset[0].protein_abundance is not None and "protein_abundance" not in without:
            params.append([self.all_homoset[homo].protein_abundance for homo in self.all_homoset.homo_namelist])
        if self.all_homoset[0].mRNA_abundance is not None and "mRNA_abundance" not in without:
            params.append([self.all_homoset[homo].mRNA_abundance for homo in self.all_homoset.homo_namelist])
        if self.all_homoset[0].decay_rate is not None and "decay_rate" not in without:
            params.append([self.all_homoset[homo].decay_rate for homo in self.all_homoset.homo_namelist])
        params.append([self.all_homoset[homo].is_secreted for homo in self.all_homoset.homo_namelist])
        if self.all_homoset[0].cys_elements is not None and "cys_elements" not in without:
            params.append([self.all_homoset[homo].cys_elements for homo in self.all_homoset.homo_namelist])
        if self.all_homoset[0].tot_volume is not None and "tot_volume" not in without:
            params.append([self.all_homoset[homo].tot_volume for homo in self.all_homoset.homo_namelist])
        if self.all_homoset[0].mean_hydrophobicity is not None and "mean_hydrophobicity" not in without:
            params.append([self.all_homoset[homo].mean_hydrophobicity for homo in self.all_homoset.homo_namelist])
        if self.all_homoset[0].glucose_cost is not None and "glucose_cost" not in without:
            params.append([self.all_homoset[homo].glucose_cost for homo in self.all_homoset.homo_namelist])
        if self.all_homoset[0].synthesis_steps is not None and "synthesis_steps" not in without:
            params.append([self.all_homoset[homo].synthesis_steps for homo in self.all_homoset.homo_namelist])
        if self.all_homoset[0].is_recent is not None and "is_recent" not in without:
            params.append([0 if self.all_homoset[homo].ishighpreserved else self.all_homoset[homo].is_recent
                           for homo in self.all_homoset.homonamelist])
        if self.all_homoset[0].meancai is not None and "meancai" not in without:
            params.append([self.all_homoset[homo].meancai for homo in self.all_homoset.homo_namelist])
        dataset = self.all_homoset.averagehomo_matrix if full else self.all_homoset.averagehomo_matrix.mean(1)
        if algo == "lasso":
            # http://scikit-learn.org/stable/modules/generated/sklearn.linear_model.LassoCV.html#sklearn.linear_model.LassoCV
            lm = LassoCV(eps=eps, n_alphas=n_alphas, fit_intercept=True, normalize=True,
                         precompute=False, max_iter=1000, verbose=2, tol=0.0001, warm_start=False,
                         positive=False, random_state=None, selection='cyclic')
        elif algo == "nn":
            # http://scikit-learn.org/stable/modules/generated/sklearn.linear_model.Perceptron.html#sklearn.linear_model.Perceptron
            lm = Perceptron(penalty=None, alpha=0.0001, fit_intercept=True, max_iter=None, tol=None,
                            shuffle=True, verbose=2, eta0=1.0, n_jobs=-1, random_state=0, class_weight=None,
                            warm_start=False)
        lm.fit(np.asarray(params), dataset)
        coeff = lm.get_params()
        coef_ = lm.coeff_
        intercept = lm.intercept_
        return coeff, coef_, intercept

    def getRelation2G3DD(self,intrachromosome="utils/meta/3Dmodel/interactions_HindIII_fdr0.01_intra_cerevisiae.csv", interchromose=[
                    "utils/meta/3Dmodel/cerevisia_inter1.csv"
                    "utils/meta/3Dmodel/cerevisia_inter2.csv"
                    "utils/meta/3Dmodel/cerevisia_inter3.csv"
                    "utils/meta/3Dmodel/cerevisia_inter4.csv"
                    "utils/meta/3Dmodel/cerevisia_inter5.csv"], n=32000, seq='cds'):
        """
        https://www.nature.com/articles/ncomms6876

        retrieve the data for the species sacharomyces cerevisiae and Schizosaccharomyces pombe
        and find if similarity distances of CUB using entropy between genes of this species is predictive
        of closeness of genes in the nucleus.
        """
        
        
        # get gene distance matrix from entropy value distance or Andres Schindelin metrics
        # compare to see how much the distance between one can explain the distance between another by
        # regression
        # retrieve the data.
        intra = pd.read_csv(intrachromosome)
        inter = pd.concat([pd.read_csv() for interchro in interchromose])
        # getting all the genes
        ftp = FTP('ftp.ensemblgenomes.org')
        kingdom = 'fungi'
        ftp.login()
        ftp.cwd('pub/release-40/' + kingdom + '/fasta/')
        data = []
        ftp.retrlines('NLST', data.append)
        species_name = 'saccharomyces_cerevisiae'
        for d in data:
            ftp.cwd(d)
            if d in species_name:
                link = []
                ftp.cwd(seq)
                ftp.retrlines('NLST', link.append)
                with open("utils/data/temp.fa.gz", "wb") as file:
                    for i in link:
                        if i[-9:] == "all.fa.gz":
                            ftp.retrbinary("RETR " + i, file.write)
                with gzip.open("utils/data/temp.fa.gz", "rt") as handle:
                    for record in SeqIO.parse(handle, "fasta"):
                        valH, _, _, _ = utils.computeyun(record.seq._data, setnans=setnans, normalized=normalized, by=by)
                        record.
                os.remove("utils/data/temp.fa.gz")
            # we associate valH to a position retrieve the positions
            # we associate compute the CUF as well
        
        # for all values
        # we take n subsets for which we compute the average CUB, then for each we compute a CUB distance value as a big distance
        # matrix of size n x n n should be 32 000
        # we create another distance matrix using an Andres distance metrics on the CUF values
        
        # for 1 intra value we create a first matrix with distance = 1 if interact or if neighboors, else 0. 
        # we then look at the inter and add for each other chrom in the rest of the matrix the same values, adding 1 for the last position
        # of n to the one of n+1
        # we then add the next intra and do the same for inter
        # we then update the matrix by adding n for positions related by a min path of n
        # we then bin the positions similarly. into the n groups defined by genomic positions. 
        # we create a new matrix n x n and for each groups, we compute the average distance to the other groups...

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
            plt.imshow(1./(1+np.array(utils.phylo_distances)))
            plt.savefig("utils/templot/evolutionarydistances.pdf")
            plt.show()
        else:
            raise UnboundLocalError("compute the phylo distance matrix first (look at the doc)")
# SPECIAL FUNCTION

    def _dictify(self, save_workspace, save_homo, add_homosets):
        """
        Used by the saving function. transform the workspace object into a dictionary that can be
        json serializable
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
                        print v
                        self.working_homoset.update({v: self.all_homoset[v]})
            elif key == "all_homoset":
                pass
            else:
                ret.update({key: val})
        return ret
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
