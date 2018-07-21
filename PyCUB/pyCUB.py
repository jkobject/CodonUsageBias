""""
Created by Jeremie KALFON
Date : 21 FEV 2018
University of Kent, ECE paris
jkobject.com

PyCUB is a Project which goal is to understand the particular dynamics of the codon usage
bias and ...

you can find more about pycub here ..

if you find you need to add any..
"""
import os
import json
import zipfile
import shutil
from ftplib import FTP
import gzip
try:
    from urllib2 import urlopen as urlopen
except:
    from urllib.request import urlopen as urlopen

from joblib import Parallel, delayed
import multiprocessing

#from rpy2.robjects.packages import importr
from ete2 import NCBITaxa
#from rpy2 import robjects
#import rpy2.robjects.packages as rpackages
#from rpy2.robjects.vectors import StrVector


import pandas as pd
import numpy as np

import espece as spe
import homoset as hset
import utils
import homology as h

from bokeh.plotting import *
from bokeh.models import *
import matplotlib.pyplot as plt
from sklearn import manifold as man
from sklearn.decomposition import PCA


import pdb


class PyCUB(object):
    """PyCUB is the main object of the project that allows the user to access most of the functions

        When using it, please follow the documentation and examples on notebooks thought you can
        still use it as you please and use some of the nice tricks provided here and in python

            Params:
            ------
            species: dictionary of Espece objects from the name of the species.
                                  (see espece.py)
            _is_saved : bool is this session saved
            _is_loaded : bool is this session loaded from somewhere.
            session : str representing the name of the current session
                    (used when saving for example)
            links : dict of all the links readily available in PyCUB.
                    for the project of Jeremie KALFON
                please use whatever datasets you may find usefull
                (you can also download from Ensembl)
            working_homoset : homoset object that store the homology that we want to work on
            all_homoset : homoset object that store the all the homology

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
        # TODO update the metadata infos and push them to the github
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

    def __init__(self, species={}, _is_saved=False,
                 _is_loaded=False, working_homoset=None, all_homoset=None, session='session1'):
        """
        will initialize the object with the different values you might have from another project

        Params:
        ------
        species: dictionary of Espece objects from the name of the species.
                                        (see espece.py)
        homoset : homoset object that store the all the homology
        session :


        _is_saved : trivial system only boolean
        _is_loaded : trivial system only boolean

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

    def get_data(self, From='yun', homonames=None, kingdom='compara=fungi', sequence='cdna',
                 additional='type=orthologues', saveonfiles=False, normalized=True, setnans=False,
                 by="entropy", using="normal", inpar=True):
        """
                Download the data from somewhere on the web (Ensembl, Yun(with links))

                you can provide a lot of different values to scrape Ensembl's datasets
                it will compute from ensembl to retrieve a similar dataset as what yun's
                data is.

                Params:
                -------
                From: flag 'yun' or 'ensembl':
                    homonames: what particular homologies you want to scrap
                    kingdom: same for kingdoms
                    sequence: the type of sequences you want to use
                    additional: additional information about the scrapping
                    saveonfiles: save the unprocessed data before populating working homoset
                normalized: bool if you want the values to be normalized by the length of the codons
                            (lengths are always saved)
                setnans: bool if you want to save the nans as metadata
                by: flag 'entropy', 'entropyLocation' (entropy location), 'frequency'
                using: flag 'random' 'normal' 'permutation' 'full'
                inpar: bool or int for parallel computing and number of core


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
                print "you can't copute codon frequency with Yun's data..."
                return False
            Parallel(n_jobs=8)(delayed(utils.getyun)(key, val) for key, val in
                               self.links['yun'].iteritems())
            self.load(All=False if homonames is not None else True, filename=homonames,
                      From=From, by=by, inpar=inpar)
        elif From == 'ensembl':
            if homonames == 'all' or homonames is None:
                with open('utils/meta/homolist.json', "r") as f:
                    homonames = json.loads(f.read())
            self.all_homoset = hset.HomoSet(homo_namelist=homonames, datatype=by)
            print "doing all " + str(len(homonames)) + " homologies"
            if bool(inpar):
                values = Parallel(n_jobs=num_cores)(delayed(utils.loadfromensembl)(
                    name, kingdom, sequence,
                    additional, saveonfiles,
                    normalized, setnans, i, by, using) for i, name in enumerate(homonames))
                for i, val in enumerate(values):
                    self.all_homoset.update({homonames[i]: val})

            else:
                for i, name in enumerate(homonames):
                    self.all_homoset.update({name:
                                             utils.loadfromensembl(name, kingdom, sequence,
                                                                   additional, saveonfiles,
                                                                   normalized, setnans, i, by, using)})
            # TODO: test full pipeline with frequency/entropy/entropylocation
            taxons, species = self.all_homoset.preprocessing(withtaxons=True)
            print "computing tRNA copy numbers"
            for i, spece in enumerate(species):
                espece_val = spe.Espece(name=spece, taxonid=taxons[i])
                espece_val.get_tRNAcopy()
                self.species.update({spece: espece_val})
            self.all_homoset.loadhashomo()
        else:
            print 'not the right From'

    def get_metadata_Ensembl(self, kingdoms):
        """
        download it and put it where it belongs in the Espece object
        parse the server https://fungi.ensembl.org/info/website/ftp/index.html
        will also get the metadata from the kingdoms that you are analysing

        Params:
        ------
        kingdoms: flag the type of kingdoms you wanna have 'fungi' 'bacteria' 'plants'
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

        Params:
        ------
        From: flag designer of the function to load metadatas
        inpar: bool for parallel processing
        """
        if not os.path.exists('utils/meta'):
            os.mkdir('utils/meta')
        if From == 'jerem':
            num_cores = -1 if inpar else 1
            Parallel(n_jobs=num_cores)(delayed(utils.mymeta)(key, val) for key, val in
                                       self.links['mymeta'].iteritems())

    def import_metadataTobias(self):
        """
        will import the metadata obtained from tobias for the fungi species affiliated to
        cerevisiae to each species for further diagnostics.
        """
        data = pd.read_csv("Yun_Species_Context.csv")
        for i, species in enumerate(data["Genome"]):
            if species in self.species:
                self.species[species].metadata["num_genes"] = data["No_Genes"][i]
                self.species[species].metadata["plant_pathogen"] = data["plant_pathogen"][i]
                self.species[species].metadata["animal_pathogen"] = data["animal_pathogen"][i]
                self.species[species].metadata["genome_size"] = data["Genome_Size"][i]
                self.species[species].metadata["plant_symbiotic"] = data["mycorrhizal"][i] or data["endophyte"][i]
                self.species[species].metadata["brown_rot"] = data["brown_rot"][i]
                self.species[species].metadata["white_rot"] = data["white_rot"][i]

# LOADINGS AND SAVINGS

    def load(self, session=None, All=False, filename='first500', From=None, by='entropyLocation', inpar=True):
        """
                Get the data that is already present on a filename

                Either load from Yun's datasets or from an already saved session.
                Is being called by get_data. But you can call it to just use one of Yun's files
                as well

                Params:
                ------
                From: if this flag is set to 'yun' it means that the filename is made of Yundata
                in which case we will create directly the homology map in the same time as the rest
                of the PyCUB object.
                    All: bool set to true if load everything from Yun
                    by: same flag as get_data (for Yun's files here).
                    filename: str the particular filename when not loading them all
                session: str if a session name is provided, then will load a zip file from
                this session's folder

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
                    self.species.update({
                        row['name']: spe.Espece(name=row['b'],
                                                link=dflink.loc[dflink['name'] == row['name'],
                                                                'b'].tolist()[0])})
                    utils.speciestable.update({i: row['name']})
                    i += 1
                self.all_homoset.species_namelist = df['name'].tolist()

                # getting the homologies now
                dou = 0
                if inpar:
                    values = Parallel(n_jobs=-1)(delayed(utils.homoyun)(
                        separation, filename,
                        homology, by=by) for homology in homo_namelist)
                    for i, val in enumerate(values):
                        self.all_homoset.update({homo_namelist[i]: h.homology(
                            # TODO
                            full=val[0].as_matrix(), names=val[1],
                            nans=np.sum(val[2].as_matrix(), axis=1),
                            lenmat=val[3].as_matrix(), doub=val[4])})
                        dou += np.count_nonzero(val[4])
                else:
                    for homo in homo_namelist:
                        val = utils.homoyun(separation, filename, homo, by=by)
                        self.all_homoset.update({homo: h.homology(
                            # TODO
                            full=val[0].as_matrix(), names=val[1],
                            nans=np.sum(val[2].as_matrix(), axis=1),
                            lenmat=val[3].as_matrix(), doub=val[4])})
                        dou += np.count_nonzero(val[4])
                # create the hashomomatrix
                self.all_homoset.loadhashomo()
                print "you had " + str(dou) + " same species homologies"
                print "reviewed " + str(len(homo_namelist)) + " homologies "

                # if we haven't change the working with processing
                self._is_saved = False
                if not note:
                    shutil.rmtree(filename)
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
            elif From is None:
                if not os.path.isfile(folder + filename + '.json'):
                    print "unzipping " + folder + filename + '.json.gz'
                    os.system("gzip -d " + folder + filename + '.json.gz')
                with open(folder + filename + ".json", "r") as f:
                    print "loading from " + filename
                    self._undictify(json.loads(f.read()))
                    print "it worked !"
                    os.system("gzip " + folder + filename + '.json')

            print "you now have " + str(np.count_nonzero(self.all_homoset.hashomo_matrix)) +\
                " genes in total"

        else:
            print "hey, it looks like this object has already loaded some things"
            print "please use loadmore or use another object"
            print "you can delete this one with 'del' "

    def loadmore(self, filename='first500', by='entropyLocation'):
        """
                Get the data that is already present on a filename when you already have data

                is usefull to load more of Yun's datasets.
                is called when load is set to All


                Params:
                ------
                From: if this flag is set to 'yun' it means that the filename is made of Yundata
                in which case we will create directly the homology map in the same time as the rest
                of the PyCUB object. Here it is the only available option.
                    filename: str the filename to additionaly load
                    by: flag same as before

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
            print "You should try load first, this object is empty"

    def save(self, name, save_workspace=True, save_homo=True, cmdlinetozip="gzip"):
        """
        call to save your work. you should call save on specific data structure if this is what you
        want to save.
        Will call other object's save, will transform all the variable into dict and save the dicts
        as json files. will save the df also as json files. PyCUB and homoset have
        their own json file.
        adding some params because else the object may be too big

        Params:
        ------
        name: str the name of the particular save on this session
        save_workspace: bool to fale not to save working_homoset
        save_homo: bool to false not to save all_homoset
        cmdlinetozip: you need to tell the platform how to zip on your system (default to MacOS zipping)

        """
        filename = "utils/save/" + self.session + '/' + name + ".json"
        print "writing in " + name
        dictify = self._dictify(save_workspace, save_homo)
        data = json.dumps(dictify, indent=4, separators=(',', ': '))
        dirname = os.path.dirname(filename)
        if not os.path.exists(dirname):
            os.makedirs(dirname)
        with open(filename, 'w') as f:
            f.write(data)
            print "it worked !"
        # now we zip to save 90% memory space
        print "only work on mac for now, please write the cmd line to zip a file HERE"
        if cmdlinetozip == 'mac':
            os.system("ditto -c -k --sequesterRsrc " + filename + ' ' + filename + '.zip')
            os.remove(filename)
        if cmdlinetozip == 'gzip':
            os.system("gzip " + filename)
        # with zipfile.ZipFile(filename + '.zip', 'w') as myzip:
        #    myzip.write(name + json)
        self._is_saved = True


# PREPROCESSINGS

    def get_working_homoset(self, clusternb=None, species=None, homologies=None, cleanhomo=None, cleanspecies=None):
        """
        create a subset of all_homoset on which you would like to do further computation

        To use once you have clustered homology groups, else takes everything.
        Can also be used just to get a subset of the all homosets.


        Params:
        ------
        clusternb: set the cluster of the group you want to get need to be between 1 and homogroupnb
        homologies: list[str] the subset as a list you want to get from all_homoset
        (can be additional to a clusternb)
        species: list[str] the subset as a list you want to get from all_homoset
        (can be additional to a clusternb)
        Returns:
        ------
        a HomoSet object (see homoset.py)
        """
        if clusternb is not None:
            # if clustering has been done
            leng = len(self.all_homoset.clusters)
            if leng > 0:
                # TODO: retrieve all other datas in case we do this function later in the pipeline
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
                        perc = self.all_homoset.hashomo_matrix.sum(1) / self.all_homoset.hashomo_matrix.shape[1]
                        homo_name = [self.all_homoset.homo_namelist[i] for i in ind if perc[i] > cleanhomo]
                    else:
                        homo_name = [self.all_homoset.homo_namelist[i] for i in ind]
                    for x in homo_name:
                        homo.homodict.update({x: self.all_homoset.homodict[x]})
                    homo.homo_namelist = homo_name
                    homo.species_namelist = self.all_homoset.species_namelist

            else:
                print "you have not clusterized you 'all_homoset' if you want to just use" +\
                    "'order_from_matrix on it. else if you want to take everything just do working = all ;)"
                return False
        else:
            homo = self.all_homoset
        if homologies is not None:
            homo.homodict = {k: homo.homodict[k] for k in homologies}
        if species is not None:
            other = [item for item in homo.species_namelist if item not in species]
            homo.remove(sepcies=other)
        homo.loadhashomo()
        if cleanspecies is not None:
            homo.clean_species(thresh=cleanspecies)
            homo.loadhashomo()
        homo.loadfullhomo()
        self.working_homoset = homo
        return homo

    def get_subset(self, homoset, clusternb=None, species=None, homologies=None):
        """
        To use once if you want to further refine a set of homologies


        Params:
        ------
        clusternb: set the cluster of the group you want to get need to be between 1 and homogroupnb
        per_species:  set to true if you want to manually select a subset of species
        per_homologies: set to true if you want to manually select a subset of homologies
        homologies: the subset as a list or a tuple of int
        species: the subset as a list, or a list of int

        Return:
        ------
        a HomoSet object (see homoset.py)
        """
        # TODO: totest
        if homologies is not None:
            homoset.homodict = {k: homoset.homodict[k] for k in homologies}
        if species is not None:
            other = [item for item in homoset.species_namelist if not (item in species)]
            homoset.remove(sepcies=other)
        homoset.loadhashomo()
        homoset.loadfullhomo()
        return homoset

    def get_full_genomes(self, kingdom='fungi', seq='cds'):
        ftp = FTP('ftp.ensemblgenomes.org')
        ftp.login()
        ftp.cwd('pub/release-40/' + kingdom + '/fasta/')
        data = []
        ftp.retrlines('NLST', data.append)
        for d in data:
            if d[-10:] == 'collection':
                ftp.cwd(d)
                subdata = []
                ftp.retrlines('NLST', subdata.append)
                for sub in subdata:
                    if sub in self.species_namelist:
                        link = []
                        ftp.cwd(sub + '/' + seq)
                        ftp.retrlines('NLST', link.append)
                        with open("data/temp.fa.gz", "wb") as file:
                            if sub[0] > 'r':
                                da = link[2]
                            elif sub[0] > 'c':
                                da = link[1]
                            else:
                                da = link[2]
                            ftp.retrbinary("RETR " + da, file.write)
                        with gzip.open("data/temp.fa.gz", "rt") as handle:
                            val, gccount = selfcomputefullloc(handle)
                            self.species[sub].fullentropy = val.mean()
                            self.species[sub].var_entropy = val.var()
                            self.species[sub].fullGCcount = gccount
                        ftp.cwd('../..')
                        os.remove("data/temp.fa.gz")
                ftp.cwd('..')
            else:
                if sub in self.species_namelist:
                    link = []
                    ftp.cwd(sub + '/' + seq)
                    ftp.retrlines('NLST', link.append)
                    with open("data/temp.fa.gz", "wb") as file:
                        if sub[0] > 'r':
                            da = link[2]
                        elif sub[0] > 'c':
                            da = link[1]
                        else:
                            da = link[2]
                        ftp.retrbinary("RETR " + da, file.write)
                    with gzip.open("data/temp.fa.gz", "rt") as handle:
                        val, gccount = selfcomputefullloc(handle)
                        self.species[sub].fullentropy = val.mean()
                        self.species[sub].var_entropy = val.var()
                        self.species[sub].fullGCcount = gccount
                    ftp.cwd('../..')
                    os.remove("data/temp.fa.gz")

    def compute_full_entropy(self, handle):
        """
        """
        GCcount = 0
        val = []
        for x, record in enumerate(SeqIO.parse(handle, "fasta")):
            GCcount += (record.seq.count('G') + record.seq.count('C')) / len(record)
            pos = 0
            valH = np.zeros(len(amino)) if by != 'frequency' else np.zeros(59)
            for k, amin in enumerate(amino):
                nbcod = len(utils.codons[amin])  # replace Cleng
                count = np.zeros(nbcod)
                X = np.zeros(nbcod)
                mn = np.ones(nbcod) / nbcod
                for i, cod in enumerate(utils.codons[amin]):
                    count[i] = record.seq.count(cod)
                lengsubseq = count.sum()  # replace subSlength
                if by == 'frequency':
                    if lengsubseq == 0:
                        valH[pos:pos + nbcod] = np.NaN if setnans else 1. / nbcod
                    else:
                        E = count / lengsubseq
                        valH[pos:pos + nbcod] = E
                    pos += nbcod
                else:
                    if lengsubseq == 0:
                        valH[k] = np.NaN if setnans else 0.5
                    else:
                        Yg = multinomial.pmf(x=count, n=lengsubseq, p=mn)
                        # efor part
                        div, i = divmod(lengsubseq, nbcod)
                        X[:int(i)] = np.ceil(div) + 1
                        X[int(i):] = np.floor(div)
                        Eg = multinomial.pmf(x=X, n=lengsubseq, p=mn)
                        # end here
                        valH[k] = -np.log(Yg / Eg) / lengsubseq if normalized else -np.log(Yg / Eg)
            val.append(valH)
        val = np.array(val)
        return val, GCcount / x

    def get_evolutionary_distance(self, display_tree=False, size=40):
        """
        uses metadata of the ancestry tree and computes a theoretical evolutionary
        distance matrix between each species

        can optionaly take any hierarchical evolutionary file between a group of species

        Returns:
        -------
        the distance matrix (numpy array)

        """
        # TODO: totest
        ncbi = NCBITaxa()
        taxons = []
        for key, val in self.species.iteritems():
            taxons.append(val.taxonid)
        tree = ncbi.get_topology(taxons)  # taxons
        # finding what this tree looks like
        if display_tree:
            print tree.get_ascii(attributes=["sci_name", "rank"])
        with open('metaphylo/temp_tree.phy', 'w') as f:  # maybe will be newick format...
            f.write(tree.write())
        """
        try:
            # https://stackoverflow.com/questions/19894365/running-r-script-from-python
            base = importr('base')
            utiles = importr('utils')
        except:
            print "you need to have R installed to compute the distance"
            return
        packnames = ('treeio')
        names_to_install = [x for x in packnames if not rpackages.isinstalled(x)]
        if len(names_to_install) > 0:
            utiles.install_packages(StrVector(names_to_install))
        robjects.r('''
            treeText <- readLines(metaphylo/temp_tree.phy)
            treeText <- paste0(treeText, collapse="")
            library(treeio)
            tree <- read.tree(text = treeText) ## load tree
            distMat <- cophenetic(tree)
            write.table(distMat,"metaphylo/phylodistMat_temp.csv")
        ''')
        """
        df = pd.read_csv("metaphylo/phylodistMat_temp.csv")
        utils.phylo_distances = df
        utils.meandist = df.sum().sum() / (len(data)**2 - len(data))

        plt.figure(figsize=(size, 200))
        plt.title('evolutionary distances')
        plt.imshow(np.array(df))
        plt.savefig("utils/evolutionarydistances.pdf")
        plt.show()

    def speciestable(self):
        return utils.speciestable

    def compute_averages(self, homoset):
        """
        compute the average entropy
        """
        # TODO: totest
        if homoset.homo_matrix is None:
            homoset.loadfullhomo()
        if homoset.hashomo_matrix is None:
            homoset.loadhashomo()
        numspecies = homoset.hashomo_matrix.sum(axis=1)
        ind = homoset.homo_matrixnames.argsort()
        pos = 0
        for i in range(len(numspecies)):
            a = homoset.homo_matrix[ind[pos:pos + numspecies[i]]]
            b = homoset.fulleng[ind[pos:pos + numspecies[i]]]
            self.species[homoset.species_namelist].average_entropy = a.mean(axis=0)
            self.species[homoset.species_namelist].average_size = b.mean(axis=0)
            self.species[homoset.species_namelist].var_entropy = a.var(axis=0)
            self.species[homoset.species_namelist].var_size = b.var(axis=0)
            pos += numspecies[i]
            self.species[homoset.species_namelist].tot_homologies = numspecies[i]
        print "average all species : " + homoset.homo_matrix.sum(axis=0)

    def compare_species(self, tophylo=False, totRNA=False):
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

        #TODO: add the pressure from the amount of tRNA the species has (j value)
        """
        # TODO: totest
        e_subspecies = []
        s_subspecies = []
        vare_subspecies = []
        vars_subspecies = []
        suff = 0
        for specie in self.species.iteritems():
            if specie.average_entropy is not None:
                e_subspecies.append[specie.average_entropy]
                s_subspecies.append[specie.average_size]
            # compare it to the phylogenetic distance matrix and the sizes
            if specie.copynumbers is not None and totRNA and specie.copynumbers["num"] > 200:
                vare_subspecies.append[specie.var_entropy]
                vars_subspecies.append[specie.var_size]
                suff += 1
        print "we have " + str(suff) + " species with sufficient statistics in their tRNA values"
        # compare it to the tRNA CN entropy values
        if alg == 'tsne':
            red = man.TSNE(n_components=n, perplexity=perplexity).fit_transform(self.all_homoset.averagehomo_matrix)
        elif alg == 'pca':
            red = PCA(n_components=n).fit_transform(self.all_homoset.averagehomo_matrix)
        colors = [rgb() for homo in self.all_homoset]
        source = ColumnDataSource(
            data=dict(x=red[:, 0], y=red[:, 1],
                      label=["species : %s" % utils.speciestable[x__] for x__ in self.names,
                             "has trna : %s" % str(homo.isrecent) for homo in self.all_homoset],
                      color=colors))
        output_notebook()
        hover = HoverTool(tooltips=[("label", "@label")])
        p = figure(title="T-sne of homologous gene X for each species",
                   tools=[hover, BoxZoomTool(), WheelZoomTool(), SaveTool(), ResetTool()])
        p.circle(x='x', y='y', source=source, color='color', size=[])
        show(p)

        # TODO: compute the variation between mean and full entropies
        # TODO:

    def compare_homologies(self, homosapiens=False, mindistance=10, preserved=True, size=10,
                           minpreserv=400, minsimi=0.9, nans=True):
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

        Params:
        ------

        """
        # TODO: to test
        if not homosapiens:
            for _, homo in self.all_homoset.homodict.iteritems():
                # if the homology is one of species that are closely related solely (meaning it is a recent one)
                newhomology = True
                mean = []
                if preserved:
                    if len(homo.full) > minpreserv and homo.similarity_scores.mean() > minsimi:
                        homo.ishighpreserved = True
                        continue
                    else:
                        homo.ishighpreserved = False
                for species in self.all_homoset.homology.names:
                    mean.append(self.all_homoset.phylo_distance[species].mean())
                    if mean > mindistance:
                        newhomology = False
                if newhomology:
                    mean = np.array(mean).mean()
                    print "found a homology with mean " + str(mean)
                    homology.isrecent = mean
                else:
                    homology.isrecent = False
        else:
            pass
        # display the differences between recent homologies and older ones
        # TODO: code this complex plotting
        if self.all_homoset.averagehomo_matrix is None:
            self.all_homoset.averagehomo_matrix = np.array([homo.mean for _, homo in self.all_homoset.homodict.iteritems()])
        if alg == 'tsne':
            red = man.TSNE(n_components=n, perplexity=perplexity).fit_transform(self.all_homoset.averagehomo_matrix)
        elif alg == 'pca':
            red = PCA(n_components=n).fit_transform(self.all_homoset.averagehomo_matrix)
        colormap = ["#1abc9c", "#3498db", "#2ecc71", "#9b59b6", '#34495e', '#f1c40f',
                    '#e67e22', '#e74c3c', '#7f8c8d', '#f39c12']
        colors = [colormap[int(bool(homo.isrecent)) + 2 * int(bool(homo.ishighpreserved))] for homo in self.all_homoset]
        labels = []
        for key, homo in self.all_homoset.homodict.iteritems():
            labels.append("homology : " + key + "\n distance : " + str(homo.isrecent) + "\n nansNumber : " + str(homo.nans.mean()))
        source = ColumnDataSource(data=dict(x=red[:, 0], y=red[:, 1], label=labels, color=colors))
        output_notebook()
        hover = HoverTool(tooltips=[("label", "@label")])
        p = figure(title="T-sne of homologous gene X for each species",
                   tools=[hover, BoxZoomTool(), WheelZoomTool(), SaveTool(), ResetTool()])
        p.circle(x='x', y='y', source=source, color='color', size=[size * homo.var.mean() for _, homo in self.all_homoset.homodict.iteritems()])
        show(p)

    def plot_distances(self):
        """
        plot the phylogenetic distance matrix
        """
        if utils.phylo_distances is not None:
            plt.imshow(np.array(utils.phylo_distances))
        else:
            print "compute the phylo distance matrix first (look at the doc)"
# SPECIAL FUNCTION

    def _dictify(self, save_workspace, save_homo):
        """
        Used by the saving function. transform the workspace object into a dictionary that can be
        json serializable
        adding some params because else the object may be too big

        Params:
        ------
        save_workspace: bool to save working_homoset
        save_homo: bool to save all_homoset
        """
        dictispecies = {}
        for key, val in self.species.iteritems():
            dictispecies.update({key: val._dictify()})
        return {"species": dictispecies,
                "all_homoset": self.all_homoset._dictify() if save_homo and
                (self.all_homoset is not None) else None,
                "working_homoset": self.working_homoset._dictify() if save_workspace and
                (self.working_homoset is not None) else None
                }

    def _undictify(self, data):
        """
        same function but to retransform everything
        Here we don't use other classes undictify functions but we just recreate them by passing it
        to their init methods which is clearer.

        Params:
        ------
        data: dict to undictify into the workspace object
        """
        species = {}
        for key, val in data["species"].iteritems():
            species.update({key: spe.Espece(data=val)})
        self.species = species
        self._is_saved = False
        self._is_loaded = True
        self.all_homoset = hset.HomoSet(data=data["all_homoset"]) if data["all_homoset"] is not None else None
        self.working_homoset = hset.HomoSet(data=data["working_homoset"]) if \
            data["all_homoset"] is not None else None


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
