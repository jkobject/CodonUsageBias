""""
Created by Jeremie KALFON
Date : 21 FEV 2018
University of Kent, ECE paris
jkobject.com

"""
import os
import glob
import json
import zipfile
import shutil
try:
    from urllib2 import urlopen as urlopen
except:
    from urllib.request import urlopen as urlopen

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import espece as spe
import homoset as hset
import utils
import homology as h

import pdb


class PyCUB(object):
    """docstring for PyCUB

            Params:
            ------
            species: dictionary of Espece objects from the name of the species.
                                            (see espece.py)
            _is_saved : trivial system only boolean
            _is_loaded : trivial system only boolean
            network :
            working_homoset : homoset object that store the homology that we want to work on
            all_homoset : homoset object that store the all the homology

    """
    # this dictionary shows which codons encode the same AA
    codons = {
        'CYS': ['TGT', 'TGC'],
        'ASP': ['GAT', 'GAC'],
        'SER': ['TCT', 'TCG', 'TCA', 'TCC', 'AGC', 'AGT'],
        'GLN': ['CAA', 'CAG'],
        #'MET': ['ATG'],
        'ASN': ['AAC', 'AAT'],
        'PRO': ['CCT', 'CCG', 'CCA', 'CCC'],
        'LYS': ['AAG', 'AAA'],
        #'STOP': ['TAG', 'TGA', 'TAA'],
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

    amino = ['Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Gln', 'Glu', 'Gly', 'His',
             'Ile', 'Leu', 'Lys', 'Phe', 'Pro', 'Ser', 'Thr', 'Tyr', 'Val']

    links = {'yun': {
        'homology1t500.zip': 'https://www.dropbox.com/s/mqi5jdgvst21a0c/homology1t500.zip?dl=1',
        'homology601t1000.zip': 'https://www.dropbox.com/s/mqi5jdgvst21a0c/homology601t1000.zip?dl=1',
        'homology1001t1500.zip': 'https://www.dropbox.com/s/mqi5jdgvst21a0c/homology1001t1500.zip?dl=1',
        'homology1501t2000.zip': 'https://www.dropbox.com/s/mqi5jdgvst21a0c/homology1501t2000.zip?dl=1',
        'homology2001t2500.zip': 'https://www.dropbox.com/s/mqi5jdgvst21a0c/homology2001t2500.zip?dl=1',
        'homology2501t3000.zip': 'https://www.dropbox.com/s/mqi5jdgvst21a0c/homology2501t3000.zip?dl=1',
        'homology3001t3500.zip': 'https://www.dropbox.com/s/mqi5jdgvst21a0c/homology3001t3500.zip?dl=1',
        'homology3501t4000.tgz': 'https://www.dropbox.com/s/mqi5jdgvst21a0c/homology3501t4000.tgz?dl=1',
        'homology4001t4500.tgz': 'https://www.dropbox.com/s/mqi5jdgvst21a0c/homology4001t4500.tgz?dl=1',
        'homology4501t5000.tgz': 'https://www.dropbox.com/s/mqi5jdgvst21a0c/homology4501t5000.tgz?dl=1'},
        'ensembl': {}}

    def __init__(self, species={}, _is_preprocessed=False, _is_saved=False,
                 _is_loaded=False, working_homoset=None, all_homoset=None, session='session1'):
        """
        will..

        Params:
        ------
        species: dictionary of Espece objects from the name of the species.
                                        (see espece.py)
        homoset : homoset object that store the all the homology


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

    def get_data(self, From='yun', all=False, names=None, families=None, kingdom=None):
        """
                Download the data from somewhere on the web

                http://rest.ensemblgenomes.org/
        """

        if From == 'yun':
            for key, val in self.links['yun'].iteritems():
                url = val
                print "downloading " + key + " with urllib"
                if not os.path.exists('utils/data/' + key):
                    f = urlopen(url)
                    data = f.read()
                    with open('utils/data/' + key, "wb") as code:
                        code.write(data)

        # elif From == 'ensembl':

    def get_metadata(self, kingdoms):
        """
        download it and put it where it belongs in the Espece object
        parse the server https://fungi.ensembl.org/info/website/ftp/index.html
        will also get the metadata from the kingdoms that you are analysing

        """

    def get_mymetadata(self):
        """
        for Yun's data to get some more, go ahead and design your own metadata retrieval here.
        obviously you woud need to change some other functions.

        for me it is mean protein abundances in cerevisiae cells.
        """

# LOADINGS AND SAVINGS

    def load(self, session=None, filename='first500', From=None, separation="homology"):
        """
                Get the data that is already present on a filename

                Params:
                ------
                from_Yun: if this flag is set to True it means that the filename is made of Yundata
                in which case we will create directly the homology map in the same time as the rest
                of the PyCUB object.

        """
        folder = "utils/data/" if session is None else "utils/save/" + session + '/'
        if not self._is_loaded:
            if From == 'yun' and session is None:
                # then we process it according to how Yun Displays its data, trying to fill in
                # as much as we can
                homo_namelist = []
                homodict = {}
                speciesdict = {}
                nameB = ""
                self.all_homoset = hset.HomoSet()
                if not os.path.exists(folder + filename):
                    zip_ref = zipfile.ZipFile(folder + filename + '.zip', 'r')
                    zip_ref.extractall(folder)
                    zip_ref.close()
                    note = False
                else:
                    note = True
                filename = folder + filename
                print "Reviewing all the " + str(len(os.listdir(filename))) + " files"
                for f in sorted(os.listdir(filename)):
                    if f.endswith(".txt"):
                        nameA = f.split(separation)[0]
                        if(nameA != nameB):
                            nameB = nameA
                            homo_namelist.append(nameB)

                self.all_homoset.homo_namelist = homo_namelist
                df, dflink = utils.retrievelist()
                i = 0
                df = df.sort_values(by='name')
                for _, row in df.iterrows():
                    self.species.update({
                        row['name']: spe.Espece(name=row['b'],
                                                link=dflink.loc[dflink['name'] == row['name'],
                                                                'b'].tolist()[0])})
                    utils.speciestable.update({i: row['name']})
                    i += 1
                self.all_homoset.species_namelist = df['name'].tolist()  # we create an ordered full list of species
                matrix = np.zeros((len(self.all_homoset.species_namelist),
                                   len(self.all_homoset.homo_namelist)), dtype=np.bool)
                i = 0
                dou = 0
                for homology in homo_namelist:
                    try:
                        # Here we call a function looking for the 18 files and creating
                        # a list of names of species presents and a matrix with the
                        # 18 values from the 18 files
                        gentab, specieslist = utils.readcods_homology(separation, filename, homology)
                    except OSError:
                        print "you do not have the files here"
                    except:
                        print "UNKNOWN problem with homology " + homology
                    espe = ''
                    j = 0
                    names = []
                    mask = np.ones(gentab.shape[0], dtype=bool)
                    print "at homology " + homology + " we have " + str(len(specieslist)) + \
                        " species possessing it"
                    for n, species in enumerate(specieslist):  # we can go faster as it is sorted
                        # creating a list of number linked to our speciestable

                        # removing doublon
                        if species == espe:
                            dou += 1
                            mask[n] = False
                        else:
                            while species != self.all_homoset.species_namelist[j]:
                                j += 1
                            names.append(j)
                            matrix[j, i] = True
                            espe = species
                            j += 1
                    homodict.update({homology: h.homology(full=gentab[mask, :], names=names)})
                    i += 1
                self.all_homoset.hashomo_matrix = matrix
                print "you had " + str(dou) + " same species homologies (it can't be processed! for now)"
                print "reviewed " + str(len(homo_namelist)) + " homologies "
                # if we haven't change the working with processing
                self.all_homoset.homodict = homodict
                self._is_saved = False
                if not note:
                    shutil.rmtree(filename)
            elif From is None:
                if not os.path.isfile(folder + filename):
                    zip_ref = zipfile.ZipFile(folder + filename + '.json.zip', 'r')
                    zip_ref.extractall(folder)
                    zip_ref.close()
                    note = False
                else:
                    note = True

                with open(folder + filename + ".json", "r") as f:
                    jsoni = f.read()
                    print "loading from " + filename
                    self._undictify(json.loads(jsoni))
                    print "it worked !"
                if not note:
                    os.remove(folder + filename + ".json")

            self._is_loaded = True
            print "you now have " + str(np.count_nonzero(self.all_homoset.hashomo_matrix)) + " genes in total"

        else:
            print "hey, it looks like this object has already loaded some things"
            print "please use loadmore or use another object"
            print "you can delete this one with 'del' "

    def loadmore(self, filename='first500', From='yun', separation="homology"):
        """
                Get the data that is already present on a filename

                Params:
                ------
                from_Yun: if this flag is set to True it means that the filename is made of Yundata
                in which case we will create directly the homology map in the same time as the rest
                of the PyCUB object.

        """
        folder = "utils/data/"
        if self._is_loaded:
            if From == 'yun':
                # then we process it according to how Yun Displays its data, trying to fill in
                # as much as we can
                homo_namelist = []
                speciesdict = {}
                nameB = ""
                if not os.path.exists(folder + filename):
                    zip_ref = zipfile.ZipFile(folder + filename + '.zip', 'r')
                    zip_ref.extractall(folder)
                    zip_ref.close()
                    note = False
                else:
                    note = True
                filename = folder + filename
                print "Reviewing all the " + str(len(os.listdir(filename))) + " files"
                for f in sorted(os.listdir(filename)):
                    if f.endswith(".txt"):
                        nameA = f.split(separation)[0]
                        if(nameA != nameB):
                            nameB = nameA
                            homo_namelist.append(nameB)

                # comparing two lists
                dup = [item for item in homo_namelist if not (item in self.all_homoset.homo_namelist)]
                size = len(homo_namelist) - len(dup)
                if size != 0:
                    print "there is " + str(size) + " duplicate from previous loads.. not cool"
                    homo_namelist = dup
                # update homonamelist
                self.all_homoset.homo_namelist.extend(homo_namelist)
                matrix = np.zeros(((len(self.all_homoset.species_namelist), len(homo_namelist))), dtype=np.bool)
                matrix = np.hstack((self.all_homoset.hashomo_matrix, matrix))
                i = self.all_homoset.hashomo_matrix.shape[1]
                print "start the iteration process, I hope you haven't clusterized your data yet..else it won't work (for now)"
                # we have a dict with num,name
                # we have a list of all species
                # we have a list of the species in our homologies
                # we want to tag the matrix for them and to
                # update the species object and to create
                # the name list in our homology object
                dou = 0
                for homology in homo_namelist:
                    try:
                        gentab, specieslist = utils.readcods_homology(separation, filename, homology)
                        # replace preprocess yun
                    except OSError:
                        print "you do not have the files here"
                    except:
                        print "UNKNOWN problem with homology " + homology
                    espe = ''
                    j = 0
                    names = []
                    print "at homology " + homology + " we have " + str(len(specieslist)) + \
                        " species possessing it"
                    for species in specieslist:  # we can go faster as it is sorted
                        # creating a list of number linked to our speciestable

                        # removing doublon
                        if species == espe:
                            dou += 1
                        else:
                            while species != self.all_homoset.species_namelist[j]:
                                j += 1
                            names.append(j)
                            # pdb.set_trace()
                            matrix[j, i] = True
                            espe = species
                            j += 1
                    self.all_homoset.update({homology: h.homology(full=gentab, names=names)})
                    i += 1
                self.all_homoset.hashomo_matrix = matrix
                print "you had " + str(dou) + " same species homologies (it can't be processed! for now)"
                print "reviewed " + str(len(homo_namelist)) + " homologies "
                self._is_saved = False
                if not note:
                    shutil.rmtree(filename)
                self._is_loaded = True
                print "you now have " + str(np.count_nonzero(self.all_homoset.hashomo_matrix)) + " genes in total"
        else:
            print "You should try load first, this object is empty"

    def save(self, name, save_workspace=True, save_homo=True, cmdlinetozip="ditto -c -k --sequesterRsrc "):
        """
        call to save your work. you should call save on specific data structure if this is what you
        want to save.
        Will call other object's save, will transform all the variable into dict and save the dicts
        as json files. will save the df also as json files. PyCUB and homoset have
        their own json file.
        adding some params because else the object may be too big

        Params:
        ------
        save_workspace: don't save working
        save_homo: don't save all

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
        os.system(cmdlinetozip + filename + ' ' + filename + '.zip')
        # with zipfile.ZipFile(filename + '.zip', 'w') as myzip:
        #    myzip.write(name + json)
        os.remove(filename)
        self._is_saved = True


# PREPROCESSINGS

    def preprocess(self):
        """
        apply a preprocessing to the loaded DNA string files to have them resemble YUN's data
        (use another special file with all the computations detailed there.)

        Code inspired from YUN's PhD.
        """

    def get_working_homoset(self, per_cluster=True, clusternb=1, species=None, homologies=None):
        """
        To use once you have clustered homology groups, else takes everything.


        Params:
        ------
        clusternb: set the cluster of the group you want to get need to be between 1 and homogroupnb
        per_species:  set to true if you want to manually select a subset of species
        per_homologies: set to true if you want to manually select a subset of homologies
        homologies: the subset as a list
        species: the subset as a list

        Return:
        ------
        a HomoSet object (see homoset.py)
        """
        isbyhomo = False
        if per_cluster:
            # if clustering has been done
            leng = len(self.all_homoset.clusters)
            if leng > 0:
                if leng == len(self.all_homoset.species_namelist):
                    # version by species
                    homo = hset.HomoSet()
                    # pdb.set_trace()
                    ind = np.argwhere(np.asarray(self.all_homoset.clusters) == clusternb - 1)[:, 0]
                    species_name = [self.all_homoset.species_namelist[i] for i in ind]
                    homo.homodict = self.all_homoset
                    homo.homodict.remove(species_name)
                    homo.homo_namelist = self.all_homoset.homo_namelist
                    homo.species_namelist = species_name

                elif leng == len(self.all_homoset.homodict):
                    # version by homologies
                    isbyhomo = True
                    homo = hset.HomoSet()
                    ind = np.argwhere(np.asarray(self.all_homoset.clusters) == clusternb - 1)[:, 0]
                    homo_name = [self.all_homoset.homo_namelist[i] for i in ind]
                    homo.homodict = dict((x, self.all_homoset[x]) for x in homo_name)
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
            other = [item for item in homo.species_namelist if not (item in species)]
            homo.remove(sepcies=other)
        homo.loadhashomo(isbyhomo)
        homo.loadfullhomo()
        self.working_homoset = homo
        return homo


# SPECIAL FUNCTION

    def _dictify(self, save_workspace, save_homo):
        """
        Used by the saving function. transform the object into a dictionary that can be
        json serializable
        adding some params because else the object may be too big

        Params:
        ------
        save_workspace: don't save working
        save_homo: don't save all
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
        """
        species = {}
        for key, val in data["species"].iteritems():
            species.update({key: spe.Espece(data=val)})
        self.species = species
        self._is_saved = False
        self._is_loaded = True
        self.all_homoset = hset.HomoSet(data=data["all_homoset"]) if data["all_homoset"] is not None else None
        self.working_homoset = hset.HomoSet(data=data["working_homoset"]) if data["all_homoset"] is not None else None


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
