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

    def __init__(self, species={}, _is_preprocessed=False, _is_saved=False,
                 _is_loaded=False, _have_homo=False,
                 working_homoset=False, all_homoset=False, session='session1'):
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
        if os.path.isdir('utils/save/' + session):
            print 'you already have a session here (just a warning)'
        else:
            os.mkdir('utils/save/' + session)

    def get_data(self, all=False, names=None, families=None, kingdom=None):
        """
                Download the data from somewhere on the web

                http://rest.ensemblgenomes.org/
        """

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

    def load(self, session=None, filename='first500', fromYun=True, separation="homology"):
        """
                Get the data that is already present on a filename

                Params:
                ------
                from_Yun: if this flag is set to True it means that the filename is made of Yundata
                in which case we will create directly the homology map in the same time as the rest
                of the PyCUB object.

        """
        print "working on session: " + self.session
        folder = "utils/data/" if session is None else "utils/save"
        if not _is_loaded:

            if fromYun:
                # then we process it according to how Yun Displays its data, trying to fill in
                # as much as we can
                homo_namelist = []
                homodict = self.all_homoset.homodict
                speciesdict = {}
                nameB = ""
                self.all_homoset = hset.HomoSet()
                zip_ref = zipfile.ZipFile(folder + filename + '.zip', 'r')
                zip_ref.extractall(folder)
                zip_ref.close()
                filename = folder + filename
                print "Reviewing all the " + str(len(os.listdir(filename))) + " files"
                for f in sorted(os.listdir(filename)):
                    if f.endswith(".txt"):
                        nameA = f.split(separation)[0]
                        if(nameA != nameB):
                            nameB = nameA
                            homo_namelist.append(nameB)

                self.all_homoset.homo_namelist = homo_namelist
                for homology in homo_namelist:
                    try:
                        # Here we call a function looking for the 18 files and creating
                        # a list of names of species presents and a matrix with the
                        # 18 values from the 18 files
                        genDF, specieslist = utils.readcods_homology(separation, filename, homology)
                        speciesdict.update({homology: specieslist})
                        homodict.update({homology: h.homology(full=genDF)})
                    except OSError:
                        print "you do not have the files here"
                    except:
                        print "UNKNOWN problem with homology " + homology
                # if we haven't change the working with processing
                self.all_homoset.homodict = homodict
                self._is_loaded = True
                self._preprocess_yun(speciesdict)
                self._is_saved = False
                os.remove(filename)  # we are saving space
            else:
                zip_ref = zipfile.ZipFile(folder + filename + 'json.zip', 'r')
                zip_ref.extractall(folder)
                zip_ref.close()
                with open(folder + filename + ".json", "r") as f:
                    jsoni = f.read()
                    print "loading from " + filename
                    self._undictify(json.loads(jsoni))
                    print "it worked !"
                os.remove(folder + filename + ".json")  # we are saving space
                self._is_loaded = True

        else:
            print "hey, it looks like this object has already loaded some things"
            print "please use loadmore or use another object"
            print "you can delete this one with 'del' "

    def loadmore(self, filename='first500', fromYun=True, separation="homology"):
        """
                Get the data that is already present on a filename

                Params:
                ------
                from_Yun: if this flag is set to True it means that the filename is made of Yundata
                in which case we will create directly the homology map in the same time as the rest
                of the PyCUB object.

        """
        folder = "utils/data/"
        if _is_loaded:
            if fromYun:
                # then we process it according to how Yun Displays its data, trying to fill in
                # as much as we can
                homo_namelist = []
                homodict = self.all_homoset.homodict
                speciesdict = {}
                nameB = ""
                zip_ref = zipfile.ZipFile(folder + filename + '.zip', 'r')
                zip_ref.extractall(folder)
                zip_ref.close()
                filename = folder + filename
                print "Reviewing all the " + str(len(os.listdir(filename))) + " files"
                for f in sorted(os.listdir(filename)):
                    if f.endswith(".txt"):
                        nameA = f.split(separation)[0]
                        if(nameA != nameB):
                            nameB = nameA
                            homo_namelist.append(nameB)
                # comparing two lists
                dup = set(homo_namelist).intersection(self.all_homoset.homo_namelist)
                size = len(dup)
                if size != 0:
                    print "there is " + size + " duplicate from previous loads"
                    return 0
                # update homonamelist
                self.all_homoset.homo_namelist.append(homo_namelist)
                for homology in homo_namelist:
                    try:
                        genDF, specieslist = utils.readcods_homology(separation, filename, homology)
                        speciesdict.update({homology: specieslist})
                        homodict.update({homology: h.homology(full=genDF)})
                    except OSError:
                        print "you do not have the files here"
                    except:
                        print "UNKNOWN problem with homology " + homology

                self.all_homoset.homodict = homodict
                self._is_loaded = True
                self._preprocess_yun(speciesdict)
                self._is_saved = False
                os.remove(filename)  # we are saving space
            else:
                print "make sure to provide the session name as well"
                zip_ref = zipfile.ZipFile(folder + filename + 'json.zip', 'r')
                zip_ref.extractall(folder)
                zip_ref.close()
                with open(folder + filename + ".json", "r") as f:
                    jsoni = f.read()
                    print "loading from " + filename
                    self._undictify(json.loads(jsoni))
                    print "it worked !"
        else:
            print "You should try load first, this object is empty"

    def preprocess(self):
        """
        apply a preprocessing to the loaded DNA string files to have them resemble YUN's data
        (use another special file with all the computations detailed there.)

        Code inspired from YUN's PhD. 
        """

    def find_clusters(self):
        """
        Finds, for each homologies in the working homoset, groups that are part of compact clusters
        it will be using gaussian mixture clustering or DBSCAN and order them according 
        to the density of each cluster (we are interested in the densest ones) and assess 
        the quality using 3 criterion:BIC, , . 
        """

    def get_working_homoset(self, per_cluster=True, clusternb=1, per_species=False,
                            per_homologies=False, species=[], homologies):
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
        if per_cluster and self.clusters:
            dicta =
            self.all_homoset.homodict[self.all_homoset.species_namelist[self.all_homoset.clusters]]
        else:
            dicta = self.all_homoset.homodict
        if per_homologies:
            dicta[]
            ----- do we select a homology or a species ?? - --
        if per_species:

        else:
            print "you have not clusterized you 'all_homoset' if you want to just use 'order_from\
            _matrix on it. else if you want to take everything just do working = all ;)"
            return False

    def save(self, name, save_workspace=True, save_homo=True):
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
        if self._is_preprocessed:
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
            with zipfile.ZipFile(filename + '.zip', 'w') as myzip:
                myzip.write(filename)
            os.remove(filename)

    def _preprocess_yun(self, speciesdict):
        """
        a private function to preprocess specifically Yun's data

        Params:
        -------
        speciesdict: dict of homologies and their relevant dataframe.
        """

        df, dflink = utils.retrievelist()
        for i, row in df.iterrows():
            self.species.update({row['a']: spe.Espece(name=row['b'],
                                                      link=dflink.full.loc[[row['a']], ['b']])})
        specieslist = dflink['a'].tolist()  # we create an ordered full list of species
        matrix = np.zeros((len(specieslist), len(speciesdict)), dtype=np.bool)
        i = 0
        dou = 0
        # COMPLEXITY
        for key, species in speciesdict.iteritems():  # we go throught the list of lists
            j = 0
            print "at " + key + " we have " + str(len(species)) + " homologies"
            species.sort()
            nameb = ''
            for name in species:  # TODO :can be parallelized
                # the assumption here is that they are ordered the same so we should only test for
                # the next ones
                if name == nameb:
                    dou += 1
                else:
                    # we go until the occurence of the name in specieslist
                    while name != specieslist[j]:
                        j += 1
                    matrix[j, i] = True
                    # part where we update the Espece object
                    # pdb.set_trace()
                    d = self.all_homoset.homodict[key].full.loc[[name], :]
                    d = d.fillna(0.5)   # we make sure to keep only one value. (nameA nameB problem)
                    d = d.iloc[[0]]
                    d = d.rename(index={name: key})  # need that weird stuff because of the transpose
                    self.species[name].genes = self.species[name].genes.append(d)
                    j += 1
                    nameb = name

            i += 1
        # create DF here with specieslist and hset.homo_namelist - ---
        self.all_homoset.species_namelist = specieslist
        self.all_homoset.hashomo_matrix = matrix
        print "you had " + str(dou) + " same species homologies (it can't be processed! for now)"
        print "reviewed " + str(i) + " homologies and " + str(j) + "species"

    def _preprocess_yun_more(self, speciesdict):
        """
        a private function to preprocess specifically Yun's data

        Params:
        -------
        speciesdict: dict of homologies and their relevant dataframe.
        """
        matrix = np.zeros((len(self.all_homoset.species_namelist), len(speciesdict)), dtype=np.bool)
        matrix = np.hstack((matrix, self.all_homoset.hashomo_matrix))
        i = self.all_homoset.hashomo_matrix.shape()[0]
        dou = 0
        # COMPLEXITY
        print "start the iteration process, I hope you haven't clusterized your data yet..\
        else it won't work (for now)"
        for key, species in speciesdict.iteritems():  # we go throught the list of lists
            j = 0
            print "at homology " + key + " we have " + str(len(species)) + " species possessing it"
            species.sort()
            nameb = ''
            for name in species:  # TODO :can be parallelized
                # the assumption here is that they are ordered the same so we should only test for
                # next ones
                if name == nameb:
                    dou += 1
                else:
                    # we go until the occurence of the name in self.all_homoset.species_namelist
                    while name != self.all_homoset.species_namelist[j]:
                        j += 1
                    matrix[j, i] = True
                    # part where we update the Espece object
                    d = self.all_homoset.homodict[key].full.loc[[name], :]
                    d = d.fillna(0.5)   # we make sure to keep only one value. (nameA nameB problem)
                    d = d.iloc[[0]]
                    d = d.rename(index={name: key})  # need that weird stuff because of the transpose
                    self.species[name].genes = self.species[name].genes.append(d)
                    j += 1
                    nameb = name

            i += 1
        self.all_homoset.hashomo_matrix = matrix
        print "you had " + str(dou) + " same species homologies (it can't be processed! for now)"
        print "reviewed " + str(i) + " homologies and " + str(j) + "species"

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
                "all_homoset": self.all_homoset._dictify() if save_homo else False,
                "working_homoset": self.working_homoset._dictify() if save_workspace else False
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
        self.all_homoset = hset.HomoSet(data=data["all_homoset"])
        self.working_homoset = hset.HomoSet(data=data["working_homoset"])


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
