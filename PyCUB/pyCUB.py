""""
Created by Jeremie KALFON
Date : 21 FEV 2018
University of Kent, ECE paris
jkobject.com

"""
import os
import glob
import json

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import espece as spe
import homoset as hset
import utils
import homology as h


class PyCUB(object):
    """docstring for PyCUB

            Params:
            ------
            species: dictionary of Espece objects from the name of the species.
                                            (see espece.py)
            _is_preprocessed : trivial system only boolean
            _is_saved : trivial system only boolean
            _is_loaded : trivial system only boolean
            _have_homo : trivial system only boolean
            full_homo_matrix : a numpy boolean array that store the matrix of gene presence in species
            network :
            homoset : homoset object that store the all the homology
            homo_namelist : list of all the homology names

    """

    def __init__(self, species={}, _is_preprocessed=False, _is_saved=False,
                 _is_loaded=False, _have_homo=False,
                 working_homoset=False, all_homoset=False):
        """
        will..

        Params:
        ------
        species: dictionary of Espece objects from the name of the species.
                                        (see espece.py)
        homoset : homoset object that store the all the homology


        _is_preprocessed : trivial system only boolean
        _is_saved : trivial system only boolean
        _is_loaded : trivial system only boolean
        _have_homo : trivial system only boolean

        """
        self.species = species
        self.working_homoset = working_homoset
        self.all_homoset = all_homoset
        self._is_saved = _is_saved
        self._is_loaded = _is_loaded

    def get_data(self):
        """
                Download the data from somewhere on the web
        """

    def load(self, filename='first500', fromYun=True, separation="homology"):
        """
                Get the data that is already present on a filename

                Params:
                ------
                from_Yun: if this flag is set to True it means that the filename is made of Yun's data
                in which case we will create directly the homology map in the same time as the rest
                of the PyCUB object.

        """
        if fromYun:
            # then we process it according to how Yun Displays its data, trying to fill in
            # as much as we can
            homo_namelist = []
            homodict = {}
            speciesdict = {}
            nameB = ""
            self.all_homoset = hset.HomoSet()
            self.working_homoset = hset.HomoSet()
            filename = "data/" + filename
            print "Reviewing all the " + str(len(os.listdir(filename))) + " files"
            for f in sorted(os.listdir(filename)):
                if f.endswith(".txt"):
                    nameA = f.split(separation)[0]
                    if(nameA != nameB):
                        nameB = nameA
                        homo_namelist.append(nameB)

            self.working_homoset.homo_namelist = homo_namelist
            self.all_homoset.homo_namelist = homo_namelist
            for homology in homo_namelist:
                # try:
                genDF, specieslist = utils.readcods_homology(separation, filename, homology)
                speciesdict.update({homology: specieslist})
                homodict.update({homology: h.homology(full=genDF)})
                # except OSError:
                #    print "you do not have the files here"
               # except:
               #     print "problem with homology " + homology

            self.working_homoset.homodict = homodict
            self.all_homoset.homodict = homodict
            self._have_homo = True
            self._is_loaded = True
            self._preprocess_yun(speciesdict)
            self._is_preprocessed = True
        else:
            filename = "/utils/save/" + filename + ".json"
            with open(filename, "r") as f:
                jsoni = f.read()
                print "loading from " + name
                self._undictify(json.loads(jsoni))
                print "it worked !"

    def preprocess(self):
        """
        will ...
        """

    def save(self, name):
        """
        call to save your work. you should call save on specific data structure if this is what you
        want to save.
        Will call other object's save, will transform all the variable into dict and save the dicts as
        json files. will save the df also as json files. PyCUB and homoset have their own json file.

        Params:
        ------

        """
        if self._is_preprocessed:
            filename = "/utils/save/" + name + ".json"
            print "writing in " + name
            dictify = self._dictify()
            data = json.dump(dictify, indent=4, separators=(',', ': '))
            with open(filename, "w") as f:
                f.write(data)
                print "it worked !"

    def _preprocess_yun(self, speciesdict):
        """
        a private function to preprocess specifically Yun's data

        Params:
        -------
        speciesdict: dict of homologies and their relevant dataframe.
        """
        df = utils.retrievelist()
        for i, row in df.iterrows():
            self.species.update({row['a']: spe.Espece(row['b'])})
        specieslist = df['a'].tolist()  # we create an ordered full list of species
        matrix = np.zeros((len(specieslist), len(speciesdict)), dtype=np.bool)
        i = 0
        dou = 0
        nameb = ''
        specieslist.sort()

        for key, species in speciesdict.iteritems():  # we go throught the list of lists
            j = 0
            print "at " + key + " we have " + str(len(species)) + " homologies"
            species.sort()
            for name in species:  # TODO :can be parallelized
                # the assumption here is that they are ordered the same so we should only test for the
                # next ones
                if name == nameb:
                    dou += 1
                   # print name + " has multiple occurences of the gene in " + key + "cannot work with that now"
                else:
                    while name != specieslist[j]:  # we go until the occurence of the name in specieslist
                        j += 1
                    matrix[j, i] = True
                    # part where we update the Espece object
                    Serie = self.all_homoset.homodict[key].full.loc[name, :]
                    df2 = pd.DataFrame(Serie).transpose()
                    print df2
                    df2 = df2.rename(index={1: ' '})  # need that weird stuff because of the transpose...
                    print df2
                    self.species[name].genes.append(df2)
                    j += 1
                    nameb = name

            i += 1
        # create DF here with specieslist and hset.homo_namelist - ---
        dftemp = pd.DataFrame(matrix, specieslist, self.all_homoset.homo_namelist)
        self.all_homoset.has_homo_matrix = dftemp
        self.working_homoset.has_homo_matrix = dftemp
        print "you had " + str(dou) + " same species homologies (it can't be processed!)"
        print "reviewed " + str(i) + " homologies and " + str(j) + "species"

    def _dictify(self):
        dictispecies = {}
        for key, val in self.species.iteritems():
            dictispecies.update({key: val._dictify()})
        return {"species": dictispecies,
                "all_homoset": self.all_homoset._dictify(),
                "working_homoset": self.working_homoset._dictify(),
                }

    def _undictify(self, data):
        """
        same function but to retransform everything 
        Here we don't use other classes undictify functions but we just recreate them by passing it to 
        their init methods which is clearer. 
        """
        species = {}
        for key, val in data["species"].iteritems():
            species.update({key: spe.Espece(data=val)})
        self.species = species
        self._is_saved = False
        self._is_loaded = True
        self.all_homoset = hset.HomoSet(data=data["all_homoset"])
        self.working_homoset = hset.HomoSet(data=data["working_homoset"])
