""""
Created by Jeremie KALFON
Date : 21 FEV 2018
University of Kent, ECE paris
jkobject.com

"""
import utils
import homoset as hset
import os
import numpy as np
import glob
import pandas as pd
import espece as spe


class PyCUB(object):
    """docstring for PyCUB

            Params:
            ------
            species: dictionary of Espece objects from the name of the species.
                                            (see espece.py)

    """

    species = {}
    is_preprocessed = False
    is_saved = False
    is_loaded = False
    have_homo = False
    matrix = None
    network = None
    homoset = None
    homo_namelist = None

    def __init__(self):
        """
        will..
        """

    def get_data(self):
        """
                Download the data from somewhere on the web
        """

    def load(self, folder='first500', fromYun=True, separation="homology"):
        """
                Get the data that is already present on a folder

                Params:
                ------
                from_Yun: if this flag is set to True it means that the folder is made of Yun's data
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
            self.homoset = hset.HomoSet()
            folder = "data/" + folder
            print "Reviewing all the " + str(len(os.listdir(folder))) + " files"
            for f in sorted(os.listdir(folder)):
                if f.endswith(".txt"):
                    nameA = f.split(separation)[0]
                    if(nameA != nameB):
                        nameB = nameA
                        homo_namelist.append(nameB)

            self.homo_namelist = homo_namelist
            for homology in homo_namelist:
                try:
                    genDF, species = utils.readcods_homology(separation, folder, homology)
                    speciesdict.update({homology: species})
                    homodict.update({homology: genDF})
                except OSError:
                    print "you do not have the files here"
                except:
                    print "problem with homology " + homology

            self.homoset.homodict = homodict
            self.have_homo = True
            self.is_loaded = True
            self._preprocess_yun(speciesdict)

    def preprocess(self):
        """
        will ...
        """

    def save(self):
        """
        call to save your work. you should call save on specific data structure if this is what you 
        want to save. 
        Will call other object's save, will transform all the variable into dict and save the dicts as 
        json files. will save the df also as json files. PyCUB and homoset have their own json file.

        Params:
        ------

        """

    def homologize(self, as_mat=True, plot=True):
        """
        Compute an homology group :
        from matrix computation 
        or from network computation:

        """

    def get_homoset(self, size, simi_align=True, max_clique=False, per_specie=True,
                    per_gene=True, name=[]):
        """
        Will compute either on the matrix or the network with different size and flags

        Params:
        -----
        size :  a value which will determine the size (max size of the retrived homology)
        50 means the first 50 in the sim
        if similarity alignment

        simi_align: flag for computing the matrix on similarity alignement

        max_clique: flag for computing on the max clique 

        Return:
        ------
        a HomoSet object (see homoset.py)
        """

    def _preprocess_yun(self, homomat):
        """
        a private function to preprocess specifically Yun's data

        Params:
        -------
        homomat: dict of homologies and their relevant dataframe.
        """
        df = utils.retrievelist()
        for i, row in df.iterrows():
            self.species.update({row['a']: spe.Espece(row['b'])})
        specieslist = df['a'].tolist()  # we create an ordered full list of species
        matrix = np.zeros((len(specieslist), len(homomat)), dtype=np.bool)
        i = 0
        dou = 0
        nameb = ''
        specieslist.sort()
        for key, species in homomat.iteritems():  # we go throught the list of lists
            j = 0
            print len(species)
            species.sort()
            for name in species:  # TODO :can be parallelized
                # the assumption here is that they are ordered the same so we should only test for the
                # next ones
                if name == nameb:
                    dou += 1
                else:
                    while name != specieslist[j]:
                        j += 1
                    matrix[j, i] = True
                    # part where we update the Espece object
                    Serie = self.homoset.homodict[key].loc[name, :]
                    df2 = pd.DataFrame(Serie).transpose()
                    df2 = df2.rename(index={1: 'Ahaha'})
                    self.species[name].genes.append(df2)
                    j += 1
                    nameb = name

            i += 1
        self.matrix = matrix
        print "you had " + str(dou) + "double homologies for the same species (it can't be processed)"
        print "num of homologies :" + str(i) + " per number of species : " + str(j)
