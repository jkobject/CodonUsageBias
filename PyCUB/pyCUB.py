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
from sklearn.cluster import SpectralClustering
import json


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

    def __init__(self, species = {}, _is_preprocessed = False, _is_saved = False, 
        _is_loaded = False, _have_homo = False, full_homo_matrix = None, network = None,
        homoset = None, homo_namelist = []):
        """
        will..

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
        self.species = species
        self._is_preprocessed = _is_preprocessed
        self._is_saved = _is_saved
        self._is_loaded = _is_loaded
        self._have_homo = _have_homo
        self.full_homo_matrix = full_homo_matrix
        self.network = network
        self.homoset = homoset
        self.homo_namelist = homo_namelist

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
            self.homoset = hset.HomoSet()
            filename = "data/" + filename
            print "Reviewing all the " + str(len(os.listdir(filename))) + " files"
            for f in sorted(os.listdir(filename)):
                if f.endswith(".txt"):
                    nameA = f.split(separation)[0]
                    if(nameA != nameB):
                        nameB = nameA
                        homo_namelist.append(nameB)

            self.homo_namelist = homo_namelist
            for homology in homo_namelist:
                try:
                    genDF, species = utils.readcods_homology(separation, filename, homology)
                    speciesdict.update({homology: species})
                    homodict.update({homology: genDF})
                except OSError:
                    print "you do not have the files here"
                except:
                    print "problem with homology " + homology

            self.homoset.homodict = homodict
            self._have_homo = True
            self._is_loaded = True
            self._preprocess_yun(speciesdict)
            self._is_preprocessed = True
        else:
            filename = "/utils/save/" + filename + ".json"
            with open(filename, "r") as f:
                jsoni = f.read()
                print "loading from " + name
                _undictify(self, json.loads(jsoni))
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
        filename = "/utils/save/" + name + ".json"
        print "writing in " + name
        dictify = _dictify(self)
        data = json.dump(dictify, indent=4, separators=(',', ': '))
        with open(filename, "w") as f:
            f.write(data)
            print "it worked !"

    def homologize_from_matrix(self, clustering='kmeans', plot=True, homogroupnb=2):
        """
        Compute an homology group :
        from matrix computation using the full_homo_matrix
        (or from network computation in homologize_from_network)

        Can be computed many times and will updata homoset with the most recent homoset found
        if homoset exists, it will save it. 

        Params:
        -------
        clustering: flags to 'kmeans', 'spectral', ... to use different sk-learn algorithms

        plot: flags to true for the function to output ploting of the affinity matrix with and without the 
        clusters
        
        homogroupnb: nb of groups you want to extract

        """
        if as_mat:
            mat = self.full_homo_matrix
            spect = SpectralClustering(n_clusters=homogroupnb, n_jobs=-1)
            spect.fit(mat)
            clust = SpectralClustering.labels_
            orderedmat = np.zeros(mat.shape)
            begin = 0
            for i in range(clusternb):
                a = clust.find(i)
                orderedmat[begin:begin + len(a)] = mat[a]
                begin += len(a)
            _have_homo = True
            self.homogroupnb = homogroupnb

            if plot:


            return affinity_matrix_

    def get_group(self, number = 1):
        """
        To use once you have clustered homology groups, else takes everything. 
        

        Params:
        ------
        number: set the number of the group you want to get need to be between 1 and homogroupnb
        """
        if :

    def get_homoset(self, size, clustering=True, clusternb=2, max_clique=False, per_specie=True,
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
        self.full_homo_matrix = matrix
        print "you had " + str(dou) + "double homologies for the same species (it can't be processed)"
        print "num of homologies :" + str(i) + " per number of species : " + rstr(j)

    def _dictify(self):
        dictispecies = {}
        dictispecies.update({key: val._dictify}) for key, val in self.species.iteritems()
        return {"species": dictispecies,
                "_is_preprocessed": self._is_preprocessed,
                "_is_saved": self._is_saved,
                "_is_loaded": self._is_loaded,
                "_have_homo": self._have_homo,
                "full_homo_matrix": self.full_homo_matrix.tolist(),
                # network
                "homoset": self.homoset.dictify(),
                "homo_namelist": homo_namelist
                }


    def _undictify(self, data):
        species = {}
        species.update({key: spe.Espece(val["name"], val["genes"])}) for key, val in data["species"].iteritems()
        self.species = species
        self._is_preprocessed = data["_is_preprocessed"]
        self._is_saved = data["_is_saved"]
        self._is_loaded = data["_is_loaded"]
        self._have_homo = data["_have_homo"]
        self.full_homo_matrix = np.array(data["full_homo_matrix"], dtype=np.bool)
        #data[" network"]
        self.homoset = hset.HomoSet(data=data["homoset"])
        self.homo_namelist = data["homo_namelist"]
