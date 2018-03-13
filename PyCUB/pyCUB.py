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


from sklearn.cluster import SpectralClustering
from sklearn.cluster import KMeans
from sklearn.metrics.pairwise import cosine_similarity
from scipy import sparse


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
                 _is_loaded=False, _have_homo=False, full_homo_matrix=None, network=None,
                 homoset=None, homo_namelist=[]):
        """
        will..

        Params:
        ------
        species: dictionary of Espece objects from the name of the species.
                                        (see espece.py)
        full_homo_matrix : a numpy boolean array that store the matrix of gene presence in species
        network :
        homoset : homoset object that store the all the homology
        homo_namelist : list of all the homology names

        _is_preprocessed : trivial system only boolean
        _is_saved : trivial system only boolean
        _is_loaded : trivial system only boolean
        _have_homo : trivial system only boolean

        """
        self.species = species
        self.full_homo_matrix = full_homo_matrix
        self.network = network
        self.homoset = homoset
        self.homo_namelist = homo_namelist

        self._is_preprocessed = _is_preprocessed
        self._is_saved = _is_saved
        self._is_loaded = _is_loaded
        self._have_homo = _have_homo

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
        if self._is_preprocessed:
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

        def compute(self, mat, homogroupnb, clust, species):
            orderedmat = np.zeros(mat.shape)
            ltemp = [0] * len(self.homo_namelist)
            begin = 0
            for i in range(homogroupnb):
                a = np.argwhere(clust == i)[:, 0]
                orderedmat[begin:begin + len(a)] = mat[a]
                ltemp[begin:begin + len(a)] = self.homo_namelist[a]
                begin += len(a)
            self.homo_namelist = ltemp
            # affinity_matrix_ =
            if plot:
                # the regular matrix
                plt.figure(figsize=(40, 40))
                plt.title('the regular matrix')
                # plt.imshow(np.matmul(cub.matrix.T ,cub.matrix))
                plt.imshow(mat)
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
                plt.imshow(orderedmat)
                plt.show()

                # affinity of the ordered matrix
                mat_sparse = sparse.csr_matrix(orderedmat)
                similarities = cosine_similarity(mat_sparse)
                plt.figure(figsize=(40, 40))
                plt.title('the affinity of the ordered matrix')
                plt.imshow(similarities)
                plt.show()

            return orderedmat

        if self._is_preprocessed:
            mat = self.full_homo_matrix
            if clustering == "spectral":
                spect = SpectralClustering(n_clusters=homogroupnb, n_jobs=-1)
                spect.fit(mat)
                clust = spect.labels_
                self.full_homo_matrix = compute(self, mat, homogroupnb, clust)
                self._have_homo = True
                self.homogroupnb = homogroupnb

            elif clustering == "kmeans":
                kmn = KMeans(n_clusters=homogroupnb, n_jobs=-1)
                kmn.fit(mat)
                clust = kmn.labels_
                self.full_homo_matrix = compute(self, mat, homogroupnb, clust)
                self._have_homo = True
                self.homogroupnb = homogroupnb
            else:
                print "you entered a wrong clustering algorithm"
                return False

    def get_group(self, number=1):
        """
        To use once you have clustered homology groups, else takes everything. 


        Params:
        ------
        number: set the number of the group you want to get need to be between 1 and homogroupnb
        """

        # if:

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

    def nb_homo_per_species(self):
        """
        """
        if self._is_preprocessed:
            sumed = np.sum(self.full_homo_matrix, axis=1)
            plt.figure(figsize=(40, 10))
            plt.title('number of homologies per species')
            plt.bar(range(len(sumed)), sumed)
            print "you can always look at a particular range of species with 'homo_namelist' "

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
                    Serie = self.homoset.homodict[key].loc[name, :]
                    df2 = pd.DataFrame(Serie).transpose()
                    df2 = df2.rename(index={1: ' '})  # need that weird stuff because of the transpose...
                    self.species[name].genes.append(df2)
                    j += 1
                    nameb = name

            i += 1
        self.homo_namelist = specieslist
        self.full_homo_matrix = matrix
        print "you had " + str(dou) + " same species homologies (it can't be processed!)"
        print "reviewed " + str(i) + " homologies and " + str(j) + "species"

    def _dictify(self):
        dictispecies = {}
        for key, val in self.species.iteritems():
            dictispecies.update({key: val._dictify})
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
        for key, val in data["species"].iteritems():
            species.update({key: spe.Espece(val["name"], val["genes"])})
        self.species = species
        self._is_preprocessed = data["_is_preprocessed"]
        self._is_saved = data["_is_saved"]
        self._is_loaded = data["_is_loaded"]
        self._have_homo = data["_have_homo"]
        self.full_homo_matrix = np.array(data["full_homo_matrix"], dtype=np.bool)
        # data[" network"]
        self.homoset = hset.HomoSet(data=data["homoset"])
        self.homo_namelist = data["homo_namelist"]
