""""
Created by Jérémie KALFON
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
    specieslist = []
    homoset = None

    def __init__(self, arg):
        self.arg = arg

    def get_data():
        """
                Download the data from somewhere on the web
        """

    def load(folder='first500', fromYun=True, separation="homology"):
        """
                Get the data that is already present on a folder

               	Params:
               	------
               	from_Yun: if this flag is set to True it means that the folder is made of Yun's data
               	in which case we will create directly the homology map in the same time as the rest
               	of the PyCUB object.

        """
        if fromYun:  # then we process it according to how Yun Displays its data, trying to fill in
        			# as much as we can
        	homo_namelist = []
        	homodict = {}
        	self.homoset = homoset.HomoSet()
        	print "Reviewing all the " + len(os.listdir(folder)) + " files"
            for f in os.listdir(folder):
                nameA = f.split(separation)[0]
                if(nameA != nameB):
                    nameB = nameA
                    homo_namelist.append(nameB)

            for homology in homo_namelist:
	            try:
	                genDF, species = utils.readcods_gene(separation,folder,homology)
	                self.specieslist.append(species) 
	               	homolist{homology: genDF}
	            except OSError:
	                print "you do not have the files here"

	        self.homoset.homodict = homodict 
            self.have_homo = True
            self.is_loaded = True 
            _preprocess_yun()

    def preprocess():

    def save():

    def homologize(as_mat=True, plot=True):
        """
        Compute an homology group :
        from matrix computation 
        or from network computation:

        """

    def get_homoset(size, simi_align=True, max_clique=False, per_specie=True,
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





    def _preprocessYun(self):
    	speciesmat = self.specieslist
    	df = utils.retrievelist()
    	for row in df.rows
    		species.update({row['b'] : spe.Espece(row['a']) })
    	self.specieslist = df['a'].tolist()
    	matrix = np.zeros((self.specieslist.shape, speciesmat.shape[1]), dtype=np.bool)
