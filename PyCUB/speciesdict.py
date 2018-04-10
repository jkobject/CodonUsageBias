""""
Created by Jeremie KALFON
Date : 9 Apr 2018
University of Kent, ECE paris
jkobject.com

"""


import json
import os
import glob

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import espece as spe
import utils
import homology as h

import collections


class Speciesdict(collections.MutableMapping):


----- transfer to pyCUB - ---

"""docstring for Speciesdict

                Object that acts as a dictionary but is not as fast

                params:
                ------
                has_homo_matrix : a numpy boolean array that store the matrix of gene presence in species
                full_homo_matrix : a numpy array similar to has_homo_matrix
                                    but containing the codon entropy vectors instead
                homodict = dictionnary of dataframes of codon usage per species from homology names
                homo_namelist : list of all the homology names
    """
    path = 'utils/save/'
    is_loaded = False

    @staticmethod
    def __init__(self, session, *args, **kwargs):
        self.store = dict()
        self.update(dict(*args, **kwargs))
        self.path += session
        self.is_loaded = True

    def __getitem__(self, key):
        if key in self.store:
            return self.store[key]
        else:
            if os.

    def __setitem__(self, key, value):
        self.store[key] = value

    def __delitem__(self, key):
        del self.store[key]

    def __iter__(self):
        return iter(self.store)

    def __len__(self):
        return len(self.store)

    def _offload(self):

    def _load(self):

    def _dictify(self):
