""""
Created by Jeremie KALFON
Date : 21 FEV 2018
University of Kent, ECE paris
jkobject.com

"""
import pandas as pd
import glob
import numpy as np


class homology(object):
    """in homology we store an homology with its dataframe, 
    it reduced reduced df with dim reduction and its clusters
    it is supposed to be store in a homoogy dictionary 

    Params:
    ------
    full : df of one homology with entropy value vector per species
    reduced :  df of one homology dimensionality reduced
    clusters : list of cluster labels
    """

    reduced = False
    clusters = False
    full = False

    def __init__(self, data=False, full=False):
        """
        will..
        """
        if not (type(data) is bool):
            self.reduced = pd.from_dict(data["reduced"]) if data["reduced"] else False,
            self.clusters = data["clusters"]
            self.full = pd.from_dict(data["full"] if data["full"] else False)

        elif not (type(full) is bool):
            self.full = full

    def clusterize_():
        """
        will..
        """

    def plot():
        """
        will..
        """

    def _dictify(self):
        """
        will..
        """
        return {"reduced": self.reduced.to_dict() if self.reduced else False,
                "clusters": self.clusters,
                "full": self.full.to_dict() if self.full else False}
