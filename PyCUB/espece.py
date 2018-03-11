""""
Created by Jeremie KALFON 
Date : 21 FEV 2018
University of Kent, ECE paris
jkobject.com

"""
import pandas as pd


class Espece(object):
    """docstring for Espece


            Params:
            ------
            code: a dict from gene_name to dna_seq string
            metadata: a dict containing different metadata information
            genes:  a dataframe of genename * codons containing the entropy_location for each of them
    """

    code = {}
    metadata = []

    def __init__(self, name):
        self.name = name
        self.genes = pd.DataFrame()

    def plot():
        """
        will...
        """

    def compute_entropy_loc():
        """
        will..
        """
        pass

    def compute_all_genome_entropy():
        """
        will..
        """
        pass

    def save():
        """
        will..
        """
