""""
Created by Jérémie KALFON 
Date : 21 FEV 2018
University of Kent, ECE paris
jkobject.com

"""


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

    def plot():

    def compute_entropy_loc():
        pass

    def compute_all_genome_entropy():
        pass

    def save():
