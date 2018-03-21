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
            name : 
    """

    code = {}
    metadata = []

    def __init__(self, name=False, genes=False, data=False):

        if not (type(data) is bool):
            self.name = data["name"]
            self.genes = pd.DataFrame.from_dict(data["genes"])
        elif not (type(genes) is bool):
            self.genes = pd.DataFrame.from_dict(genes)
        else:
            self.genes = pd.DataFrame()
            self.name = name

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

    def _dictify(self):
        return {"name": self.name,
                "genes": self.genes.to_dict() if not (type(self.genes) is bool) else False,
                "code": self.code,
                "metadata": self.metadata}
