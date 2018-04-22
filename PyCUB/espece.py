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
            name : 
            is_stored
    """

    code = {}
    metadata = {}
    is_stored = False
    name = None

    def __init__(self, name=None, genes=None, data=None, link=None):

        if data is not None:
            self.name = data["name"]
            self.metadata = data["metadata"]
            self.code = data["code"]
        else:
            self.name = name
            self.metadata.update({'link': link})

    def plot():
        """
        will...
        """

    def get_metadata():
        """
        get additional information on the species such as :
        its metadata from ENSEMBL
        its gff dict 

        """

    def get_data():
        """
        get its dna from ENSEMBL
        """

    def get_gff_info():
        """
        get its gff dict 
        """

    def compute_entropy_loc():
        """
        use the DNA string to compute the entropy loc using Yun's pipeline
        """
        pass

    def compute_all_genome_entropy():
        """
        use all the entropy locations and use Yun's way to find the full genome entropy
        """
        pass

    def _dictify(self):
        """
        Used by the saving function. transform the object into a dictionary that can be 
        json serializable
        """
        return {"name": self.name,
                "code": self.code,
                "metadata": self.metadata}
