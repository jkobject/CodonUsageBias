""""
Created by Jérémie KALFON
Date : 21 FEV 2018
University of Kent, ECE paris
jkobject.com

"""
import pandas as pd
import glob
import numpy as np


def readcods_homology(separation, folder="first500", homo_name="YAL019W", aminonb=18):
    """
    read the all the files for one homology and goves uou

    :param folder: the folder you wanan look onto
    :param homoname: the type of gene you are looking for
    :return two dictionaries : one being a set of dictionaries looking like the file we have
                                the other a dictionnary of homologous genes and their
                                18 components values for each species
            a list of species



    """
    gendict = {}
    # We want to get the basic information from the files
    first_file = glob.glob(folder + "/" + gene + separation + "His.*")[0]
    print "looking at gene :" + gene
    meta = pd.read_csv(first_file).dropna()
    rows = meta.shape[0]
    gentab = np.zeros((aminonb, rows))
    amino = []
    i = 0
    for file in glob.glob(folder + "/" + gene + "*.*"):
        if file[-7:-4] != 'ror':  # TODO: write if it belongs to a list of amino
                                    # acid (given by the user)
            amino.append(file[-7:-4])
            struct = pd.read_csv(file).reset_index()
            gentab[:, i] = struct['entropylocation'].values.tolist()
            species = struct['species'].values.tolist()
            i += 1
    genDF = pd.DataFrame(data=gentab, index=species, column=amino)
    return genDF, species

    def retrievelist():
        return pd.read_csv("first500/order_name461", header=None, names=['a', 'b'])
