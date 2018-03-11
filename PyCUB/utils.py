""""
Created by Jeremie KALFON
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
    :param homo_name: the type of homology
     you are looking for
    :return two dictionaries : one being a set of dictionaries looking like the file we have
                                the other a dictionnary of homologous genes and their
                                18 components values for each species
            a list of species



    """
    gendict = {}
    # We want to get the basic information from the files
    # so we read one and extract them
    first_file = glob.glob(folder + "/" + homo_name + separation + "His.*")[0]
    # print "looking at homology :" + homo_name
    meta = pd.read_csv(first_file)
    rows = meta.shape[0]
    gentab = np.zeros((rows, aminonb))
    amino = []
    i = 0
    for file in glob.glob(folder + "/" + homo_name + separation + "*.*"):
        if file[-7:-4] != 'ror':  # TODO: write if it belongs to a list of amino
                                    # acid (given by the user)
            amino.append(file[-7:-4])
            # we change the nan values as .5 as it is the mean of our distribution
            # and we don't want to bias it
            struct = pd.read_csv(file).reset_index().fillna(0.5)
            gentab[:, i] = struct['entropyLocation'].values.tolist()
            species = struct['species'].values.tolist()
            i += 1
    genDF = pd.DataFrame(data=gentab, index=species, columns=amino)
    return genDF, species


def retrievelist():
    return pd.read_csv("data/first500/order_name461.csv", header=None, names=['a', 'b'])
