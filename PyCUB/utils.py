""""
Created by Jeremie KALFON
Date : 21 FEV 2018
University of Kent, ECE paris
jkobject.com

"""
import pandas as pd
import glob
import numpy as np
import requests
import sys
import os
import pdb
try:
    from urllib2 import urlopen as urlopen
except:
    from urllib.request import urlopen as urlopen

import homology as h

speciestable = {}


def homoyun(i, homology, separation, filename, species_namelist, by):
    matrix = np.zeros(len(species_namelist), dtype=np.bool)
    try:
        gentab, specieslist, st, lenmat = readcods_homology(separation, filename, homology, by=by)
        # replace preprocess yun
    except OSError:
        print "you do not have the files here"
    espe = ''
    j = 0
    doub = 0
    names = []
    print "at homology " + homology + " we have " + str(len(specieslist)) + \
        " species possessing it"
    for species in specieslist:  # we can go faster as it is sorted
        # creating a list of number linked to our speciestable

        # removing doublon
        if species == espe:
            names.append(names[-1])
            doub += 1
        else:
            while species != species_namelist[j]:
                j += 1
            names.append(j)
            matrix[j] = True
            espe = species
            j += 1
    return [{homology: h.homology(full=gentab, names=names, nans=st, lenmat=lenmat)}, matrix, doub]


def getyun(key, val):
    url = val
    print "downloading " + key + " with urllib"
    if not os.path.exists('utils/data/' + key):
        f = urlopen(url)
        data = f.read()
        with open('utils/data/' + key, "wb") as code:
            code.write(data)
    else:
        print "file is already there"


def mymeta(key, val):
    url = val
    print "downloading " + key + " with urllib"
    if not os.path.exists('utils/meta' + key):
        f = urlopen(url)
        data = f.read()
        with open('utils/meta' + key, "wb") as code:
            code.write(data)


def readcods_homology(separation, folder="first500", homo_name="YAL019W", by='entropyLocation', aminonb=18):
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
    # We want to get the basic information from the files
    # so we read one and extract them
    first_file = glob.glob(folder + "/" + homo_name + separation + "His.*")[0]
    # print "looking at homology :" + homo_name
    meta = pd.read_csv(first_file)
    rows = meta.shape[0]
    gentab = np.zeros((rows, aminonb))
    lenmat = np.zeros((rows, aminonb))
    st = np.zeros(rows, dtype=np.bool)
    i = 0
    for file in glob.glob(folder + "/" + homo_name + separation + "*.*"):
        if file[-7:-4] != 'ror':  # TODO: write if it belongs to a list of amino
            # acid (given by the user)
            # we change the nan values as .5 as it is the mean of our distribution
            # and we don't want to bias it

            struct = pd.read_csv(file).reset_index().sort_values(by='species')
            st += struct[by].isna().values
            struct = struct.fillna(0.5)
            lenmat[:, i] = struct['length'].values.tolist()
            gentab[:, i] = struct[by].values.tolist()
            species = struct['species'].values.tolist()
            i += 1
    return gentab, species, st, lenmat


def loadfromensembl(homology=None, sequence=None):
    """


    Params:
    ------
    homology: list of String of homology you want to retrieve
    """
    server = "http://rest.ensemblgenomes.org"
    if homology is not None:
        ext = "/homology/id/"
        for homo in homology:
            ext += homo + '?'
            if sequence is not None:
                # dna   cdna    cds ncrna   Protein EMBL    GENBANK MySQL   TSV GTF GFF3
                ext += 'sequence=' + sequence
            r = requests.get(server + ext, headers={"Content-Type": "application/json"})
            if not r.ok:
                r.raise_for_status()
                sys.exit()
            data = r.json()
            with open('utils/data/' + homo + '.json', "wb") as code:
                code.write(data)
    # TODO: find what should be interesting to get on the ensembl websiste here
    # http://rest.ensemblgenomes.org/
    # TODO: find where the GFF3 (gene annotation) files are


def retrievelist():
    return pd.read_csv("utils/meta/order_name461.csv", header=None, names=['name', 'b']),\
        pd.read_csv("utils/meta/names_with_links.csv", header=None, names=['name', 'b'])
