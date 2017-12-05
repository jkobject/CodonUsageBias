import numpy as np
import scipy as sp
from sklearn import preprocessing as prep
from sklearn import manifold as man
from sklearn import cluster
import pandas as pd
import glob
import matplotlib.pyplot as plt
from bokeh.plotting import *
from bokeh.models import *





def readcods(folder ='kalfonDTA', gene="YAL019W", aminonb = 18, ):
    """
    read the things and gives you back a nice dict

    :param folder: the folder you wanan look onto
    :param type: the type of gene you are looking for
    :return: two dictionaries : one being a set of dictionaries looking like the file we have
                                the other a dictionnary of homologous genes and their
                                18 components values for each species
            a list of species

    dict
        gene
            amino
                species
                gene
                entropyvalue
                entropylocation
                length

    """
    gendict = {}
    meta = pd.read_csv(folder + "/" + gene+"homologyAla.txt").dropna()
    rows = meta.shape[0]
    species = meta['species'].values.tolist()
    gentab = np.zeros((aminonb,rows))
    i =0
    for file in glob.glob(folder + "/"+gene+"*.txt"):
        amino = file[-7:-4]
        if amino != 'ror':
            struct = pd.read_csv(file).dropna().reset_index()
            gentab[i,:] = struct['entropyLocation'].values
            i+=1
            gendict.update({amino:struct})
    return gentab,gendict,species


def sortrank (dict):
    "sort first according to the entropy value"
    dictsort = dict
    for type in dict:
        for codon in dict[type]:
            dictsort[type][codon] = dict[type][codon].sort_values(by=['entropyValue', 'entropyLocation'])
    return dictsort

def prepro(gene):
    """
    :param gene: array of homologous gene
    :return: normalized gene, tsne 30, tsne 40, scaled matrix, scaling matrix
    """
    gene = np.swapaxes(gene, 0, 1)
    stdscal = prep.StandardScaler().fit(gene)
    scaled = stdscal.transform(gene)
    tsnedstrong = man.TSNE(perplexity =40).fit_transform(scaled)
    tsned = man.TSNE(perplexity=30).fit_transform(scaled)
    return gene,tsned,tsnedstrong,scaled,stdscal

def clustgenes(gene):
    "cluster by all numerical values "
    kmean = cluster.MiniBatchKMeans(n_clusters=10)
    centroids = kmean.fit(gene)
    return centroids


def clustK(dict):
    "cluster by all numerical values "
    dictclust = dict
    for type in dict:
        for codon in dict[type]:
            kmean = cluster.MiniBatchKMeans(n_clusters=5)
            dictclust[type][codon] = kmean.fit(dict[type][codon].as_matrix()[:,2:5])
    return dictclust

def sortspecies(dict):
    """
        "sort according to species instead of codons"

    species
        gene+codon
            values


    """
    species = {}
    array = np.ndarray
    nspe = dict.values()[0].values()[0].count()[1]
    for it in range(nspe):
        countb = 0
        for type in dict:
            counta = 0
            for codon in dict[type]:
                a = dict[type][codon]
                b = a.iloc[[it]]
                array[countb][counta]= b.as_matrix()[2:5]
                ++counta
        ++countb
        species.update({dict[type][codon]['species'][it]: arraytot})

    return species

def bokeplot(tsneval,species, getimage=True):
    source = ColumnDataSource(
        data=dict(
            x=tsneval[:,0],
            y=tsneval[:,1],
            label=["species : %s" % (x_) for x_ in species]
        )
    )
    hover = HoverTool(tooltips=[
        ("label", "@label"),
    ])
    p = figure(title="T-sne of homologous gene X for each species",
               tools=[hover,BoxZoomTool(),WheelZoomTool(),SaveTool(),ResetTool()])
    p.circle('x', 'y', size= 10, source=source, color="#2222aa")
    if(getimage):
        output_file("bokeplotlast.html")
    show(p)
    return p