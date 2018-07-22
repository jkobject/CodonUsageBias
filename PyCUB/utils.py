""""
Created by Jeremie KALFON
Date : 21 FEV 2018
University of Kent, ECE paris
jkobject.com

"""
import pandas as pd
from scipy.stats import multinomial
from scipy.stats import multivariate_normal
from numpy.random import randint
from math import log
import glob
import numpy as np
import requests
import os
import pdb
import json
try:
    from urllib2 import urlopen as urlopen
except:
    from urllib.request import urlopen as urlopen

import homology as h

"""
    the utils is where utilitary functions are stored (in particular preprocessing functions here)

    Values:
    ------
    speciestable: dict[int] will be used by all other objects to make numbers corresponds to species
    codons: dictionaary[list[str]] making codons corresponds to amino acids,
    amino: list[str] amino acids in their alphabetical ordering
    aminosort: list[str] amino acids sorted by the number of their codons and the alphabetical orderings
    of their codons.
    MAXITR: for computational purposes when you reach a certain number of iteration in computing partition functions
            in entropy locations, you stop as you have a good enough partition function
    LMAX: int the maximal length of amino acid in some functions, it will do a full partition function computation
            if lower than this value. all the values here are choosen for a max amount of computation in full partition
            function of MAXITR

    smax: list[int]
            other values to reach 1 million iterations for each codon number
    smax2: list[int]
            other values to reach 1 million iterations for each codon number when similarity is not counted when
            computing the partition function.

"""
speciestable = {}
phylo_distances = None  # pandas dataframe
codons = {
    'ALA': ['GCA', 'GCC', 'GCG', 'GCT'],
    'ARG': ['CGA', 'CGC', 'CGG', 'CGT', 'AGG', 'AGA'],
    'ASN': ['AAC', 'AAT'],
    'ASP': ['GAT', 'GAC'],
    'CYS': ['TGT', 'TGC'],
    'GLN': ['CAA', 'CAG'],
    'GLU': ['GAG', 'GAA'],
    'GLY': ['GGT', 'GGG', 'GGA', 'GGC'],
    'HIS': ['CAT', 'CAC'],
    'ILE': ['ATC', 'ATA', 'ATT'],
    'LEU': ['TTA', 'TTG', 'CTC', 'CTT', 'CTG', 'CTA'],
    'LYS': ['AAA', 'AAG'],
    # 'MET': ['ATG'],
    'PHE': ['TTT', 'TTC'],
    'PRO': ['CCT', 'CCG', 'CCA', 'CCC'],
    'SER': ['TCT', 'TCG', 'TCA', 'TCC', 'AGC', 'AGT'],
    # 'STOP': ['TAG', 'TGA', 'TAA'],
    'THR': ['ACC', 'ACA', 'ACG', 'ACT'],
    # 'TRP': ['TGG'],
    'TYR': ['TAT', 'TAC'],
    'VAL': ['GTA', 'GTC', 'GTG', 'GTT']}
anticodons = {
    'ALA': ["CGA", "CGG", "CGT", "CGC"],
    'ARG': ["GCA", "GCG", "GCT", "GCC", "TCT", "TCC"],
    'ASN': ["TTA", "TTG"],
    'ASP': ["CTA", "CTG"],
    'CYS': ["ACA", "ACG"],
    'GLN': ["GTT", "GTC"],
    'GLU': ["CTT", "CTC"],
    'GLY': ["CCA", "CCG", "CCT", "CCC"],
    'HIS': ["GTA", "GTG"],
    'ILE': ["TAA", "TAG", "TAT"],
    'LEU': ["AAT", "AAC", "GAA", "GAG", "GAT", "GAC"],
    'LYS': ["TTT", "TTC"],
    # 'MET': ['ATG'],
    'PHE': ["AAA", "AAG"],
    'PRO': ["GGA", "GGG", "GGT", "GGC"],
    'SER': ["AGA", "AGG", "AGT", "AGC", "TCA", "TCG"],
    # 'STOP': ['TAG', 'TGA', 'TAA'],
    'THR': ["TGA", "TGG", "TGT", "TGC"],
    # 'TRP': ['TGG'],
    'TYR': ["ATA", "ATG"],
    'VAL': ["CAA", "CAG", "CAT", "CAC"]}

amino2reduce = {"A": 'ALA',
                "R": 'ARG',
                "N": 'ASN',
                "D": 'ASP',
                "C": 'CYS',
                "Q": 'GLN',
                "E": 'GLU',
                "G": 'GLY',
                "H": 'HIS',
                "I": 'ILE',
                "L": 'LEU',
                "K": 'LYS',
                "M": 'MET',
                "F": 'PHE',
                "P": 'PRO',
                "S": 'SER',
                "T": 'THR',
                "W": 'TRP',
                "Y": 'TYR',
                "V": 'VAL'}

amino = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS',
         'ILE', 'LEU', 'LYS', 'PHE', 'PRO', 'SER', 'THR', 'TYR', 'VAL']
aminosort = ['LYS', 'ANS', 'THR', 'ILE', 'ARG', 'GLN', 'HIS', 'PRO', 'GLU', 'ASP',
             'ALA', 'GLY', 'VAL', 'TYR', 'SER', 'LEU', 'CYS', 'PHE']
LMAX = 2500
MAXITR = 1000000
smax = [0, 0, LMAX, 1414, 180, 0, 40]
smax2 = [0, 0, LMAX, LMAX, 520, 0, 145]
MAX = 20  # max value of leng after which we encounter
# weird results in compute jerem's computation nbcod**leng


def homoyun(separation, folder="first500", homo_name="YAL019W",
            by='entropyLocation', aminonb=18, aminolist=None):
    """
    read the all the files for one homology and returns everything inside (Yun)

    Params:
    ------
    folder: the folder you wanan look onto
    homo_name: the type of homology
     you are looking for
    by:
    aminonb
    Returns:
    --------
    list[df,df,df,df,list[bool]]: the list of df of 'by' values,  df of species, df of places where there is
    nan values, df of lengths of genes, location where there is doublons is yun's data (copy number of the gene)

    """
    # We want to get the basic information from the files
    # so we read one and extract them
    # first_file = glob.glob(folder + "/" + homo_name + separation + "His.*")[0]
    # print "looking at homology :" + homo_name
    aminolist = amino if aminolist is None else aminolist
    struct = pd.DataFrame({'': []})
    for i, file in enumerate(glob.glob(folder + "/" + homo_name + separation + "*.*")):
        if file[-7:-4] != 'ror' and file[-7:-4].upper() in aminolist:
            if i == 0:
                struct = pd.read_csv(file).sort_values(by='species').reset_index()
            else:
                struct = pd.concat([struct, pd.read_csv(file).sort_values(by='species').reset_index()], axis=1)
    lenmat = struct['length'].fillna(0).values.astype(int)
    nanvalues = struct[by].isna().values
    gentab = struct[by].fillna(0.5).values
    species = struct['species'].values.transpose()[0]  # maybe we are missing some by only taking one
    doub = [False]
    doub.extend([False if species[i] is not species[i - 1] else True for i in range(1, len(species))])
    print "at homology " + homo_name
    return [gentab, species, nanvalues, lenmat, doub]


def getyun(key, val):
    """
    function to get data from an url 'val' and save it under a file 'key'
    """
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
    """
    same as getyun but saves under meta
    """
    url = val
    print "downloading " + key + " with urllib"
    if not os.path.exists('utils/meta/' + key):
        f = urlopen(url)
        data = f.read()
        with open('utils/meta/' + key, "wb") as code:
            code.write(data)
            print "downloaded"
    else:
        print key + " already present"


def retrievenames():
    """
    returns the dfs for yun's data (names of the species alphabetically ordered and their links)
    """
    return pd.read_csv("utils/meta/order_name461.csv", header=None, names=['name', 'b']),\
        pd.read_csv("utils/meta/names_with_links.csv", header=None, names=['name', 'b'])


def loadfromensembl(homology, kingdom='compara=fungi', sequence='cdna',
                    additional='type=orthologues', saveonfiles=False, normalized=False,
                    setnans=False, number=0, by="entropy", using="normal"):
    """
    Load from ensembl the datas required in parameters ( look at PyCUB.get_data for more information)
    returns a fully populated homology object.
    """
    server = "http://rest.ensemblgenomes.org"
    print 'homology: ' + homology + ' : ' + str(number)
    ext = "/homology/id/" + homology + '?'
    if sequence is not None:
        # dna   cdna    cds ncrna   Protein EMBL    GENBANK MySQL   TSV GTF GFF3
        ext += 'sequence=' + sequence
    if kingdom is not None:
        ext += ';' + kingdom
    if additional is not None:
        ext += ';' + additional
    try:
        r = requests.get(server + ext, headers={"Content-Type": "application/json"})
    except ConnectionError:
        print "problem at " + homology
        return loadfromensembl(homology, kingdom=kingdom, sequence=sequence,
                               additional=additional, saveonfiles=saveonfiles, normalized=normalized,
                               setnans=setnans, number=number, by=by, using=using)
    if not r.ok:
        r.raise_for_status()
    data = r.json()['data'][0]['homologies']
    if saveonfiles:
        with open('utils/data/' + homology + '.json', "wb") as code:
            code.write(json.dump(data))
    species, GCcount, lenmat, H, nans, similarities, KaKs_Scores, taxons, proteinids = process(
        data, normalized=normalized, setnans=setnans, by=by)
    if by == 'entropyLocation':
        H = getloc(H, np.array(lenmat), using=using)
    # here we add two things into names but only as a temporary saving measures removed by the
    # application fo preprocessing in homoset.
    homo = h.homology(names=[species, taxons], full=H, lenmat=lenmat,
                      nans=nans, KaKs_Scores=KaKs_Scores, similarity_scores=similarities,
                      proteinids=proteinids, GCcount=GCcount)
    homo.order(withtaxons=True)  # a first ordering of the data, usefull afterward in the preprocessing
    return homo


# ################# YUN Processing Pipeline  ####################################


def process(data, normalized=False, setnans=False, by='entropy'):
    """
    function used by loadfromensembl() to process the retrieved data
    """
    species = []
    H = []
    nans = []
    lenmat = []
    KaKs_Scores = []
    similarities = []
    taxons = []
    proteinids = []
    GCcounts = []
    for dat in data:  # each species
        # https://en.wikipedia.org/wiki/Ka/Ks_ratio
        if dat["dn_ds"] is not None:
            KaKs_Scores.append(dat["dn_ds"])
        dat = dat['target']
        if dat["perc_id"] is not None:
            similarities.append(dat["perc_id"])
        if dat["taxon_id"] is not None:
            taxons.append(dat["taxon_id"])
        if dat["protein_id"] is not None:
            proteinids.append(dat["protein_id"])
        species.append(dat['species'])
        valH, len_i, nan, GCcount = computeyun(dat['align_seq'].encode('ascii', 'ignore').replace("-", ""),
                                               normalized=normalized, setnans=setnans, by=by)
        H.append(valH)
        nans.append(nan)
        lenmat.append(len_i)
        GCcounts.append(GCcount)
    return species, np.array(GCcounts, dtype=int), np.array(lenmat, dtype=int), np.array(H),\
        np.array(nans), KaKs_Scores if len(KaKs_Scores) != 0 else None,\
        similarities if len(similarities) != 0 else None,\
        taxons if len(taxons) != 0 else None,\
        proteinids if len(proteinids) != 0 else None


def computeyun(data, setnans=False, normalized=False, by='entropy'):
    """
    function used by process to compte the 'entropy' or 'frequency' of the processed data
    """
    c = [data[i:i + 3] for i in range(0, len(data), 3)]
    valH = np.zeros(len(amino)) if by != 'frequency' else np.zeros(59)  # the number of codons usefull
    len_i = []
    nans = 0
    pos = 0
    GCcount = (data.count('G') + data.count('C')) / len(data)
    for k, amin in enumerate(amino):
        nbcod = len(codons[amin])  # replace Cleng
        count = np.zeros(nbcod)
        X = np.zeros(nbcod)
        mn = np.ones(nbcod) / nbcod
        for i, cod in enumerate(codons[amin]):
            for j, val in enumerate(c):
                if val == cod:
                    count[i] += 1
                    c.pop(j)
        lengsubseq = count.sum()  # replace subSlength
        if by == 'frequency':
            if lengsubseq == 0:
                valH[pos:pos + nbcod] = np.NaN if setnans else 1. / nbcod
                nans += 1
            else:
                E = count / lengsubseq
                valH[pos:pos + nbcod] = E
            pos += nbcod
        else:
            if lengsubseq == 0:
                valH[k] = np.NaN if setnans else 0.5
                nans += 1
            else:
                Yg = multinomial.pmf(x=count, n=lengsubseq, p=mn)
                # efor part
                div, i = divmod(lengsubseq, nbcod)
                X[:int(i)] = np.ceil(div) + 1
                X[int(i):] = np.floor(div)
                Eg = multinomial.pmf(x=X, n=lengsubseq, p=mn)
                # end here
                valH[k] = -np.log(Yg / Eg) / lengsubseq if normalized else -np.log(Yg / Eg)
        len_i.append(lengsubseq)
    return valH, len_i, nans, GCcount


def getloc(valH, geneleng, using='normal'):
    """
    the function to compute the entropy location (adpated from Yun Deng's code University of Kent 2018)
    need the entropy values from the genes an homology and the length values as well.

    Params:
    ------

    valH
    """
    # TODO: if you concatenate more data together, it should be more efficient
    # then you do for x in valH and you get the size of each and concatenate, order, do the computation
    # then you retrieve the list back with the known ordering and sizes
    numberindiv = len(valH)  # works like a shape[0]
    acodon = [[2, 3, 4, 5, 6, 8, 11, 12, 16],
              [9], [0, 7, 13, 15, 17], [1, 10, 14]]
    val = np.zeros((18, numberindiv))
    for ind, nbcod in enumerate([2, 3, 4, 6]):
        lengs = []
        for i in acodon[ind]:
            lengs.extend(geneleng[:, i])
        # we get a series of length for each species
        # for each amino acid of this codon length
        # we have all 2 then all 3 then all ..
        lengs = np.array(lengs)
        # we sort by length
        argu = lengs.argsort()
        lengs = lengs[argu]
        lt = -1
        # will possibly find similar lengths values and thus
        # reduce computation time
        valHloc = np.zeros(len(lengs))
        for y, leng in enumerate(lengs):
            posx, posy = divmod(argu[i], numberindiv)
            # we access the position knowing that we need to convert row number
            # into the possible position on the entropy matrix
            E = valH[posy, acodon[ind][posx]]
            # if H is also superior or inferior to the threshold, according to the
            # computation of the entropy location, both their values should be the
            # same.
            if lt != leng:
                lt = leng
                mn = np.ones(nbcod) / nbcod
                X = np.zeros(nbcod)
                i, div = divmod(leng, nbcod)
                X[:div] = np.ceil(i) + 1
                X[div:] = np.floor(i)
                Eg = multinomial.pmf(x=X, n=leng, p=mn)
                if Eg == 0:
                    valHloc[y] = 0
                    continue
                ref = computepartition(nbcod, leng, 1000000, using=using)
                try:
                    ref = np.divide(np.log(np.divide(Eg, ref)), leng)
                    hist, edges = np.histogram(ref, int((ref.max() - ref.min()) * 10000))
                except:
                    pdb.set_trace()
                    valHloc[y] = 0
                    continue

                lhist = len(edges)
                if edges[0] < 0:
                    nhist = np.zeros(lhist - 2)
                    nhist[0] = np.sum(hist[:2])
                    nhist[1:lhist - 2] = hist[2:lhist - 1]
                    edgesnew = np.zeros(lhist - 2)
                    edgesnew[:lhist - 2] = edges[1:lhist - 1]
                    hist = nhist
                else:
                    edgesnew = np.zeros(lhist - 1)
                    edgesnew[:] = edges[:lhist - 1]

                hist = np.cumsum(np.multiply(hist, np.divide(
                    Eg, np.exp(np.multiply(edgesnew, leng)))))
            if E > edgesnew[-1]:
                location = last_match_index(hist, 1)
            else:
                location = np.argmax(edgesnew >= E)
                if location != 1:
                    location_value = hist[location]
                    if location_value > 0.95 and location_value < 1:
                        location -= 1
                    if location_value >= 1:
                        location = last_match_index(hist, 1)
            valHloc[y] = hist[location]
        # we reset the right ordering
        valHloc[argu] = valHloc[:]
        for y, e in enumerate(acodon[ind]):
            # we reset the position of the aminoacid numbers.
            val[e, :] = valHloc[numberindiv * y:numberindiv * (y + 1)]
    return val.T


def computepartition(nbcod, leng, MAXITR, using='normal'):
    """
    according to a 'using' parameter, will select the partition function algorithm to use
        jeremcompute: bool to true if use a jeremie twist to this function
        all function are made and found by jeremie kalfon, derived from earlier
        work by Yun Deng Phd student at the University of Kent 2O18

    Basically works by computing every possible codon presence for a defined amino acid (and it number of
    repetitions)
    done by computing the pdf from a multinomial according to those values

    """
    if using == 'computejerem':
        return computejerem(nbcod, leng)
    elif using == 'permutation':
        return computepartition_without_permutation(nbcod, leng)
    elif using == 'full':
        return computepartition_sorted_full(nbcod, leng)
    elif using == 'normal':
        return computepartition_normal_approximation(nbcod, leng)
    else:
        return NameError(" give a Using from [permutation, random, normal, full]")


def randomdraw(nbcod, leng):
    """
    one of the four partition functions algorithms
    this one either use the regular Yun version or try to sort the values to limit the number
    of sampling to do and do a limited dynamic programming version.

    Basically works by computing a random subset of every possible codon presence
    for a defined amino acid (and it number of repetitions)
    done by computing the pdf from a multinomial according to those values
    """
    ind = 0
    listvect = []
    # we iterate this way to have only a max value or everything thanks to a
    # sampling without replacement
    while ind != MAXITR:
        prevect = [0] * nbcod
        prevect[1:] = list(np.sort(randint(0, leng + 1, nbcod - 1)))
        for i in range(nbcod - 1):
            prevect[i] = prevect[i + 1] - prevect[i]
        prevect[-1] = leng - prevect[-1]
        listvect.append(list(prevect))
        ind += 1
    listvect.sort()
    return multinomial(leng, np.ones(nbcod) / nbcod).pmf(listvect)


def computepartition_sorted_full(nbcod, leng):
    """
    one of the four partition functions algorithms
    Basically works by computing every possible codon presence for a defined amino acid (and it number of
    repetitions)
    done by computing the pdf from a multinomial according to those values
    """
    mn = np.ones(nbcod) / nbcod
    if mlen(leng, nbcod) == 'full':
        # if we are ok to do full method
        a = 0
        val = []
        for i in range(leng + 1):
            a = leng - i
            # complexity O(leng)
            if nbcod > 2:
                for j in range(a + 1):
                    b = a - j
                    # complexity O(x) x = ((leng**2)/2) - (leng/2)
                    if nbcod > 3:
                        # parallelize here
                        for k in range(b + 1):
                            c = b - k
                            # complexity O(x) x = ((leng**3)*(5/12)) - ((leng**2)*(3/4)) - (2*leng/3)
                            if nbcod > 4:
                                for l in range(c + 1):
                                    d = c - l
                                    for m in range(d + 1):
                                        e = d - m
                                        val.append([i, j, k, l, m, e])
                                        # unknown complexity
                            else:
                                val.append([i, j, k, c])
                    else:
                        val.append([i, j, b])
            else:
                val.append([i, a])
        return multinomial(leng, mn).pmf(val)
    else:
        return randomdraw(nbcod, leng)


def computepartition_without_permutation(nbcod, leng):
    """
    one of the four partition functions algorithms
    this one sorts the values and don't compute duplicates (that will end up with the same
    binomial pdf) to limit the number of sampling to do and do a limited dynamic programming version.

    Basically works by computing every possible codon presence for a defined amino acid (and it number of
    repetitions)
    done by computing the pdf from a multinomial according to those values
    """
    mn = np.ones(nbcod) / nbcod
    if mlen(leng, nbcod) == 'full':
        # if we are ok to do full method
        a = 0
        val = []
        for i in range((leng / 2) + 1):
            a = leng - i
            if nbcod > 2:
                for j in range(i, a / 2 + 1):
                    b = a - j
                    if nbcod > 3:
                        for k in range(j, b / 2 + 1):
                            c = b - k
                            if nbcod > 4:
                                for l in range(k, c / 2 + 1):
                                    d = c - l
                                    for m in range(l, d / 2 + 1):
                                        e = d - m
                                        val.append([i, j, k, l, m, e])
                                        # unknwon complexity
                            else:
                                val.append([i, j, k, c])
                    else:
                        val.append([i, j, b])
            else:
                val.append([i, a])
        densities = multinomial(leng, mn).pmf(val)
        values = list(val)
        for i, v in enumerate(val):
            positions = list(permute(v, 0, len(v) - 1))
            temp = []
            for position in positions:
                pos = list(position)
                if pos != v:
                    temp.append(pos)
            values.extend(temp)
            densities = np.append(densities, np.repeat(densities[i], len(temp)))
        ind = sorted(range(len(values)), key=values.__getitem__)
        densities[:] = densities[ind]

        return densities
    else:
        # we draw from a multinomial function
        return randomdraw(nbcod, leng)


def computejerem(nbcod, leng):
    """
    one of the four partition functions algorithms
    this ones uses the repetitions in the multinomial equation to provide the fastest way to compute every possible values
    in this multinomial equation
    moreover, it uses the symmetry given by the probability vector to not compute permutations and just
    copy and paste permutation values

    Basically works by computing every possible codon presence for a defined amino acid (and it number of
    repetitions)
    done by computing the pdf from a multinomial according to those values
    """
    if mlen(leng, nbcod) == 'full':
        # if we are ok to do full method

        # TODO: check if it is faster to replace by log computations
        # Here we have a problem if we compute directly with a range starting with zero as we will
        # end up dividing by zero, the best idea is to reuse computationas for 2/3/4 codons to compute
        # the first part at where one the the values [A,B,nC,D,E,.] is at zero
        def comp2(densities, comp, initialval, leng, val):
            densities.append(comp)
            val.append(list(initialval))
            for i in range(1, leng + 1):
                # Here a is the last value in [A,B,C,D,E,.] which can be computed as the total length
                # minus the other values. it is propagated into the other for loops as one can see for the
                # computations where nbcod > 2
                a = leng - i
                comp = (comp / i) * (a + 1.)
                densities.append(comp)
                initialval[-1] = a
                initialval[-2] = i
                val.append(list(initialval))

            return densities, val

        def comp3(densities, comp, initialval, leng, val):
            for i in range(1, leng + 1):
                a = leng - i
                comp = (comp / i) * (a + 1.)
                # here comp is the memory of what part of the multinomial equation is at that point
                # To go from the multinomial of [A,B,C,D,E,.] to the one of [A,B,C,D,E+1,.-1] it is
                # just a division and a multiplication so we iteratively divide and multiply
                initialval[-1] = a
                initialval[-2] = 0
                initialval[-3] = i
                dens, va = comp2([], comp, initialval, a, [])
                val.extend(va)
                densities.extend(dens)
            return densities, val

        def comp4(densities, comp, initialval, leng, val):
            for i in range(1, leng + 1):
                a = leng - i
                comp = (comp / i) * (a + 1.)
                initialval[-1] = a
                initialval[-2] = 0
                initialval[-3] = 0
                initialval[-4] = i
                dens, va = comp2([], comp, initialval, a, [])
                val.extend(va)
                densities.extend(dens)
                dens, va = comp3([], comp, initialval, a, [])
                val.extend(va)
                densities.extend(dens)
            return densities, val

        def comp6(densities, comp, initialval, leng, val):
            for j in range(1, leng + 1):
                a = leng - j
                comp = (comp / j) * (a + 1.)
                dens, va = comp2([], comp, [0, j, 0, 0, 0, a], a, [])
                densities.extend(dens)
                val.extend(va)
                dens, va = comp3([], comp, [0, j, 0, 0, 0, a], a, [])
                densities.extend(dens)
                val.extend(va)
                dens, va = comp4([], comp, [0, j, 0, 0, 0, a], a, [])
                densities.extend(dens)
                val.extend(va)
            for i in range(1, leng + 1):
                a = leng - i
                comp = (comp / i) * (a + 1.)
                dens, va = comp2([], comp, [i, 0, 0, 0, 0, a], a, [])
                densities.extend(dens)
                val.extend(va)
                dens, va = comp3([], comp, [i, 0, 0, 0, 0, a], a, [])
                densities.extend(dens)
                val.extend(va)
                dens, va = comp4([], comp, [i, 0, 0, 0, 0, a], a, [])
                densities.extend(dens)
                val.extend(va)
                comp_i = comp
                for j in range(1, a + 1):
                    b = a - j
                    comp_i = (comp_i / j) * (b + 1.)
                    dens, va = comp2([], comp_i, [i, j, 0, 0, 0, b], b, [])
                    densities.extend(dens)
                    val.extend(va)
                    dens, va = comp3([], comp_i, [i, j, 0, 0, 0, b], b, [])
                    densities.extend(dens)
                    val.extend(va)
                    dens, va = comp4([], comp_i, [i, j, 0, 0, 0, b], b, [])
                    densities.extend(dens)
                    val.extend(va)
            return densities, val

        if nbcod > 1:
            val = []
            densities = []
            initialval = [0] * (nbcod)
            initialval[-1] = leng
            dens, va = comp2([], 1., initialval, leng, [])
            val.extend(va)
            densities.extend(dens)
        if nbcod > 2:
            dens, va = comp3([], 1., [0] * (nbcod), leng, [])
            val.extend(va)
            densities.extend(dens)

        if nbcod > 3:
            dens, va = comp4([], 1., [0] * (nbcod), leng, [])
            val.extend(va)
            densities.extend(dens)

        if nbcod == 6:
            dens, va = comp6([], 1., [0] * (nbcod), leng, [])
            val.extend(va)
            densities.extend(dens)

        # in order not to have zeros that would block the iterations (keeping all values to zero)
        # we keep it high enough by not dividing it with the inverse probability which is very high
        # around 10^20
        densities = np.array(densities, dtype=np.float)
        if leng > MAX:
            val = leng % MAX
            for i in range(leng / MAX):
                densities = np.divide(densities, nbcod**MAX)
            if val != 0:
                densities = np.divide(densities, nbcod**val)
        else:
            densities = np.divide(densities, nbcod**leng)
        return densities
    else:
        # we draw from a multinomial function
        return randomdraw(nbcod, leng)


def computepartition_normal_approximation(nbcod, leng, probavector=None):
    """
    one of the 4 partition function algorithm, works by
    approximates the binomial distribution with a multidimensional gaussian function
    as a gaussian approximates a binomial under many trials and this is basically what we do.
    """
    # TODO: debug the singular matrix problem
    pdb.set_trace()
    proba = np.ones(nbcod) / nbcod if probavector is None else probavector
    mean = proba * leng
    ker = np.kron(proba, proba).reshape(nbcod, nbcod) * (-leng)
    np.fill_diagonal(ker, (proba * (1 - proba)) * leng)
    if leng**nbcod < MAXITR:
        if nbcod == 2:
            x, y = np.mgrid[0:leng, 0:leng]
            pos = np.empty(x.shape + (nbcod,))
            pos[:, :, 0] = x
            pos[:, :, 1] = y
        elif nbcod == 3:
            x, y, z = np.mgrid[0:leng, 0:leng, 0:leng]
            pos = np.empty(x.shape + (nbcod,))
            pos[:, :, :, 0] = x
            pos[:, :, :, 1] = y
            pos[:, :, :, 2] = z
        elif nbcod == 4:
            x, y, z, v = np.mgrid[0:leng, 0:leng, 0:leng, 0:leng]
            pos = np.empty(x.shape + (nbcod,))
            pos[:, :, :, :, 0] = x
            pos[:, :, :, :, 1] = y
            pos[:, :, :, :, 2] = z
            pos[:, :, :, :, 3] = v
        elif nbcod == 6:
            x, y, z, v, w, o = np.mgrid[0:leng, 0:leng, 0:leng, 0:leng, 0:leng, 0:leng]
            pos = np.empty(x.shape + (nbcod,))
            pos[:, :, :, :, :, :, 0] = x
            pos[:, :, :, :, :, :, 1] = y
            pos[:, :, :, :, :, :, 2] = z
            pos[:, :, :, :, :, :, 3] = v
            pos[:, :, :, :, :, :, 4] = w
            pos[:, :, :, :, :, :, 5] = o
    else:
        pos = np.random.randint(leng, size=(nbcod, MAXITR))
    F = multivariate_normal(mean, ker)
    Z = F.pdf(pos).flatten()
    return Z


def mlen(leng, nb):
    """
    compute wether we need to do the full version of the algorithm or the
    sampling version (see documentation of utils for the threshold values)
    """
    if nb == 6:
        return 'samp' if leng > smax[nb] else 'full'
    elif nb == 4:
        return 'samp' if leng > smax[nb] else 'full'
    elif nb == 3:
        return 'samp' if leng > smax[nb] else 'full'
    elif nb == 2:
        return 'samp' if leng > smax[nb] else 'full'

# This function takes three parameters:
# 1. String
# 2. Starting index of the string
# 3. Ending index of the string.


def permute(vect, start, end):
    """
    compute all possible permutations of the set vect
    """
    if start == end:
        return [tuple(vect)]
    else:
        perm = set([])
        val = []
        for i in range(start, end + 1):
            vect[start], vect[i] = vect[i], vect[start]
            val = permute(list(vect), start + 1, end)
            for j in val:
                perm.add(j)
            vect[start], vect[i] = vect[i], vect[start]  # backtrack
        return perm


def last_match_index(a, value):
    idx = np.array(np.unravel_index(((a < value)[::-1]).argmax(), a.shape))
    return a.shape - idx - 1
