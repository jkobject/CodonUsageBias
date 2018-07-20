""""
Created by Jeremie KALFON
Date : 21 FEV 2018
University of Kent, ECE paris
jkobject.com

"""
import requests
import pdb
import numpy as np
import utils
from scipy.stats import multinomial
try:
    from urllib2 import urlopen as urlopen
except:
    from urllib.request import urlopen as urlopen


class Espece(object):
    """docstring for Espece


            Params:
            ------
            code: a dict from gene_name to dna_seq string
            metadata: a dict containing different metadata information
            name : the full scientific name of the species
            is_stored
    """

    code = None
    metadata = {
        "link": None,
        "num_genes": 0,
        "plant_pathogen": None,
        "animal_pathogen": None,
        "genome_size": 0,
        "plant_symbiotic": None,  # endophyte or mycorrhizal
        "brown_rot": None,
        "white_rot": None
    }
    is_stored = False
    name = ''
    taxonid = ''
    copynumbers = None
    average_entropy = None
    average_size = None
    var_entropy = None
    var_size = None
    fullentropy = None
    fullGCcount = None
    tRNAentropy = None

    def __init__(self, **kwargs):
        """
        """
        data = kwargs.get('data', None)
        if data is not None:
            self.name = data.get(["name"], None)
            self.metadata = data.get(["name"], None)["metadata"]
            self.code = data.get(["code"], None)
            self.taxonid = data.get(["taxonid"], None)
            self.copynumbers = data.get(["copynumbers"], None)
            self.code = data.get('code', None)
            self.is_stored = data.get('is_stored', None)
            self.average_entropy = data.get("average_entropy", None)
            self.average_size = data.get("average_size", None)
            self.var_entropy = data.get("var_entropy", None)
            self.var_size = data.get("var_size", None)
            self.fullentropy = data.get("fullentropy", None)
            self.fullGCcount = data.get("fullGCcount", None)
            self.tRNAentropy = data.get("tRNAentropy", None)
        else:
            self.code = kwargs.get('code', None)
            self.is_stored = kwargs.get('is_stored', None)
            self.taxonid = kwargs.get('taxonid', '')
            self.name = kwargs.get('name', '')
            self.copynumbers = kwargs.get('copynumbers', None)
            self.average_entropy = kwargs.get("average_entropy", None)
            self.average_size = kwargs.get("average_size", None)
            self.var_entropy = kwargs.get("var_entropy", None)
            self.var_size = kwargs.get("var_size", None)
            self.fullentropy = kwargs.get("fullentropy", None)
            self.fullGCcount = kwargs.get("fullGCcount", None)
            self.tRNAentropy = kwargs.get("tRNAentropy", None)

    def plot():
        """
        will present some interesting info about this species.
        """

    def get_metadata():
        """
        will call all functions to retrieve all possible metadatas for this species.
        """

    def get_tRNAcopy(self, getentropy=True, setnans=False):
        """
        Retrieves tRNA copy numbers from ensembl DB
        will print the number of tRNAs and the number of tRNAs with
        a knwon codons ( the usefull ones)

        will stop and set a trace for the user to inspect the data
        to do so: please write "dat" in the console. if you see something that
        should be corrected please do so from the console directly or from the code
        if there seems to be an error in the code

        if it is an error in the db that you can't do anything, like a mismatched codon
        and amino acid, you can't do much. resume the process by typing "c" in the console.

        Params:
        -------
        species: string, the species from which you want the Trna copy number
        Returns:
        ---------
        dict[dict] : the tRNA copy number for the species given.
        """
        server = "http://rest.ensemblgenomes.org"
        print 'species: ' + self.name
        ext = "/lookup/genome/" + self.name + '?'
        add = "biotypes=tRNA;level=transcript"
        r = requests.get(server + ext + add, headers={"Content-Type": "application/json"})
        if not r.ok:
            r.raise_for_status()
        data = r.json()
        copynumber = {}
        for key, val in utils.anticodons.iteritems():
            copynumber.update({key: {}})
            for v in val:
                v.replace('T', 'U')
                copynumber[key].update({v.replace('T', 'U'): 0})
        num = 0
        j = 0
        for j, dat in enumerate(data):
            if dat["name"] is not None and len(dat["name"]) != 0:
                if dat["name"][0:4] == 'tRNA':
                    try:
                        if dat["description"] is None:
                            if dat["name"][10:13] != '':
                                if len(dat["name"]) == 14:
                                    codn = dat["name"][10:13]
                                if len(dat["name"]) == 13:
                                    codn = dat["name"][9:12]
                                if 'T' in codn:
                                    codn = codn.replace('T', 'U')
                            else:
                                continue
                        else:
                            codn = dat["description"][23:26]
                        ami = dat["name"][5:8].upper()
                        if ami == 'SEC':
                            ami = 'SER'
                            copynumber[ami][codn] += 1
                            num += 1
                        elif ami in ['TRP', 'MET', 'UND', 'SUP', 'UNK']:
                            continue
                        elif ami == 'PSE':
                            codn = dat["name"][12:15].upper()
                            for key, val in copynumber.iteritems():
                                if type(val) is dict:
                                    for k, v in val.iteritems():
                                        if k == codn:
                                            copynumber[key][k] += 1
                        else:
                            copynumber[ami][codn[::-1]] += 1
                            num += 1
                    except KeyError:
                        print "KeyError"
                        pdb.set_trace()
                elif dat["name"][0:3] == 'trn':
                    try:
                        codn = dat["name"][5:8].upper()
                        if 'T' in codn:
                            codn = codn.replace('T', 'U')
                        ami = utils.amino2reduce[dat["name"][3]]
                        if ami in ['TRP', 'MET', 'UND', 'SUP', 'UNK']:
                            continue
                        else:
                            copynumber[ami][codn[::-1]] += 1
                            num += 1
                    except KeyError:
                        print "KeyError"
                        pdb.set_trace()
            elif dat["description"] is not None and len(dat["description"]) > 10:
                if dat["description"][0:4] == 'tRNA':
                    try:
                        codn = dat["description"][23:26]
                        ami = dat["description"][5:8].upper()
                        if ami == 'SEC':
                            ami = 'SER'
                            copynumber[ami][codn] += 1
                            num += 1
                        elif ami in ['TRP', 'MET', 'UND', 'SUP', 'UNK']:
                            continue
                        else:
                            copynumber[ami][codn[::-1]] += 1
                            num += 1
                    except KeyError:
                        print "KeyError"
                        pdb.set_trace()
        if num == 0:
            print "empty data"
        print "we got " + str(j) + " datapoints and managed to extract " + str(num)
        # we find probabilities of tRNA
        k = 0
        tRNAentropy = np.zeros(18) if getentropy else None
        for _, v in copynumber.iteritems():
            n = np.array(v.values()).sum()
            if n > 0:
                for _, val in v.iteritems():
                    val = val / n
            if getentropy:
                nbcod = len(v)  # replace Cleng
                count = v.values()
                X = np.zeros(nbcod)
                mn = np.ones(nbcod) / nbcod
                if n == 0:
                    tRNAentropy[k] = np.NaN if setnans else 0.5
                else:
                    Yg = multinomial.pmf(x=count, n=n, p=mn)
                    # efor part
                    div, i = divmod(n, nbcod)
                    X[:i] = np.ceil(div) + 1
                    X[i:] = np.floor(div)
                    Eg = multinomial.pmf(x=X, n=n, p=mn)
                    # end here
                    tRNAentropy[k] = -np.log(Yg / Eg)
                k += 1
        # Here we can compute as well the entropy of the tRNA CNs when there is suficient number of val
        # else we can set it to zero (to NaN) this allows us to directly compare two values
        copynumber.update({'num': num})  # total number
        copynumber.update({'datapoints': j})  # possible number of datapoints
        self.copynumbers = copynumber
        self.tRNAentropy = tRNAentropy if getentropy else None

    def gettaxons(self):
        """
        """
        species_namelist = set([])
        taxons = []
        species = []  # we need to have another species list for the ordering
        # to be kept as it should (a set is ordered to be accessed in log time)
        utils.speciestable = {}
        if self.homodict is not None:
            i = 0
            helper = {}
            for key, val in self.homodict.iteritems():
                for j, name in enumerate(val.names):
                    species_namelist.add(utils.speciestable[name])
            self.species_namelist = species_namelist
            return taxons, species

    def getfullgenome(self):
        """
        get the full genome from ensembl.
        """
        # preprocess the data from ensembl
        # TODO: find how to get full genome from ensembl
        self.fullentropy, _, _, self.fullGCcount = utils.computeyun(data=fullcdna, normalized=False, setnans=False, by='entropy')
        pass

    def get_epigenomes(self):
        """
        get from ensembl all the data about the epigenome that could help asking interesting questions
        about the CUB

        Params:
        ------
        """
        pass

    def relationtomymeta(self):
        """
        you can add you own metadata and then think about ways to look at wether or not they
        influence the codon usage bias values.
        """
        pass
        # compute the influence of their relationship to plants (symbiotic or rot)
        # compute the influence of their nocivity
        # compute the influence of the size of their prot cod genes, the genome size, the number of genes
        # compute the influence of the GC count

    def _dictify(self):
        """
        Used by the saving function. transform the object into a dictionary that can be
        json serializable
        """
        return {"name": self.name,
                "code": self.code,
                "taxonid": self.taxonid,
                "copynumbers": self.copynumbers,
                "metadata": self.metadata,
                "is_stored": self.is_stored,
                "average_entropy": self.average_entropy,
                "average_size": self.average_size,
                "var_entropy": self.var_entropy,
                "var_size": self.var_size,
                "fullentropy": self.fullentropy,
                "fullGCcount": self.fullGCcount,
                "tRNAentropy": self.tRNAentropy}
