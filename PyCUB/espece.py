""""
Created by Jeremie KALFON
Date : 21 FEV 2018
University of Kent, ECE paris
jkobject.com

"""
import requests
import numpy as np
import utils
from scipy.stats import multinomial
try:
    from urllib2 import urlopen as urlopen
except:
    from urllib.request import urlopen as urlopen


class Espece(object):
    """docstring for Espece

    This is an object that contains all required information of a species for PyCUB and some nice functions
    to interact for each species


    Args:
        code: a dict from gene_name to dna_seq string (deprecated)
        metadata: a dict containing different metadata information that one can gather, preferentially boolean
            flags to classify the species for plottings and comparings
        name: the full scientific name of the species
        is_stored, state if the data is stored in HD (deprecated)
        link: the link to the ensembl genome
        num_genes: the number of gene of this species
        genome_size: the bp size of the coding genome
        name: the name of the species
        taxonid: the number assoxiated to this taxon
        copynumbers: the approx. copynumbers if any of each tRNA known of this species
        average_entropy: the mean CUB value for each amino acids (CUBD dimension)
        average_size: the mean size of each homologies
        var_entropy: the mean of fullvarentropy
        fullentropy: the array containing all CUB values from the homologies of this species
        fullvarentropy: the variance for each amino acids of the full CUB values of this species
        fullGCcount: the GC content of the full coding genome
        varGCcount: the variance of the GC content of the full coding genome
        tRNAentropy: the entropy values of the copynumbers of the tRNAs if sufficient tRNAs exist
        tot_homologies: the total number of homologies to cerevisiae

    """

    code = None  # dict
    num_genes = 0  # int
    genome_size = 0  # int
    link = None  # str
    metadata = {
        "isplant_pathogen": False,
        "isanimal_pathogen": False,
        "isplant_symbiotic": False,  # endophyte or mycorrhizal
        "isbrown_rot": False,
        "iswhite_rot": False
    }
    is_stored = False  # bool
    name = ''  # str
    taxonid = None  # str
    copynumbers = None  # dict
    average_entropy = None  # array float
    average_size = None  # float
    var_entropy = None  # float
    fullentropy = None  # array float
    fullvarentropy = None  # array float
    fullGCcount = None  # int
    varGCcount = None  # float
    meanGChomo = None  # float
    tRNAentropy = None  # array float
    tot_homologies = None  # int
    meanecai = None  # float

    def __init__(self, **kwargs):
        """
        can intialize the file from kwargs as a raw dictionnary for json format (output of dictify) or from regular args.
        """
        data = kwargs.get("data", None)
        if data is not None:
            self.name = data.get("name", None)
            if data.get("metadata", None) is None:
                self.metadata = {
                    "isplant_pathogen": False,
                    "isanimal_pathogen": False,
                    "isplant_symbiotic": False,  # endophyte or mycorrhizal
                    "isbrown_rot": False,
                    "iswhite_rot": False
                }
            else:
                self.metadata = data.get("metadata", None)
            self.code = data.get("code", None)
            self.taxonid = data.get("taxonid", None)
            self.copynumbers = data.get("copynumbers", None)
            self.is_stored = data.get('is_stored', None)
            self.average_entropy = np.asarray(data["average_entropy"]) if data.get("average_entropy", False) else None
            self.average_size = data.get("average_size", None)
            self.var_entropy = data.get("var_entropy", None)
            self.fullentropy = np.asarray(data["fullentropy"]) if data.get("fullentropy", False) else None
            self.fullGCcount = data.get("fullGCcount", None)
            self.tRNAentropy = np.asarray(data["tRNAentropy"]) if data.get("tRNAentropy", False) else None
            self.num_genes = data.get("num_genes", 0)
            self.genome_size = data.get("genome_size", 0)
            self.link = data.get("link", None)
            self.fullvarentropy = np.asarray(data["fullvarentropy"]) if data.get("fullvarentropy", False) else None
            self.varGCcount = data.get("varGCcount", None)
            self.tot_homologies = data.get("tot_homologies", None)
            self.meanGChomo = data.get("meanGChomo", None)
            self.meanecai = data.get("meanecai", None)

        else:
            self.code = kwargs.get('code', None)
            self.is_stored = kwargs.get('is_stored', None)
            self.taxonid = kwargs.get('taxonid', '')
            self.name = kwargs.get('name', '')
            self.copynumbers = kwargs.get('copynumbers', None)
            self.average_entropy = kwargs.get("average_entropy", None)
            self.average_size = kwargs.get("average_size", None)
            self.var_entropy = kwargs.get("var_entropy", None)
            self.fullentropy = kwargs.get("fullentropy", None)
            self.fullGCcount = kwargs.get("fullGCcount", None)
            self.tRNAentropy = kwargs.get("tRNAentropy", None)
            self.num_genes = kwargs.get("num_genes", 0)
            self.genome_size = kwargs.get("genome_size", 0)
            self.link = kwargs.get("link", None)
            self.fullvarentropy = kwargs.get("fullvarentropy", None)
            self.varGCcount = kwargs.get("varGCcount", None)
            self.tot_homologies = kwargs.get("tot_homologies", None)
            self.metadata = kwargs.get("metadata", None)
            self.meanGChomo = kwargs.get("meanGChomo", None)
            self.meanecai = kwargs.get("meanecai", None)

    def __str__(self):
        """
        will present some interesting info about this species.
        """
        if self.name:
            print "species: " + self.name
        print "----------------------------------"
        if self.taxonid:
            print "taxonid" + str(self.taxonid)
        print "metadata" + str(self.metadata)
        print "----------------------------------"
        if self.copynumbers is not None:
            print "copynumbers of tRNA: " + str(self.copynumbers)
        if self.average_size is not None:
            print "average size: " + str(self.average_size)
        if self.tRNAentropy is not None:
            print "tRNA entropy: " + str(self.tRNAentropy)
        if self.num_genes:
            print "number of genes: " + str(self.num_genes)
        if self.genome_size:
            print "genome size: " + str(self.genome_size)
        if self.tot_homologies is not None:
            print "total number of homologies to cerevisiae: " + str(self.tot_homologies)
        print "----------------------------------"
        if self.average_entropy is not None:
            print "average entropy: " + str(self.average_entropy)
        if self.var_entropy is not None:
            print "variance of entropy: " + str(self.var_entropy)
        if self.fullvarentropy is not None:
            print "full variance of entropy: " + str(self.fullvarentropy)
        if self.varGCcount is not None:
            print "variance of the GC content: " + str(self.varGCcount)
        if self.meanecai is not None:
            print "mean ECAI: " + str(self.meanecai)

# CpCpGpApApTpApTpApTpTpCpCpGpApApTpApTpApTpTpCpCpGpApApTpApTpApTpTpCpCpGpApApTpApTpApTpTpTpTpCpCpGpApApTpApTpApTpTp
# GbGbCbTbTbAbTbAbTbAbAbGbGbCbTbTbAbTbAbTbAbAbGbGbCbTbTbAbTbAbTbAbAbGbGbCbTbTbAbTbAbTbAbAbAbGbGbCbTbTbAbTbAbTbAbAbAb

    def get_tRNAcopy(self, by="entropy", setnans=False, kingdom='fungi', baseCNvalue=2):
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

        Args:
            species: string, the species from which you want the Trna copy number

        Returns:
            Will populate copynumbers. And tRNAentropy if by="entropy"
            Or will not do anything if the species is unavailable and will print it

        Raises:
            AttributeError: this is a wrong argument try frequency or entropy
        """
        server = "http://rest.ensemblgenomes.org" if kingdom != 'vertebrate' else "http://rest.ensembl.org"
        print 'species: ' + self.name
        ext = "/lookup/genome/" + self.name + '?'
        add = "biotypes=tRNA;level=transcript"
        r = requests.get(server + ext + add, headers={"Content-Type": "application/json"})
        if not r.ok:
            print " !! ---> unavailable species"
            return
        data = r.json()
        copynumber = {}
        for key, val in utils.anticodons.iteritems():
            copynumber.update({key: {}})
            for v in val:
                v.replace('T', 'U')
                copynumber[key].update({v.replace('T', 'U'): baseCNvalue})
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
        if num == 0:
            print "empty data"
        print "we got " + str(j) + " datapoints and managed to extract " + str(num)
        # we find probabilities of tRNA
        k = 0
        if num > 100:
            tRNAentropy = np.zeros(18) if by == "entropy" else None
            for _, v in copynumber.iteritems():
                n = np.array(v.values()).sum()
                if n > 0:
                    for _, val in v.iteritems():
                        val = val / n
                # Else we keep the raw frequency values
                if by == "entropy":
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
                elif by != "frequency":
                    raise AttributeError("this is a wrong argument try frequency or entropy")

        # Here we can compute as well the entropy of the tRNA CNs when there is suficient number of val
        # else we can set it to zero (to NaN) this allows us to directly compare two values
        copynumber.update({'num': num})  # total number
        copynumber.update({'datapoints': j})  # possible number of datapoints
        self.copynumbers = copynumber
        if by == "entropy" and num > 100:
            self.tRNAentropy = tRNAentropy

    def gettaxons(self, kingdom='fungi'):
        """
        Pars the ensemblgenomes REST API to retrieve the taxons id 

        for the species from which we would not have any (downloaded via Yun for example)

        Raises:
            HTTPrequestError: not able to connect to the server
        """
        # http: // rest.ensemblgenomes.org / info / genomes / arabidopsis_thaliana?
        server = "http://rest.ensemblgenomes.org" if kingdom != 'vertebrate' else "http://rest.ensembl.org"
        print 'species: ' + self.name
        ext = "/info/genomes/" + self.name + '?'
        r = requests.get(server + ext, headers={"Content-Type": "application/json"})
        if not r.ok:
            r.raise_for_status()
        data = r.json()
        self.taxonid = data["species_taxonomy_id"]

    def get_epigenomes(self):
        """
        get from ensembl all the data about the epigenome that could help asking interesting questions about the CUB
        """
        # curl 'http://rest.ensemblgenomes.org/overlap/id/AT3G52430?feature=array_probe' - H 'Content-type:application/json'
        # curl 'http://rest.ensemblgenomes.org/overlap/id/AT3G52430?feature=repeat' - H 'Content-type:application/json'
        pass

    def _dictify(self):
        """
        Used by the saving function. transform the object into a dictionary that can be json serializable

        Returns:
            A dict holding every element to be jsonized
        """
        return {"name": self.name,
                "code": self.code,
                "num_genes": self.num_genes,
                "genome_size": self.genome_size,
                "link": self.link,
                "fullvarentropy": self.fullvarentropy.tolist() if self.fullentropy is not None else None,
                "varGCcount": self.varGCcount,
                "tot_homologies": self.tot_homologies,
                "taxonid": self.taxonid,
                "copynumbers": self.copynumbers,
                "metadata": self.metadata,
                "meanGChomo": self.meanGChomo,
                "meanecai": self.meanecai,
                "is_stored": self.is_stored,
                "average_entropy": self.average_entropy.tolist() if self.fullentropy is not None else None,
                "average_size": self.average_size,
                "var_entropy": self.var_entropy,
                "fullentropy": self.fullentropy.tolist() if self.fullentropy is not None else None,
                "fullGCcount": self.fullGCcount.tolist() if self.fullGCcount is not None else None,
                "tRNAentropy": self.tRNAentropy.tolist() if self.tRNAentropy is not None else None}
