""""
Created by Jeremie KALFON
Date : 21 FEV 2018
University of Kent, ECE paris
jkobject.com

"""
import requests
import pdb
import utils
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

    code = {}
    metadata = {}
    is_stored = False
    name = ''
    taxonid = ''
    copynumbers = None

    def __init__(self, name=None, code=None, data=None, link=None, taxonid='', copynumbers=None):
        """
        """
        if data is not None:
            self.name = data["name"]
            self.metadata = data["metadata"]
            self.code = data["code"]
            self.taxonid = data["taxonid"]
            self.copynumbers = data["copynumbers"]
        if taxonid is not '':
            self.taxonid = taxonid
        if name is not None:
            self.name = name
        if link is not None:
            self.metadata.update({'link': link})
        if copynumbers is not None:
            self.copynumbers = copynumbers

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

    def get_tRNAcopy(self):
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
        else:
            print num
        print j
        copynumber.update({'num': num})
        copynumber.update({'tot_trna': j})
        self.copynumbers = copynumber

    def compute_all_genome_entropy():
        """
        use all the entropy locations and use Yun's way to find the full genome entropy
        """
        pass

    def getfullgenome():
        """
        """
        pass

    def _dictify(self):
        """
        Used by the saving function. transform the object into a dictionary that can be
        json serializable
        """
        return {"name": self.name,
                "code": self.code,
                "taxonid": self.taxonid,
                "copynumbers": self.copynumbers,
                "metadata": self.metadata}
