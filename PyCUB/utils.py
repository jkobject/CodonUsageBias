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
import math
import glob
import numpy as np
import requests
import os
from math import sqrt
import json
try:
    from urllib2 import urlopen as urlopen
except:
    from urllib.request import urlopen as urlopen

import homology as h
import pdb
"""
    the utils is where utilitary functions are stored (in particular preprocessing functions here)

    Values:
    ------
    speciestable: dict[int] will be used by all other objects to make numbers corresponds to species
    indexcai: a dict[amino: value] to known what is the frequency of this amino acid in a gene X of the reference species for the
        reference_index() and computeCAI() function
    indexecai: a dict[amino: value] to known what is the frequency of this amino acid in a gene X of the reference species for the
        reference_index() and computeECAI() function it is different to the CAI ones as it does not contains
        the references for a set of highly expressed genes but the ones for each homologies
    listvect: used by the fast partition compute function
    phylo_distances: df where the phylodistance df is stored (see PyCUB.pyCUB.get_evolutionnary_distance())
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

    anticodons: dict the anticodons (usefull sometimes see tRNA_copynumbers)
    amino2reduce: dict the mapping of amino acid values to their correponding char
    codamino: dict the map from codons to aminoacids
    CUBD: int the dimension of the CUB values (either 18 or 49)
    MAX: int a max value of leng for the computepartition function above which the size is so big that it overflows python
        arbitrary sized values are not present in python funtions that are optimized
    amino2meta: dict mapping of amino acids to some of their chemical characteristics displayed below
    synthsteps: list the number of step to synthetize this amino acid
    glucosecost: list the number of glucode to synth this amino acid
    hydrophob: list a value of how hydrophobic this amino acid is
    volume: list the volume of this molecule ( in Ang^3)
    isoelectricpoint: list the Pi value of this amino
    conservation: list how much this amino acid is conserved (deprecated)
    colormap: list the colormap used in PyCUB
    callback: str the callback used in the homology interactive plot (used to change colors of dots according to what data
        is selected)
    callback_allgenes: 
    callback_allhomo: str the callback used in the PyCUB interactive plot (used to change colors of dots according to what data
        is selected)
"""
speciestable = {}
indexcai = {}
indexecai = {}

listvect = None

phylo_distances = None  # pandas dataframe

codamino = {
    'GCA': 'ALA',
    'GCC': 'ALA',
    'GCG': 'ALA',
    'GCT': 'ALA',

    'CGA': 'ARG',
    'CGC': 'ARG',
    'CGG': 'ARG',
    'CGT': 'ARG',
    'AGG': 'ARG',
    'AGA': 'ARG',

    'AAC': 'ASN',
    'AAT': 'ASN',

    'GAT': 'ASP',
    'GAC': 'ASP',

    'TGT': 'CYS',
    'TGC': 'CYS',

    'CAA': 'GLN',
    'CAG': 'GLN',

    'GAG': 'GLU',
    'GAA': 'GLU',

    'GGT': 'GLY',
    'GGG': 'GLY',
    'GGA': 'GLY',
    'GGC': 'GLY',

    'CAT': 'HIS',
    'CAC': 'HIS',

    'ATC': 'ILE',
    'ATA': 'ILE',
    'ATT': 'ILE',

    'TTA': 'LEU',
    'TTG': 'LEU',
    'CTC': 'LEU',
    'CTT': 'LEU',
    'CTG': 'LEU',
    'CTA': 'LEU',

    'AAA': 'LYS',
    'AAG': 'LYS',

    'TTT': 'PHE',
    'TTC': 'PHE',

    'CCT': 'PRO',
    'CCG': 'PRO',
    'CCA': 'PRO',
    'CCC': 'PRO',

    'TCT': 'SER',
    'TCG': 'SER',
    'TCA': 'SER',
    'TCC': 'SER',
    'AGC': 'SER',
    'AGT': 'SER',

    'ACC': 'THR',
    'ACA': 'THR',
    'ACG': 'THR',
    'ACT': 'THR',

    'TAT': 'TYR',
    'TAC': 'TYR',

    'GTA': 'VAL',
    'GTC': 'VAL',
    'GTG': 'VAL',
    'GTT': 'VAL'
}
codons = {
    'ALA': ['GCA', 'GCC', 'GCG', 'GCT'],  # GC
    'ARG': ['CGA', 'CGC', 'CGG', 'CGT', 'AGG', 'AGA'],  # CG
    'ASN': ['AAC', 'AAT'],  # ASN - LYS
    'ASP': ['GAT', 'GAC'],  # ASP - GLU
    'CYS': ['TGT', 'TGC'],  # ~
    'GLN': ['CAA', 'CAG'],  # GLN - HIS
    'GLU': ['GAG', 'GAA'],
    'GLY': ['GGT', 'GGG', 'GGA', 'GGC'],  # GG
    'HIS': ['CAT', 'CAC'],
    'ILE': ['ATC', 'ATA', 'ATT'],  # ILE -- MET
    'LEU': ['TTA', 'TTG', 'CTC', 'CTT', 'CTG', 'CTA'],  # CT
    'LYS': ['AAA', 'AAG'],
    # 'MET': ['ATG'],
    'PHE': ['TTT', 'TTC'],  # PHE - LEU
    'PRO': ['CCT', 'CCG', 'CCA', 'CCC'],  # CC
    'SER': ['TCT', 'TCG', 'TCA', 'TCC', 'AGC', 'AGT'],  # TC
    # 'STOP': ['TAG', 'TGA', 'TAA'],
    'THR': ['ACC', 'ACA', 'ACG', 'ACT'],  # AC
    # 'TRP': ['TGG'],
    'TYR': ['TAT', 'TAC'],
    'VAL': ['GTA', 'GTC', 'GTG', 'GTT']}  # GT
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
MAXITR = 5000000
smax = [0, 0, LMAX, 1414, 100, 0, 35]
smax2 = [0, 0, LMAX, LMAX, 520, 0, 145]
CUBD = 18
MAX = 20  # max value of leng after which we encounter
# weird results in compute jerem's computation nbcod**leng
#
amino2meta = {'ALA': [1, 0.5, 1.8, 88.6, 6.11, 5],  # synthsteps, glucosecost, hydrophob, volume, isoelectricpoint
              'ARG': [10, 1.39, -4.5, 173.4, 10.76, 5.2],
              'ASN': [1, 0.79, -3.5, 114.1, 10.76, 4],
              'ASP': [1, 0.61, -3.5, 111.1, 2.98, 5.3],
              'CYS': [9, 0.75, 2.5, 108.5, 5.02, 4.4],
              'GLN': [2, 0.92, -3.5, 143.8, 5.65, 3.4],
              'GLU': [1, 0.86, -3.5, 138.4, 3.08, 4.4],
              'GLY': [4, 0.31, -0.4, 60.1, 6.06, 6.6],
              'HIS': [1, 1.46, -3.2, 153.2, 7.64, 3.1],
              'ILE': [11, 1.21, 4.5, 166.7, 6.04, 7.8],
              'LEU': [7, 1.21, 3.8, 166.7, 6.04, 7.3],
              'LYS': [10, 1.31, -3.9, 168.6, 9.47, 4.4],
              'PHE': [9, 1.84, 2.8, 162.9, 5.91, 4.2],
              'PRO': [4, 0.99, -1.6, 189.9, 6.3, 5.6],
              'SER': [3, 0.49, -0.8, 112.7, 5.68, 3.8],
              'THR': [6, 0.69, -0.7, 89, 5.6, 3.8],
              'TYR': [9, 1.77, -1.3, 193.6, 5.63, 3.2],
              'VAL': [4, 0.96, 4.2, 140, 6.02, 6.5]
              }
synthsteps = [1, 10, 1, 1, 9, 2, 1, 4, 1, 11,
              7, 10, 9, 4, 3, 6, 9, 4]

glucosecost = [0.5, 1.39, 0.79, 0.61, 0.75, 0.92, 0.86,
               0.31, 1.46, 1.21, 1.21, 1.31, 1.84, 0.99, 0.49, 0.69, 1.77, 0.96]

hydrophob = [1.8, -4.5, -3.5, -3.5, 2.5, -3.5, -3.5, -0.4, -3.2, 4.5,
             3.8, -3.9, 2.8, -1.6, -0.8, -0.7, -1.3, 4.2]

volume = [88.6, 173.4, 114.1, 111.1, 108.5, 143.8, 138.4, 60.1, 153.2,
          166.7, 166.7, 168.6, 162.9, 189.9, 112.7, 89, 193.6, 140]

isoelectricpoint = [6.11, 10.76, 10.76, 2.98, 5.02, 5.65, 3.08, 6.06, 7.64,
                    6.04, 6.04, 9.47, 5.91, 6.3, 5.68, 5.6, 5.63, 6.02]

conservation = [5, 5.2, 4, 5.3, 4.4, 3.4, 4.4, 6.6, 3.1, 7.8, 7.3, 4.4, 4.2, 5.6, 3.8, 3.8, 3.2, 6.5]

colormap = ['#f39c12', "#1abc9c", "#3498db", "#2ecc71", "#9b59b6", '#34495e', '#f1c40f', '#e67e22', '#e74c3c', '#7f8c8d']
callback = """
            // JavaScript code goes here
            console.log("callback-------callback");
            var colors = ['#f39c12', "#1abc9c", "#3498db", "#2ecc71",
                        "#9b59b6", '#34495e', '#f1c40f','#e67e22', '#e74c3c', '#7f8c8d'];
            // the model that triggered the callback is cb_obj:
            var b = cb_obj.get("active");
            var rgbToHex = function(rgb){
              var hex = Number(rgb).toString(16);
              if (hex.length < 2) {
                   hex = "0" + hex;
              }
              return hex;
            };
            var fullColorHex = function(r,g,b) {
                var red = rgbToHex(r);
                var green = rgbToHex(g);
                var blue = rgbToHex(b);
                return "#"+red+green+blue;
            };
            // models passed as args are automagically available
            var data = source.data;
            var len = data.color.length;
            var col = [];
            var othercol = '#1abc9c';
            if(b === 1){ //Doublon
                var j = 0;
                for(i=0;i<len;++i){
                    if(data.doub[i]){
                        col.push(colors[j]);
                        col[i-1] = colors[j];
                        j++;
                    }else{col.push("#999999");}
                }
                source.data.color = col;
            }
            if(b ===2){ //Nans
                var temp = 0
                max = Math.max(...data.nans)
                for(i=0;i<len;++i){
                    temp =  Math.floor(255 * (1 - ((data.nans[i] + 1) / 19)));
                    col.push(fullColorHex(temp, temp, temp));
                    }
                source.data.color = col;
            }
            if(b === 0) { //Cluster
                var j = 0;
                if(data.clusters){
                    for(i=0;i<len;++i){col.push(colors[1+data.clusters[i]]);}
                }else{
                    for(i=0;i<len;++i){col.push(othercol);}
                }
                source.data.color = col;
            }
            if(b === 3) { //KaKs_Scores
                max = Math.max(...data.KaKs_Scores)
                for(i=0;i<len;++i){
                    // normalizing
                    temp =  Math.floor(255 * (1 - (data.KaKs_Scores[i] / max)));
                    col.push(fullColorHex(temp,temp,76))
                }
                source.data.color = col;
            }
            if(b === 4){ //similarity_scores
                max = Math.max(...data.similarity_scores)
                for(i=0;i<len;++i) {
                    temp =  Math.floor(255 * (1 - (data.similarity_scores[i] / max)));
                    col.push(fullColorHex(temp,temp,76))
                }
                data.color = col;
            }
            if(b === 5){ //Lengths
                max = Math.max(...data.lengths)
                for(i=0;i<len;++i) {
                    temp =  Math.floor(255 * (1 - (data.lengths[i] / max)));
                    col.push(fullColorHex(temp, temp, temp));
                }
                data.color = col;
            }
            if(b===6){ //GC
                var temp = 0
                for(i=0;i<len;++i){
                    temp =  Math.floor(255 * data.gc[i]);
                    col.push(fullColorHex(52, temp, temp));
                    }
                source.data.color = col;
            }
            if(b===7){ //ECAI
                var temp = 0
                for(i=0;i<len;++i){
                    temp =  Math.floor(255 * data.ecai[i]);
                    col.push(fullColorHex(52, temp, temp));
                    }
                source.data.color = col;
            }
            if(b===9){ //CAI
                var temp = 0
                for(i=0;i<len;++i){
                    temp =  Math.floor(255 * data.cai[i]);
                    col.push(fullColorHex(52, temp, temp));
                    }
                source.data.color = col;
            }
            source.change.emit();
            """
# TODO: maybe use a nonlinearity to increase the differences amongst color when the datapoints
# are super close to each others (gc values, kaks scores..)
callback_allgenes = """
            // JavaScript code goes here
            console.log("callback-------callback");
            var colors = ['#f39c12', "#1abc9c", "#3498db", "#2ecc71",
                        "#9b59b6", '#34495e', '#f1c40f','#e67e22', '#e74c3c', '#7f8c8d'];
            // the model that triggered the callback is cb_obj:
            var b = cb_obj.get("active");
            var rgbToHex = function(rgb){
              var hex = Number(rgb).toString(16);
              if (hex.length < 2) {
                   hex = "0" + hex;
              }
              return hex;
            };
            var fullColorHex = function(r,g,b) {
                var red = rgbToHex(r);
                var green = rgbToHex(g);
                var blue = rgbToHex(b);
                return "#"+red+green+blue;
            };
            // models passed as args are automagically available
            var data = source.data;
            var len = data.color.length;
            var col = [];
            var othercol = '#1abc9c';
            
            if(b === 0){ //show clusters
                var j = 0;
                if(data.clusters){
                    for(i=0;i<len;++i){col.push(colors[1+data.clusters[i]]);}
                }else{
                    for(i=0;i<len;++i){col.push(othercol);}
                }
                source.data.color = col;
            }
            if(b === 1){//show num_genes
                var temp = 0
                max = Math.max(...data.num_genes)
                for(i=0;i<len;++i){
                    temp =  Math.floor(92 * data.num_genes[i]/max);
                    col.push(fullColorHex(temp, 152, 219));
                    }
                source.data.color = col;
            }
            if(b === 2) { //show genome_size
                var temp = 0
                max = Math.max(...data.genome_size)
                for(i=0;i<len;++i){
                    temp =  Math.floor(252 * data.genome_size[i]/max);
                    col.push(fullColorHex(temp, 152, 219));
                    }
                source.data.color = col;
            }
            if(b === 3) {  //show show full/mean CUB diff
                var temp = 0
                for(i=0;i<len;++i){
                    temp =  Math.floor(200 * data.efulldiff[i]);
                    col.push(fullColorHex(temp, 204, 113));
                    }
                source.data.color = col;
            }
            if(b === 4){ //Show GC count
                var temp = 0
                for(i=0;i<len;++i) {
                    temp =  Math.floor(255 * data.gccount[i]);
                    col.push(fullColorHex(temp, temp, 113));
                }
                data.color = col;
            }
            if(b === 5){ //show fullmean GC diff
                var temp = 0
                for(i=0;i<len;++i) {
                    temp =  Math.floor(255 * data.gcfulldiff[i]);
                    col.push(fullColorHex(temp, 204, temp));
                }
                data.color = col;
            }
            if(b === 6){ //show GC variance
                var temp = 0
                for(i=0;i<len;++i) {
                    temp =  Math.floor(455 * data.varGCcount[i]);
                    col.push(fullColorHex(temp, temp, 113));
                }
                data.color = col;
            }
            if(b === 7){ //show distance to tRNA UB"
                source.data.color = source.data.recent;
            }
            if(b === 8){ //  show tRNA number
                var temp = 0
                max = Math.max(...data.tRNA_number)
                for(i=0;i<len;++i){
                    temp =  Math.floor(54 + (200 * data.tRNA_number[i]/max));
                    col.push(fullColorHex(52, 73, temp));
                    }
                source.data.color = col;
            }
            if(b === 9){ //  show avgphilodistance
                var temp = 0
                max = Math.max(...data.distances)
                min = Math.min(...data.distances)
                for(i=0;i<len;++i){
                    temp =  Math.floor(54 + (254 * data.distances[i]/max));
                    col.push(fullColorHex(52, 73, temp));
                    }
                source.data.color = col;
            }
            /*
            if(b === 10){ //  show full phylodistance
                var temp = 0
                for(i=0;i<len;++i){
                    temp =  Math.floor(64 + 100 * data.gc[i]));
                    col.push(fullColorHex(52, 73, temp));
                    }
                source.data.color = col;
            }*/
            if(b > 10){ //  show avgphilodistance
                for(i=0;i<len;++i){
                    if(data[String(b)][i]){
                        col.push(colors[2])
                    }else{
                        col.push(colors[3])
                    }
                    }
                source.data.color = col;
            }
            source.change.emit();
            """
callback_allhomo = """
            // JavaScript code goes here
        console.log("callback-------callback");
        var colors = ['#f39c12', "#1abc9c", "#3498db", "#2ecc71",
                    "#9b59b6", '#34495e', '#f1c40f','#e67e22', '#e74c3c', '#7f8c8d'];
        // the model that triggered the callback is cb_obj:
        var b = cb_obj.get("active");
        var rgbToHex = function(rgb){
          var hex = Number(rgb).toString(16);
          if (hex.length < 2) {
               hex = "0" + hex;
          }
          return hex;
        };
        var fullColorHex = function(r,g,b) {
            var red = rgbToHex(r);
            var green = rgbToHex(g);
            var blue = rgbToHex(b);
            return "#"+red+green+blue;
        };
        // models passed as args are automagically available
        var data = source.data;
        var len = data.color.length;
        var col = [];
        var othercol = '#1abc9c';
        
        if(b === 0){ //"show Recent/preserved"
            source.data.color = data.recent;
        }
        if(b === 1){//"showclusters"
            var j = 0;
            if(data.clusters){
                for(i=0;i<len;++i){col.push(colors[1+data.clusters[i]]);}
            }else{
                for(i=0;i<len;++i){col.push(othercol);}
            }
            source.data.color = col;
        }
        if(b === 2) { //"show avg similarity_scores"
            var temp = 0
            max = Math.max(...data.similarity_scores)
            for(i=0;i<len;++i){
                temp =  Math.floor(252 * data.similarity_scores[i]/max);
                col.push(fullColorHex(temp, 152, 219));
                }
            source.data.color = col;
        }
        if(b === 3) {  //"show avg KaKs_Scores"
            var temp = 0
            max = Math.max(...data.KaKs_Scores)
            for(i=0;i<len;++i){
                temp =  Math.floor(200 * data.KaKs_Scores[i]/max);
                col.push(fullColorHex(temp, 204, 113));
                }
            source.data.color = col;
        }
        if(b === 4){ //"show Nans avg"
            var temp = 0
            max = Math.max(...data.nans)
            for(i=0;i<len;++i){
                temp =  Math.floor(200 * data.nans[i]/max);
                col.push(fullColorHex(temp, temp, 113));
            }
            source.data.color = col;
        }
        if(b === 5){ //lenmat
            var temp = 0
            max = Math.max(...data.lenmat)
            for(i=0;i<len;++i){
                temp =  Math.floor(200 * data.lenmat[i]/max);
                col.push(fullColorHex(temp, temp, 113));
            }
            source.data.color = col;
        }
        if(b === 6){ //GCcount
            var temp = 0
            for(i=0;i<len;++i) {
                temp =  Math.floor(255 * data.GCcount[i]);
                col.push(fullColorHex(temp, 204, temp));
            }
            source.data.color = col;
        }
        if(b === 7){ //weight
            var temp = 0
            max = Math.max(...data.weight)
            for(i=0;i<len;++i){
                temp =  Math.floor(200 * data.weight[i]/max);
                col.push(fullColorHex(temp, temp, 113));
            }
            source.data.color = col;
        }
        if(b === 8){ //protein_abundance
            var temp = 0
            max = Math.max(...data.protein_abundance)
            for(i=0;i<len;++i){
                temp =  Math.floor(200 * data.protein_abundance[i]/max);
                col.push(fullColorHex(temp, 204, 113));
            }
            source.data.color = col;
        }
        if(b === 9){ //mRNA_abundance
            var temp = 0
            max = Math.max(...data.mRNA_abundance)
            for(i=0;i<len;++i){
                temp =  Math.floor(200 * data.mRNA_abundance[i]/max);
                col.push(fullColorHex(113, temp, temp));
            }
            source.data.color = col;
        }
        if(b === 10){ //decay_rate
            var temp = 0
            max = Math.max(...data.decay_rate)
            for(i=0;i<len;++i){
                temp =  Math.floor(200 * data.decay_rate[i]/max);
                col.push(fullColorHex(204, temp, 113));
            }
            source.data.color = col;
        }
        if(b === 11){ //is_secreted
            var temp = 0
            max = Math.max(...data.is_secreted)
            for(i=0;i<len;++i){
                temp =  Math.floor(200 * data.is_secreted[i]/max);
                col.push(fullColorHex(temp, 113, temp));
            }
            source.data.color = col;
        }
        if(b === 12){ //cys_elements
            var temp = 0
            max = Math.max(...data.cys_elements)
            for(i=0;i<len;++i){
                temp =  Math.floor(200 * data.cys_elements[i]/max);
                col.push(fullColorHex(temp,temp,76));
            }
            source.data.color = col;
        }
        if(b === 13){ //tot_volume
            var temp = 0
            max = Math.max(...data.tot_volume)
            for(i=0;i<len;++i){
                temp =  Math.floor(200 * data.tot_volume[i]/max);
                col.push(fullColorHex(52, 76, temp));
            }
            source.data.color = col;
        }
        if(b === 14){ //mean_hydrophobicity
            var temp = 0
            max = Math.max(...data.mean_hydrophobicity)
            for(i=0;i<len;++i){
                temp =  Math.floor(200 * data.mean_hydrophobicity[i]/max);
                col.push(fullColorHex(52, temp, 113));
            }
            source.data.color = col;
        }
        if(b === 15){ //glucose_cost
            var temp = 0
            max = Math.max(...data.glucose_cost)
            for(i=0;i<len;++i){
                temp =  Math.floor(200 * data.glucose_cost[i]/max);
                col.push(fullColorHex(temp, 52, temp));
            }
            source.data.color = col;
        }
        if(b === 16){ //synthesis_steps
            var temp = 0
            max = Math.max(...data.synthesis_steps)
            for(i=0;i<len;++i){
                temp =  Math.floor(200 * data.synthesis_steps[i]/max);
                col.push(fullColorHex(temp, temp, 76));
            }
            source.data.color = col;
        }
        if(b === 17){ //meanecai
            var temp = 0
            for(i=0;i<len;++i) {
                temp =  Math.floor(455 * data.meanecai[i]);
                col.push(fullColorHex(temp, temp, 113));
            }
            source.data.color = col;
        }
        if(b === 18){ //meancai
            var temp = 0
            for(i=0;i<len;++i) {
                temp =  Math.floor(455 * data.meancai[i]);
                col.push(fullColorHex(temp, temp, 113));
            }
            source.data.color = col;
        }
        if(b === 19){ //conservation
            var temp = 0
            max = Math.max(...data.conservation)
            for(i=0;i<len;++i){
                temp =  Math.floor(200 * data.conservation[i]/max);
                col.push(fullColorHex(temp, temp, 76));
            }
            source.data.color = col;
        }            
        source.change.emit();
            """
callback_plotall = """
            // JavaScript code goes here
            console.log("callback-------callback");
            var colors = ['#f39c12', "#1abc9c", "#3498db", "#2ecc71",
                        "#9b59b6", '#34495e', '#f1c40f','#e67e22', '#e74c3c', '#7f8c8d'];
            // the model that triggered the callback is cb_obj:
            var b = cb_obj.get("active");
            var rgbToHex = function(rgb){
              var hex = Number(rgb).toString(16);
              if (hex.length < 2) {
                   hex = "0" + hex;
              }
              return hex;
            };
            var fullColorHex = function(r,g,b) {
                var red = rgbToHex(r);
                var green = rgbToHex(g);
                var blue = rgbToHex(b);
                return "#"+red+green+blue;
            };
            // models passed as args are automagically available
            var data = source.data;
            var len = data.color.length;
            var col = [];
            var othercol = '#1abc9c';
            if(b === 0){
                source.data.color = source.data.homologies;
            }
            if(b ===2){
                var temp = 0
                for(i=0;i<len;++i){
                    temp =  Math.floor(50 + (156 * data.meanentropy[i]));
                    col.push(fullColorHex(52, 73, temp));
                    }
                source.data.color = col;
            }
            if(b === 1) {
                var temp = 0
                max = Math.max(...data.lengths)
                console.log(max)
                for(i=0;i<len;++i){
                    temp =  Math.floor(192 * (data.lengths[i]/max));
                    col.push(fullColorHex(temp, 152, 219));
                    }
                source.data.color = col;
            }
            source.change.emit();
            """

# CpCpGpApApTpApTpApTpTpCpCpGpApApTpApTpApTpTpCpCpGpApApTpApTpApTpTpCpCpGpApApTpApTpApTpTpTpTpCpCpGpApApTpApTpApTpTp
# GbGbCbTbTbAbTbAbTbAbAbGbGbCbTbTbAbTbAbTbAbAbGbGbCbTbTbAbTbAbTbAbAbGbGbCbTbTbAbTbAbTbAbAbAbGbGbCbTbTbAbTbAbTbAbAbAb


def homoyun(separation, folder="first500", homo_name="YAL019W",
            by='entropyLocation', aminolist=None):
    """
    read the all the files for one homology and returns everything inside (Yun)

    Args:
        folder: str the folder you wanan look onto
        homo_name: str the type of homology you are looking for
        by: str flag which of entropy value or Avalue (entropyLocation) you want
        aminolist: list[str] the amino acids you are interested in. None to say All

    Returns:
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
    for i, allnan in enumerate(nanvalues.all(1)):
        if allnan:
            nanvalues = np.vstack(nanvalues[0:i], nanvalues[i + 1:])
            lenmat = np.vstack(lenmat[0:i], lenmat[i + 1:])
            gentab = np.vstack(gentab[0:i], gentab[i + 1:])
            species = np.vstack(species[0:i], species[i + 1:])
            doub = np.vstack(doub[0:i], doub[i + 1:])
    nanvalues = np.count_nonzero(nanvalues, 1)
    print "at homology " + homo_name
    return [gentab, species, nanvalues, lenmat, doub]


def getyun(key, val):
    """
    function to get data from an url 'val' and save it under a file 'key'

    Args:
        key: str name of file
        val: str url of file
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

    Args:
        key: str name of file
        val: str url of file
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

    Returns:
        pd dfs of the ordered names and the names with links
    """
    return pd.read_csv("utils/meta/order_name461.csv", header=None, names=['name', 'b']),\
        pd.read_csv("utils/meta/names_with_links.csv", header=None, names=['name', 'b'])


def loadfromensembl(homology, kingdom='fungi', sequence='cdna',
                    additional='type=orthologues', saveonfiles=False, normalized=False,
                    setnans=False, number=0, by="entropy", using="normal", getCAI=False):
    """
    Load from ensembl the datas required in parameters ( look at PyCUB.get_data for more information)
    returns a fully populated homology object.

    Args:
        homology: str the homology code
        additional: str additional information on the retrieved sequence
        kingdom: str flags the relevant kingdom of you current session [fungi,plants,bacteria, animals]
        sequence: str flags the type of sequence you consider the full genome is (coding or non coding or full) [cds, all, cda]
        by: str flags what type of computation should be done [entropy,frequency, entropylocation]
        normalized: bool to true if should we normalize the entorpy by length
        saveonfiles: bool to true if the retrieved data should be saved on a file
        setnans: bool to true if nans should be set to NaN instead of an avg value
        using: the method to compute the partition function if using entropy location
        getCAI: wether or not to compute CAI !! need to have called the corresponding function on 
            Pycub before hand !!

    Returns:
        a populated PyCUB.homology of the homology object by [names, taxons, full, lenmat, homocode, nans,
        KaKs_Scores, similarity_scores, proteinids, GCcount, geneids, refs, ecai, refgene, refprot,
        tot_volume, mean_hydrophobicity, glucose_cost, synthesis_steps, isoelectricpoint,cai,
        conservation, uncounted]
        OR None if the homology is empty

    Raises:
        ConnectionError: "tried 50 times but still not able to connect"

    """
    server = "http://rest.ensemblgenomes.org"
    print 'homology: ' + homology + ' : ' + str(number)
    ext = "/homology/id/" + homology + '?'
    if sequence is not None:
        # dna   cdna    cds ncrna   Protein EMBL    GENBANK MySQL   TSV GTF GFF3
        ext += 'sequence=' + sequence
    if kingdom is not None:
        ext += ';compara=' + kingdom
    if additional is not None:
        ext += ';' + additional
    try:
        r = requests.get(server + ext, headers={"Content-Type": "application/json"})
    except ConnectionError:
        print "problem at " + homology
        if number > 50:
            raise IOError("tried 50 times but still not able to connect")
        return loadfromensembl(homology, kingdom=kingdom, sequence=sequence,
                               additional=additional, saveonfiles=saveonfiles, normalized=normalized,
                               setnans=setnans, number=number + 1, by=by, using=using)
    if not r.ok:
        r.raise_for_status()
    data = r.json()['data'][0]['homologies']
    if not data:
        return None
    if saveonfiles:
        with open('utils/data/' + homology + '.json', "wb") as code:
            code.write(json.dump(data))
    species, GCcount, lenmat, H, nans, similarities, KaKs_Scores, taxons, proteinids,\
        geneid, ref, ecai, cai, refgene, refprot, vol, cost, hydrophob, synthcost, isoepoint, conservation, others = process(
            data, normalized=normalized, setnans=setnans, by=by, getCAI=getCAI)
    if by == 'entropyLocation':
        H = getloc(H, np.array(lenmat), using=using)
    # here we add two things into names but only as a temporary saving measures removed by the
    # application fo preprocessing in homoset.
    homo = h.homology(names=[species, taxons], full=H, lenmat=lenmat, homocode=homology,
                      nans=nans, KaKs_Scores=KaKs_Scores, similarity_scores=similarities,
                      proteinids=proteinids, GCcount=GCcount, geneids=geneid, refs=ref, ecai=ecai, cai=cai, refgene=refgene,
                      refprot=refprot, tot_volume=vol, mean_hydrophobicity=hydrophob, glucose_cost=cost,
                      synthesis_steps=synthcost, isoelectricpoint=isoepoint, conservation=conservation, othercods=others)
    homo.order(withtaxons=True)  # a first ordering of the data, usefull afterward in the preprocessing
    return homo


# ################# YUN Processing Pipeline  ####################################
# CpCpGpApApTpApTpApTpTpCpCpGpApApTpApTpApTpTpCpCpGpApApTpApTpApTpTpCpCpGpApApTpApTpApTpTpTpTpCpCpGpApApTpApTpApTpTp
# GbGbCbTbTbAbTbAbTbAbAbGbGbCbTbTbAbTbAbTbAbAbGbGbCbTbTbAbTbAbTbAbAbGbGbCbTbTbAbTbAbTbAbAbAbGbGbCbTbTbAbTbAbTbAbAbAb

def process(data, normalized=False, setnans=False, by='entropy', getCAI=False):
    """
    function used by loadfromensembl() to process the retrieved data

    Args:
        data: list[str] the list of codons in this sequence
        normalized: to normalize of not the data (entropy)
        setnans: to show Nans or replace them
        by: 'entropy' the CUB metric to use
        getCAI: wether or not to compute CAI !! need to have called the corresponding function on 
            Pycub before hand !!

    Returns:
        species: list[str] the name of each homologous species
        GCcounts: np.array[float] the GC count for each homologous genes
        lenmat: np.array[int] the number of codons encoding each amino acids
        H: np.array[float] the CUB value for each homologous gene
        nans: np.array[int] the number of unexistant amino acid in the sequence
        KaKs_Scores: np.array[float] a proxy of the similarity of the homologous genes to the reference one
        similarities: np.array[float] the similarity of the homologous genes to the reference one
        taxons: list[str] the taxons of each species
        proteinids: list[str] the id of each homologous protein
        geneid: list[str] the id of each homologous gene
        ref: np.array[float] the reference gene's CUB value
        ecai: np.array[float] the ecai of each homologous gene
        cai: np.arrau[float] the cai of each homologous genes
        conservation: str the total amino acid conservation of the protein
        uncounted: float frequency of uncounted elements
        id: str the source gene id
        protid: str the source protein id
        vol: int, the volume of the corresponding protein
        cost: float the glucose cost to synthetize the corresponding protein
        hydrophob: float the hydrophobicity of the corresponding protein
        synthcost: int the timesteps to synthetize the corresponding protein
        isoepoint: float the Pi value of the corresponding protein
        if any of them does not exist, will return None instead
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
    geneid = []
    ecai = []
    cailist = []
    others = 0
    uncounted = 0
    # TODO retrieve also the reference genome
    # and compute its CAI
    referencespecies = data[0]["source"]['align_seq'].encode('ascii', 'ignore').replace("-", "")
    uncounted = len(referencespecies) - (referencespecies.count('A') + referencespecies.count('T') +
                                         referencespecies.count('C') + referencespecies.count('G'))
    if uncounted:
        print "ref uncounted = " + str(uncounted)
        referencespecies = referencespecies.replace("Y", "T").replace("R", "G").replace("K", "G")\
            .replace("M", "A").replace('S', 'C').replace("W", "A").replace("B", "C").replace("D", "T")\
            .replace("H", "T").replace("V", "G").replace("N", "C")
    referencespecies = [referencespecies[i:i + 3] for i in range(0, len(referencespecies), 3)]
    ref, _, _ = computeyun(list(referencespecies), normalized=normalized, setnans=setnans, by=by)
    reference_index(list(referencespecies))
    vol, cost, hydrophob, synthcost, isoepoint, conservation = None, None, None, None, None, None
    vol, cost, hydrophob, synthcost, isoepoint, conservation = compute_meta(list(referencespecies))
    for n, dat in enumerate(data):  # each species
        # https://en.wikipedia.org/wiki/Ka/Ks_ratio

        if dat["dn_ds"] is not None:
            KaKs_Scores.append(dat["dn_ds"])
        dat = dat['target']
        if dat["perc_id"] is not None and dat["perc_pos"] is not None:
            similarities.append(dat["perc_id"] / dat["perc_pos"] if dat["perc_pos"] != 0 else 0)
        if dat["taxon_id"] is not None:
            taxons.append(dat["taxon_id"])
        if dat["protein_id"] is not None:
            proteinids.append(dat["protein_id"])
        species.append(dat['species'])
        if dat["id"] is not None:
            geneid.append(dat["id"])
        codseq = dat['align_seq'].encode('ascii', 'ignore').replace("-", "")
        GCcount = float(codseq.count('G') + codseq.count('C'))
        uncounted = (len(codseq) - GCcount) - (codseq.count('A') + codseq.count('T'))
        if uncounted:
            print "uncounted = " + str(uncounted)
            codseq = codseq.replace("Y", "T").replace("R", "G").replace("K", "G").replace("M", "A").replace('S', 'C')\
                .replace("W", "A").replace("B", "C").replace("D", "T").replace("H", "T").replace("V", "G").replace("N", "C")
            GCcount = float(codseq.count('G') + codseq.count('C'))
        GCcount = GCcount / len(codseq)
        codseq = [codseq[i:i + 3] for i in range(0, len(codseq), 3)]
        if getCAI:
            cai = computeCAI(list(codseq))
            cailist.append(cai)
        c, other = computeECAI(list(codseq))
        ecai.append(c)
        valH, len_i, nan = computeyun(codseq, normalized=normalized, setnans=setnans, by=by)
        if valH is None:  # in case the species was in fact empty (db problem sometimes)
            del KaKs_Scores[-1]
            del similarities[-1]
            del taxons[-1]
            del proteinids[-1]
            del species[-1]
            continue
        H.append(valH)
        nans.append(nan)
        lenmat.append(len_i)
        GCcounts.append(GCcount)
        others += other
    return species, np.array(GCcounts, dtype=int), np.array(lenmat, dtype=int), np.array(H),\
        np.array(nans), np.array(similarities) if len(similarities) != 0 else None,\
        np.array(KaKs_Scores) if len(KaKs_Scores) != 0 else None, taxons if len(taxons) != 0 else None,\
        proteinids if len(proteinids) != 0 else None, geneid, ref, np.array(ecai), np.array(cailist) if len(cailist) > 0 else None,\
        data[0]["source"]["id"], data[0]["source"]["protein_id"], vol, cost, hydrophob, synthcost, isoepoint,\
        conservation, float(others) / (n + 1)


def compute_meta(data):
    """
    gets all the metavalue from this gene sequence according ot the chemical informations
    of the amino acids composing it

    Args:
        data: list[str] the list of codons in this sequence

    Returns:
        vol: int, the volume of the corresponding protein
        cost: float the glucose cost to synthetize the corresponding protein
        hydrophob: float the hydrophobicity of the corresponding protein
        synthcost: int the timesteps to synthetize the corresponding protein
        isoepoint: float the Pi value of the corresponding protein
    """
    vol, cost, hydrophob, synthcost, conservation = 0, 0, 0, 0, 0
    iso = []
    for cod in data:
        if cod not in ['ATG', 'TGG', 'TAA', 'TAG', 'TGA']:
            synthsteps, glucosecost, hydrophoby, volume, isoepoint, conservat = amino2meta[codamino[cod]]
            vol += volume
            cost += glucosecost
            hydrophob += hydrophoby
            iso.append(isoepoint)
            synthcost += synthsteps
            conservation += conservat
    iso = np.array(iso)
    mean = iso.mean()
    tempsup = 10000
    tempmin = -10000
    for a in isoelectricpoint:
        if a < tempsup and a > mean:
            tempsup = a
        elif a > tempmin and a < mean:
            tempmin = a
    isoepoint = (tempmin + tempsup) / 2
    return vol, cost, hydrophob / len(data), synthcost, isoepoint, conservation


def reference_index(data, forcai=False):
    """
    compute he RCSU value of this gene

    this gene will be the reference gene from which we will compute the CAI
    of other genes can be used both for CAI and ECAI

    Args:
        data: list[str] the list of codons in this sequence
        forcai: set to true if it is to compute the indexcai which is used to compute the cai values
    """
    # RCSU values are CodonCount/((1/num of synonymous codons) * sum of
    # all synonymous codons)
    global indexcai
    global indexecai
    for k, amin in enumerate(amino):
        subcodons = codons[amin]
        nbcod = len(subcodons)  # replace Cleng
        count = np.zeros(nbcod)
        for i, cod in enumerate(subcodons):
            count[i] = data.count(cod)
        lengsubseq = count.sum()  # replace subSlength
        rcsu = []
        if nbcod != 1 and lengsubseq != 0:
            denominator = float(lengsubseq) / nbcod
            # calculate the RSCU value for each of the codons
            for i in range(len(subcodons)):
                rcsu.append(count[i] / denominator)
            # now generate the index W=RCSUi/RCSUmax:
            rcsu_max = max(rcsu)
            for i, codon in enumerate(subcodons):
                if forcai:
                    indexcai[codon] = rcsu[i] / rcsu_max
                else:
                    indexecai[codon] = rcsu[i] / rcsu_max


def computeCAI(data):
    """
    compute the CAI according to the reference index WHICH NEEDS TO BE PREVIOUSLY COMPUTED

    Args:
        data: list[str] the list of codons in this sequence

    Returns:
        the CAI
    """
    cai_value, cai_length = 0, 0
    # if no index is set or generated, the default SharpEcoliIndex will
    # be used.
    for codon in data:
        if codon not in ['ATG', 'TGG', 'TAA', 'TAG', 'TGA']:
            # these two codons are always one, exclude them:
            val = indexcai.get(codon, 1)
            if val != 0:
                cai_value += math.log(val)
                cai_length += 1

    return math.exp(cai_value / (cai_length - 1.0))


def computeECAI(data):
    """
    compute the Evolutionary CAI according to the reference index WHICH NEEDS TO BE PREVIOUSLY COMPUTED

    Args:
        data: list[str] the list of codons in this sequence

    Returns:
        the ECAI, the amount of codons that were in ['ATG', 'TGG', 'TAA', 'TAG', 'TGA']
    """
    cai_value, cai_length, other = 0, 0, 0
    # if no index is set or generated, the default SharpEcoliIndex will
    # be used.
    for codon in data:
        if codon not in ['ATG', 'TGG', 'TAA', 'TAG', 'TGA']:
            # these two codons are always one, exclude them:
            val = indexecai.get(codon, 1)
            if val != 0:
                cai_value += math.log(val)
                cai_length += 1
        else:
            other += 1

    return math.exp(cai_value / (cai_length - 1.0)), other


def computeyun(data, setnans=False, normalized=False, by='entropy'):
    """
    function used by process to compte the 'entropy' or 'frequency' of the processed data

    Args:
        data: list[str] the list of codons in this sequence
        by: str flags what type of computation should be done [entropy,frequency]
        normalized: bool to true if should we normalize the entorpy by length
        setnans: bool to true if nans should be set to NaN instead of an avg value

    Returns:
        the valH (CUB values array)list[float], the length of
        each amino acid codon set len_i (list[int]), the number of nans(int), the GCcount (float)

    """
    global CUBD
    if by != 'frequency':
        valH = np.zeros(len(amino))
    if by != 'entropy':  # the number of codons usefull
        CuF = np.zeros(59)
    CUBD = len(amino) if by != 'frequency' else 59
    len_i = []
    nans = 0
    pos = 0
    for k, amin in enumerate(amino):
        subcodons = codons[amin]
        nbcod = len(subcodons)  # replace Cleng
        count = np.zeros(nbcod)
        X = np.zeros(nbcod)
        mn = np.ones(nbcod) / nbcod
        for i, cod in enumerate(subcodons):
            count[i] = data.count(cod)
        """
        for val in data:
            for i, cod in enumerate(subcodons):
                print val, cod
                print val == cod
                if val == cod:
                    count[i] += 1
                    pad.pop(j)
                    j -= 1
                    break
            j += 1
        """
        lengsubseq = count.sum()  # replace subSlength
        len_i.append(lengsubseq)
        if by == 'frequency'or by == "entropy" + "frequency":
            if lengsubseq == 0:
                CuF[pos:pos + nbcod] = np.NaN if setnans else 1. / nbcod
                nans += 1
            else:
                E = count / lengsubseq
                CuF[pos:pos + nbcod] = E
            pos += nbcod
        if by == "entropy" or by == "entropy" + "frequency":
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
    if by == 'frequency':
        valH = CuF
    elif by == 'entropy' + 'frequency':
        if nans == k:
            return None, None, None, None
        return valH, CuF, len_i, nans / 2
    else:
        if nans == k:
            return None, None, None
        return valH, len_i, nans


def getloc(valH, geneleng, using='computejerem'):
    """
    the function to compute the entropy location (adpated from Yun Deng's code University of Kent 2018)
    need the entropy values from the genes an homology and the length values as well.

    Args:
        valH: the ordered entropy values
        geneleng: the corresponding length values
        using: str flags the function to use from computejerem, permutation, full, normal

    Returns:
        the corresponding computed Avalue
    """
    numberindiv = len(valH)  # works like a shape[0]
    acodon = [[2, 3, 4, 5, 6, 8, 11, 12, 16],
              [9], [0, 7, 13, 15, 17], [1, 10, 14]]
    val = np.zeros((18, numberindiv))
    for ind, nbcod in enumerate([2, 3, 4, 6]):
        lengs = []
        for y in acodon[ind]:
            lengs.extend(geneleng[:, y])
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
            if leng < nbcod:
                valHloc[y] = 0
                continue
            posx, posy = divmod(argu[y], numberindiv)
            # we access the position knowing that we need to convert row number
            # into the possible position on the entropy matrix
            E = valH[posy, acodon[ind][posx]]
            # if H is also superior or inferior to the threshold, according to the
            # computation of the entropy location, both their values should be the
            # same.
            # print '{0:.2f}%\r'.format(float(((y + (ind * len(lengs))) * 100)) / (len(lengs) * 4)),
            if lt != leng:
                print str(nbcod) + ',' + str(leng), str(leng**nbcod),
                lt = leng
                mn = np.ones(nbcod) / nbcod
                X = np.zeros(nbcod)
                # efor part
                i, div = divmod(leng, nbcod)
                X[:div] = np.ceil(i) + 1
                X[div:] = np.floor(i)
                Eg = multinomial.pmf(x=X, n=leng, p=mn)
                # end
                if Eg == 0:
                    valHloc[y] = 0
                    continue
                ref = computepartition(nbcod, leng, using=using)
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
                bdist = np.cumsum(np.multiply(hist, np.divide(Eg, np.exp(np.multiply(edgesnew, leng)))))
                del hist
            if E > edgesnew[-1]:
                location = last_match_index(bdist, 1)
            else:
                location = np.argmax(edgesnew >= E)
                if location != 1:
                    location_value = bdist[location]
                    if location_value > 0.95 and location_value < 1:
                        location -= 1
                    if location_value >= 1:
                        location = last_match_index(bdist, 1)
            valHloc[y] = bdist[location]
        # we reset the right ordering
        valHloc[argu] = valHloc[:]
        for y, e in enumerate(acodon[ind]):
            # we reset the position of the aminoacid numbers.
            val[e, :] = valHloc[numberindiv * y:numberindiv * (y + 1)]
    return val.T

# CpCpGpApApTpApTpApTpTpCpCpGpApApTpApTpApTpTpCpCpGpApApTpApTpApTpTpCpCpGpApApTpApTpApTpTpTpTpCpCpGpApApTpApTpApTpTp
# GbGbCbTbTbAbTbAbTbAbAbGbGbCbTbTbAbTbAbTbAbAbGbGbCbTbTbAbTbAbTbAbAbGbGbCbTbTbAbTbAbTbAbAbAbGbGbCbTbTbAbTbAbTbAbAbAb


def computepartition(nbcod, leng, using='computejerem'):
    """
    according to a 'using' parameter, will select the partition function algorithm to use
        jeremcompute: bool to true if use a jeremie twist to this function
        all function are made and found by jeremie kalfon, derived from earlier
        work by Yun Deng Phd student at the University of Kent 2O18

    Basically works by computing every possible codon presence for a defined amino acid (and it number of
    repetitions)
    done by computing the pdf from a multinomial according to those values

    Params:
        leng: int length of sequence
        nbcod: int number of codons
        using: str flags the function to use from computejerem, permutation, full, normal

    Returns:
        the partition array from the selected partition function

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
        raise AttributeError(" give a Using from [permutation, random, normal, full]")


def randomdraw(nbcod, leng):
    """
    one of the four partition functions algorithms
    this one either use the regular Yun version or try to sort the values to limit the number
    of sampling to do and do a limited dynamic programming version.

    Basically works by computing a random subset of every possible codon presence
    for a defined amino acid (and it number of repetitions)
    done by computing the pdf from a multinomial according to those values
    It will use a specific dynamic programming trick to increase the computation by more than 200*

    Params:
        leng: int length of sequence
        nbcod: int number of codons

    Returns:
        the partition array
    """
    global listvect
    # we iterate this way to have only a max value or everything thanks to a
    # sampling without replacement
    print 'randomdraw\r',
    if listvect is None:
        listvect = np.zeros((MAXITR, nbcod))
        for ind in xrange(MAXITR):
            prevect = np.zeros(nbcod)
            prevect[1:] = np.sort(randint(0, leng + 1, nbcod - 1))
            for i in range(nbcod - 1):
                prevect[i] = prevect[i + 1] - prevect[i]
            prevect[-1] = leng - prevect[-1]
            listvect[ind] = prevect
    elif nbcod != listvect.shape[1]:
        listvect = np.zeros((MAXITR, nbcod))
        for ind in xrange(MAXITR):
            prevect = np.zeros(nbcod)
            prevect[1:] = np.sort(randint(0, leng + 1, nbcod - 1))
            for i in range(nbcod - 1):
                prevect[i] = prevect[i + 1] - prevect[i]
            prevect[-1] = leng - prevect[-1]
            listvect[ind] = prevect
    else:
        add = leng - listvect[0].sum()
        posx = np.arange(0, MAXITR)
        posy = randint(0, nbcod - 1, MAXITR)
        listvect[posx, posy] += add

    return multinomial(leng, np.ones(nbcod) / nbcod).pmf(listvect)


def computepartition_sorted_full(nbcod, leng):
    """
    one of the four partition functions algorithms
    Basically works by computing every possible codon presence for a defined amino acid (and it number of
    repetitions)
    done by computing the pdf from a multinomial according to those values

    Params:
        leng: int length of sequence
        nbcod: int number of codons

    Returns:
        the partition array
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

    Params:
        leng: int length of sequence
        nbcod: int number of codons

    Returns:
        the partition array
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

    Params:
        leng: int length of sequence
        nbcod: int number of codons

    Returns:
        the partition array

    """
    if mlen(leng, nbcod) == 'full':
        # if we are ok to do full method
        print "computejerem\r",
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


def computepartition_normal_approximation(nbcod, leng, probavector):
    """
    one of the 4 partition function algorithm, works by
    approximates the binomial distribution with a multidimensional gaussian function
    as a gaussian approximates a binomial under many trials and this is basically what we do.

    Params:
        leng: int length of sequence
        nbcod: int number of codons
        probavector: the probability vector of the codons for this amino acid

    Returns:
        the partition array
    """
    # TODO: debug the singular matrix problem
    pdb.set_trace()
    mean = probavector * leng
    ker = np.kron(probavector, probavector).reshape(nbcod, nbcod) * (-leng)
    np.fill_diagonal(ker, (probavector * (1 - probavector)) * leng)
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

# CpCpGpApApTpApTpApTpTpCpCpGpApApTpApTpApTpTpCpCpGpApApTpApTpApTpTpCpCpGpApApTpApTpApTpTpTpTpCpCpGpApApTpApTpApTpTp
# GbGbCbTbTbAbTbAbTbAbAbGbGbCbTbTbAbTbAbTbAbAbGbGbCbTbTbAbTbAbTbAbAbGbGbCbTbTbAbTbAbTbAbAbAbGbGbCbTbTbAbTbAbTbAbAbAb


def mlen(leng, nb):
    """
    compute wether we need to do the full version of the algorithm or the
    sampling version (see documentation of utils for the threshold values)
    This function takes three parameters:
    1. String
    2. Starting index of the string
    3. Ending index of the string.

    Args:
        leng: int length of sequence
        nbcod: int number of codons

    Returns:
        'samp' or 'full' according to the max size
    """
    if nb == 6:
        return 'samp' if leng > smax[nb] else 'full'
    elif nb == 4:
        return 'samp' if leng > smax[nb] else 'full'
    elif nb == 3:
        return 'samp' if leng > smax[nb] else 'full'
    elif nb == 2:
        return 'samp' if leng > smax[nb] else 'full'


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
    """
    the last place in the array where this value is
    above all others

    Args:
        value: float, the value from which to find the last matching index
        a: the array

    Returns:
        the index (int)
    """
    idx = np.array(np.unravel_index(((a < value)[::-1]).argmax(), a.shape))
    return a.shape - idx - 1


def rgb2hex(rgb):
    """
    compute the hex str from an rgd triplet

    Args:
        rgb: triplet of ints between 0-255

    Returns:
        the corresponding hex str
    """
    return '#%02x%02x%02x' % rgb


def endresdistance(a, b):
    """
    a statistical distance measure between two random vectors
    """
    me = (a + b) / 2
    return sqrt(kl(a, me) + kl(b, me))


def kl(a, b):
    """
    kulback-leiber distance computation
    """
    return (np.log(a / b) * a).sum()


class dotdict(dict):
    """
    dot.notation access to dictionary attributes

    May be used to pass a dict within it and access it by ".value" instead of "[value]"
    """
    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__
