""""
Created by Jeremie KALFON
Date : 9 Apr 2018
University of Kent, ECE paris
jkobject.com

"""


import json
import os
import glob

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import espece as spe
import utils
import homology as h

import collections


class Homodict(collections.MutableMapping):


"""docstring for Homodict

    Object that acts as a dictionary but is not as fast

    params:
    ------

    """
    path = 'utils/save/'
    is_loaded = False
    keys = []

###########################################

    @staticmethod
    def __init__(self, session, fast=False, val={}):
    	"""
		Params:
		------
		session: String, the name of the current session
		fast: Bool, if we want to have everything in the ram or in the disk 
		(fast means in the ram)
		val: a dictionnary if we want to initialize it with that.

    	"""
		self.is_loaded = fast
		self.path += session
    	self.store = dict()
    	self.update(val)
        

    def __getitem__(self, key):
        if self.is_loaded:
            return self.store[key]
        else:
            return _load(self, key)

    def __setitem__(self, key, value):
    	if not self.is_loaded:
    		if key in self.keys:
            	self._save(key, value)
        else:
        	self.store[key] = value

    def __delitem__(self, key):
    	if self.is_loaded:
        	del self.store[key]
        else:
        	path = self.path + key + ".json"
			os.remove(path)

    def __del__(self):
    	if not is_loaded:
			for key in self.keys:
				del self[key]
		del self.path
		del self.store 
		del self.keys
		del self.is_loaded


    def __iter__(self):
    	if self.is_loaded:
        	return iter(self.store)
        else:
        	return self.keys

    def __len__(self):
        return len(self.keys)

    def update(self, val):
    	if self.is_loaded:
        	self.store.update(val)
        	for key, _ in val.iteritems():
   		     	self.keys.append(key)
        else: 
        	for key, value in val.iteritems():
        		self.save(key,val)
        		self.keys.append(key)

#######################################################

     def compute_fast(self):
    	if not self.is_loaded:
		    print "will compute faster but take a lot of memory.. this might take time"
		    for key in self.keys:
		        self.update(dict(key, self._load(key)))
		    self.is_loaded = True
		else:
			print 'already done'

    def offload(self):
    	if self.is_loaded:
    		print "will free up memory but access to this dict will take time (this takes time as well"
    		for key in keys:
    			self._save(self,key)
    		self.is_loaded = False

    def _save(self, key, val):
    	path = self.path + key + ".json"
        with open(path "w") as f:
            jsoni = f.write(val)

    def _load(self, key):
        path = self.path + key + ".json"
        with open(path "r") as f:
	        jsoni = f.read()
	        os.remove(path)
	        return h.homology(data=jsoni)


    def _dictify(self):
    	if not self.is_loaded:
    		self.offload()
    	return self.path
        	
