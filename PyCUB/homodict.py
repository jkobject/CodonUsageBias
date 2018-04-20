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
    isload = False
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
		self.isload = fast
		self.path += session
    	self.store = dict()
    	self.update(val)

    def __getitem__(self, key):
        if self.isload:
            return self.store[key]
        else:
            return self._load(key)

    def __setitem__(self, key, value):
    	if self.isload:
    		self.store[key] = value
        else:
        	if key in self.keys:
    			# writing on a file clears it.
            	self._save(key, value)

    def __delitem__(self, key):
    	if self.isload:
        	del self.store[key]
        else:
        	path = self.path + key + ".json"
			os.remove(path)
		keys.remove(key)

    def __del__(self):
    	if not isload:
			for key in self.keys:
				del self[key]
		del self.path
		del self.store 
		del self.keys
		del self.isload


    def __iter__(self):
    	if self.isload:
        	return iter(self.store)
        else:
        	return iter(self.keys)

    def iteritems(self):
    	if self.isload:
    		return self.store.iteritems()
    	else:
    		return iter([(key,self._load(key)) for key in keys])

    def __len__(self):
        return len(self.keys)

    def update(self, val):
    	if self.isload:
        	self.store.update(val)
        	for key, _ in val.iteritems():
   		     	self.keys.append(key)
        else: 
        	for key, value in val.iteritems():
        		self.save(key,val)
        		self.keys.append(key)

#######################################################

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


###########################################################

     def compute_fast(self):
    	if not self.isload:
		    print "will compute faster but take a lot of memory.. this might take time"
		    for key in self.keys:
		        self.update(dict(key, self._load(key)))
		    self.isload = True
		else:
			print 'already done'

    def offload(self):
    	if self.isload:
    		print "will free up memory but access to this dict will take time (this takes time as well"
    		for key in keys:
    			self._save(self,key)
    		self.isload = False


    def _dictify(self):
    	if not self.isload:
    		self.offload()
    	return self.path

    def retdict(self):
    	print "you are retrieving the full dict"
    	if not self.isload:
    		self.compute_fast()
    	return self.store
        	
