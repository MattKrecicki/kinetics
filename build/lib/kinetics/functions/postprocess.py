# -*- coding: utf-8 -*-
"""
Created on Sun May 19 08:33:30 2024

@author: matt krecicki
@email: matthewkrecicki@gmail.com

functions that enable the user to storage transeint results and recover the
results for later post-processing


"""

import h5py
import numpy as np
from kinetics.errors.checkerrors import _isstr, _inlist

ALLOWED_TYPES = ["invpke"]


class container:

    
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)


class resultscontainer:
    
    
    def __init__(self, results, typ="pke"):
        
        _inlist(typ, "results type: {}".format(typ), ALLOWED_TYPES)
        self.typ = typ
        
        if type(results) is str:
            self.__recover(results)
        else:
            self.inputs = results.inputs
            self.outputs = results.outputs
    
    
    def export(self, filename):
        
        _isstr(filename, "h5 export filename {}".format(filename))
        
        if self.typ == "invpke":
            a=1
            
        
    
    def __recover(self, filename):
        """function recovers result from previous simulation
        

        Parameters
        ----------
        filename : str
            path to h5 file containing results.

        Returns
        -------
        None.

        """
        _isstr(filename, "h5 results file name {}".format(filename))