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

ALLOWED_TYPES = ["invpke", "pke"]

class container:
    
    def __init__(self):
        pass


class resultscontainer:
    
    
    def __init__(self, results, typ="pke"):
        
        _inlist(typ, "results type: {}".format(typ), ALLOWED_TYPES)
        self.typ = typ
        
        if type(results) is str:
            self.__recover(results)
        else:
            self.__dict__ = results.__dict__.copy()
            self.typ = typ
    
    
    def export(self, filename):
        
        _isstr(filename, "h5 export filename {}".format(filename))
        
        if self.typ == "pke":
            self.__exportpke(filename)
        elif self.typ == "invpke":
            self.__exportinvpke(filename)
    
    
    def __exportGroup(self, group, ATTR, obj=None):
        """function exports a groups results data to hdf5 file"""
        
        for i in ATTR:
            try:
                _isstr(i, "Attribute")
                if obj is not None:
                    data = getattr(obj, i)
                else:
                    data = getattr(self, i)
                if type(data) in [np.ndarray, list]:                 
                    if type(data) is list: data = np.asarray(data)
                    group.create_dataset(i, data=data, dtype=str(data.dtype))
                elif type(data) is str:
                    group.create_dataset(i, data=data.encode("ascii",
                                                             "ignore"))
                else:
                    group.create_dataset(i, data=data, dtype=type(data))
            except:
                print("WARNING: {} could not be exported!".format(i))
    
    
    def __buildGroup(self, f, key, attrs):
        """function reconstructs a single groups results data from hdf5 file"""
        
        con = container()
        
        for i in attrs:
            try:
                data = f[key][i][()]
                if type(data) is bytes:
                    data = str(data, "utf-8")
                setattr(con, i, data)
            except:
                print("{} not found in results".format(i))
        
        setattr(self, key, con)
    
    
    def __exportpke(self, filename):
        
        from kinetics.functions.pointkinetics import PKE_INPUTS_DICT
        PKE_OUT_LIST = ["nt", "dnt", "power", "rho", "flux", "timepoints"]
        
        with h5py.File(filename, "w") as f:
            
            #export inputs to hdf5 file
            self.__exportGroup(f.create_group("inputs"),
                               list(PKE_INPUTS_DICT.keys()),
                               obj=getattr(self, "inputs"))
            #output outputs to hdf5 file
            self.__exportGroup(f.create_group("outputs"), PKE_OUT_LIST,
                               obj=getattr(self, "outputs"))
    
    
    def __exportinvpke(self, filename):
        
        from kinetics.functions.inversepointkinetics import INV_PKE_DICT
        INV_OUT_LIST = ["nt", "dnt", "power", "rhototal", "flux", "timepoints"]
        
        with h5py.File(filename, "w") as f:
            
            #export inputs to hdf5 file
            self.__exportGroup(f.create_group("inputs"),
                               list(INV_PKE_DICT.keys()),
                               obj=getattr(self, "inputs"))
            #output outputs to hdf5 file
            self.__exportGroup(f.create_group("outputs"), INV_OUT_LIST,
                               obj=getattr(self, "outputs"))
    
    
    def __recoverpke(self, filename):
        
        from kinetics.functions.pointkinetics import PKE_INPUTS_DICT
        PKE_OUT_LIST = ["nt", "dnt", "power", "rho", "flux", "timepoints"]
        
        with h5py.File(filename, "r+") as f:
            self.__buildGroup(f, "inputs", list(PKE_INPUTS_DICT.keys()))
            self.__buildGroup(f, "outputs", PKE_OUT_LIST)
    
    
    def __recoverinvpke(self, filename):
        
        from kinetics.functions.inversepointkinetics import INV_PKE_DICT
        INV_OUT_LIST = ["nt", "power", "rhototal", "flux", "timepoints"]
        
        with h5py.File(filename, "r+") as f:
            self.__buildGroup(f, "inputs", list(INV_PKE_DICT.keys()))
            self.__buildGroup(f, "outputs", INV_OUT_LIST)
    
    
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
        
        if self.typ == "pke":
            self.__recoverpke(filename)
        elif self.typ == "invpke":
            self.__recoverinvpke(filename)
            
        