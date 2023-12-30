# -*- coding: utf-8 -*-
"""outputs.py

output containers for various kinetic solvers 

@author: matt krecicki
@email: matthewkrecicki@gmail.com

"""

import h5py
from kinetics.errors.checkerrors import _isstr

# -----------------------------------------------------------------------------
# ----- point kinetics solver outputs
# -----------------------------------------------------------------------------


class pointkineticscontainer:
    """class contains outputs from point kinetics simulation and allows user to
    export and recover outputs"""
    
    PKE_OUTPUTS = ["nt", "power", "dnt", "timepoints", "rho"]
    SRC_PKE_OUTPUTS = []
    
    
    def __init__(self, **kwargs):
        """function initalizes point kinetics container"""
        
        self.__dict__.update(kwargs)

    
    def export(self, filename):
        """function exports solution to hdf5 file for storage"""
        
        #make sure given file name is a string
        _isstr(filename, "hdf5 output file name")
        
        #determine which set of values need to be exported
        if self.typ == "PKE":
            keys = self.PKE_OUTPUTS
        elif self.typ == "SPKE":
            keys = self.SRC_PKE_OUTPUTS

        
        #export keys to h5 file
        with h5py.File(filename, "w") as f:
            for i in range(len(keys)):
                f.create_dataset(keys[i], data=getattr(self, keys[i]))
            
                
    
    def recover(self, filename):
        """function recovers solution from hdf5 file"""
        
        _isstr(filename, "hdf5 output file name")
        
        with h5py.File(filename, "r+") as f:
            
            keys = list(f.keys())
            
            for key in keys:
                setattr(self, key, f[key][()])
            
    
        
