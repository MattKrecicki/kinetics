# -*- coding: utf-8 -*-
"""outputs.py

output containers for various kinetic solvers 

@author: matt krecicki
@email: matthewkrecicki@gmail.com

"""

import numpy as np
import h5py


class pointkineticsoutputs:
    
    
    def __init__(self, **kwargs):
        """function initalizes point kinetics container"""
        
        self.__dict__.update(kwargs)
    
    
    def export(self, filename):
        pass
    
    
    def recover(self, filename):
        pass
    
        
