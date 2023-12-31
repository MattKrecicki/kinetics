# -*- coding: utf-8 -*-
"""customerrors.py

A set of custom error raises to identify when the error occurs.

"""

import numpy as np
from kinetics.errors.checkerrors import _ispositive, _isobject, _isequallength

# -----------------------------------------------------------------------------
# Run error checking for inputs container
# -----------------------------------------------------------------------------

#dict containing inputs descriptions and error checking info
#format is the following:
#   idx 0 is the first data type
#   idx 1 is the second data type, which is required if first data type is
#   a np.ndarray or a list
#   idx 2 is the input key description
#   idx 3 is the units of the input
#   idx 4 is positive value required error checking flag


def _checkdict(Dict, inputObj):
    """function runs input object error checking based on input dict"""
    
    for key in list(Dict.keys()):
        var = getattr(inputObj, key)
        _isobject(var, Dict[key][0], Dict[key][2])
        
        if Dict[key][0] in [list, np.ndarray]:
            for val in var:
                _isobject(val, Dict[key][1], Dict[key][2])
                if Dict[key][4]:
                    _ispositive(val, Dict[key][2])
        elif Dict[key][0] is not object:
            _ispositive(var, Dict[key][2])


def _pkecheck(inputObj):
    
    #make sure beta and lamda are the same length
    _isequallength(getattr(inputObj, "lamda"), len(getattr(inputObj, "beta")),
                   "decay neutron precusor decay constant and yield fractions")
    
    



