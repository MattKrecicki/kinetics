# -*- coding: utf-8 -*-
"""
Created on Mon Jan 17 13:08:44 2022

@author: matt krecicki
@email: mkrecicki@gatech.edu

"""

import h5py
from ntptransient.utilities.checkerrors import _isstr, _inlist


INT_SOL_ARR_KEYS = ["Power", "PrecursorDensity1", "PrecursorDensity2",
                    "PrecursorDensity3", "PrecursorDensity4", "PrecursorDensity5",
                    "PrecursorDensity6", "TFuel", "I135", "Xe135", "Pm149",
                    "Sm149", "U235", "U238", "RhoFuelTemp", "RhoFuelGasDensity",
                    "RhoMESGasDensity",  "RhoMERGasDensity",  "RhoXe135",
                    "RhoSm149", "DrumAngle", "RhoConDrum", "J4AEnthalpy",
                    "J5AEnthalpy",  "J4BEnthalpy",  "J5BEnthalpy", "J6Enthalpy", 
                    "J7Enthalpy", "J8Enthalpy", "MixEnthalpy", "J10Enthalpy",
                    "J11Enthalpy", "J12Enthalpy", "DemandPressure",
                    "DemandTemperature"]


NONINT_SOL_ARR_KEYS = ["Time", "MassFlow", "TurbineMassFlow", "BypassMassFlow",
                       "J3APressure", "J3ATemperature", "J4APressure",
                       "J4ATemperature", "J5APressure", "J5ATemperature",
                       "J3BPressure", "J3BTemperature", "J4BPressure",
                       "J4BTemperature", "J5BPressure", "J5BTemperature",
                       "J6Pressure", "J6Temperature", "J7Pressure",
                       "J7Temperature", "J8Pressure", "J8Temperature",
                       "J9Pressure", "J9Temperature", "J10Pressure",
                       "J10Temperature", "J11Pressure", "J11Temperature",
                       "J12Pressure", "J12Temperature", "PumpSpeed",
                       "TBCVPosition", "PumpPower", "TurbinePower"]


class export:
    """this class exports the solution arrays to an hdf5 file"""
    
    def __init__(self, filename, intSolArr, nonintSolArr, intSolTime):
        
        _isstr(filename, "hdf5 output file name")
        self.filename = filename
        
        with h5py.File(self.filename, "w") as f:
            
            #write time integrated quanities to hdf5 file
            intSol = f.create_group("time integrated solution")
            self._writesoluionArrays(intSol, intSolArr, INT_SOL_ARR_KEYS,\
                timeArr=intSolTime)
            
            #write time non-integrated qunaities to hdf5 file
            nonIntSol = f.create_group("time non-integrated solution")
            self._writesoluionArrays(nonIntSol,  nonintSolArr,\
                NONINT_SOL_ARR_KEYS)


    def _writesoluionArrays(self, group, solArr, keys, timeArr=False):
        """function writes dictionary to hdf5 file subset"""
        
        if type(timeArr) is not bool:
            group.create_dataset("time", data=timeArr)
        for i in range(len(keys)):
            group.create_dataset(str(keys[i]), data=solArr[:,i])
    

class recover:
    """class recovers analysis from hdf5 output file"""
    
    def _recoverSolutions(self, results, dictKey):
        """function recovers solution"""
        
        res = results[dictKey]
        for i in res.keys():
            setattr(self, i, res[i][()])
    
    
    def __init__(self, h5file, intKeys=NONINT_SOL_ARR_KEYS,
                 nonintKeys=NONINT_SOL_ARR_KEYS):
        """initalization of results object"""
        
        _isstr(h5file, "output file name")
        self.filename = h5file
        self._intKeys = intKeys
        self._nonintKeys = nonintKeys
        
        with h5py.File(self.filename, "r+") as f:
            keys = list(f.keys())
            _inlist('time integrated solution', "hdf5 results keys", keys)
            _inlist('time non-integrated solution', "hdf5 results keys", keys)
            self._recoverSolutions(f, 'time integrated solution')
            self._recoverSolutions(f, 'time non-integrated solution')

    
    
    
    
    
        