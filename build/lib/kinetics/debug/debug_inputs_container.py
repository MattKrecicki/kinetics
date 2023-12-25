# -*- coding: utf-8 -*-
"""
Created on Tue Jan 10 10:20:49 2023

@author: matt krecicki
"""

import numpy as np
from ntpSystem.containers.inputs import inputsContainer


inp = inputsContainer(verbose=True, typ="steadystate")


inp.addFuelGeom(noFuelElems = 127, fuelL = 1.20, fuelFlowA = 2.0253e-04,
                fuelHper = 3.0687e-01)

