# -*- coding: utf-8 -*-
"""

Generate saturated properties using CoolProp

"""


from ntpSystem.functions.generatesaturationprops import SaturationProperties
from ntpSystem.functions.propertytable import SAT_HYDROGEN_FILE



# -----------------------------------------------------------------------------
# Uncomment if you want to generate saturated properties
# -----------------------------------------------------------------------------


SaturationProperties(nvals=1000, coolant='Hydrogen',
                     Pmin=7357.83, Pmax=1.2964e+06, 
                     outputFile=SAT_HYDROGEN_FILE)
    



