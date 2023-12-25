# -*- coding: utf-8 -*-
"""add_materials

Create a new database with hydrogen properties as a function of 
pressure and enthalpy.

"""

from ntpSystem import setDataPath
from ntpSystem.functions.propertytable import PropertyTable,\
    ChangeDataDependencies
    
from ntpSystem.functions.materials import Material
from ntpSystem.functions.addmaterial import addMaterial

# There is a capability to add multiple materials at once
# ---------------------------------------------------------
# from ntpSystem.functions.materials import Materials
# from ntpSystem.functions.addmaterial import addMaterials


# -----------------------------------------------------------------------------
# Read the original file
# -----------------------------------------------------------------------------

H5_ORIG_FILE = "ThermoPhysicalProperties.h5"
# Read the table with hydrogen properties
table = PropertyTable(H5_ORIG_FILE)

# Read the original table to obtain pressure- and temperature-dependent pty
H2_PT = table.read("H2")    

# Change dependnecies
H2_PH = ChangeDataDependencies(H2_PT, xattr="P", yattr="h", 
                               gridNy=[[1E+2, 1E+3, 1E+4, 1E+5, 1E+06, 1E+07],
                                       [200, 200, 200, 200, 300, 300, 300]])



# -----------------------------------------------------------------------------
# Define a new H2 material
# -----------------------------------------------------------------------------

# H2_PH.T  -  represents the enthalpy in J/kg

mat_H2 =\
    Material("H2_PH", H2_PH.T, H2_PH.P, reference='same as H2')

# properties to be defined for the new material
ptynames = [ 'cp', 'cv', 'mol', 'my', 'r', 's', 'tc', 'v', 'T']

for pty in ptynames:
    vals = H2_PH.getpty(pty)
    if pty == 'T':
        pty = 'tempK'
    mat_H2.addproperty(pty, vals)


# -----------------------------------------------------------------------------
# Add the new material to the database
# -----------------------------------------------------------------------------
H5_NEW_FILE = "ThermoPhysicalProperties_04Apr23.h5"
h5pathOrig = setDataPath(H5_ORIG_FILE)
h5pathNew = setDataPath(H5_NEW_FILE)

# add the new H2_PH data into the database
addMaterial(mat_H2, newFile=h5pathNew, origFile=h5pathOrig, overwrtFile=False,
                copyMat=False, copyName=None)
