.. _newmaterialcard:

====================
Adding New Materials
====================


The user has the option of defining and adding new materials into a data library red by ``ntpSystem.``
The process includes:

- :ref:`Defining a material object <matobject>`
- :ref:`Adding the new object to an existing (or new) library <addmatobject>`
- :ref:`Existing library <exitinglibrary>` with allowed properties.
- :ref:`Manipulate data before adding to library <matmanipulate>` with allowed properties.

:ref:`A full example is included here <matfullexample>`

.. _matobject:

Material
---------


Define a new material ``Material`` with all its properties to
the material object. This function acts as an interface between the ``material`` and
``propertytable``.


Start by importing the following items:

.. code::

	from ntpSystem.functions.materials import Material


The **syntax** to define a new material is:

.. code::

	newmat = Material(matId, temperatures, pressures=None, reference=None)


where:

- ``matId`` name of the material (str)
- ``temperatures`` temperature points (array)
- ``pressures``  Pressure points (array). Optional card.
- ``reference`` reference string. Optional card.
	
The **object** includes an ``addproperty``:

.. code::

	newmat.addproperty(pty, vals)
	
where:

- ``pty`` name of the property (str)
- ``vals`` values for the specific property 1-D or 2-D array depending on whether only temperature or both temperatures and pressures are provided as dependencies. In the 2-D array the columns represent the temperatures while the rows represent the pressure values.


**Example-1**: pressure and temperature dependent properties

.. code::

	newmat = Material("newMat", np.array([300, 900, 1800]), np.array([10E+6, 11E+6]))
	newmat.addproperty("my", np.array([[15.0, 13.5, 9.0], [14.9, 13.4, 8.9]]))
	newmat.addproperty("tc", np.array([[15.0, 13.5, 9.0], [14.9, 13.4, 8.9]])) 


**Example-2**: only temperature dependent properties

.. code::

	UC = Material("UC", np.array([300, 900, 1800]))
	UC.addproperty("tc", np.array([15.0, 13.5, 9.0]))



.. _addmatobject:

Adding Material Object to Library
---------------------------------

The ``ntpSystem`` allows to add a new material set, defined by the ``Material`` class, to an existing or a new h5 file/database.
The user has the flexibility to modify or replace the original h5 file. In addition, the user can remove or store the original material sets.


Start by importing the following items:

.. code::

	from ntpSystem.functions.addmaterial import addMaterial


The **syntax** to add the new material to the library is:

.. code::

	addMaterial(material, newFile, origFile, overwrtFile,
              copyMat, copyName)


where:

- ``material`` : Material object that contains the information about a certain material. Includes temperature-pressure dependencies, properties (e.g. tc), and material Id
- ``newFile`` : str name of the new h5 data file to store all the materials
- ``origFile`` : str, optional name of the original h5 data file
- ``overwrtFile`` : bool, optional a flag to indicate whether the original file must be overwritten. False as default. This parameter acts as a layer of protection in case newFile already exists. In such a case, overwrtFile must be set to True, otherwise an error will be raised.
- ``copyMat`` : bool, optional a flag to indicate whether an existing material should be saved under a different name. For example if UO2 already exists and the user defines a new UO2. In such a case, copyMat=False deletes the original UO2, otherwise it stores it under a different name. True as default
- ``copyName`` : str, optional the name under which an existing material should be saved. The default relies on ORIG_COPY_<material-name>, e.g. ORIG_COPY_UO2.


	
**Example**:

.. code::

	UC = Material("UC", np.array([300, 900, 1800]))
  UC.addproperty("tc", np.array([15.0, 13.5, 9.0]))
  addMaterial(material=UC,
               newFile="../tests/ThermoPhysicalPropertiesMod.h5",
               origFile="../database/ThermoPhysicalProperties.h5",
               copyMat=True)


.. _exitinglibrary:

Existing Data Library
---------------------

The package comes with a built in library that contains the
following materials in a source denoted here as the **table**:

	* H2
	* H2_PH
	* Molybdenum
	* Tungsten
	* UC
	* UN
	* UO2
	*	Zircaloy
	*	ZrC	
	

The *currently* allowed list is given in the table below:

============= ==============================================================
Property			Description
============= ==============================================================
cp	          Heat capacity (constant pressure), J/kg/K
------------- --------------------------------------------------------------
cv	      		Heat capacity (constant volume), J/kg/K
------------- --------------------------------------------------------------
g		  				Gamma=Cp/Cv, dimensionless
------------- --------------------------------------------------------------
h     				Enthalpy, J/kg
------------- --------------------------------------------------------------
my		    		Viscosity, kg/m/s
------------- --------------------------------------------------------------
pr  					Prandtl Number, dimensionless
------------- --------------------------------------------------------------
r  						Density, kg/m^3
------------- --------------------------------------------------------------
s  						Entropy, J/kg/K
------------- --------------------------------------------------------------
tc  					Thermal Conductivity, W/m/K
------------- --------------------------------------------------------------
v  						Sonic Velocity, m/s
------------- --------------------------------------------------------------
tempK  				Temperature, Kelvin
------------- --------------------------------------------------------------
P  						Pressure, Pascal
------------- --------------------------------------------------------------
mol						Mole fraction, dimensionless
------------- --------------------------------------------------------------
nu						Poisson ratio, dimensionless
------------- --------------------------------------------------------------
alpha					Coefficient of thermal expansion, m/m/K
------------- --------------------------------------------------------------
alphaT				Zero stress temperature, K
------------- --------------------------------------------------------------
E							Modulus of elasticity, Pa
============= ==============================================================


.. _matmanipulate:
	
Manipulating the Data Library
-----------------------------

Replace the original pressure and/or temperature dependnecies with new x- and y-attributes.


Start by importing the following items:

.. code::

	from ntpSystem.functions.propertytable import ChangeDataDependencies


The **syntax** is:

.. code::

	ChangeDataDependencies(matobject, xattr, yattr, ny, gridNy)


where:

- matobject PropertyData object
- xattr : (string). name of the x-attribute/property, e.g. "P" - pressure
- yattr : (string). name of the y-attribute/property, e.g. "h" - enthalpy. Does not have to be provided
- ny : (int). Max. number of points for the y-grid. If set to None, the original grid will be used. However, this may result in prohibitively large number of points.
- gridNy : (2-dim list) the first list defined the cutoffs of yattrs and the second list provides the number of grid points. The 2nd list must be n+1 in length where the first list should be n in length. e.g., [[1E+3, 1E+4, 1E+5, 1E+06], enthalpy values [50, 50, 50, 50, 200]] - the 1st number is the number of pts below 1E+03 J/kg and the last above 1E+06 J/kg


**Example**:

.. code::

	H5_ORIG_FILE = "ThermoPhysicalProperties.h5"
	# Read the table with hydrogen properties
	table = PropertyTable(H5_ORIG_FILE)
	
	# Read the original table to obtain pressure- and temperature-dependent pty
	H2_PT = table.read("H2")    
	
	# Change dependnecies
	H2_PH = ChangeDataDependencies(H2_PT, xattr="P", yattr="h", 
	                               gridNy=[[1E+3, 1E+4, 1E+5, 1E+06],
	                                       [20, 20, 20, 20, 20]])


.. _matfullexample:

Full Example
------------


.. code::

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
	
	ORIG_FILE = "../ignoretests/ThermoPhysicalProperties.h5"
	H5_ORIG_FILE = "ThermoPhysicalProperties.h5"
	# Read the table with hydrogen properties
	table = PropertyTable(H5_ORIG_FILE)
	
	# Read the original table to obtain pressure- and temperature-dependent pty
	H2_PT = table.read("H2")    
	
	# Change dependnecies
	H2_PH = ChangeDataDependencies(H2_PT, xattr="P", yattr="h", 
	                               gridNy=[[1E+3, 1E+4, 1E+5, 1E+06],
	                                       [20, 20, 20, 20, 20]])
	
	
	
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