.. _propertytable:

Interface to the Property Table and Data
========================================

Description
-----------

Interface into pre-generated **mono**- or **multi-dimensional** property
tables.

Provides a very simple look-up table like interface into standard
thermophysical tables. Assumes a specific structure of the data stored
in a hierarchical data format (HDF) file. Is designated to be used for
the purpose of NTP thermal analyses.

PropertyTable
-------------

The following section provides the **inputs, returns, and the methods**
of the ``PropertyTable`` class

Inputs:
^^^^^^^

-  ``h5path`` : Full path for the hdf5 file

Methods:
^^^^^^^^

-  ``read`` : read the property data set of a specific material

-  ``getref`` : obtain the reference for the material database

-  ``materials``: obtain a list of all the materials in the database


Execute
^^^^^^^

.. code:: 

    from ntpThermo.functions.propertytable import PropertyTable

.. code:: 

    # Read the database and store in table
    table = PropertyTable("../database/ThermoPhysicalProperties.h5")

The ``materials`` method
------------------------

**Obtains all the materials included in the datafile**

This is a useful class if the user is not aware of the materials
included in the datafile.

Returns:
^^^^^^^^

-  ``list``: all the materials in the database

Execute
^^^^^^^

.. code:: 

    table.materials()


.. parsed-literal::

    ['H2', 'Molybdenum', 'Tungsten', 'UC', 'UN', 'UO2', 'Zircaloy', 'ZrC']



The ``getref`` method
---------------------

**Return the reference for a specific material**

Inputs:
^^^^^^^

-  ``mat`` : string, name of the material in the databse, e.g. “H2”


Example:
^^^^^^^^

.. code:: python

           >>> table.getref("H2")

Execute
^^^^^^^

.. code:: 

    table.getref("H2")


.. parsed-literal::

    'NASA, 1993. Computer Program for Thermal and Transport Properties of Para Hydrogen from 20 to 10000 K, NASA Lewis Research Center. NASA-TP-3378. NASA, Space Technology Mission Directorate Game Changing Development Program, FY19 Annual Review https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20190031810.pdf'



The ``read`` method
-------------------

**Obtain a complete data set for a specific material.**

Inputs:
^^^^^^^

-  ``mat`` : string, name of the material in the databse, e.g. “H2”


Examples:
^^^^^^^^^

.. code:: python

           >>> table = PropertyTable("ThermoPhysicalProperties.h5)
           >>> H2 = table.read("H2")
           >>> UO2 = table.read("UO2")

Execute:
^^^^^^^^

.. code:: 

    H2 = table.read("H2")

PropertyData
------------

The following section provides the **inputs, returns, and the methods**
of the ``PropertyData`` class

Inputs:
^^^^^^^

-  ``matset`` : hdf5 data field

Methods:
^^^^^^^^

-  ``evaluate`` : obtain a specific property value

-  ``plot`` : plot a property value and the data points from the table

-  ``whatis`` : description of the parameter with its units

-  ``properties``: obtain a list of all the available properties



Execute
^^^^^^^

.. code:: 

    from ntpThermo.functions.propertytable import PropertyData
    import matplotlib.pyplot as plt

.. code:: 

    H2 = table.read("H2")

The ``propertynames`` method
----------------------------

**Print all the available properties for a specific material**


Execute
^^^^^^^

.. code:: 

    H2.propertynames()


.. parsed-literal::

    ['cp', 'cv', 'g', 'h', 'my', 'pr', 'r', 's', 'tc', 'v']



The ``whatis`` method
---------------------

**Return information on a specific property**

Inputs:
^^^^^^^

``pty`` : string, name of the thermal property, e.g. “tc” thermal
conductivity


Execute
^^^^^^^

.. code:: 

    H2.whatis("tc")


.. parsed-literal::

    Units(name='Thermal Conductivity', units='W/m/K')


.. code:: 

    h = H2.whatis("tc")
    print(h.name)
    print(h.units)


.. parsed-literal::

    Thermal Conductivity
    W/m/K
    

The ``evaluate`` method
-----------------------

**Evaluate a specific property for given temperature and/or pressure.**

Pressure and/or temperatures can be provided as arguments, or by name.
If just the temperature is used, either directly pass a ``None``
pressure, e.g. ``evaluate("tc", None, 600)`` or use named arguments with
``evaluate("tc", temperature=600)``. Similarly for just pressure, but
the option also exists to not pass anything as well, e.g.
``evalute("tc", 20)``

Inputs:
^^^^^^^

-  ``pty`` : string, name of the thermal property, e.g. “tc” thermal
   conductivity
-  ``pressure`` : float, optional, pressure in MPa
-  ``temperature`` : float, optional, temperature in Kelvin


Notes:
^^^^^^

-  2-D interpolation is allowed for temperature and pressure.
-  1-D interpolation is allowed only for temperature.


Execute:
^^^^^^^^

.. code:: 

    table = PropertyTable("../database/ThermoPhysicalProperties.h5")
    H2 = table.read("H2")

.. code:: 

    print(H2.evaluate("tc", pressure=6, temperature=1754))
    print(H2.evaluate("cp", pressure=6, temperature=1754))


.. parsed-literal::

    0.7120571
    16535.7
    

The ``plot`` method
-------------------

**Comparative plot for the value of a specific property.**

Inputs:
^^^^^^^

-  ``pty`` : string, name of the thermal property, e.g. “tc” thermal
   conductivity
-  ``pressure`` : float, optional, pressure in MPa
-  ``temperature`` : float, optional, temperature in Kelvin
-  ``desc`` : str, optional, description of the property. Default value
   is read through the whatis method.
-  ``units`` : str, optional, description of the property’s units.
   Default value is read through the whatis method.

Execute:
^^^^^^^^

.. code:: 

    H2.plot("tc", pressure=18, temperature=557)



.. image:: propertytable_files//propertytable_108_0.png


The ``mixbinary`` method
------------------------

**Evaluate a weighted value for a binary mixture.**

The method evaluates the weighted property given two values that
represent the same property for two materials given a certain weight
fraction at a specific temperature and/or pressure .

Methodology: weighted properties for binary systems
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

| There are several options to weight the properties of binary mixtures.
| **1. Serial:**

:math:`\begin{align}
\begin{split}
z = w\times z_1 + (1-w)\times z_2 
\end{split}
\end{align}`

**2. Parallel:**

:math:`\begin{align}
\begin{split}
z = \left( \frac{w}{z_1} + \frac{1-w}{z_2}  \right)^{-1} 
\end{split}
\end{align}`

where :math:`w` is the weight fraction for material of type 1 having a
property (e.g. conductivity) value of :math:`z_1`. The property value
for material 2 is :math:`z_2`.

There is a third option that is used to weight the properties when
particles (e.g. UO:math:`_2`) are embedded within a matrix material,
such as tungsten.

**3. Bruggeman–Fricke methodology**

:math:`\begin{align}
\begin{split}
z = z_1 + (1-w)(z_2-z_1)\left(\frac{z_1}{z_2} \right)^{1/3}
\end{split}
\end{align}`

where :math:`w` is the volume fraction for the particles’ material
having a property value of :math:`z_1`. The property value for the
matrix material is :math:`z_2`.

Inputs:
^^^^^^^

-  ``val1`` : float, value for material of type 1, e.g. thermal
   conductivity for Tungsten
-  ``val2`` : float, value for material of type 2, e.g. thermal
   conductivity for Moly
-  ``w1`` : float, weight fraction of material-1
-  ``method`` : string, weighting method, e.g. “parallel”


Notes:
^^^^^^

The Bruggeman-Fricke weighting method assumes that ``val1`` represent
the particles, ``val2`` represent the matrix, and ``w1`` is the volume
fraction of the particles in the matrix.


Execute:
^^^^^^^^

.. code:: 

    from functions.propertytable import mixbinary

.. code:: 

    mixbinary(val1= 0.5, val2= 4.6, w1= 0.164, method= "Parallel")


.. parsed-literal::

    1.9617877857386554



Summary example
---------------

Obtain the **weighted conductivity** of the UO2 matrix (40% by volume)
embedded within a tungsten materix

**Read the data**

.. code:: 

    data = PropertyTable("../database/ThermoPhysicalProperties.h5")

**Read the thermal conductivity for UO2 and Tungsten**

.. code:: 

    UO2 = data.read("UO2")
    Tungsten = data.read("Tungsten")

**Obtain the conductivity values for UO2 and Tungsten**

.. code:: 

    valUO2 = UO2.evaluate("tc", temperature= 1585)
    valW = Tungsten.evaluate("tc", temperature= 1585)

**Obtain the weighted value**

.. code:: 

    val = mixbinary(val1= valUO2, val2= valW, w1= 0.4, method= "Bruggeman-Fricke")

.. code:: 

    print("Conductivity for UO2 is {} W/mK".format(round(valUO2, 2)))
    print("Conductivity for Tungsten is {} W/mK".format(round(valW, 2)))
    print("Weighted conductivity is {} W/mK".format(round(val, 2)))


.. parsed-literal::

    Conductivity for UO2 is 1.78 W/mK
    Conductivity for Tungsten is 109.14 W/mK
    Weighted conductivity is 18.1 W/mK
    

**Plot results**

.. code:: 

    UO2.plot("tc", temperature=1585, desc ='Thermal conductivity', units='W/m/K')



.. image:: propertytable_files//propertytable_138_0.png


.. code:: 

    Tungsten.plot("tc", temperature=1585, desc ='Thermal conductivity', units='W/m/K')



.. image:: propertytable_files//propertytable_139_0.png


