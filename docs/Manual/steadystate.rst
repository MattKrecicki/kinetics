.. _ss_simulation:

=======================
Steady state simulation
=======================


In order to define the cards needed to run a steady state problem, the user is required to
define all the required input cards and only then execute the steady state problem.
The inputs are defined using the ``SteadyStateInputs`` container.


.. code::

	from ntpSystem.functions.steadystatesolver import SteadyStateSystem
	ssOutps = SteadyStateSystem(ssInps)
	ssOutps.Solve()  # solve steady state


An integral part of this system code is to define the input files using the ``ntpThermo`` code.


.. _ssinpscards:

Inputs
------

Import the container: 

.. code::

	from ntpSystem.containers.steadyinputs import SteadyStateInputs

Initialize the container: 

.. code::

	ssInps = SteadyStateInputs(prntMsgs)
	
``prntMsgs`` is a flag to indicate if messages should be printed during execution of the code. The default is True.


Feed the container with data using the ``AddData`` method:

.. code::

	ssInps.AddData(component, **kwargs)
	
- ``component`` is the name of the component for which the data should be provided for. The different types of components are provided in the following table.
- ``kwargs`` are variable names and their values. The variables names are pre-defined for each component.


Once all the data is provided **verify that the data is properly fed** using the ``ValidateInputs`` method:

.. code::

	ssInps.ValidateInputs()
	

.. _component_cards:

Components
----------


The different *components* that must be populated with data:

============================= ==============================================================
Component											Description
============================= ==============================================================
:ref:`component_react`		    Basic definitions for core, nozzle, reflector, and piping components
----------------------------- --------------------------------------------------------------
:ref:`component_operat`	      operational parameters
----------------------------- --------------------------------------------------------------
:ref:`component_turbopump`	  descriptive data of the turbopump components
============================= ==============================================================

:ref:`A full example of how to pupulate input cards is included. <inputcardsfullexample>`

Each component has assigned keys some which are madatory to be populated and some optional.
The keys for all the components are **case sensitive**!


.. _component_react:

reactor
^^^^^^^

=================== ============================================================== =============  ================
Key                 Description                                                     Required       Range
=================== ============================================================== =============  ================
fuel_inpfile        input file for ntpThermo									                       yes            string
------------------- -------------------------------------------------------------- -------------  ----------------
fuel_chId           string name for the channel						                           yes            string
------------------- -------------------------------------------------------------- -------------  ----------------
fuel_n              total number of coolant channels                                 yes              >=1
------------------- -------------------------------------------------------------- -------------  ----------------
fuel_layer0         index of the very first layer in the channel                     yes              >=0
------------------- -------------------------------------------------------------- -------------  ----------------
fuel_nlayers        number of layers belonging to this element type                  yes              >=1
------------------- -------------------------------------------------------------- -------------  ----------------
mesupply_inpfile    input file for ntpThermo									                       yes            string
------------------- -------------------------------------------------------------- -------------  ----------------
mesupply_chId       string name for the channel						                           yes            string
------------------- -------------------------------------------------------------- -------------  ----------------
mesupply_n          total number of coolant channels                                 yes              >=1
------------------- -------------------------------------------------------------- -------------  ----------------
mesupply_layer0     index of the very first layer in the channel                     yes              >=0
------------------- -------------------------------------------------------------- -------------  ----------------
mesupply_nlayers     number of layers belonging to this element type                  yes              >=1
------------------- -------------------------------------------------------------- -------------  ----------------
mereturn_inpfile    input file for ntpThermo									                       yes            string
------------------- -------------------------------------------------------------- -------------  ----------------
mereturn_chId       string name for the channel						                           yes            string
------------------- -------------------------------------------------------------- -------------  ----------------
mereturn_n          total number of coolant channels                                 yes              >=1
------------------- -------------------------------------------------------------- -------------  ----------------
mereturn_layer0     index of the very first layer in the channel                     yes              >=0
------------------- -------------------------------------------------------------- -------------  ----------------
mereturn_nlayers    number of layers belonging to this element type                  yes              >=1
------------------- -------------------------------------------------------------- -------------  ----------------
nozzle_inpfile    	input file for ntpThermo									                       yes            string
------------------- -------------------------------------------------------------- -------------  ----------------
nozzle_chId         string name for the channel						                           yes            string
------------------- -------------------------------------------------------------- -------------  ----------------
nozzle_n            total number of coolant channels                                 yes              >=1
------------------- -------------------------------------------------------------- -------------  ----------------
nozzle_layer0       index of the very first layer in the channel                     yes              >=0
------------------- -------------------------------------------------------------- -------------  ----------------
nozzle_nlayers      number of layers belonging to this element type                  yes              >=1
------------------- -------------------------------------------------------------- -------------  ----------------
reflector_inpfile   input file for ntpThermo									                       yes            string
------------------- -------------------------------------------------------------- -------------  ----------------
reflector_chId      string name for the channel						                           yes            string
------------------- -------------------------------------------------------------- -------------  ----------------
reflector_n         total number of coolant channels                                 yes              >=1
------------------- -------------------------------------------------------------- -------------  ----------------
reflector_layer0    index of the very first layer in the channel                     yes              >=0
------------------- -------------------------------------------------------------- -------------  ----------------
reflector_nlayers   number of layers belonging to this element type                  yes              >=1
------------------- -------------------------------------------------------------- -------------  ----------------
pipeFuel_inpfile   	input file for ntpThermo									                       yes            string
------------------- -------------------------------------------------------------- -------------  ----------------
pipeFuel_chId      	string name for the channel						                           yes            string
------------------- -------------------------------------------------------------- -------------  ----------------
pipeFuel_n         	total number of coolant channels                                 yes              >=1
------------------- -------------------------------------------------------------- -------------  ----------------
pipeFuel_layer0    	index of the very first layer in the channel                     yes              >=0
------------------- -------------------------------------------------------------- -------------  ----------------
pipeFuel_nlayers    number of layers belonging to this element type                  yes              >=1
------------------- -------------------------------------------------------------- -------------  ----------------
pipeMod_inpfile   	input file for ntpThermo									                       yes            string
------------------- -------------------------------------------------------------- -------------  ----------------
pipeMod_chId      	string name for the channel						                           yes            string
------------------- -------------------------------------------------------------- -------------  ----------------
pipeMod_n         	total number of coolant channels                                 yes              >=1
------------------- -------------------------------------------------------------- -------------  ----------------
pipeMod_layer0    	index of the very first layer in the channel                     yes              >=0
------------------- -------------------------------------------------------------- -------------  ----------------
pipeMod_nlayers    	number of layers belonging to this element type                  yes              >=1
=================== ============================================================== =============  ================


**Example**

.. code::

	ssInps.AddData("reactor",
	    fuel_inpfile = '.\\inputs_largeNerva\\FuelCircuit.txt', 
	    fuel_chId = 'fe', fuel_n = 564, fuel_layer0 = 0, fuel_nlayers = 20,
	    mesupply_inpfile = '.\\inputs_largeNerva\\ReflMEcircuit.txt',  
	    mesupply_chId = 'mod', mesupply_n = 241, mesupply_layer0 = 0, mesupply_nlayers = 10,
	    mereturn_inpfile = '.\\inputs_largeNerva\\ReflMEcircuit.txt',  
	    mereturn_chId = 'mod', mereturn_n = 241, mereturn_layer0 = 10, mereturn_nlayers = 10,
	    nozzle_inpfile = '.\\inputs_largeNerva\\ReflMEcircuit.txt',
	    nozzle_chId = 'ref', nozzle_n = 400, nozzle_layer0 = 0, nozzle_nlayers = 10,
	    reflector_inpfile = '.\\inputs_largeNerva\\ReflMEcircuit.txt',
	    reflector_chId = 'ref', reflector_n = 400, reflector_layer0 = 10, reflector_nlayers = 10,
	    pipeFuel_inpfile = '.\\inputs_largeNerva\\PipeFuel.txt',
	    pipeFuel_chId = 'pipeFuel', pipeFuel_n = 1, pipeFuel_layer0 = 0, pipeFuel_nlayers = 5,
	    pipeMod_inpfile = '.\\inputs_largeNerva\\PipePump.txt',
	    pipeMod_chId = 'pipePump', pipeMod_n = 1, pipeMod_layer0 = 0, pipeMod_nlayers = 5)


Notes:
- All the components must be defined. At the moment, the components must include the fuel, nozzle, reflector, supply, return, fuel-pipe, and moderator-pipe.
- Each component can be defined using a single and unique file. However, a single file can also be used to define multiple channels. In addition, the solution of each file can be obtained using the various convergence schemes that exist in ``ntpThermo``.
- 







.. _component_operat:

operation
^^^^^^^^^

=================== ============================================================== =============  ================
Key                 Description                                                     Required       Range
=================== ============================================================== =============  ================
Pt                   hydrogen storage tank pressure [Pa]                             yes            >=0
------------------- -------------------------------------------------------------- -------------  ----------------
Tt                   hydorgen storage tank temperature [K]                           yes            >=0
------------------- -------------------------------------------------------------- -------------  ----------------
Pc                   chamber pressure [Pa]                                           yes            >=0
------------------- -------------------------------------------------------------- -------------  ----------------
Tc                   chamber temperature [K]                                         yes            >=0
------------------- -------------------------------------------------------------- -------------  ----------------
power                total power [watts]                                             no             >=0
------------------- -------------------------------------------------------------- -------------  ----------------
mdot                 total mass flow rate [kg/s]                                     no             >=0
------------------- -------------------------------------------------------------- -------------  ----------------
At                   nozzle throat area [m^2]                                        yes            >=0
------------------- -------------------------------------------------------------- -------------  ----------------
nozzleExpR           nozzle expansion ratio                                          yes            >=0
=================== ============================================================== =============  ================

**Example**

.. code::

	ssInps.AddData("operation",
	    Pc = 6.89E+06,              # desired chamber pressure, Pa
	    Tc = 2794,                  # desired chamber temperature, K
	    Pt = 0.22e+6,               # tank pressure, Pa
	    Tt = 21.0,                  # tank temperature, K
	    power = 581640222.4,        # absolute total power in watts
	    mdot = None,                # system mass flow rate kg/s
	    At = 0.005662964739541421,  # nozzle throat area, m2
	    nozzleExpR = 300.,          # nozzle expansion ratio 
	    )  # operation


The following options to calculate power and mass flow rate exist:
- User cannot define ``mdot`` and ``power`` together.
- Either mass flow rate or power can be defined.
- If none of the mass flow rate or power is defined, isentropic relation is used to calculate the mass flow rate in the system.



.. _component_turbopump:

turbopump
^^^^^^^^^

=================== ============================================================== =============  ================
Key                 Description                                                     Required       Range
=================== ============================================================== =============  ================
noPumps             number of turbopumps in engine system                           yes            >=1
------------------- -------------------------------------------------------------- -------------  ----------------
pumpPratio          pump outlet-to-inlet pressure                                   yes            >=1.0
------------------- -------------------------------------------------------------- -------------  ----------------
pumpSs              target specific suction (20 is max.), unitless                  5.0           [0.1, 20]
------------------- -------------------------------------------------------------- -------------  ----------------
pumpMarginPin       margin inlet pump pressure needed for caviation                 0.1           [0, 1]
------------------- -------------------------------------------------------------- -------------  ----------------
pumpEfficiency      constant pump efficiency                                        no            [0, 1]
------------------- -------------------------------------------------------------- -------------  ----------------
turbType            turbine type                                                    axial         axial / radial
------------------- -------------------------------------------------------------- -------------  ----------------
turbMat             material name                                                   no            see below
------------------- -------------------------------------------------------------- -------------  ----------------
turbDensity         material density [kg/m^3]                                       no             >=0.0
------------------- -------------------------------------------------------------- -------------  ----------------
turbStress          material stress [Pa]                                            no             >=0.0
------------------- -------------------------------------------------------------- -------------  ----------------
turbD               turbine diameter [m]                                            no             >=0.0
------------------- -------------------------------------------------------------- -------------  ----------------
turbRefPoint        turbine inlet/outlet reference point for density                outlet        outlet / inlet / average
------------------- -------------------------------------------------------------- -------------  ----------------
turbEfficiency      constant turbine efficiency                                     no            >=0
------------------- -------------------------------------------------------------- -------------  ----------------
gearR               gear ratio                                                      1.0            >=1.0
------------------- -------------------------------------------------------------- -------------  ----------------
strsF               stress factor                                                   no            >=0
------------------- -------------------------------------------------------------- -------------  ----------------
turbh2D             turbine h/D ratio                                               no            >=0
------------------- -------------------------------------------------------------- -------------  ----------------
safetyF             safety factor                                                   1.2            >=0
------------------- -------------------------------------------------------------- -------------  ----------------
turbAxB02B          inner-to-outer radius of the disk ratio                         4.0            [3, 20]
------------------- -------------------------------------------------------------- -------------  ----------------
turbAxTaper         Blade taper (fig.10.8 Emrich's book.)                           0.5            [0.4, 1.0]
------------------- -------------------------------------------------------------- -------------  ----------------
turbRadr2R          Blade-to-disc (r/R) radii ratio                                 0.4            [0.0, 1.0]
------------------- -------------------------------------------------------------- -------------  ----------------
cMesh               coarse mesh division for finding optimum efficiency             30             [5, 10000]
------------------- -------------------------------------------------------------- -------------  ----------------
fMesh               fine mesh division for finding optimum efficiency               30             [5, 10000]
=================== ============================================================== =============  ================



Allowed **turbopump** materials:

.. code::

	['stainless steel', 'aluminum', 'brass', 'bronze', 'inconcel', 'titanium']

The user can also define their own material by defining: ``turbDensity`` and ``turbStress``.


The following cards from the table above are optional:

.. code::

	[
	'turbMat', 'turbDensity', 'turbStress', 'pumpEfficiency',
	'turbEfficiency', 'gearR', 'turbD', 'turbAxB02B', 'turbAxTaper',
	'turbRadr2R', 'turbh2D', 'strsF', 'pumpD']
	

There are two modes of calculation in the case of a turbine:
1. ``turbD`` is directly provided and the max shaft speed is calculated together with the gear ratio.
2. ``gearR``gear ratio is provided and ``turbD`` is calculated on-the-fly.



.. _inputcardsfullexample:


Full Example
------------


.. code::

	#------------------------------------------------------------------------------
	#----- Data for a dual NERVA pump
	#------------------------------------------------------------------------------
	noFuels = 564
	noMods = 241
	totPow = 581640222.4  #  555E+06
	# The real total power should be 581640222.4
	
	mdot = 12.68  # kg/s
	pumpPratio = 15.65 / 0.22  # pump pressure ratio
	
	#------------------------------------------------------------------------------
	#----- define inputs
	#------------------------------------------------------------------------------
	
	
	#initalize input container-----------------------------------------------------
	ssInps = SteadyStateInputs()
	
	
	# reactor geometry-------------------------------------------------------------
	 
	ssInps.AddData("reactor",
	    fuel_inpfile = '.\\inputs_largeNerva\\FuelCircuit.txt', 
	    fuel_chId = 'fe', fuel_n = 564, fuel_layer0 = 0, fuel_nlayers = 20,
	    mesupply_inpfile = '.\\inputs_largeNerva\\ReflMEcircuit.txt',  
	    mesupply_chId = 'mod', mesupply_n = 241, mesupply_layer0 = 0, mesupply_nlayers = 10,
	    mereturn_inpfile = '.\\inputs_largeNerva\\ReflMEcircuit.txt',  
	    mereturn_chId = 'mod', mereturn_n = 241, mereturn_layer0 = 10, mereturn_nlayers = 10,
	    nozzle_inpfile = '.\\inputs_largeNerva\\ReflMEcircuit.txt',
	    nozzle_chId = 'ref', nozzle_n = 400, nozzle_layer0 = 0, nozzle_nlayers = 10,
	    reflector_inpfile = '.\\inputs_largeNerva\\ReflMEcircuit.txt',
	    reflector_chId = 'ref', reflector_n = 400, reflector_layer0 = 10, reflector_nlayers = 10,
	    pipeFuel_inpfile = '.\\inputs_largeNerva\\PipeFuel.txt',
	    pipeFuel_chId = 'pipeFuel', pipeFuel_n = 1, pipeFuel_layer0 = 0, pipeFuel_nlayers = 5,
	    pipeMod_inpfile = '.\\inputs_largeNerva\\PipePump.txt',
	    pipeMod_chId = 'pipePump', pipeMod_n = 1, pipeMod_layer0 = 0, pipeMod_nlayers = 5)
	
	
	# operational conditions-------------------------------------------------------
	ssInps.AddData("operation",
	    Pc = 6.89E+06,              # desired chamber pressure, Pa
	    Tc = 2794,                  # desired chamber temperature, K
	    Pt = 0.22e+6,               # tank pressure, Pa
	    Tt = 21.0,                  # tank temperature, K
	    power = 581640222.4,        # absolute total power in watts
	    mdot = None,                # system mass flow rate kg/s
	    At = 0.005662964739541421,  # nozzle throat area, m2
	    nozzleExpR = 300.,          # nozzle expansion ratio 
	    )  # operation
	
	
	
	# turbopumppump data-----------------------------------------------------------
	ssInps.AddData("turbopump",
	    noPumps = 2,               # number of pumps in expander cycle
	    pumpSs = 20.0,              # target pump specific suction
	    pumpPratio = pumpPratio,         # guess for outlet-to-inlet pump pressure ratio
	    pumpMarginPin = 0.15,       # margin inlet pump pressure for caviation
	    pumpEfficiency = None,     # pump efficiency
	    pumpD = None,              # the user can specify the pump impeller diameter 
	    turbMat = None,            # turbine material
	    turbDensity = 1500,  # 2000.0,      # turbine density [kg/m^3]
	    turbStress = 300E+06,     # turbine density [Pa]
	    turbType = 'axial',        # turbine type ['axial', 'radial']
	    turbD = None,              # turbine rotor diameter [m]
	    turbRefPoint = 'outlet',   # ref. point [inlt/outlet] for density calc.
	    turbEfficiency = None,     # turbine efficiency
	    turbAxTaper = 0.5,         # Blade taper
	    turbAxB02B = 4.0,          # inner-to-outer radius of the disk ratio
	    # strsF = 0.2,               # stress factor for turbine
	    safetyF = 1.2,             # safety factor for turbine
	    gearR = 1.0,               # gear ratio turbine-to-pump speed
	    cMesh = 100,                # coarse mesh division for finding optimum eff.
	    fMesh = 100,                # fine mesh division for finding optimum eff.
	    )  # turbopump
	
	# Add power in the pump and turbine
	
	
	# validate that all attributes are defined-------------------------------------
	ssInps.ValidateInputs()