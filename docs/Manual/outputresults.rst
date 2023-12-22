.. _printresults:

==============
Print results
==============


The results and inputs can be printed to a txt file.
The following syntax should be used.

**Import**:

.. code::

	from ntpSystem.functions.printoutput import PrintTxtOutput
	ssOutps = SteadyStateSystem(ssInps)
	ssOutps.Solve()  # solve steady state

**Create the inputs container and execute the steady-state simulation**:

	
.. code::

	ssOutps = SteadyStateSystem(ssInps)
	ssOutps.Solve() 
	
	
**Print results and inputs to a text file**:

	
.. code::

	PrintTxtOutput(results, inps=None, dirpath="./", outputName="output.txt", pmode="w")
                   
where,

- ``results`` : CoreResults
    an object that contains all the attributes associated to all the
    results. Must be provided
- ``inps`` : SteadyStateInputss object
    an object that contains all the input attributes. Default is None - in which case no inputs are printed.
- ``dirpath`` : str
    full directory path where the file will be saved. Default is: "./" - current working directory.
- ``outputName`` : str
    name of the outputfile. Default is "output.txt".
- ``pmode`` : str
    print mode, e.g. ``w`` stands for write and ``a`` for append.


**Example-1**: with inputs.

.. code::

	PrintTxtOutput(ssOutps, ssInps, outputName="outputNerva")


**Example-2**: without inputs.

.. code::

	PrintTxtOutput(ssOutps, outputName="outputNerva")