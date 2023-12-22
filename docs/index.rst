NTP System Package (ntpSystem's) documentation!
===============================================

.. only:: html

The framework includes a steady-state solver used to                                                                   
converge on pumping requirements and conserve enthalpy prior to the start of a transient. 
The novel approach           
here relies on the steady-state component within the framework to generate operational maps. The latter indicate       
the feasibility of control based on nozzle chamber conditions for a given engine design in lieu of publicly unavailableeferred citation for attribution
pump curves. The framework is applied to investigate uprate in power maneuvers in conjunction with                     
sensitivity studies to show the feasibility of an engine startup   
                                                    
    
		V. Manickam, D. Kotlyar, 2022. 
		Implementation of a comprehensive reduced order methodology for transient analysis of nuclear thermal propulsion engines, 
		Nuclear Engineering and Design, 395, 111841.  https://doi.org/10.1016/j.nucengdes.2022.111841


.. toctree::
   :maxdepth: 1
   :caption: Contents:
   
   overview
   contributors
   Manual/index
   Examples/index
   Methodology/index
   develop/index
   install
   license
   glossary



.. code::

                      ----------------
                      \              /
                       \   H2 Tank  /
                        ------------     pump         turbine
                             |                 shaft   _
                             v           /|___________/ |
                              --------->| |___________| |
                                         \|           \_|
                                          |       -----^-------
                                          |      |             |      
                             o------------v      |             |
                             |                   |             |
                             v                   |-----|><|----
                             |                   |    bypass   |
                             |                   |             |
       -------------<--------                    |             |
       |                     |    Mix Tee        |             |
       |     REFLECTOR------~~---->o-------------              |
       |         ^           |     ^                           |
       |         |           |     |          -----------------       
       |         |           |     |          |
       -----> NOZZLE         |     |          |
                             |     ^          v
                             |     |          |
                             |   Return       |
                             v     |        Fuel
                           Supply  |          |
                             |     |          |
                              -----           v
                                           Chamber
                                              |
                                              v
                                         Nozzle throat
                                              |
                                              v
                                              
                                       
                                         
     


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
