# -*- coding: utf-8 -*-
"""
Created on Fri Feb  2 17:17:59 2024

@author: matt krecicki
@email: matthewkrecicki@gmail.com

References:
    1. Verification of iterative matrix solutions for multipoint kinetics equations
    2. Investigating decoupling effects in CFV type of reactors using Avery’s coupled reactors theory

Need to move the dependencies to the region object and out of the container object

container object should have the interpolation functions within it and the ablity to generate the matrix


"""


import numpy as np
from kinetics.containers.inputs import regionKineticsData
from kinetics.containers.inputs import multiPointKineticsInputsContainer as \
    inputsContainer


# -----------------------------------------------------------------------------
# ----- Setup Multipoint Kinetics Data Container for Transient Analysis -------
# -----------------------------------------------------------------------------

mpkdata = inputsContainer(nregions=3, dependencies=["hrod"], typ="avery")


# ----- initialize region 1 data ----------------------------------------------
r1 = \
    regionKineticsData(Id="1", typ="avery", x=0.0, y=-1.0, z=0.0, 
                       dependencies=["hrod"], volume=1.0)

#add region 1's 0cm rod height function
r1.add(#average energy released per fission
       Q=200.0,
       #one-group neutron velocity
       v=1.86E+04,        
       #delayed neutron group decay constant
       lamdaki = np.array([0.40529]),
       #delayed neutron fraction groups
       Bki = np.array([0.00376]),
       #detail how couling is generated
       coupling = ["1", "2", "3"],
       # coupling to:  1      2       3
       Kjk = np.array([0.92804, 0.03756, 0.03754]),
       # delayed coupling, for this case assumed to be the same as the prompt population   
       Kjkd = np.array([0.92804, 0.03756, 0.03754]),
       # coupling to:  1      2       3
       Ljk =  np.array([0.3613, 0.6154, 0.6126])*1e-6,
       Ljkd = np.array([0.3613, 0.6154, 0.6126])*1e-6,
       #specify dependency variables
       hrod=0.0)

#add region 1's 35cm rod height function
r1.add(Q=200.0,
       #one-group neutron velocity
       v=1.86E+04,
       #volume of fissionable region
       #delayed neutron group decay constant
       lamdaki = np.array([0.40529]),
       #delayed neutron fraction groups
       Bki = np.array([0.00376]),
       #detail how couling is generated
       coupling = ["1", "2", "3"],
       # coupling to:  1      2       3
       Kjk = np.array([0.92804, 0.03756, 0.03754]),
       # delayed coupling, for this case assumed to be the same as the prompt population   
       Kjkd = np.array([0.92804, 0.03756, 0.03754]),
       # coupling to:  1      2       3
       Ljk =  np.array([0.3613, 0.6154, 0.6126])*1e-6,
       Ljkd = np.array([0.3613, 0.6154, 0.6126])*1e-6,
       #specify dependency variables
       hrod=35.0)

#add 1st region to data container
mpkdata.add(r1)

# ----- initialize region 2 data ----------------------------------------------
r2 = \
    regionKineticsData(Id="2", typ="avery", x=1.0, y=0.0, z=0.0, 
                       dependencies=["hrod"], volume=1.0)


#add region 2's 0cm rod height function
r2.add(#average energy released per fission
       Q=200.0,
       #one-group neutron velocity
       v=1.86E+04,
       #volume of fissionable region
       volume=1.0,
       #delayed neutron group decay constant
       lamdaki = np.array([0.40529]),
       #delayed neutron fraction groups
       Bki = np.array([0.00376]),
       #detail how couling is generated
       coupling = ["1", "2", "3"],
       # coupling to:  1      2       3
       Kjk = np.array([0.03754, 0.92804, 0.03756]),
       # delayed coupling, for this case assumed to be the same as the prompt population   
       Kjkd = np.array([0.03754, 0.92804, 0.03756]),
       # coupling to:  1      2       3
       Ljk =  np.array([0.6126, 0.3613, 0.6154])*1e-6,
       Ljkd = np.array([0.6126, 0.3613, 0.6154])*1e-6,
       #specify dependency variables
       hrod=0.0)


#add region 2's 35cm rod height function
r2.add(#average energy released per fission
       Q=200.0,
       #one-group neutron velocity
       v=1.86E+04,
       #volume of fissionable region
       volume=1.0,
       #delayed neutron group decay constant
       lamdaki = np.array([0.40529]),
       #delayed neutron fraction groups
       Bki = np.array([0.00376]),
       #detail how couling is generated
       coupling = ["1", "2", "3"],
       # coupling to:  1      2       3
       Kjk = np.array([0.03754, 0.92804, 0.03756]),
       # delayed coupling, for this case assumed to be the same as the prompt population   
       Kjkd = np.array([0.03754, 0.92804, 0.03756]),
       # coupling to:  1      2       3
       Ljk =  np.array([0.6126, 0.3613, 0.6154])*1e-6,
       Ljkd = np.array([0.6126, 0.3613, 0.6154])*1e-6,
       #specify dependency variables
       hrod=35.0)

#add 2nd region to data container
mpkdata.add(r2)

# ----- initialize region 3 data ----------------------------------------------

r3 = \
    regionKineticsData(Id="3", typ="avery", x=-1.0, y=0.0, z=0.0, 
                       dependencies=["hrod"], volume=1.0)
    

#add region 3's 0cm rod height function
r3.add(#average energy released per fission
       Q=200.0,
       #one-group neutron velocity
       v=1.86E+04,
       #delayed neutron group decay constant
       lamdaki = np.array([0.40529]),
       #delayed neutron fraction groups
       Bki = np.array([0.003657]),
       #detail how couling is generated
       coupling = ["1", "2", "3"],
       # coupling to:  1      2       3
       Kjk = np.array([0.03756, 0.03754, 0.92804]),
       # delayed coupling, for this case assumed to be the same as the prompt population   
       Kjkd = np.array([0.03756, 0.03754, 0.92804]),
       # coupling to:  1      2       3
       Ljk =  np.array([0.6154, 0.6126, 0.3613])*1e-6,
       Ljkd = np.array([0.6154, 0.6126, 0.3613])*1e-6,
       #specify dependency variables
       hrod=0.0)


#add region 3's 35cm rod height function
r3.add(#average energy released per fission
       Q=200.0,
       #one-group neutron velocity
       v=1.86E+04,
       #volume of fissionable region
       volume=1.0,
       #delayed neutron group decay constant
       lamdaki = np.array([0.40529]),
       #delayed neutron fraction groups
       Bki = np.array([0.003657]),
       #detail how couling is generated
       coupling = ["1", "2", "3"],
       # coupling to:  1      2       3
       Kjk = np.array([0.03756, 0.03754, 0.92804]),
       # delayed coupling, for this case assumed to be the same as the prompt population   
       Kjkd = np.array([0.03756, 0.03754, 0.92804]),
       # coupling to:  1      2       3
       Ljk =  np.array([0.6154, 0.6126, 0.3613])*1e-6,
       Ljkd = np.array([0.6154, 0.6126, 0.3613])*1e-6,
       #specify dependency variables
       hrod=35.0)

#add 3rd region to data container
mpkdata.add(r3)

# -----------------------------------------------------------------------------
# ----- Assemble final container ----------------------------------------------
# -----------------------------------------------------------------------------

mpkdata.validate()

# use export function to save inputs for future use

