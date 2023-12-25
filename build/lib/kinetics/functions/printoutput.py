# -*- coding: utf-8 -*-
"""printoutput.py

Function that prints a txt based output.

All the output and input parameters are printed by default. 
The index keys for all the parameters are included at the end of the output
file.

"""

import os

from datetime import datetime


from ntpSystem.errors.checkerrors import _isstr, _inlist
from ntpSystem.functions.steadystatesolver import COMPONENTS_DICT
from ntpSystem.containers.steadyinputs import COMPONENTS_DICT as INPS_DICT



def PrintTxtOutput(results, inps=None, dirpath="./", outputName="output.txt",
                   pmode="w"):
    """Print txt output file

    Parameters
    ----------
    results : SteadyState Results
        an object that contains all the attributes associated to all the
        results
    inps : SteadyStateInputs
        an object that contains all the attributes for the input cards
    dirpath : str
        full directory path where the file will be saved
    outputName : str
        name of the outputfile
    pmode : str
        print mode, e.g. ``w`` stands for write and ``a`` for append.


    Raises
    ------
    OSError
        If the path ``dirpath`` is entered, but does not exist.
    TypeError
        If ``dirpath``, ``outputName``, ``frmt``, ``pmode`` not of correct type

    """

    # error checking
    _checkerrors(dirpath, outputName, pmode)
    
    # define the full path to the file
    fullpath = dirpath + outputName

    # print the output file
    with open(fullpath, "w") as txtFile:
        if pmode == "w":
            # print logo
            _printLogo(txtFile)
            # print date and time
            _printDateTime(txtFile)  # print only as a header

        _printLabel(txtFile, "OUTPUT RESULTS")
        for cmpntName, cmpntDict in COMPONENTS_DICT.items():
            # print index key for all the outputs
            res = getattr(results, cmpntName)
            _printOutputs(txtFile, res, cmpntName, cmpntDict)
            
        if inps is not None:
            _printLabel(txtFile, "PROVIDED INPUTS")
            for inpName, inpDict in INPS_DICT.items():
                _printInputs(txtFile, inps, inpName, inpDict)
        

def _printDateTime(txtFile):
    """Print date and time of the current execution"""
    # strPrint = "----- INDEX KEY OUTPUTS -----"
    timeNow = datetime.now()
    dt_string = timeNow.strftime("%d/%m/%Y %H:%M:%S")
    txtFile.write("\n\n")
    nperc = 150
    txtFile.write("*"*nperc+"\n")
    txtFile.write("\t"*20+"Georgia Institute of Technology"+"\n")
    txtFile.write("\t"*20+"CoRE Group"+"\n")
    txtFile.write("\t"*20+"Dan Kotlyar"+"\n")
    txtFile.write("\t"*20+"Execution date and time "+dt_string+"\n")
    txtFile.write("*"*nperc+"\n\n\n")


def _printLabel(txtFile, strLabel):
    """print title label"""
    txtFile.write("\n\n")
    nperc = 150
    txtFile.write("%"*nperc+"\n")
    txtFile.write("\t"*20+strLabel+"\n")
    txtFile.write("%"*nperc+"\n\n\n")


def _printOutputs(txtFile, results, strPrint, dictParam):
    """print output values for a certain component"""

    txtFile.write("****************** "+strPrint+"\n")

    lenkey = 10  # number of characters for the longest key
    lendesc = 10  # number of characters for the longest description
    lenvals = 40 
    # find the number of required spacing in the first and second columns
    for key, desc in dictParam.items():
        if len(key) > lenkey:
            lenkey = len(key)  # first column
        if len(desc) > lendesc:
            lendesc = len(desc)  # second column

    # formatting
    frmt1 = "{:"+str(lenkey)+"s}"  # format for column 1
    frmt2 = "{:"+str(lendesc)+"s}"  # format for column 2
    frmt3 = "{:"+str(lenvals)+"s}"  # format for column 3
    
    # print title
    txtFile.write(frmt1.format('-'*lenkey)+"|"+frmt2.format('-'*lendesc)+"|"+
                  frmt3.format('-'*lenvals)+"\n")
    txtFile.write(
        frmt1.format("Parameter")+"|"+frmt2.format("Description")+"|"+
        frmt3.format("Value")+"\n")
    txtFile.write(frmt1.format('-'*lenkey)+"|"+frmt2.format('-'*lendesc)+"|"+
                  frmt3.format('-'*lenvals)+"\n")
    for key, desc in dictParam.items():
        value = str(results[key])
        txtFile.write(frmt1.format(key)+"|"+frmt2.format(desc)+
                      "|"+frmt3.format(value)+"\n")
    # underline under the entire table
    txtFile.write(frmt1.format('-'*lenkey)+"|"+frmt2.format('-'*lendesc)+
                  frmt3.format('-'*lenvals)+"\n")
    txtFile.write("\n\n")


def _printInputs(txtFile, inputs, strPrint, dictParam):
    """print input values for a certain component"""

    txtFile.write("****************** "+strPrint+"\n")

    lenkey = 10  # number of characters for the longest key
    lendesc = 10  # number of characters for the longest description
    lenvals = 40 
    # find the number of required spacing in the first and second columns
    for key, descL in dictParam.items():
        desc = descL[0]
        if len(key) > lenkey:
            lenkey = len(key)  # first column
        if len(desc) > lendesc:
            lendesc = len(desc)  # second column

    # formatting
    frmt1 = "{:"+str(lenkey)+"s}"  # format for column 1
    frmt2 = "{:"+str(lendesc)+"s}"  # format for column 2
    frmt3 = "{:"+str(lenvals)+"s}"  # format for column 3
    
    # print title
    txtFile.write(frmt1.format('-'*lenkey)+"|"+frmt2.format('-'*lendesc)+"|"+
                  frmt3.format('-'*lenvals)+"\n")
    txtFile.write(
        frmt1.format("Parameter")+"|"+frmt2.format("Description")+"|"+
        frmt3.format("Value")+"\n")
    txtFile.write(frmt1.format('-'*lenkey)+"|"+frmt2.format('-'*lendesc)+"|"+
                  frmt3.format('-'*lenvals)+"\n")
    for key, descL in dictParam.items():
        desc = descL[0]
        value = getattr(inputs, key)
        value = str(value)
        txtFile.write(frmt1.format(key)+"|"+frmt2.format(desc)+
                      "|"+frmt3.format(value)+"\n")
    # underline under the entire table
    txtFile.write(frmt1.format('-'*lenkey)+"|"+frmt2.format('-'*lendesc)+
                  frmt3.format('-'*lenvals)+"\n")
    txtFile.write("\n\n")



def _checkerrors(dirpath, outputName, pmode):
    """check that all variables are of correct type"""
    _isstr(dirpath, "Directory path")
    _isstr(outputName, "File output name")
    _isstr(pmode, "Print mode")
    _inlist(pmode, "Print mode", ["a", "w"])


    # check if directory exist
    if not os.path.isdir(dirpath):
        raise OSError("The directory {} does not exist.".format(dirpath))


def _printLogo(txtFile):
    " print logo"
    txtFile.write("/////////////////////////////////////////////////////////\n") 
    txtFile.write("     NUCLEAR THERMAL PROPULSION SYSTEM SOLVER            \n") 
    txtFile.write("                 EXPANDER CYCLE                          \n") 
    txtFile.write(" STEADY-STATE and TRANSIENT COMPUTATIONAL FRAMEWORK      \n") 
    txtFile.write("                                                         \n") 
    txtFile.write("/////////////////////////////////////////////////////////\n") 
    txtFile.write("                                                         \n")
    txtFile.write("                                                         \n")
    txtFile.write("               ----------------                          \n")
    txtFile.write("               \              /                          \n")
    txtFile.write("                \   H2 Tank  /                           \n")
    txtFile.write("                 ------------     pump         turbine   \n")
    txtFile.write("                      |                 shaft   _        \n")
    txtFile.write("                      v           /|___________/ |       \n")
    txtFile.write("                       --------->| |___________| |       \n")
    txtFile.write("                                  \|           \_|       \n")
    txtFile.write("                                   |       -----^------- \n")
    txtFile.write("                                   |      |             |\n")
    txtFile.write("                      o------------v      |             |\n")
    txtFile.write("                      |                   |             |\n")
    txtFile.write("                      v                   |-----|><|---- \n")
    txtFile.write("                      |                   |    bypass   |\n")
    txtFile.write("                      |                   |             |\n")
    txtFile.write("-------------<--------                    |             |\n")
    txtFile.write("|                     |    Mix Tee        |             |\n")
    txtFile.write("|     REFLECTOR------~~---->o-------------              |\n")
    txtFile.write("|         ^           |     ^                           |\n")
    txtFile.write("|         |           |     |          ----------------- \n")
    txtFile.write("|         |           |     |          |                 \n")
    txtFile.write("-----> NOZZLE         |     |          |                 \n")
    txtFile.write("                      |     ^          v                 \n")
    txtFile.write("                      |     |          |                 \n")
    txtFile.write("                      |   Return       |                 \n")
    txtFile.write("                      v     |        Fuel                \n")
    txtFile.write("                    Supply  |          |                 \n")
    txtFile.write("                      |     |          |                 \n")
    txtFile.write("                       -----           v                 \n")
    txtFile.write("                                    Chamber              \n")
    txtFile.write("                                       |                 \n")
    txtFile.write("                                       v                 \n")
    txtFile.write("                                  Nozzle throat          \n")
    txtFile.write("                                       |                 \n")
    txtFile.write("                                       v                 \n")
    txtFile.write("                                                         \n")


