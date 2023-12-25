# -*- coding: utf-8 -*-
"""customerrors

A set of custom error raises to identify when the error occurs.

"""

# -----------------------------------------------------------------------------
# Raise Errors Section
# -----------------------------------------------------------------------------


class ssInputsError(Exception):
    def __init__(self, message):
        str1 = "\n!!!An issue while assigning steady state input values:\n"
        str1 += "-------------------------------------------------------- "
        super().__init__(str1+"\n"+message)

class ssOutputsError(Exception):
    def __init__(self, message):
        str1 = "\n!!!An issue while assigning steady state output values:\n"
        str1 += "-------------------------------------------------------- "
        super().__init__(str1+"\n"+message)