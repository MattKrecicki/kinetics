# -*- coding: utf-8 -*-
"""
Created on Fri Feb  2 10:34:06 2024

@author: matt krecicki


"""



class feedback:
    
    def __checkinputs(self):
        """function runs basic error checking on user defined inputs"""
        pass
    
    def __init__(self, **kwargs):
        """function initalizes feedback class"""
        
        self.__dict__.update(kwargs)
        self.__checkinputs()
        
        