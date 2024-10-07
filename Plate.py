# -*- coding: utf-8 -*-
"""
Created on Tue Sep  4 17:35:26 2012

@author: root
"""
# Class plate is being developed to be used generally.  Stiffened panels and 
# more complex geometries become easier ot analyze with its implementation.  
# It locates a plates physical dimensions and spacial location as defined by 
# the repective endpoints
class Plate:
    def __init__(self, breadth, thickness, x1, y1, x2, y2, E, panel=None, ptype=None):
        self._breadth = breadth
        self._thickness = thickness
        self._x1 = x1
        self._y1 = y1
        self._x2 = x2
        self._y2 = y2
        self._E = E
        self._ptype = ptype
        self.panel = panel
        
    def getPlatebreadth(self):
        ''' Returns the length of the plate'''
        return self._breadth
    
    def getPlatethickness(self):
        ''' Returns the thickness ofa plate'''
        return self._thickness
        
    def getPlatex1(self):
        ''' Returns the x1 coordinate of the plate'''
        return self._x1
        
    def getPlatey1(self):
        ''' Returns the y1 coordinate of the plate'''
        return self._y1
        
    def getPlatex2(self):
        ''' Returns the x2 coordinate of the plate'''
        return self._x2
        
    def getPlatey2(self):
        ''' Returns the y2 coordinate of the plate'''
        return self._y2
    
    def getE(self):
        '''Returns the modulus of elasticicty'''
        return self._E
    
    def getPType(self):
        '''Returns the plate type'''
        return self._ptype
    
    def getPanel(self):
        '''Returns a reference to an associated panel object'''
        return self.panel
        
    