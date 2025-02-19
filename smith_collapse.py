# -*- coding: utf-8 -*-
#Simplified implementation of the Smith collapse method for ultimate strength
# (c) 2025 Regents of the Univesity of Michigan

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.integrate import PchipInterpolator

class Element:
    '''
    Basic class for any structural element within Smith's method
    '''

    def __init__(self, idString, stress, strain, loc_x, loc_y, area, interpType='cubic'):
        '''
        Initialize the element with its properties

        Parameters:
        -----------
        idString:   String identifying the element for output
        stress:     Numpy array of stress values, positive for tension, negative for compression, MPa
        strain:     Numpy array of strain values, unitless, should cover negative (compression) and positive (tension) values.  Must be monotonically increasing
        loc_x:      Location of element in the cross-section, in the horizontal (x) direction, mm
        loc_y:      Location of element in the cross-section, in the vertical (y) direction, mm
        area:       Area of the element, mm^2
        interpType  Interpolation type for the stress-strain curve, default is 'cubic' other options is 'Pchip', scipy's piecewise cubic hermite interpolating polynomial that does not overshoot the data for data with sharp changes.
        '''

        #Check that the stress and strain arrays are the same length
        if len(stress) != len(strain):
            raise ValueError('Stress and strain arrays must be the same length')
        
        #Check that a valid option was passed for the interpolation type
        if interpType == 'cubic':
            self.interp = CubicSpline(strain, stress)
        elif interpType == 'Pchip':
            self.interp = PchipInterpolator(strain, stress)
        else:
            raise ValueError('Interpolation type must be cubic or Pchip')
        

        self.idString = idString
        self.area = area
        self.loc_x = loc_x
        self.loc_y = loc_y
        self.stress = stress
        self.strain = strain

        self.min_strain = strain[0]
        self.max_strain = strain[-1]

        self.stress_at_min_strain = stress[0]
        self.stress_at_max_strain = stress[-1]


    
    def force(self, strain):
        '''
        Calculate the force in the element based on the curvature of the cross-section

        Parameters:
        -----------
        strain:  local strain at the element, unitless

        Returns:
        --------
        Force in the element, N
        '''

        stress = 0

        if strain < self.min_strain:
            stress = self.stress_at_min_strain
        elif strain > self.max_strain:
            stress = self.stress_at_max_strain
        else:
            stress = self.interp(strain)

        #Calculate the force
        force = stress * self.area

        return force


class SmithCollapse:
    '''
    Class to implement the Smith collapse method for ultimate strength
    '''

    def __init__():
        pass

    def addElement(self, element):
        '''Allow us to build up the cross section slowly'''
        pass

    def ApplyCuravture(self, curvature, NAy):
        '''Takes a curvature and a proposed NA position, calcs resulting moment and force'''
        pass

    def plotCollapse(self):
        '''Iterates through curvature values and plots the collapse curve
        finding the equalibrium point on each'''
        pass
    
