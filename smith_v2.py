# -*- coding: utf-8 -*-
#Simplified implementation of the Smith collapse method for ultimate strength
#Rewrite of the original code to be more modular and easier to read
#Based on code from: Matt Lankowksi and Richard Thill
# (c) 2025 Regents of the Univesity of Michigan

import matplotlib.pyplot as plt
import logging
import TPanel_trans
import math
import HansenC as HC
import numpy as np
from scipy.interpolate import CubicSpline

logger = logging.getLogger(__name__)

class ElementBase:
    '''
    Base class for all individual elements
    '''
    def __init__(self, name, area, x_loc, y_loc):
        '''
        Initialize the element with its name, area, and location.
        Parameters  
        ----------
        name : str
            Name of the element.
        area : float
            Area of the element, m^2
        x_loc : float
            X location of the element, m, positive to starboard.
        y_loc : float
                Y location of the element, m, postive is up.
        '''
        self.name = name
        self.area = area
        self.x_loc = x_loc
        self.y_loc = y_loc
    
    def getForceAndStress(self, strain):
        '''
        Calculate the force and stress in the element based on the strain.
        Parameters
        ----------
        strain : float
            Strain in the element.
        Returns
        -------
        tuple
            Force and stress in the element, N, and MPA.
        '''
        raise NotImplementedError("This method should be implemented by subclasses.")

class ElasticElement:
    '''
    Class for a purely elastic element
    '''
    def __init__(self, name, area, x_loc, y_loc, E):
        '''
        Initialize the elastic element with its name, area, location, and elastic modulus.
        Parameters
        ----------
        name : str
            Name of the element.
        area : float
            Area of the element, m^2
        x_loc : float
            X location of the element, m, positive to starboard.
        y_loc : float
                Y location of the element, postive is up.
        E : float
            Elastic modulus of the element, MPa.
        '''
        super().__init__(name, area, x_loc, y_loc)
        self.E = E
    
    def getForceAndStress(self, strain):
        '''
        Calculate the force and stress in the elastic element based on the strain.
        Parameters
        ----------
        strain : float
            Strain in the element.
        Returns
        -------
        tuple
            Force and stress in the element, N, and MPA.
        '''
        force = self.E * self.area * strain*10e6
        stress = self.E * strain
        return force, stress

class EP_compression_element(ElementBase):
    '''
    Class for an element that is elastic-pefectly plastic in tension, and
    follows a load-shortening curve in compression.
    '''
    def __init__(self, name, area, x_loc, y_loc, sigma_yield, E, interp_compression):
        '''
        Initialize the elastic element with its name, area, location, and elastic modulus.
        Parameters
        ----------
        name : str
            Name of the element.
        area : float
            Area of the element, m^2
        x_loc : float
            X location of the element, m, positive to starboard.
        y_loc : float
                Y location of the element, postive is up.
        sigma_yield : float
            Yield stress of the element, MPa.
        strain_yield : float
            Yield strain of the element, MPa.
        interp_compression : function
            Function to interpolate the load-shortening curve in compression.
        '''
        super().__init__(name, area, x_loc, y_loc)
        self.sigma_yield = sigma_yield
        self.E = E
        self.strain_yield = self.sigma_yield / self.E #Yield strain
        self.interp_compression = interp_compression
        
    
    def getForceAndStress(self, strain):
        '''
        Calculate the force and stress in the elastic element based on the strain.
        Parameters
        ----------
        strain : float
            Strain in the element.
        Returns
        -------
        tuple
            Force and stress in the element, N, and MPA.
        '''
        stress = 0

        if (strain > self.strain_yield):
            #Plastic region
            stress = self.sigma_yield
        elif (strain > 0):
            #Elastic region
            stress = self.E * strain
        elif (strain < 0):
            #Compression region
            stress = self.interp_compression(strain) #in MPa
        #might need to add elif about last data point in compression

        #if np.isnan(stress):
            #print ("NaN in stress")
            #print ("provided strain was: ", strain)

        force = stress*self.area*10e6 #in N
        return force, stress

class SmithMethod:
    '''
    Class for the Smith method of ultimate strength.
    '''
    def __init__(self):
        '''
        Initialize the Smith method with an empty list of elements.
        '''
        self._elements = []
        self._need_update = True
        
    def discretize(self, t_panel_list):
        
        count = -1

        for panel in t_panel_list:
            data_i = HC.HansenC(panel)
            data_i._strn_flip = -np.flip(data_i._strn)
            data_i._strss_flip = -np.flip(data_i._strss)
            cs = CubicSpline(data_i._strn_flip, data_i._strss_flip, bc_type='natural')
            count += 1
            stiff_spacing = panel.getb()
            numElements = panel.getnstiff()+2
            numStiff = panel.getnstiff()
            for j in range(numElements):
                if j < numStiff:
                    y_i = panel.get_bot()+((j+1)*(stiff_spacing)*math.sin(math.radians(-panel.getOrnt())))
                    x_i = 0
                    Ys = panel.getsmatl().getYld() #ATTN! currently just uses yield strength of first element
                    E = panel.getsmatl().getE() #ATTN! currently just uses E of first element
                    yield_strn = Ys/E
                    area = panel.gettp() * panel.getb() + panel.gettw() * panel.gethw() + panel.gettf() * panel.getbf()
                    element = EP_compression_element("test", area, x_i, y_i, Ys, E, cs)
                    self.addElement(element)
                elif j == numStiff:
                    y_0 = panel.get_bot()
                    y_last = panel.get_bot()+((numStiff)*(stiff_spacing)*math.sin(math.radians(-panel.getOrnt())))
                    x_0 = 0
                    x_last = 0
                    span_0 = stiff_spacing
                    span_last = panel.getB()-(stiff_spacing*(numStiff))-span_0
                    Area_0 = panel.gettp() * span_0
                    Area_last = panel.gettp() * span_last
                    Ys = panel.getsmatl().getYld() #ATTN! currently just uses yield strength of first element
                    E = panel.getsmatl().getE() #ATTN! currently just uses E of first element
                    yield_strn = Ys/E
                    element = EP_compression_element("test", Area_0, x_0, y_0, Ys, E, cs)
                    self.addElement(element)
                    element = EP_compression_element("test", Area_last, x_last, y_last, Ys, E, cs)
                    self.addElement(element)
        return 

    def addElement(self, element):
        '''
        Adds a member of the element class to the list of elements in the cross section'
        '''
        self._elements.append(element)
        self._need_update = True
        logger.debug(f"Added element: {element.name} to Smith method.")

    def getOverallProperties(self,  mirror = True):
        '''
        Calculate the overall properties of the cross section.
        Returns
        -------
        tuple
            Area, x location, and y location of the centroid of the cross section.
        '''
        
        total_area = 0
        total_mx = 0
        total_my = 0
        total_ixx = 0
        total_iyy = 0
        total_shift_denom = 0
        yloc_list = []
        xloc_list = []

        if mirror == True:
            factor = 2 
        else:
            factor = 1

        # Calculate areas and moments about origin
        for element in self._elements:
            total_area += element.area * factor
            total_mx += element.area * element.x_loc 
            total_my += element.area * element.y_loc 
            total_ixx += element.area * (element.y_loc**2) 
            total_iyy += element.area * (element.x_loc**2) 
            total_shift_denom += element.area * element.E * 10e6 #in N
            yloc_list.append(element.y_loc)
            xloc_list.append(element.x_loc)
        
        #Catch case of now panels ready
        if total_area == 0:
            return 0, 0, 0, 0, 0, 0, 0, 0
        else:
            x_na = total_mx / total_area
            y_na = total_my / total_area
            ixx = total_ixx - total_area * y_na**2
            iyy = total_iyy - total_area * x_na**2
            return total_area, x_na, y_na, ixx, iyy, total_shift_denom, yloc_list, xloc_list
        
    def setup(self):
        '''Determine the curvature bounds for the ultimate strength calculator'''

        total_area, x_na, y_na, ixx, iyy, total_shift_denom, yloc_list, xloc_list = self.getOverallProperties()
        c1 = abs(y_na - max(yloc_list))
        c2 = abs(y_na - min(yloc_list))
        c = max(c1,c2)
        #print ("c=", c)
        Ys = self._elements[0].sigma_yield #ATTN! currently just uses yield strength of first element
        E = self._elements[0].E #ATTN! currently just uses E of first element
        YldCrv = Ys/(c*E) #yield curvature in 1/m
        #print ("Yield curvature=", YldCrv, "1/m")
        
        num_crv_inc = 20

        #for sagging condition
        Crv_max_sag = -2*YldCrv
        Crv_min_sag = Crv_max_sag/10
        Crv_step_sag = (Crv_max_sag - Crv_min_sag)/num_crv_inc

        #for hogging condition
        Crv_max_hog = 2*YldCrv
        Crv_min_hog = Crv_max_hog/10
        Crv_step_hog = (Crv_max_hog - Crv_min_hog)/num_crv_inc
        
        return Crv_min_sag, Crv_step_sag, Crv_min_hog, Crv_step_hog, num_crv_inc
    
    def plotSection(self):
        #Make a list of all X/Y locations for plotting
        x = []
        y = []
        for element in self._elements:
            x.append(element.x_loc)
            y.append(element.y_loc)
        #Plot the elements
        plt.figure()
        plt.scatter(x, y)
        plt.title("Smith Collapse Cross Section at Centers of Elements")
        plt.xlabel("X Location (m)")    
        plt.ylabel("Y Location (m)")
        plt.grid()
        plt.axis('equal')
        plt.show()

    
    def applyCurvature(self, curvature, NAGuess = None,  mirror = True):
        '''
        Apply a curvature to the cross section and calculate the forces and stresses in each element.
        Parameters
        ----------
        curvature : float
            Curvature of the cross section, 1/m.
        Returns
        -------
        tuple
            force overall in the cross section N, and moment in the cross section kN*m.
        '''
        
        #don't forget to add factor for mirroring

        if mirror == True:
            factor = 2 
        else:
            factor = 1

        #Calculate the overall properties of the cross section
        if NAGuess is None:
            area, x_na, y_na, ixx, iyy, yloc_list, xloc_list = self.getOverallProperties()
        else:
            y_na = NAGuess
            
        print ("NAGuess is: ", y_na)

        #Calculate the strain in each element
        moment = 0
        force_overall = 0 
        for element in self._elements:
            strain = curvature * -(element.y_loc - y_na)
            force, stress = element.getForceAndStress(strain)
            moment += force * (element.y_loc - y_na)*1./1000 * factor #in kN*m
            force_overall += force * factor #in N
        
        return force_overall, moment
    
    
    def solveCurvature(self, curvature, NAGuess = None):
        '''
            Solve for the bending moment of the cross section using the Smith method.
            Parameters: 
            ----------
            curvature : float
                Curvature of the cross section, 1/m.
            Returns
            -------     
            tuple
                Resisting bending moment kN*m, and NA position, meters
        '''
        area, x_na, y_na, ixx, iyy, total_shift_denom, yloc_list, xloc_list = self.getOverallProperties()

        count = 0

        if NAGuess is not None:
            y_na = NAGuess
        
        force_tol = 10 #in N
        force = 2 * force_tol #just to start
        while abs(force) > force_tol:
            force, moment = self.applyCurvature(curvature, y_na)
            if abs(force) > force_tol:
                shift = force/(curvature*total_shift_denom) #in m
                y_na -=  shift #might need to add abs(force) to get this to converge more quickly
            #print ("force is =", force)
            count += 1
            #print ("NA0 is updated to", y_na)
            if count > 100:
                print ("NA0 not converging")
                break
            
        return moment, y_na

    def getCollapseCurve(self):
        '''
            Increments through curvatures from setup function.
            Returns
            -------     
            tuple
                Resisting bending moment kN*m, and curvature 1/m
        '''
        crv_array = []
        M_array = []

        Crv_min_sag, Crv_step_sag, Crv_min_hog, Crv_step_hog, num_crv_inc = self.setup()

        #for sagging condition first
        for i in range(0,num_crv_inc+1):
            curv_applied_sag = Crv_min_sag + i*Crv_step_sag
            moment, y_na = self.solveCurvature(curv_applied_sag)
            crv_array.append(curv_applied_sag)
            M_array.append(moment)
        #for hogging condition next
        for i in range(0,num_crv_inc+1):
            curv_applied_hog = Crv_min_hog + i*Crv_step_hog
            moment, y_na = self.solveCurvature(curv_applied_hog)
            crv_array.append(curv_applied_hog)
            M_array.append(moment)

        plt.plot(crv_array,M_array, 'ro')
        plt.xlabel('Curvature [1/m]')
        plt.ylabel('VBM [kN*m]')
        plt.title('Collapse Curve for Section')
        plt.grid(True)
        plt.show()

        print ("here's crv array: ", crv_array)
        print ("here's M array: ", M_array)
        
        return M_array, crv_array
    
    def getUltimateMoment(self):
        '''
            Finds the ultimate collapse moment from the collapse curve.
            Returns
            -------     
            float
                Mult, ultimate moment kN*m
        '''
        M_array, crv_array = self.getCollapseCurve()
        #Find the max moment in the array
        Mult = max(M_array, key=abs)
        print ("here's the ultimate moment: ", Mult)
        
        return Mult
