# -*- coding: utf-8 -*-
#Simplified implementation of the Smith collapse method for ultimate strength
#Rewrite of the original code to be more modular and easier to read
#Based on code from: Matt Lankowksi and Richard Thill
# (c) 2025 Regents of the Univesity of Michigan

import matplotlib.pyplot as plt

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
    def __init__(self, name, area, x_loc, y_loc, sigma_yield, strain_yield, interp_compression):
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
        self.strain_yield = strain_yield
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
            stress = self.interp_compression(strain)

        force = stress*self.area*10e6
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
        
    
    def addElement(self, element):
        '''
        Adds a member of the element class to the list of elements in the cross section'
        '''
        self._elements.append(element)
        self._need_update = True
        logger.debug(f"Added element: {element.name} to Smith method.")

    def getOverallProperties(self):
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

        # Calculate areas and moments about origin
        for element in self._elements:
            total_area += element.area
            total_mx += element.area * element.x_loc
            total_my += element.area * element.y_loc
            total_ixx += element.area * (element.y_loc**2)
            total_iyy += element.area * (element.x_loc**2)
        
        #Catch case of now panels ready
        if total_area == 0:
            return 0, 0, 0, 0, 0
        else:
            x_na = total_mx / total_area
            y_na = total_my / total_area
            ixx = total_ixx - total_area * y_na**2  
            iyy = total_iyy - total_area * x_na**2
            return total_area, x_na, y_na, ixx, iyy

    
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

    
    def applyCurvature(self, curvature, NAGuess = None):
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
        
        #Calculate the overall properties of the cross section
        if NAGuess is None:
            area, x_na, y_na, ixx, iyy = self.getOverallProperties()
        else:
            y_na = NAGuess

        #Calculate the strain in each element
        moment = 0
        force_overall = 0 
        for element in self._elements:
            strain = curvature * (element.y_loc - y_na)
            force, stress = element.getForceAndStress(strain)
            moment += force * (element.y_loc - y_na)*1./1000.
            force_overall += force
        
        return force_overall, moment
    
    
    def solveCurvature(self, curvature):
        '''
            Solve for the curvature of the cross section using the Smith method.
            Parameters: 
            ----------
            curvature : float
                Curvature of the cross section, 1/m.
            Returns
            -------     
            tuple
                Resisting bending moment kN*m, and NA position, meters
        '''

