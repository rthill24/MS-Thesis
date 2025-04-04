# Class for getting the stress and effective area of an element at any strain.
# Takes the strain, stress, and effective area data of a panel and creates an
# Element class object that stores this data and allows for interpolation of
# stress or effective area based on any strain within the data range
#
# Author: Matthew Lankowski

import numpy as np
import sys
from scipy import interpolate

import matplotlib.pyplot as plt

class Element:
    '''Class for getting the stress and effective area of a structural element
    at any strain.'''
    
    def __init__(self, panel, strain_unsort, stress_unsort, AE_unsort, xloc=0, yloc=0, qloc=['NA','NA','NA'], el_type='normal', output=0):
        '''
        Initialize an Element object, requiring a TPanel object, stress,
        strain, and effective area data, and optional attributes as follows
        
        Parameters
        ----------
        panel - a TPanel class object
        strain_unsort - strain values
        stress_unsort - stress values with indices corresponding to strain_unsort values
        AE_unsort - effective area values with indices corresponding to strain_unsort values
        xloc (optional) - x-location of neutral axis of panel
        yloc (optional) - y-location of neutral axis of panel
        qloc (optional) - qualitative location of panel components (used for PaikCorrosion method)
        el_type (optional) - used to distinguish normal panels from hard corners ('normal' or 'HC', respectively)
        
        Returns
        -------
        No return values, stores variables
        '''
        
        self._output = output        
        
        if len(strain_unsort) != len(stress_unsort) or len(strain_unsort) != len(AE_unsort) or len(stress_unsort) != len(AE_unsort):
            if self._output == 1:
                print ("ERROR - strain, stress, and effective area array sizes must match!")
            sys.exit()
        
        self._panel = panel
        self._yloc = yloc
        self._xloc = xloc
        self._qlocp = qloc[0]
        self._qlocw = qloc[1]
        self._qlocf = qloc[2]
        self._el_type = el_type
        
        if self._el_type == 'HC':
            self._AE = AE_unsort
            self._strn = np.zeros(len(self._AE))
            self._strs = np.zeros(len(self._AE))
            return
        
        #Routine for sorting data to achieve monotonically increasing strain
        #(required for interpolator)
        straindict = dict()
        for i in range(0,len(strain_unsort)):
            straindict[strain_unsort[i]] = stress_unsort[i], AE_unsort[i]
        self._strn = sorted(straindict)
        self._strs = np.zeros(len(self._strn))
        self._AEsorted = np.zeros(len(self._strn))
        count = 0
        for i in self._strn:
            self._strs[count] = straindict[i][0]
            self._AEsorted[count] = straindict[i][1]
            count += 1

        #routine for storing data in four regions: 1.plastic tension, 2.elastic
        #tension, 3.elastic compression, 4.plastic compression
        strnR1 = []
        strnR2 = []
        strnR3 = []
        strnR4 = []     
        strsR1 = []
        strsR2 = []
        strsR3 = []
        strsR4 = []
        
        AER1 = []
        AER2 = []
        AER3 = []
        AER4 = []
        
        #find the intersection of regions 3 and 4 (where max stress occurs)        
        idx = self._strs.argmax() 
        strnatmaxstress = self._strn[idx]
        minstress = min(self._strs)
        for i in range(0,len(self._strs)):
            #Define the the perfectly plastic tension Region 1 as the region of
            #maximum negative stress
            if self._strs[i] == minstress:
                strnR1.append(self._strn[i])
                strsR1.append(self._strs[i])
                AER1.append(self._AEsorted[i])
            #Define linearly elastic tension Region 2 as the region between
            #the maximum negative stress and zero stress
            elif self._strs[i] > minstress and self._strs[i] < 0:
                strnR2.append(self._strn[i])
                strsR2.append(self._strs[i])
                AER2.append(self._AEsorted[i])
            #Define nonlinear elastic compression Region 3 as the region between
            #zero and maximum positive stress
            elif self._strs[i] >= 0:                
                if self._strn[i] <= strnatmaxstress:
                    strnR3.append(self._strn[i])
                    strsR3.append(self._strs[i])
                    AER3.append(self._AEsorted[i])
                #Define plastic collapse Region 4 as the region betwen maximum
                #positive stress and maximum positive strain
                else:
                    strnR4.append(self._strn[i])
                    strsR4.append(self._strs[i])
                    AER4.append(self._AEsorted[i])
        
        #Convert lists to numpy arrays
        self._strn = np.array(self._strn)
        self._strnR1 = np.array(strnR1)
        self._strnR2 = np.array(strnR2)
        self._strnR3 = np.array(strnR3)
        self._strnR4 = np.array(strnR4)
        self._strsR1 = np.array(strsR1)
        self._strsR2 = np.array(strsR2)
        self._strsR3 = np.array(strsR3)
        self._strsR4 = np.array(strsR4)
        
        self._AER1 = np.array(AER1)
        self._AER2 = np.array(AER2)
        self._AER3 = np.array(AER3)
        self._AER4 = np.array(AER4)
        
        self._AE = np.array(self._AEsorted)
        
        #Define interpolation schemes for different regions
        if len(self._strnR1) == 0:
            self._strnR1 = np.zeros(2)
            self._strsR1 = np.zeros(2)
            self._intrpstrsR1 = 0
        else:
            self._intrpstrsR1 = min(self._strsR1)
        self._intrpAER1 = self._panel.getArea()
        if len(self._strnR2) == 0:
            self._strnR2 = np.zeros(2)
            self._strsR2 = np.zeros(2)
            self._intrpstrsR2 = 0
        else:
            self._intrpstrsR2 = interpolate.InterpolatedUnivariateSpline(self._strnR2, self._strsR2)
        self._intrpAER2 = self._panel.getArea()
        if len(self._strnR3) == 0:
            self._strnR3 = np.zeros(2)
            self._strsR3 = np.zeros(2)
            self._AER3 = np.ones(2)*self._panel.getArea()
            self._intrpstrsR3 = 0
            self._intrpAER3 = self._panel.getArea()
        else:
            self._intrpstrsR3 = interpolate.InterpolatedUnivariateSpline(self._strnR3, self._strsR3)
            self._intrpAER3 = interpolate.InterpolatedUnivariateSpline(self._strnR3, self._AER3)
        if len(self._strnR4) == 0:
            self._strnR4 = np.zeros(2)
            self._strsR4 = np.zeros(2)
            self._AER4 = np.ones(2)*self._panel.getArea()
            self._intrpstrsR4 = 0
            self._intrpAER4 = self._panel.getArea()
        else:
            self._intrpstrsR4 = interpolate.InterpolatedUnivariateSpline(self._strnR4, self._strsR4)
            self._intrpAER4 = self._AEsorted[-1] #doesn't change in region 4 so just use last
        
        '''store min's and max's in each region for getStress()'''
        self.strnR1_min = min(self._strnR1)
        self.strnR1_max = max(self._strnR1)
        self.strnR2_min = min(self._strnR2)
        self.strnR2_max = max(self._strnR2)
        self.strnR3_min = min(self._strnR3)
        self.strnR3_max = max(self._strnR3)
        self.strnR4_min = min(self._strnR4)
        self.strnR4_max = max(self._strnR4)
        self.AE_min = min(self._AE)
        self.AE_max = max(self._AE)
        

    def getStress(self, straini):
        '''return stress at a given strain'''

        if self._el_type == 'HC':
            HC_elastic_stress = straini*self._panel._smatl.getE()
            if abs(HC_elastic_stress) <= self._panel.getYsavg():
                return HC_elastic_stress
            else:
                if straini > 0:
                    return self._panel.getYsavg()
                else:
                    return -self._panel.getYsavg()
     
        if straini < self.strnR1_min:
            if self._output == 1:
                print ("WARNING - strain below interpolation range. Stress w/ min strain in data set used")        
#            return self._strsR1[0]
            return -1
        elif straini >= self.strnR1_min and straini <= self.strnR1_max:
            return self._intrpstrsR1
        elif straini > self.strnR1_max and straini < self.strnR2_min:
            return (self._strsR1[len(self._strsR1)-1] + self._strsR2[0]) / 2
        elif straini >= self.strnR2_min and straini <= self.strnR2_max:
            return self._intrpstrsR2(straini)
        elif straini > self.strnR2_max and straini < self.strnR3_min:
            return (self._strsR2[len(self._strsR2)-1] + self._strsR3[0]) / 2
        elif straini >= self.strnR3_min and straini <= self.strnR3_max:
            return self._intrpstrsR3(straini)
        elif straini > self.strnR3_max and straini < self.strnR4_min:
            return (self._strsR3[len(self._strsR3)-1] + self._strsR4[0]) / 2
        elif straini >= self.strnR4_min and straini <= self.strnR4_max:
            return self._intrpstrsR4(straini)
        elif straini > self.strnR4_max:
            if self._output == 1:
                print ("WARNING - strain above interpolation range. Stress w/ max strain in data set used")    
#            return self._strsR4[-1]
            return -1
                
    def getAE(self, straini):                
        '''return effective area at a given strain'''
        if self._el_type == 'HC':
            return self._panel.getArea()
            
        if straini <= 0:
            return self._panel.getArea()
        elif straini > 0 and straini <= self.strnR3_max:
            return self._intrpAER3(straini)
        elif straini > self.strnR3_max:
            return self._intrpAER4
        
    def getForce(self, straini):
        '''return force at a given strain'''
        stress = float(self.getStress(straini))
        area = self.getPanel().getArea()
        
        return stress * area
    
    def getPanel(self):
        '''return t_panel object'''
        return self._panel
                
    def getYloc(self):
        '''return y-location (vertical), yloc'''
        return self._yloc
        
    def getXloc(self):
        '''return x-location (horizontal), xloc'''
        return self._xloc
        
    def Copy(self):
        return self
        
    def Update(self, panel, strain_unsort, stress_unsort, AE_unsort):
        '''the purpose of this function is to be able to update element without
        deleting other custom attributes not defined in init.'''
        self.__init__(panel, strain_unsort, stress_unsort, AE_unsort, self._xloc, self._yloc, [self._qlocp, self._qlocw, self._qlocf], self._el_type)