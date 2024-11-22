# Devine, Thomas 
# The University of Michigan, Marine Structures Design Labratory
# June 2012
"""
Section class

Uses input of a series of T-panel trans objects and takes their properties to produce 
relevant quantities a user would require.

Currently
-Section Area
-Section VerticalCentroid
-Section Moment of inertia about the horizonatal centroidal axis
"""
import sys, os, fileinput, math, string, numpy, scipy, random
import Structures, TPanel_trans, Plate
import matplotlib.pyplot as plt

# Append_Panels()
# this function appends panels to a list of panel 
# objects.
class section(object):
    def __init__(self):
        self._panels = []
    def Append_Panels(self, panel):
        self._panels.append(panel)

# Create explode command ennumerating all plates stored in the section.  Will
# be used to effectively test the moment of inertia of the stiuffened panels 
# and then the section example self.plates
# [b t, x1, y1, x2, y2]

    def Explode(self):
        # Create list for plate storage
        self._plates = []
        # Methodology to explode a Grillage.  Currently doesnot handle the 
        # transverse portion of the grillage
        for panel in self._panels:
            # First plate defined for each panels is the t_panel base
            self._plates.append(Plate.Plate(panel.getB(),\
                panel.gettp(),\
                panel.getXloc(),\
                panel.getYloc(),\
                panel.getXloc()+panel.getB()*math.cos(math.radians(panel.getOrnt())),\
                panel.getYloc()-panel.getB()*math.sin(math.radians(panel.getOrnt())),\
                panel.getmatlP().getE(),\
                ptype='plate',panel=panel))
            for N in range(int(panel.getnstiff())):
                # Append plates in a web then flange process
                # Web
                self._plates.append(Plate.Plate(panel.gethw(),\
                    panel.gettw(),\
                    panel.getXloc()+panel.getb()*(N+1)*math.cos(math.radians(panel.getOrnt())),\
                    panel.getYloc()+panel.getb()*(N+1)*-1*math.sin(math.radians(panel.getOrnt())),\
                    panel.getXloc()+panel.getb()*(N+1)*math.cos(math.radians(panel.getOrnt()))+panel.gethw()*math.cos(math.radians(panel.getOrnt())-math.pi/2),\
                    panel.getYloc()+panel.getb()*(N+1)*-1*math.sin(math.radians(panel.getOrnt()))+panel.gethw()*-1*math.sin(math.radians(panel.getOrnt())-math.pi/2),\
                    panel.getmatlS().getE(),\
                    ptype='web',panel=panel))
                #Flange
                self._plates.append(Plate.Plate(panel.getbf(),\
                    panel.gettf(),\
                    panel.getXloc()+panel.getb()*(N+1)*math.cos(math.radians(panel.getOrnt()))+panel.gethw()*math.cos(math.radians(panel.getOrnt())-math.pi/2)-panel.getbf()/2*math.cos(math.radians(panel.getOrnt())-math.pi),\
                    panel.getYloc()+panel.getb()*(N+1)*-1*math.sin(math.radians(panel.getOrnt()))+panel.gethw()*-1*math.sin(math.radians(panel.getOrnt())-math.pi/2)-panel.getbf()/2*-1*math.sin(math.radians(panel.getOrnt())-math.pi),\
                    panel.getXloc()+panel.getb()*(N+1)*math.cos(math.radians(panel.getOrnt()))+panel.gethw()*math.cos(math.radians(panel.getOrnt())-math.pi/2)+panel.getbf()/2*math.cos(math.radians(panel.getOrnt())+math.pi),\
                    panel.getYloc()+panel.getb()*(N+1)*-1*math.sin(math.radians(panel.getOrnt()))+panel.gethw()*-1*math.sin(math.radians(panel.getOrnt())-math.pi/2)+panel.getbf()/2*-1*math.sin(math.radians(panel.getOrnt())+math.pi),\
                    panel.getmatlS().getE(),\
                    ptype='flange',panel=panel))
        return self._plates
        
    # Function to perform all necessary calculations for the Section    
    def _upCalcs(self):
        self._Area = 0.0
        self._YFirstMoment = 0.0
        self._YSecondMoment = 0.0
        self._EI = 0.0
        for i in self._plates:
            # Plate Thickness * Plate Length
            self._Area += i.getPlatebreadth()*i.getPlatethickness()
            # Plate Thickness * Plate Length * Y Plate center
            self._YFirstMoment += (i.getPlatebreadth()*i.getPlatethickness()*(i.getPlatey2()+i.getPlatey1())/2)
        # First moment Arm/Area
        self._YCentroid = self._YFirstMoment/self._Area    
        for j in self._plates:
            #I_NA = Area*y_height^2/12 + Area*distance_to_centroid^2
            MoI = j.getPlatebreadth()*j.getPlatethickness()*math.pow((j.getPlatey2()-j.getPlatey1()),2)/12+j.getPlatebreadth()*j.getPlatethickness()*math.pow((((j.getPlatey2()+j.getPlatey1())/2)-self._YCentroid),2) 
            self._YSecondMoment += MoI
            self._EI += MoI * j.getE()
            
    def create_section_plot(self, show=True, mirror=False, color='k', ax_obj=None):
        
        if not ax_obj:
            plt.figure(1)
    #        plt.subplot(1,1,1)
    #        plt.axis([0,11000,-4500,6500])
            ax = plt.gca()
        else:
            ax = ax_obj
        ax.set_aspect('equal')
#        ax.set_autoscale_on(False)
        for plate in self._plates:
            plt.plot([plate.getPlatex1(), plate.getPlatex2()],[plate.getPlatey1(), plate.getPlatey2()],color,linewidth=1.5)
        if mirror:
            for plate in self._plates:
                plt.plot([-plate.getPlatex1(), -plate.getPlatex2()],[plate.getPlatey1(), plate.getPlatey2()],color,linewidth=1.5)
        if show:            
            plt.show()
        
    def getSectionArea(self):
        ''' Returns the cross sectional area of a miship's section'''
        return self._Area
        
    def getYCentroid(self):
        '''Returns the Y coordinate of the section centroid'''
        return self._YCentroid
        
    def getSectionYMOI(self):
        ''' Returns the Moment of inertia of the section about the X-X centroidal axis'''
        return self._YSecondMoment    
        
    def getEI(self):
        '''Returns EI of the entire section'''
        return self._EI            