
import TPanel_trans
import copy
import math
import Structures
import logging
from numpy import zeros
import midship_section
import Plate
import sys
import shutil
import os
import Section
import scipy
import sys, os, fileinput, math, string, numpy, scipy, random
import matplotlib.pyplot as plt
from operator import itemgetter
from scipy import integrate
import deterioratingStructure
import SmithCollapse
import HansenC
import Element    


def set_up_smith_collapse(self, extra_outputs=False):
        '''
        Set up smith collapse analysis for the section
        
        This method explodes the current midship section into individual stiffened
        panels and sets up a Smith-type progressive analysis for the panels.  The method
        will return the collapse analysis object that can be used to generate moment-curvature
        curves, find ultimate moments, etc...
        
        Parameters
        ----------
        extra_outputs:  Flag to return outputs of area, position, and momements for the elements
                            added to the Smith collapse analysis.  May be necessary for certain 
                            calculations (default=False)
        
        Returns
        -------
        collapse:       A reference to a SmithCollapseC object based on the current midship
                            section that can be used to perform a Smith type progressive
                            collapse analysis.
        miny:           The smallest y location within the exploded stiffened panels.  Note: this
                            value is only returned if the 'extra_outputs' flag is set to true.
        maxy:           The largest y location within the exploded stiffened panels.  Note: this
                            value is only returned if the 'extra_outputs' flag is set to true.   
        total_moment:   The total moment of the explded stiffened panels.  Note: this
                            value is only returned if the 'extra_outputs' flag is set to true. 
        total_area:     The total area of the explded stiffened panels.  Note: this
                            value is only returned if the 'extra_outputs' flag is set to true. 
        '''
              
        #collapse = SmithCollapse.SmithCollapseC('SmithCollapse', forcetol=200000.)
        #Discretize the elements
        XSection = []
        total_moment = 0.0
        total_area = 0.0
        miny = 1e20
        maxy = 0.0
        nPanels = 0
        for i in range(len(self.grillages)):
            
            current_panel = self.grillages[i].getTTPanRef()
            if current_panel.type != 'HC':
                
                #Calculate stress, strain, and area using Hansen C
                HCRelations = HansenC.HansenC(current_panel)
                strn = HCRelations._strn
                strss = HCRelations._strss                
                AE = numpy.max(HCRelations._AE)
                
                #Calculate panel spacing
                panel_xspacing = current_panel.getb()*math.cos(math.radians(current_panel.ornt))
                panel_yspacing = -current_panel.getb()*math.sin(math.radians(current_panel.ornt))
                
                #Explode grillage into individual stiffened panels and add to Smith collapse analysis
                numElements = current_panel.getnstiff()+2
                for j in range(numElements-1):
                    if j < numElements - 2:
                        #Get plate properties and add to smith collapse analysis
                        elx = current_panel.sloc[0] + (j+1)*panel_xspacing - current_panel.getNA()*math.sin(math.radians(current_panel.ornt + 180))
                        ely = current_panel.sloc[1] + (j+1)*panel_yspacing + current_panel.getNA()*math.cos(math.radians(current_panel.ornt))
                        if ely > maxy:
                            maxy = ely
                        if ely < miny:
                            miny = ely
                        E = current_panel.getmatlP().getE()
                        yld = current_panel.getmatlP().getYld()
                        total_moment += ely * AE
                        total_area += AE
                        XSection.append(Element.Element(XSection, strn, strss, AE, ely))
                        nPanels += 1
                    else:
                        #Create elements corresponding to end plates
                        
                        #Get coordinates of 1st plate element
                        elx = current_panel.sloc[0] + 0.25*panel_xspacing 
                        ely = current_panel.sloc[1] + 0.25*panel_yspacing + 0.5*current_panel.gettp()*math.cos(math.radians(current_panel.ornt))
                        #Get coordinates of 2nd plate element
                        elx_2 = elx + panel_xspacing * (current_panel.getnstiff() + 0.5)
                        ely_2 = ely + panel_yspacing * (current_panel.getnstiff() + 0.5)
                        E = current_panel.getmatlP().getE()
                        yld = current_panel.getmatlP().getYld()
                        AE = 0.5 * current_panel.getb() * current_panel.gettp()
                        XSection.append(Element.Element(XSection, strn, strss, AE, ely))
                        XSection.append(Element.Element(XSection, strn, strss, AE, ely_2))
                        
        if extra_outputs:
            return XSection, maxy, miny, total_moment, total_area
        else:
            return XSection
        
        
