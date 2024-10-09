import numpy as np
import copy
import math
import sys
import shutil
import os
import Plate
import Cost
from operator import itemgetter
import Section
import matplotlib.pyplot as plt
import deterioratingStructure
import Corrosion as corrode
import Structures

class Midship_Repair_Cost(deterioratingStructure.TPanel_Repair):
    '''
    A repair class for the midship section class - can be extended to consider
    whatever overhaul, drydocking, or mainteance costs are neccessary
    '''
    
    def __init__(self, cost, recoat_cost):
        '''
        Constructor takes the renewal cost, crack repair cost, and recoat cost
        '''   
        deterioratingStructure.TPanel_Repair.__init__(self, cost, recoat_cost)
    
    def dry_docked_repair(self, multiplier):
        '''
        Repair cost given that the repair must be done at dry dock
        
        This returns the renewal cost of the midship section piece with a multiplier
        based on the cost to dry dock the vessel in order to make the repair.
        
        Parameters
        -----------
        No parameters
        
        Returns
        --------
        dry_dock_cost:  Cost value
        '''
        
        renew_cost = self.cost
        dry_dock_cost = renew_cost * multiplier
        
        return dry_dock_cost
    
    def underway_repair_cost(self, multipler):
        '''
        Repair cost given that the repair must be done underway
        
        This returns the renewal cost of the midship section piece with a multiplier
        based on the cost to do an underway repair for the vessel.
        
        Parameters
        ----------
        No parameters
        
        Returns
        --------
        underway_cost:  Cost value
        '''
        
        renew_cost = self.cost
        underway_cost = renew_cost * multipler
        
        return underway_cost

class Midship_Section(object):
    '''
    A class that represents a vessels midship section.
    
    This class uses a list of TPanel_Trans instances to define a midship section and can
    age, repair, and renew the grillages.  It also can store a maintenance schedule to be used
    when calculating total lifetime mainteance costs.
    
    Parameters
    -----------
    grillage_list:      List of TPanel_Trans instances
    coatLife:           The coat life of the structure to be used with Paik's corrosion model
    limit_tp:           Limiting percentage of plate thickness before repair is needed (default=.25)
    limit_tws:          Limiting percentage stiffener web thickness before repair needed (default=.25)
    limit_tfs:          Limiting percentage stiffener flange thick. before repair needed (default=.25)
    limit_twf:          Limiting percentage frame web thickness before repair needed (default=.25)
    limit_tff:          Limiting percentage frame flange thick. before repair needed (default=.25)
    age:                Initial  age for both structure and coating, if not zero (default=0)
    emergency_dry_dock: The cost of undergoing an emergency dry dock to perform a repair for
                            the vessel (default=5e5)
    planned_drydock:    The cost of undergoing a planned dry dock to perform scheduled maintenance
                            on the vessel (default=2.5e5)
    underway:           The cost of performing underway maintenance for the vessel (default=5e4)
    dry_dock_dict:      A dictionay of functional locations that are associated with instances of TPanel_trans
                            set either to 1 if they require a dry dock to replace or 0 if they can be repaired
                            underway.  If no dictionary is supplied functional locations based on Paik's corrosion
                            model will be used.  (default=None)
    '''
    def __init__(self, grillage_list, recoat_cost, coatLife=7,  limit_tp = .75, limit_tws = .75, limit_tfs = .75,
                 limit_twf = .75, limit_tff = .75, age = 0.,emergency_dry_dock=5e5,planned_dry_dock=2.5e5,underway=5e4, dry_dock_dict=None):
        
        #if not dry_dock_dict:
            #self.dry_dock_dict = {'BSH':1, 'BSLBW':1, 'BSLBF':1, 'BSV':1, 'SSLBW':1, 'SSLBF':1, 'ABV':0, 'ABH':0, 'DLBW':0, 'AOH':0, 'NA':0}
        #else:
            #self.dry_dock_dict = dry_dock_dict
        
        #self.emergency_dry_dock = emergency_dry_dock
        #self.planned_dry_dock = planned_dry_dock
        #self.underway = underway
        #self.coatLife = coatLife
        
        #self.section_age = age
        #self.schedule = []

        ### Take list of T-Panels and convert them into TransTPanelDet objects - 
        ### deteriorating T-Panel objects which can be aged, repaired, and fatigued
        self.grillages = []
        for i in grillage_list:  
            coster = Cost.Cost_trans(i)
            renewal_cost = coster.Total_Cost_()
            corrosion_model = corrode.paikCorrosion(coatLife, locationCodePlate=i.qloc[0], locationCodeWeb=i.qloc[1], locationCodeFlange=i.qloc[2])
            self.grillages.append(deterioratingStructure.TransTPanelDet(Midship_Repair_Cost(renewal_cost, recoat_cost),i,corrosion_model,limit_tp = limit_tp*i.gettp(), \
                                                 limit_tws = limit_tws*i.gettw(), limit_tfs = limit_tfs*i.gettf(),limit_twf = limit_twf*i.gettwt(), limit_tff = 1))
        
        #Get basic data about midship section: Neutral axis, EI
        section_analysis = Section.section()
        for grill in self.grillages:
            panel = grill.getTTPanRef()
            section_analysis.Append_Panels(panel)
        section_analysis.Explode()
        section_analysis._upCalcs()
        
        #compression = PC.PaikCReg()
        #self.EIT, self._NAy, self.area, self.volume = self.section_data()
        #self.EIT /= 10**9
        #self.fatigue_loaders = []
        #self.initial_SM = self.vessel_SM()
        #self.initial_fUCS = []
        #for panel in grillage_list:
            #self.initial_fUCS.append( (compression.ucs(panel)/panel.getYsavg()) )
    
    def section_modulii(self):
        '''
        Calculates a list of section modulii for the mid ship section.
        
        Parameters
        ----------
        No parameters
        
        Returns
        --------
        SM: Array of section modulii
        '''
        
        SM = np.zeros(len(self.grillages))
        for i in range(len(self.grillages)):
            panel = self.grillages[i].getTTPanRef()
            SM[i] = panel.getINA() / panel.gety_max()
        
        return SM
            
    
    def section_data(self):
        '''
        Calculates the neutral axis and EI of the midship section.
        
        This method will use the Section.py analysis module to explode each of the
        TPanels into individual plates and calculate the moment, centroid, total area,
        and volume values.
        
        Parameters
        -----------
        No parameters
        
        Returns
        --------
        EI:     Total EI of the section
        NAy:    Y location of the neutral axis
        area:   Total cross-sectional area of the section
        volume: Volume of the structure
        '''

        section_analysis = Section.section()
        volume = 0.0
        for grill in self.grillages:
            panel = grill.getTTPanRef()
            volume += panel.getTotalVolume()
            section_analysis.Append_Panels(panel)
        section_analysis.Explode()
        section_analysis._upCalcs()

        EI = section_analysis.getEI()
        NAy = section_analysis.getYCentroid()
        area = section_analysis.getSectionArea()
        
        return EI, NAy, area, volume
    
    def production_cost(self):
        '''
        Calculates the production cost of the midship section using Caldwell's
        costing method
        
        Parameters
        ----------
        No parameters
        
        Returns
        --------
        production_cost:    The total production cost for the section in its __CURRENT__ state
        '''
        
        production_cost = 0.0
        for panel in self.grillages:
            tpanel = panel.getTTPanRef()
            cost_object = Cost.Cost_trans(tpanel)
            production_cost += cost_object.Total_Cost_()
        
        
        return production_cost
    
    
    def vessel_SM(self):
        '''
        Calculates the section modulus of the entitre vessel.
        
        Parameters
        ----------
        None
        
        Returns
        --------
        SM: Section modulus
        '''
        
        sm_section = Section.section()
        for dpanel in self.grillages:
            tpanel = dpanel.getTTPanRef()
            sm_section.Append_Panels(tpanel)
        
        sm_section.Explode()
        sm_section._upCalcs()
        
        MOI = sm_section.getSectionYMOI()
        Y = sm_section.getYCentroid()
        d = []
        for dpanel in self.grillages:
            tpanel = dpanel.getTTPanRef()
            y = tpanel.sloc[1]
            d.append(math.fabs(y-Y))
        
        c = max(d)
        SM = MOI / c
        
        return SM
    
    def plot_section(self, show=False, mirror=True, ax_obj=None):
        '''Plots the section'''
        
        section_plot = Section.section()
        for grill in self.grillages:
            panel = grill.getTTPanRef()
            section_plot.Append_Panels(panel)
        section_plot.Explode()
        section_plot.create_section_plot(show=False, mirror=mirror, ax_obj=None)
        plt.xlabel('Longitudinal Position from Centerline (m)', fontsize=14)
        plt.ylabel('Vertical Position from Keel (m)', fontsize=14)
        if show:
            plt.show()
        