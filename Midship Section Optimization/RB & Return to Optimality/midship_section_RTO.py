### Module to represent a vessel's midships section with the ability to age, repair, and calculate lifetime costs
### based on a given structural definition  input as a list of TPanel instances.

### Based on the code by Matthew Lankowski and extended by Dylan Temple

import sys
sys.path.insert(0, 'C:/Users/rthill/Documents/MS-Thesis')

import pystra as ra
import numpy as np
import copy
import math
import sys
import shutil
import os
import Plate
import TPanel_trans
import Cost
import Section
import SmithCollapse
import HansenC
import Corrosion as corrode
import CrackDetailVL
import deterioratingStructure
from operator import itemgetter
from scipy import integrate
import PaikCompression as PC
import matplotlib as plt

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
    
    def underway_repair_cost(self, multiplier):
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
        underway_cost = renew_cost * multiplier
        
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
        
        if not dry_dock_dict:
            self.dry_dock_dict = {'BSH':1, 'BSLBW':1, 'BSLBF':1, 'BSV':1, 'SSLBW':1, 'SSLBF':1, 'ABV':0, 'ABH':0, 'DLBW':0, 'AOH':0, 'NA':0}
        else:
            self.dry_dock_dict = dry_dock_dict
        
        self.emergency_dry_dock = emergency_dry_dock
        self.planned_dry_dock = planned_dry_dock
        self.underway = underway
        self.coatLife = coatLife
        
        self.section_age = age
        self.schedule = []
        
        ### Take list of T-Panels and convert them into TransTPanelDet objects - 
        ### deteriorating T-Panel objects which can be aged, repaired, and fatigued
        self.grillages = []
        for i in grillage_list:  
            coster = Cost.Cost_trans(i)
            renewal_cost = coster.Total_Cost_()
            corrosion_model = corrode.paikCorrosion(coatLife, locationCodePlate=i.qloc[0], locationCodeWeb=i.qloc[1], locationCodeFlange=i.qloc[2])
            self.grillages.append(deterioratingStructure.TransTPanelDet(Midship_Repair_Cost(renewal_cost, recoat_cost),i,corrosion_model,limit_tp = limit_tp*i.gettp(), \
                                                 limit_tws = limit_tws*i.gettw(), limit_tfs = limit_tfs*i.gettf(),limit_twf = limit_twf*i.gettwt(), limit_tff = limit_tff*i.gettft(), age = age))
        
        #Get basic data about midship section: Neutral axis, EI
        section_analysis = Section.section()
        for grill in self.grillages:
            panel = grill.getTTPanRef()
            section_analysis.Append_Panels(panel)
        section_analysis.Explode()
        section_analysis._upCalcs()
        
        """ compression = PC.PaikCReg()
        self.EIT, self._NAy, self.area, self.volume = self.section_data()
        self.EIT /= 10**9
        self.fatigue_loaders = []
        self.initial_SM = self.vessel_SM()
        self.initial_fUCS = []
        for panel in grillage_list:
            self.initial_fUCS.append( (compression.ucs(panel)/panel.getYsavg()) ) """
    
    """ def section_modulii(self):
        '''
        Calculates the SM_top and SM_bot for the midship section.
        
        Parameters
        ----------
        No parameters
        
        Returns
        --------
        SM_bot:    The SM value for the bottom of the midship section
        SM_top:    The SM value for the top of the midship section
        '''
        

        SM = np.zeros(len(self.grillages))
        for i in range(len(self.grillages)):
            #panel = self.grillages[i].getTTPanRef()
            #y_max = panel.gety_max()
        
        y_top = max(y_max) """
        

            
    
    def section_data(self, mirror = True): #"mirror" will double the calculations for a symmetric section
        '''
        Calculates the midship section's properties
        
        This method will use the Section.py analysis module to explode each of the
        TPanels into individual plates and calculate the moment, centroid, total area,
        weight, I_NA, and minimum section modulus.
        
        Parameters
        -----------
        No parameters
        
        Returns
        --------
        EI:     Total EI of the section
        NAy:    Y location of the neutral axis
        area:   Total cross-sectional area of the section
        weight: weight of the structure
        I_NA:   Moment of inertia about the neutral axis
        SM_min: Minimum section modulus of the section
        My:     Yield moment of the section
        Mult:   Ultimate moment of the section based on specified yield strength
        '''
        density = self.grillages[0].getTTPanRef().getmatlP().getDensity() #only uses plating density
        yield_strength = self.grillages[0].getTTPanRef().getmatlP().getYld() #only uses plating yield strength

        section_analysis = Section.section()
        volume = 0.0
        maxy = 0
        miny = 0
        Mult = 0
        if mirror == True:
            factor = 2 
        else:
            factor = 1

        for grill in self.grillages:
            panel = grill.getTTPanRef()
            volume += panel.getTotalVolume() * factor
            section_analysis.Append_Panels(panel)
        weight = volume * density * factor
        section_analysis.Explode()
        section_analysis._upCalcs()

        EI = section_analysis.getEI() * factor
        NAy = section_analysis.getYCentroid()
        area = section_analysis.getSectionArea() * factor
        I_NA = section_analysis.getSectionYMOI() * factor

        top = np.zeros(len(self.grillages))
        bot = np.zeros(len(self.grillages))
        d = np.zeros(len(self.grillages))
        A = np.zeros(len(self.grillages))
        for i in range(len(self.grillages)):
            panel = self.grillages[i].getTTPanRef()
            top[i] = panel.get_top()
            bot[i] = panel.get_bot()
            d[i] = abs(((panel.get_top() + panel.get_bot())/2)-NAy)
            A[i] = panel.getArea() * factor
            Mult += d[i] * A[i] * yield_strength * 1000

        maxy = max(top)
        miny = min(bot)
        y_top = maxy
        y_bot = miny

        c_top = y_top - NAy
        c_bot = NAy - y_bot
        SM_top = I_NA / c_top
        SM_bot = I_NA / c_bot

        SM_min = min(SM_top, SM_bot)

        My = yield_strength * SM_min * 1000 #in kN*m
        
        return EI, NAy, area, weight, I_NA, SM_min, My, Mult
    
    def production_cost(self, mirror = True): #mirror will double the calculations for a symmetric section
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
        if mirror == True:
            factor = 2
        else:
            factor = 1
        production_cost = 0.0
        for panel in self.grillages:
            tpanel = panel.getTTPanRef()
            cost_object = Cost.Cost_trans(tpanel)
            production_cost += cost_object.Total_Cost_() * factor
        
        
        return production_cost

    def HG_stress (self, M_tot = 33052 , mirror = True): #M_tot is the total moment in the midship section in kN*m
        '''
        Calculates the hull girder bending stress for each panel at most extreme fiber for each panel. Returns compressive (+) stresses for use in reliability analysis.
        
        Parameters
        -----------
        M_tot : The total moment in the midship section in kN*m
        
        Returns
        --------
        sigma_HG:   The hull girder bending stress for each panel in the midship section based on most extreme fiber of each panel (in MPa)
        '''

        if mirror == True:
            factor = 2 
        else:
            factor = 1

        section_analysis = Section.section()
        for grill in self.grillages:
            panel = grill.getTTPanRef()
            section_analysis.Append_Panels(panel)
        section_analysis.Explode()
        section_analysis._upCalcs()

        NA_y = section_analysis.getYCentroid()
        I_NA = section_analysis.getSectionYMOI() * factor

        top = np.zeros(len(self.grillages))
        bot = np.zeros(len(self.grillages))
        c_top = np.zeros(len(self.grillages))
        c_bot = np.zeros(len(self.grillages))
        c = np.zeros(len(self.grillages))
        sigma_HG = np.zeros(len(self.grillages))

        for i in range(len(self.grillages)):
            panel = self.grillages[i].getTTPanRef()
            top[i] = panel.get_top()
            bot[i] = panel.get_bot()
            c_top[i] = abs(top[i] - NA_y)
            c_bot[i] = abs(bot[i] - NA_y)
            c[i] = max(c_top[i], c_bot[i])
            sigma_HG[i] = ((M_tot*c[i])/I_NA) / 1000 #convert to MPa

        return sigma_HG
    
    def Hughes_Panel (self, waterline, rho, p_design, sigma_HG):

        '''
        Calculates the ultimate axial stress for each panel as a ratio against yield strength of the panel material by checking against Mode I, II, and III collapse
        Reference from Chapter 14 of Hughes, Owen F. Ship Structural Design: A Rationally Based, Computer-Aided Optimization Approach. Cambridge University Press, 1996.

        Parameters
        ----------
        waterline: the height of the design waterline above baseline
        rho: density of the water
        stiff_spacing : The spacing between stiffeners in the panel
        frame_spacing: The spacing between frames in the panel
        p_design: The design pressure for the plating on the vessel
        sigma_HG: The hull girder bending stress for each panel in the midship section based on most extreme fiber

        
        Returns
        --------
        sigma_a_ult:    The ultimate axial stress for each panel in the midship section based on minumum of Mode I, II, or III collapse
        '''
        E = self.grillages[0].getTTPanRef().getmatlP().getE() #only uses plating Young's modulus
        sig_ys = self.grillages[0].getTTPanRef().getsmatl().getYld() #stiffener yield strength
        sig_yp = self.grillages[0].getTTPanRef().getmatlP().getYld() #plating yield strength
        NA_y = self.section_data()[1]

        section_analysis = Section.section()
        for grill in self.grillages:
            panel = grill.getTTPanRef()
            section_analysis.Append_Panels(panel)
        section_analysis.Explode()
        section_analysis._upCalcs()

        #Mode I Initializations
        top = np.zeros(len(self.grillages))
        bot = np.zeros(len(self.grillages))
        mid = np.zeros(len(self.grillages))
        p_total = np.zeros(len(self.grillages))
        sig_a = np.zeros(len(self.grillages))
        a = np.zeros(len(self.grillages))
        b = np.zeros(len(self.grillages))
        b_stiff = np.zeros(len(self.grillages))
        M_o = np.zeros(len(self.grillages))
        t_p = np.zeros(len(self.grillages))
        t_w = np.zeros(len(self.grillages))
        h_w = np.zeros(len(self.grillages))
        t_f = np.zeros(len(self.grillages))
        A_p = np.zeros(len(self.grillages))
        A_w = np.zeros(len(self.grillages))
        A_f = np.zeros(len(self.grillages))
        A = np.zeros(len(self.grillages))
        M_p = np.zeros(len(self.grillages))
        M_w = np.zeros(len(self.grillages))
        M_f = np.zeros(len(self.grillages))
        NA = np.zeros(len(self.grillages))
        d_p = np.zeros(len(self.grillages))
        d_w = np.zeros(len(self.grillages))
        d_f = np.zeros(len(self.grillages))
        y_f = np.zeros(len(self.grillages))
        I_p_i = np.zeros(len(self.grillages))
        I_p_NA = np.zeros(len(self.grillages))
        I_w_i = np.zeros(len(self.grillages))
        I_w_NA = np.zeros(len(self.grillages))
        I_f_i = np.zeros(len(self.grillages))
        I_f_NA = np.zeros(len(self.grillages))
        I = np.zeros(len(self.grillages))
        del_o = np.zeros(len(self.grillages))
        delta = np.zeros(len(self.grillages))
        sig_e = np.zeros(len(self.grillages))
        phi = np.zeros(len(self.grillages))
        sig_f = np.zeros(len(self.grillages))
        rho_NA = np.zeros(len(self.grillages))
        lamb = np.zeros(len(self.grillages))
        eta = np.zeros(len(self.grillages))
        mu = np.zeros(len(self.grillages))
        zeta = np.zeros(len(self.grillages))
        R_I = np.zeros(len(self.grillages))
        sig_a_u_I = np.zeros(len(self.grillages))

        #Mode II Initializations
        beta_II = np.zeros(len(self.grillages))
        zeta_II = np.zeros(len(self.grillages))
        T = np.zeros(len(self.grillages))
        tau_II = np.zeros(len(self.grillages))
        sig_F_II = np.zeros(len(self.grillages))
        b_II = np.zeros(len(self.grillages))
        M_o_II = np.zeros(len(self.grillages))
        A_p_II = np.zeros(len(self.grillages))
        A_II = np.zeros(len(self.grillages))
        M_p_II = np.zeros(len(self.grillages))
        NA_II = np.zeros(len(self.grillages))
        d_p_II = np.zeros(len(self.grillages))
        d_w_II = np.zeros(len(self.grillages))
        d_f_II = np.zeros(len(self.grillages))
        y_p_II = np.zeros(len(self.grillages))
        I_p_i_II = np.zeros(len(self.grillages))
        I_p_NA_II = np.zeros(len(self.grillages))
        I_w_i_II = np.zeros(len(self.grillages))
        I_w_NA_II = np.zeros(len(self.grillages))
        I_f_i_II = np.zeros(len(self.grillages))
        I_f_NA_II = np.zeros(len(self.grillages))
        I_II = np.zeros(len(self.grillages))
        del_o_II = np.zeros(len(self.grillages))
        rho_NA_II = np.zeros(len(self.grillages))
        h_II = np.zeros(len(self.grillages))
        delta_P = np.zeros(len(self.grillages))
        y_NA_i_II = np.zeros(len(self.grillages))
        H = np.zeros(len(self.grillages))
        delta_H = np.zeros(len(self.grillages))
        delta_II = np.zeros(len(self.grillages))
        lamb_II = np.zeros(len(self.grillages))
        eta_II = np.zeros(len(self.grillages))
        eta_p_II = np.zeros(len(self.grillages))
        mu_II = np.zeros(len(self.grillages))
        zeta_R_II = np.zeros(len(self.grillages))
        R_II = np.zeros(len(self.grillages))
        sig_a_u_II = np.zeros(len(self.grillages))
        sig_a_u = np.zeros(len(self.grillages))
        Rat_a_u = np.zeros(len(self.grillages)) #remove this one after testing

        """ #Mode III Initializations
        M_o_g = np.zeros(len(self.grillages))
        del_o_g = np.zeros(len(self.grillages))
        mu_GH = np.zeros(len(self.grillages))
        y_f_II = np.zeros(len(self.grillages))
        eta_p_GH = np.zeros(len(self.grillages))
        eta_GH = np.zeros(len(self.grillages))
        lamb_GH = np.zeros(len(self.grillages))
        zeta_GH = np.zeros(len(self.grillages))
        R_GH = np.zeros(len(self.grillages))
        mu_III = np.zeros(len(self.grillages))
        zeta_III = np.zeros(len(self.grillages))
        R_III = np.zeros(len(self.grillages))
        sig_a_u_III = np.zeros(len(self.grillages)) """

        #Mode I Failure Calculations
        for i in range(len(self.grillages)):
            panel = self.grillages[i].getTTPanRef()
            top[i] = panel.get_top()
            bot[i] = panel.get_bot()
            mid[i] = (top[i] + bot[i])/2
            if mid[i] < waterline:
                avg_water_depth = waterline - mid[i]
                p_hydro = rho * 9.81/1000 * avg_water_depth
                p_tot = p_design + p_hydro
            else: 
                p_tot = 0

            p_total[i] = p_tot/1000 #convert to MPa
            sig_a[i] = sigma_HG[i] 
            a[i] = panel.geta()
            b[i] = panel.getb()
            b_stiff[i] = panel.getbf()
            M_o[i] = p_total[i] * b[i] * a[i] * (a[i]/2) #in MN*m
            t_p[i] = panel.gettp()
            t_w[i] = panel.gettw()
            h_w[i] = panel.gethw()
            t_f[i] = panel.gettf()
            A_p[i] = t_p[i] * b[i]
            A_w[i] = t_w[i] * h_w[i]
            A_f[i] = t_f[i] * b_stiff[i]
            A[i] = A_p[i] + A_w[i] + A_f[i]
            M_p[i] = A_p[i] * (t_p[i]/2)
            M_w[i] = A_w[i] * (t_p[i] + (h_w[i]/2))
            M_f[i] = A_f[i] * (t_p[i] + h_w[i] + (t_f[i]/2))
            NA[i] = (M_p[i] + M_w[i] + M_f[i]) / (A[i])
            d_p[i] = (t_p[i]/2) - NA[i]
            d_w[i] = (t_p[i] + (h_w[i]/2)) - NA[i]
            d_f[i] = (t_p[i] + h_w[i] + (t_f[i]/2)) - NA[i]
            y_f[i] = -d_f[i]
            I_p_i[i] = ((1/12)*(b[i]*t_p[i]**3))
            I_p_NA[i] = I_p_i[i] + (A_p[i] * d_p[i]**2)
            I_w_i[i] = ((1/12)*(t_w[i]*h_w[i]**3))
            I_w_NA[i] = I_w_i[i] + (A_w[i] * d_w[i]**2)
            I_f_i[i] = ((1/12)*(b_stiff[i]*t_f[i]**3))
            I_f_NA[i] = I_f_i[i] + (A_f[i] * d_f[i]**2)
            I[i] = I_p_NA[i] + I_w_NA[i] + I_f_NA[i]
            del_o[i] = (5* (p_total[i]) * b[i] * (a[i]**4))/(384 * E * I[i]) #in m
            delta[i] = -a[i]/750 #from the text, positive if towards stiffeners
            sig_e[i] = (math.pi**2 * E * I[i]) / (A[i] * a[i]**2) #in MPa
            phi[i] = 1/(1-(sig_a[i]/sig_e[i]))
            sig_f[i] = sig_a[i] + ((M_o[i] * y_f[i]) / I[i]) + ((sig_a[i]*A[i]*(del_o[i] + delta[i])*y_f[i]*phi[i])/(I[i])) #in MPa
            if sig_f[i] < 0:
                sig_a_u_I[i] = 0.01 #this is to set the attained strength very low in case of negative sig_f
                #print ("sig_f_I less than zero")
            else:
                rho_NA[i] = (I[i]/A[i])**0.5
                lamb[i] = (a[i]/(math.pi*rho_NA[i]))*((sig_f[i]/E)**0.5)
                eta[i] = ((del_o[i]+delta[i])*y_f[i])/(rho_NA[i]**2)
                mu[i] = (M_o[i] * y_f[i]) / (I[i]*sig_f[i])
                zeta[i] = 1 - mu[i] + ((1+eta[i])/lamb[i]**2)
                R_I[i] = (zeta[i]/2) - (((zeta[i]**2)/4) - ((1-mu[i])/lamb[i]**2))**0.5
                sig_a_u_I[i] = sig_ys * R_I[i] #in MPa

            #Mode II Failure Calculations

            beta_II[i] = (b[i]/t_p[i]) * ((sig_yp/E)**0.5)
            zeta_II[i] = 1 + (2.75/(beta_II[i]**2))
            T[i] = 0.25 * (2 + zeta_II[i] - ((zeta_II[i]**2)-(10.4/(beta_II[i]**2)))**0.5)
            tau_II[i] = 0 #(p_total[i] * b[i] * (a[i]/2))/A[i] #formula from text in MPa
            sig_F_II[i] = ((T[i]-0.1)/T[i]) * sig_yp * ((1-(3*(tau_II[i]/sig_yp)**2)))**0.5
            if sig_F_II[i] < 0:
                sig_F_II[i] = 0.01 #this is to set the attained strength very low in case of negative sig_f
                #print ("sig_F_II less than zero")
            else:
                b_II[i] = T[i] * b[i]
                M_o_II[i] = p_total[i] * b_II[i] * a[i] * (a[i]/2) #in MN*m
                A_p_II[i] = t_p[i] * b_II[i]
                A_II[i] = A_p_II[i] + A_w[i] + A_f[i]
                M_p_II[i] = A_p_II[i] * (t_p[i]/2)
                NA_II[i] = (M_p_II[i] + M_w[i] + M_f[i]) / (A_II[i])
                d_p_II[i] = (t_p[i]/2) - NA_II[i]
                d_w_II[i] = (t_p[i] + (h_w[i]/2)) - NA_II[i]
                d_f_II[i] = (t_p[i] + h_w[i] + (t_f[i]/2)) - NA_II[i]
                y_p_II[i] = NA_II[i]
                I_p_i_II[i] = ((1/12)*(b_II[i]*t_p[i]**3))
                I_p_NA_II[i] = I_p_i_II[i] + (A_p_II[i] * d_p_II[i]**2)
                I_w_i_II[i] = ((1/12)*(t_w[i]*h_w[i]**3))
                I_w_NA_II[i] = I_w_i[i] + (A_w[i] * d_w_II[i]**2)
                I_f_i_II[i] = ((1/12)*(b_stiff[i]*t_f[i]**3))
                I_f_NA_II[i] = I_f_i[i] + (A_f[i] * d_f_II[i]**2)
                I_II[i] = I_p_NA_II[i] + I_w_NA_II[i] + I_f_NA_II[i]
                del_o_II[i] = (5* (p_total[i]) * b_II[i] * (a[i]**4))/(384 * E * I_II[i]) #in m
                rho_NA_II[i] = (I_II[i]/A_II[i])**0.5
                h_II[i] = panel.getNAStiff() - (t_p[i]/2)
                delta_P[i] = (h_II[i]*(A_f[i]+A_w[i])) * ((1/A_II[i])-(1/A[i]))
                y_NA_i_II[i] = panel.get_NA_base()
                H[i] = y_NA_i_II[i] - NA_y
                delta_H[i] = I[i]/(A[i]*H[i])
                delta_II[i] = (delta_P[i] + delta_H[i])
                lamb_II[i] = (a[i]/(math.pi*rho_NA_II[i]))*((sig_F_II[i]/E)**0.5)
                eta_II[i] = ((del_o_II[i]+delta_II[i])*y_p_II[i])/(rho_NA_II[i]**2)
                eta_p_II[i] = (delta_P[i]*y_p_II[i])/(rho_NA_II[i]**2)
                mu_II[i] = (M_o_II[i] * y_p_II[i]) / (I_II[i]*sig_F_II[i])
                zeta_R_II[i] = 1 - mu_II[i] + ((1+eta_II[i])/lamb_II[i]**2)
                R_II[i] = (zeta_R_II[i]/2) - (((zeta_R_II[i]**2)/4) - ((1-mu_II[i])/((1+eta_p_II[i])*lamb_II[i]**2)))**0.5
                sig_a_u_II[i] = sig_F_II[i] * R_II[i] * (A_II[i]/A[i]) #in MPa

            #Mode III failure calculations

            # Initialize M_o_g with a guess value
            """ M_o_g[i] = 10
            tolerance = 1e-6
            max_iterations = 100
            iteration = 0
            error = float('inf')

            while error > tolerance and iteration < max_iterations:
                y_f_II[i] = -d_f_II[i]
                mu_GH[i] = (M_o_g[i] * y_f_II[i]) / (I_II[i]*(-sig_ys))
                eta_p_GH[i] = (delta_P[i]*y_f_II[i])/(rho_NA_II[i]**2)
                eta_GH[i] = ((del_o_II[i]+delta_II[i])*y_f_II[i])/(rho_NA_II[i]**2)
                lamb_GH[i] = (a[i]/(math.pi*rho_NA_II[i]))*((sig_ys/E)**0.5)
                zeta_GH[i] = ((1-mu_GH[i])/(1+eta_p_GH[i]))-((1+eta_p_GH[i]+eta_GH[i])/((1+eta_p_GH[i])*lamb_GH[i]**2))
                R_GH[i] = (zeta_GH[i]/2) - (((zeta_GH[i]**2)/4) + ((1-mu_GH[i])/((1+eta_p_GH[i])*lamb_GH[i]**2)))**0.5
                sig_F_III = -sig_ys
                mu_III[i] = (M_o_g[i] * y_p_II[i]) / (I_II[i]*sig_F_III)
                zeta_III[i] = ((1-mu_III[i])/(1+eta_p_II[i]))-((1+eta_p_II[i]+eta_II[i])/((1+eta_p_II[i])*lamb_II[i]**2))
                R_III[i] = (zeta_III[i]/2) - (((zeta_III[i]**2)/4) - ((1-mu_III[i])/((1+eta_p_II[i])*lamb_II[i]**2)))**0.5

                # Calculate the error
                error = abs(R_GH[i] - R_III[i])

                # Update M_o_g using Newton's method
                if error > tolerance:
                    # Calculate the derivative of the error with respect to M_o_g
                    delta_M_o_g = 1e-6
                    M_o_g_temp = M_o_g[i] + delta_M_o_g
                    mu_GH_temp = (M_o_g_temp * y_f_II[i]) / (I_II[i]*(-sig_ys))
                    zeta_GH_temp = ((1-mu_GH_temp)/(1+eta_p_GH[i]))-((1+eta_p_GH[i]+eta_GH[i])/((1+eta_p_GH[i])*lamb_GH[i]**2))
                    R_GH_temp = (zeta_GH_temp / 2) - (((zeta_GH_temp**2) / 4) + ((1 - mu_GH_temp) / ((1 + eta_p_GH[i]) * lamb_GH[i]**2)))**0.5
                    derivative = (R_GH_temp - R_GH[i]) / delta_M_o_g

                    # Update M_o_g
                    M_o_g[i] -= (R_GH[i] - R_III[i]) / derivative

                iteration += 1
            
            sig_a_u_III[i] = sig_F_III * R_III[i] * (A_II[i]/A[i]) #in MPa """

            sig_a_u[i] = np.nanmin([sig_a_u_I[i], sig_a_u_II[i]])

            Rat_a_u[i] = sig_a_u[i]/((sig_ys+sig_yp)/2)

        return sig_a_u
    
    def HG_FOS (self, My, M_tot = 33052, mirror = True): #M_tot is the total moment in the midship section in kN*m
        '''
        Calculates the FOS against hull girder collapse
        '''
        FOS_HG = (My/M_tot) #to be replaced with Mult
        return FOS_HG

    def HG_reliability(self, My_nom, Ms_nom = 3006, Mw_r = 1, Mw_cov = 0.15, Mw_nom = 20031, Md_r = 1, Md_cov = 0.25, Md_nom = 10015, My_r = 1, My_cov = 0.15):
        '''
        Calculates the reliability of the midship section's hull girder strength
        
        Parameters
        ----------
        
        My_nom:     The nominal yield moment of the section
        Ms_nom:     The nominal moment from the SDI analysis
        Mw_r:       The ratio of the mean to nominal value for the wave moment
        Mw_cov:     The coefficient of variation for the wave moment
        Mw_nom:     The nominal wave moment
        Md_r:       The ratio of the mean to nominal value for the dynamic bending moment
        Md_cov:     The coefficient of variation for the dynamic bending moment
        Md_nom:     The nominal dynamic bending moment
        My_r:       The ratio of the mean to nominal value for the yield moment
        My_cov:     The coefficient of variation for the yield moment
        
        Returns
        --------
        beta_HG:    The reliability index of the midship section's hull girder strength
        '''
        self.limit_state = ra.LimitState(lambda Ms, Mw, Md, My: 1- ((Ms+Mw+Md)/My))
        
        #initialize stochastic model
        self.stochastic_model = ra.StochasticModel()
        
        #define random variables
        
        ## Ms is a constant value determined from SDI analysis in GHS
        self.stochastic_model.addVariable(ra.Constant("Ms", Ms_nom))
        
        ## Mw is a Gumbel distributed random variable with a mean/nominal ratio of 1 and a COV of 0.15, where nominal value is 27975 from LR rules
        self.stochastic_model.addVariable(ra.Gumbel("Mw", Mw_nom*Mw_r, Mw_cov*Mw_nom*Mw_r))
        
        ## Md is a Gumbel distributed random variable with a mean/nominal ratio of 1 and a COV of 0.25, where nominal value is 13903 from LR rules
        self.stochastic_model.addVariable(ra.Gumbel("Md", Md_nom*Md_r, Md_cov*Md_nom*Md_r))
        
        ## My is a Lognormal distributed random variable with a mean/nominal ratio of 1 and a COV of 0.15, where nominal value is provided by design optimization
        self.stochastic_model.addVariable(ra.Lognormal("My", My_nom*My_r, My_cov*My_nom*My_r))
    
        #initialize reliability analysis
        options = ra.AnalysisOptions()
        options.setPrintOutput(False)
        
        Analysis = ra.Form(
            analysis_options=options,
            stochastic_model=self.stochastic_model,
            limit_state=self.limit_state
        )
        
        Analysis.run()
        
        beta_HG = Analysis.getBeta()
        P_F_HG = Analysis.getFailure()
        
        return beta_HG, P_F_HG

    def Hughes_panel_FOS (self, sig_a_u, HG_stress):
        '''
        Calculates the FOS against panel compressive collapse based on panel w/ max ratio of HG_stress to ultimate compressive axial strength
        '''
        # Find the maximum ratio of HG_stress to sig_a_u
        ratios = np.zeros(len(HG_stress))
        for i in range(len(HG_stress)):
            ratios[i] = HG_stress[i] / sig_a_u[i]
        
        # Find the index of the maximum ratio
        max_ratio_index = np.argmax(ratios)
        
        # Use the maximum ratio to find the sig_au value and HG_stress value to use in the reliability analysis
        sig_a_u_Rel = sig_a_u[max_ratio_index]

        HG_stress_Rel = HG_stress[max_ratio_index]

        FOS_Hpanel = sig_a_u_Rel/HG_stress_Rel
        return FOS_Hpanel

    def Hughes_panel_reliability (self, sig_a_u, HG_stress, sig_a_u_Rel_r = 1.05, sig_a_u_Rel_cov = 0.075, HG_stress_Rel_r = 1, HG_stress_Rel_cov = 0.20):

        # Calculate the reliability of the panels to compressive collapse based on panel w/ max ratio of HG_stress to ultimate compressive axial strength
        # Results in one beta value for the entire section based on the most critical panel

        # Find the maximum ratio of HG_stress to sig_a_u
        ratios = np.zeros(len(HG_stress))
        for i in range(len(HG_stress)):
            ratios[i] = HG_stress[i] / sig_a_u[i]
        
        # Find the index of the maximum ratio
        max_ratio_index = np.argmax(ratios)
        
        # Use the maximum ratio to find the sig_au value and HG_stress value to use in the reliability analysis
        sig_a_u_Rel = sig_a_u[max_ratio_index]

        HG_stress_Rel = HG_stress[max_ratio_index]

        self.limit_state = ra.LimitState(lambda HG_stress_Rel, sig_a_u_Rel: 1 - (HG_stress_Rel/sig_a_u_Rel))
        
        #initialize stochastic model
        self.stochastic_model = ra.StochasticModel()
        
        #define random variables
        
        ## sig_au is a Lognormal distributed random variable with a mean/nominal ratio of 1 and a COV of 0.15, where nominal value is provided by design optimization
        self.stochastic_model.addVariable(ra.Lognormal("HG_stress_Rel", HG_stress_Rel*HG_stress_Rel_r, HG_stress_Rel_cov*HG_stress_Rel*HG_stress_Rel_r))

        ## sig_y is a Lognormal distributed random variable with a mean/nominal ratio of 1 and a COV of 0.15, where nominal value is provided by design optimization
        self.stochastic_model.addVariable(ra.Lognormal("sig_a_u_Rel", sig_a_u_Rel*sig_a_u_Rel_r, sig_a_u_Rel_cov*sig_a_u_Rel*sig_a_u_Rel_r))
        
        #initialize reliability analysis
        options = ra.AnalysisOptions()
        options.setPrintOutput(False)
        
        Analysis = ra.Form(
            analysis_options=options,
            stochastic_model=self.stochastic_model,
            limit_state=self.limit_state
        )
        
        Analysis.run()
        
        beta_Hpanel = Analysis.getBeta()
        P_F_Hpanel = Analysis.getFailure()
        
        return beta_Hpanel, P_F_Hpanel
    
    def plating_FOS (self, p_allow, p_design = 38.36): #pressure in ksi
        '''
        Calculates the FOS against bottom plating collapse
        '''
        FOS_plating = p_allow/p_design
        return FOS_plating

    def plating_reliability (self, p_allow, p_design_nom = 38.36, p_design_r = 1, p_design_cov = 0.25): #pressure in ksi
        
        self.limit_state = ra.LimitState(lambda p_d, p_a: 1 - (p_d/p_a))
        
        #initialize stochastic model
        self.stochastic_model = ra.StochasticModel()
        
        #define random variables
        
        ## p_d is a Weibull distributed random variable with a mean/nominal ratio of 1 and a COV of 0.25, where nominal value is 38.36 from LR rules
        self.stochastic_model.addVariable(ra.Weibull("p_d", p_design_nom*p_design_r, p_design_cov*p_design_nom*p_design_r))

        ## p_a is the allowable permanent set pressure
        self.stochastic_model.addVariable(ra.Constant("p_a", p_allow))
        
        #initialize reliability analysis
        options = ra.AnalysisOptions()
        options.setPrintOutput(False)
        
        Analysis = ra.Form(
            analysis_options=options,
            stochastic_model=self.stochastic_model,
            limit_state=self.limit_state
        )
        
        Analysis.run()
        
        beta_plating = Analysis.getBeta()
        P_F_plating = Analysis.getFailure()
        
        return beta_plating, P_F_plating

    def plot_section(self, show=False, mirror=True, ax_obj=None): #mirror is used to mirror the plot, does not automatically mirror calculations
        '''Plots the section'''
        
        section_plot = Section.section()
        for grill in self.grillages:
            panel = grill.getTTPanRef()
            section_plot.Append_Panels(panel)
        section_plot.Explode()
        section_plot.create_section_plot(show=False, mirror=mirror, ax_obj=None)
        plt.pyplot.xlabel('Longitudinal Position from Centerline (m)', fontsize=14)
        plt.pyplot.ylabel('Vertical Position from Keel (m)', fontsize=14)
        if show:
            plt.show()
            
    """ def calculate_fatigue_loading(self, average_cycles, average_moment):
        '''
        Generates a list of FatigueLoading objects based for the fatigue details associated
        with the section.
        
        Parameters
        ----------
        average_cycles:         The average annual fatigue cycles
        average_moment:         The average fatigue moment experienced by the section over
                                    its service life
        
        Returns
        --------
        No return value
        '''        
        self.fatigue_loaders = []
        for tpaneli in self.grillages:
            if tpaneli.getTTPanRef().type != 'HC':
                loads = []
                grill = tpaneli.getTTPanRef()
                #Determine spacing between panels in both directionss
                panel_yspacing = -grill.getb()*math.sin(math.radians(grill.ornt))
                panel_xspacing = grill.getb()*math.cos(math.radians(grill.ornt))
                
                for j in range(len(tpaneli.getFatigueDetails())):
                    #Must break panels down into single stiffened panels with small plates at either end
                    if j == 0:
                        #Get y location for first plate element
                        detail_y = grill.sloc[1] + 0.25*panel_yspacing + 0.5*grill.gettp()*math.cos(math.radians(grill.ornt))
                        area = 0.5 * grill.getb() * grill.gettp()
                    elif j == len(tpaneli.getFatigueDetails())-1:
                        #Get y location for final plate
                        entry_y = grill.sloc[1] + 0.25*panel_yspacing + 0.5*grill.gettp()*math.cos(math.radians(grill.ornt))
                        detail_y = entry_y + panel_yspacing * (grill.getnstiff() + 0.5)
                        area = 0.5 * grill.getb() * grill.gettp()
                    else:
                        #Get Y location
                        detail_y = grill.sloc[1] + (j+1)*panel_yspacing + grill.getNA()*math.cos(math.radians(grill.ornt))    
                        area = grill.getArea()
                    
                    MoI = grill.getINA() + area*(self._NAy - detail_y)**2
                    panel_stress = math.fabs(average_moment * (self._NAy - detail_y) / MoI)
                    Moment = panel_stress * area
                    stress_range = ((tpaneli.getTTPanRef().getmatlP().getE()/(10**9)) * Moment * math.fabs(detail_y - self._NAy)) / (self.EIT) / (10**6)  
                    loading = deterioratingStructure.Fatigue_Loading(stress_range, average_cycles)
                    loads.append(loading)
                    
            self.fatigue_loaders.append(loads)             

    def generate_fatigue_details(self, crack_cost):
        '''
        Adds fatigue details to each transverse TPanel based on the number of stiffeners.
        
        Parameters
        -----------
        crack_cost:             The cost to repair a fatigue crack
        
        Returns
        --------
        No return value
        '''
        
        #Create fatigue details and individual fatigue loading objects for each stiffened panel in each grillage - also add 
        #plates at either end of the grillage

        self.fatigue_loaders = []
        for tpaneli in self.grillages: 
            if tpaneli.getTTPanRef().type != 'HC':
                numDetails = tpaneli.getTTPanRef().getnstiff() + 2

                for j in range(numDetails):                    
                    fatigue_detail = deterioratingStructure(deterioratingStructure.repairCost(crack_cost), CrackDetailVL.CrackDetailVL(Nrepair=0), tpaneli._age)
                    tpaneli.add_fatigue_detail(fatigue_detail)
                  
    def age_structure(self, time, loading, fatigue_loading):
        '''
        Ages the structure 'time' years
        
        This method will age every one of the midship section grillages by the time step
        indicated by 'time' and the loading specified under the given loading object.
        
        Parameters
        ----------
        time:               The timestep to age the structure by
        loading:            A loading object
        fatigue_loading:    A list of lists of Fatigue Loading objects where item[i][j] corresponds
                                to the loading for the jth fatigue detail in the ith TPanel.
        
        Returns
        -------
        No return value 
        '''
        if time:
            self.section_age += time
            for i in range(len(self.grillages)):
                self.grillages[i].age(time,loading)
                details = self.grillages[i].getFatigueDetails()
                if details and fatigue_loading:
                    for j in range(len(details)):
                        details[j].age(time, fatigue_loading[i][j])
        
    def schedule_mainteance(self, timestep):
        
        '''
        Schedules a repair for a timestep in the sections life.
        
        This method will schedule a repair for each of the grillages defined in the 
        initial grillage list at a set timestep.
        
        Parameters
        ----------
        timestep:   timestep to complete a repair in
        
        Returns
        --------
        No return value
        '''
        self.schedule.append(timestep)
        
    
    def repair_needs(self):
        '''
        Determines which grillages are due for repair
        
        This method returns a list of T panels which need repair based on the
        thickness limits supplied upon the initialization of the midship section
        class.
        
        Parameters
        ----------
        No parameters
        
        Return Values
        -------------
        needs:  A list of indices corresponding to grillages which need repair
        '''
        needs = []        
        for i in range(len(self.grillages)):
            if self.grillages[i].needsRepair():
                needs.append(i)
        
        return needs
    
    def renew_grillages(self, grills, renew_fatigue_details=True):
        '''
        Renews a set of grillages in the cross section.
        
        This method will renew grillages corresponding to the list of indices supplied
        to the method.  This will also return the total cost to replace all the grillages.
        
        Parameters
        ----------
        
        grills:                     A list of indices corresponding to grillages to replace
        renew_fatigue_details:      A flag to indicate whether or not fatigue details should be
                                        renewed
        
        Returns
        --------
        total_cost: Cost to replace all specified grillages
        '''
        total_cost = 0.0

        for i in range(len(self.grillages)):
            if i in grills:
                panel = self.grillages[i]
                corrosion_model = corrode.paikCorrosion(self.coatLife, locationCodePlate=panel.getTTPanRef().qloc[0], locationCodeWeb=panel.getTTPanRef().qloc[1], \
                    locationCodeFlange=panel.getTTPanRef().qloc[2])
                if renew_fatigue_details:
                    fatigueModel = CrackDetail
                else:
                    fatigueModel=None
                costObject = panel.renew(corrosion_model, fatigueModel)
                total_cost += costObject.getRenewalCost()
        
        return total_cost
    
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
              
        collapse = SmithCollapse.SmithCollapseC('SmithCollapse', forcetol=200000.)
        #Discretize the elements
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
                AE = np.max(HCRelations._AE)
                
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
                        collapse.addElement(elx, ely, AE, True, strss, strn, E, yld)
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
                        collapse.addElement(elx, ely, AE, True, strss, strn, E, yld)
                        collapse.addElement(elx_2, ely_2, AE, True, strss, strn, E, yld)
                        
        if extra_outputs:
            return collapse, maxy, miny, total_moment, total_area
        else:
            return collapse
        
        
    def curvature_moment_curve(self, first_curvature, last_curvature, show=False, nC=100):
        '''
        Produces a curvature-moment curve for the cross section.
        
        This method will plot the curvature-moment curve as used in the Smith-type progressive
        collapse analysis that is used to get the ultimate moment for the vessel.  The curve is based
        on the current cross section definition.
        
        Parameters
        ----------
        first_curvature:    The initial curvature for the curve
        last_curvature:     The final curvature in the curve
        show:               Plot the curvature-moment curve and display figure (default=False)
        nC:                 The number of points to discretize the curvature into
        
        Returns
        -------
        A tuple containing (curvature, moment)
        '''
        
        #Initialize Smith-Type collapse analysis
        collapse = self.set_up_smith_collapse()
                    
        #@TODO ADD DEPTH ROUTINE
        NA_guess = 6
        curvatures = np.linspace(first_curvature, last_curvature, nC)
        moments = np.zeros(nC)
        
        for i in range(nC):
            result = collapse.applyVerticalCurvature(curvatures[i], NA_guess)
            moments[i] = result[0]
            NA_guess = result[1]
        
        if show:
            fig = plt.figure(1)
            ax = fig.add_subplot(111)
            ax.plot(curvatures, moments, 'k', linewidth=2)
            ax.set_xlabel('Curvature (1/m)', fontsize=16)
            ax.set_ylabel('Moment (N-m)', fontsize=16)
            for tl in ax.get_xticklabels():
                tl.set_y(-.02)
                tl.set_fontsize('14')
            for tl in ax.get_yticklabels():
                tl.set_fontsize('14')
            
            plt.show()
        
        return (curvatures, moments)

    def get_ultimate_moment(self, mtype='sagging'):
        '''
        Returns the ultimate moment for the midship section.
        
        This method uses a Smith-type progressive collapse analysis to calculate
        the ultimate moment for the cross section.  TThe transverse T-Panels used to
        define the midship section are exploded into single stiffened plates and then
        passed to an instance of the Smith collapse object.
        
        Parameters
        ----------
        mtype:  The type of moment, either hogging or sagging, to find the ultimate 
                    value of. (default='sagging')
        
        Returns
        -------
        abs(collapsed[0]):  The absolute value of the ultimate momemtn
        '''
            
        
        ### Initialize effective moments and areas to be used in initial neutral axis guess
        total_moment = 0.0
        total_area = 0.0
        miny = 1e20
        maxy = 0.0
        collapse, maxy, miny, total_moment, total_area = self.set_up_smith_collapse(extra_outputs=True)
                        
        ### Set up initial guess for neutral axis and minimum/maximum curvature values
        NA0 = self._NAy
        c1 = math.fabs(NA0 - maxy)
        c2 = math.fabs(NA0 - miny)
        c = np.maximum(c1, c2)
        Yldcrv = self.grillages[0].getTTPanRef().getmatlP().getYld()/(c*self.grillages[0].getTTPanRef().getmatlP().getE())
        
        if mtype.lower() == 'sagging' or mtype.lower() == 's':
            crvMax = Yldcrv * 20
            crvMin = Yldcrv / 20
        elif mtype.lower() == 'hogging' or mtype.lower() == 'h':
            crvMin = -Yldcrv * 20
            crvMax = -Yldcrv / 20
        else:
            sys.exit("argument 'mtype' takes values of 'sagging', 'hogging', 's', or 'h'.  Got value '"+str(mtype)+"' instead.")
        
        ### Set up smith collapse for vertical bending and carry out analysis
        collapse.setupVerticalBending()
        collapsed = collapse.findMaxVerticalMom(crvMin, crvMax, 1e-5, NA0)
        
        return math.fabs(collapsed[0])
    
    def predictRepair(self, span):
        '''
        Predicts which panels will need to be replaced in a given span of time
        
        Parameters
        -----------
        span:   Number of years to span prediction out to
        
        Returns
        --------
        predicted_repairs:  Indices corresponding to panels which will need
                                repairs over time span
        '''
        
        #Get current age to revert to after aging process
        current_section_age = self.section_age
        current_ages = np.zeros(len(self.grillages))
        for i in range(len(self.grillages)):
            current_ages[i] = self.grillages[i]._age
        
        #Age structure through time span and make prediction
        self.age_structure(span, None, None)
        predicted_repairs = self.repair_needs()
        
        #Renew structure and re-age to current age
        for i in range(len(self.grillages)):
            self.renew_grillages([i], renew_fatigue_details=False)
        
        for i in range(len(self.grillages)):
            self.grillages[i].age(current_ages[i], None)
        
        self.section_age = current_section_age
        
        return predicted_repairs
    
    def performScheduledMaintenance(self):
        '''
        Performs any scheduled maintenance for this year of the structure's life.  If
        there is maintenance then any grillages predicted to fail in the time span betweem
        now and the next repair will be repaired.
        
        Parameters
        ----------
        None
        
        Returns
        -------
        maintenance_cost:   Cost to replace panels
        repair_cost:        The cost of the repair method, either planned dry-docking cost or 0
        '''

        #Determine if maintenance is scheduled and if so, what the interval is
        if self.section_age in self.schedule:
            schedule_index = self.schedule.index(self.section_age)
            if schedule_index != len(self.schedule) - 1:
                next_schedule = self.schedule[schedule_index+1]
                span = next_schedule - self.section_age
                if self.section_age + span > (self.end + self.service_extension):
                    span = (self.end+self.service_extension) - self.section_age
            else:
                #if there are no successive maintenance cycles predict out to end of service life + extension
                span = (self.end + self.service_extension) - self.section_age

            #Predict repairs
            predicted_repairs = self.predictRepair(span)
            maintenance_cost = self.renew_grillages(predicted_repairs)
            repair_cost = self.planned_dry_dock
            
        else:
            maintenance_cost = 0.0
            repair_cost = 0.0 
            predicted_repairs = []
        
        return maintenance_cost, repair_cost, predicted_repairs
            
            
    def replaceCorroded(self):
        '''
        Replace corroded panels, return a repair cost and a repair type (either underway/drydock 
        costs or None)
        
        Parameters
        ----------
        No Parameters
        
        Returns
        -------
        repair_cost:    The cost of the replaced panels
        rtype:          The cost of the repair method or None
        current_needs:  List of replaced panel indices
        '''

        #Replace panels that need to be replaced due to corrosion
        repair_cost = 0.0
        current_needs = self.repair_needs()
        rtype = 0.0
        need_dry_dock = False
        if self.section_age > self.coatLife:
            if current_needs:
                rtype = self.underway
            for paneli in current_needs:
                for location in self.grillages[paneli].getTTPanRef().qloc:
                    if self.dry_dock_dict[location] == 1:
                        need_dry_dock = True
                        rtype = self.emergency_dry_dock
            repair_cost = self.renew_grillages(current_needs)
        
        return repair_cost, rtype, current_needs
    
    def fatigueCost(self):
        '''
        Gets the current fatigue cost for the section based on fatigued details defined
        from the generate_fatigue_detail method.
        
        Parameters
        -----------
        None
        
        Return
        ------
        fatigue_cost:   The current probabilistic fatigue cost for the section
        '''

        #Determine the fatigue cost
        fatigue_cost = 0.0
        for paneli in self.grillages:
            if paneli.getTTPanRef().type != 'HC':
                for detail in paneli.getFatigueDetails():
                    fatigue_cost += detail.getCrackProb() * detail.getRepairObject().getRenewalCost()       
        
        return fatigue_cost
    
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
            
    def SM_regenerateSection(self, fracSM, fracUCS, already_drydocked=False):
        '''
        This method will replace panels to maintain both a global section modulus and local buckling.  It will first 
        strategically replace panels to ensure that the section modulus of the structure is above a fraction of the 
        original.  Panels are chosen based on their overall increase to the section modulus and with an emphasis on
        not dry-docking the vessel if possible.  Panels are then replaced if the ratio of their ultimate compressive 
        strength to their yield stress is below a specified fraction of the original ratio.
        
        Parameters
        ----------
        fracSM:             The fraction of initial section modulus to compare to
        fracUCS:            The fraction of initial UCS/Ys to compare to
        already_drydocked:  A flag indicating if the vessel has been dry docked to perform previous
                                repairs in the current year. (default=False)
        
        Returns
        --------
        regeneration_cost:  The cost to replace all neccessary panels
        repair_method:      The cost of the repair method used; either underway maintenance or an
                                emergency dry docking.
        '''
        
        ###Replace panels to ensure the section modulus of the vessel is not below the set threshold
        regeneration_cost = 0.0
        repair_method = 0
        if self.section_age > self.coatLife:
            ###Replace panels if SM has fell below threshold
            if self.vessel_SM() < self.initial_SM * fracSM:
                checking_sm = True
                repair_method = self.underway
                kount = 0
                while checking_sm:
                    ranks = []
                    for i in range(len(self.grillages)):
                        paneli = self.grillages[i]
                        ###First try and only replace panels that will not incur an emergency dry docking
                        if kount == 0 and not already_drydocked:
                            added = False
                            for location in paneli.getTTPanRef().qloc:
                                if self.dry_dock_dict[location] == 0 and not added:
                                    added = True
                                    ranki = paneli._baseSM - (paneli.getTTPanRef().getINA() / paneli.getTTPanRef().gety_max())
                                    ranks.append( (i, math.fabs(ranki)) )
                                    repair_method = self.underway
                        else:
                            ###Second try all panels
                            ranki = paneli._baseSM - (paneli.getTTPanRef().getINA() / paneli.getTTPanRef().gety_max())
                            ranks.append((i, math.fabs(ranki)))
                            repair_method = self.emergency_dry_dock
                
                    ###Replace panels that will have greatest impact on section modulus
                    ranks = sorted(ranks, key=itemgetter(1), reverse=True)
            
                    for ranked in ranks:
                        if checking_sm:
                            regeneration_cost += self.renew_grillages([ranked[0]])
                            if self.vessel_SM() >= self.initial_SM * fracSM:
                                checking_sm = False
                
                    kount += 1
                    if kount > 5:
                        print ("Could not get ultimate strength to meet threshold after 5 iterations! Ending cycle.")
                        checking_sm = False
        
            ###Next replace panels due to local losses in buckling strength
            compression = PC.PaikCReg()
            for i in range(len(self.grillages)):
                fUCS = compression.ucs(self.grillages[i].getTTPanRef()) / self.grillages[i].getTTPanRef().getYsavg()
                if fUCS < fracUCS * self.initial_fUCS[i]:
                    regeneration_cost += self.renew_grillages([i])
                    for location in self.grillages[i].getTTPanRef().qloc:
                        if self.dry_dock_dict[location]:
                            repair_method = self.emergency_dry_dock
        
                
        return regeneration_cost, repair_method                    
                    
        
    def UM_regenerateSection(self, mmtLimit_sagging, mmtLimit_hogging, already_drydocked=False):
        '''
        This method will strategically replace panels until the ultimate moment in both hogging
        and sagging is higher than the limit supplied in the method.
        
        Parameters
        ----------
        mmtLimit_sagging:   Minimum allowable ultimate sagging moment
        mmtLimit_hogging:   Minimum allowable ultimate hogging moment
        already_drydocked:  A flag set to true if the method should NOT try to avoid
                                incurring an emergency dry dock cost. (default=False)
                                
        Returns
        -------
        regeneration_cost:  The cost to regenerate the section to the desired strength
        repair_method:        The cost of the repair type: either the underway repair cost, the
                                dry docking cost, or None if no repairs were needed.
        '''
        
                
        #Replace panels to ensure that the ultimate moment is larger than the moment limit
        regeneration_cost = 0.0
        repair_method = 0.0
        if self.section_age > self.coatLife:
            if self.get_ultimate_moment('s') < mmtLimit_sagging or self.get_ultimate_moment('h') < mmtLimit_hogging:
                checking_strength = True
                repair_method = self.underway
            else:
                checking_strength = False
                
            kount = 0.0
            while checking_strength:
                ranks = []
                for i in range(len(self.grillages)):
                    paneli = self.grillages[i]
                    #First try and replace panels that will not require emergency dry dockinng if no 
                    #emergency dry docking was needed to replace panels due to corrosion
                    if kount == 0 and not already_drydocked:
                        added = False
                        for location in paneli.getTTPanRef().qloc:
                            if self.dry_dock_dict[location] == 0 and not added:
                                added = True
                                #Rank panels based on greatest moment of inertia to be gained from replacement
                                ranki = (paneli._baseNA + paneli._baseArea * (self._NAy - paneli.getTTPanRef().sloc[1])**2) - (paneli.getTTPanRef().getNA() + \
                                                paneli.getTTPanRef().getArea() * (self._NAy - paneli.getTTPanRef().sloc[1])**2)
                                ranks.append((i, math.fabs(ranki)))
                                repair_method = self.underway
                    else:
                        ranki = (paneli._baseNA + paneli._baseArea * (self._NAy - paneli.getTTPanRef().sloc[1])**2) - (paneli.getTTPanRef().getNA() + \
                                        paneli.getTTPanRef().getArea() * (self._NAy - paneli.getTTPanRef().sloc[1])**2)
                        ranks.append((i, ranki))
                        repair_method = self.emergency_dry_dock
                
                #Sort ranked items by their ranking and replace.  Check after each replacement to see if
                #the ultimate moment is high enough.
                ranks = sorted(ranks, key=itemgetter(1), reverse=True)
            
                for ranked in ranks:
                    if checking_strength:
                        regeneration_cost += self.renew_grillages([ranked[0]])
                        if self.get_ultimate_moment('s') >= mmtLimit_sagging and self.get_ultimate_moment('h') >= mmtLimit_hogging:
                            checking_strength = False
                
                kount += 1
                if kount > 5:
                    print ("Could not get ultimate strength to meet threshold after 5 iterations! Ending cycle.")
                    checking_strength = False
                
        return regeneration_cost, repair_method """
        

    """ def maintenance_cost(self, end, refresh_age=0, mmtLimit_sagging=0, mmtLimit_hogging=0,
                            output=1,discount_rate=0., average_cycles=0, average_fatigue_moment=0,
                            service_extension=0, service_distribution=None,mult=1.,fracSM=0.8, fracUCS=0.8):
        '''
        Calculate the maintenance cost from the current age to the timestep specified in end.
        
        This method will age the structure from the current age to the age specified in end and estimate
        the total maintenance cost the vessel will incur in that time period.  This will include the cost
        of replacing panels due to corrosion and fatigue.  If there are maintenance items scheduled in any
        of the time steps this will also incur a cost.  
        
        Parameters
        ----------
        end:                    The timestep to end the calculation in
        refresh_age:            If this flag is set to true then the aging that occurs during
                                    this calculation will not be applied to the structure 
                                    permanently (default = 0)
        mmtLimit_sagging:       The smallest ultimate sagging moment the vessel is allowed to have before 
                                    repairs are undertaken to increase it.  If no limit is supplied
                                    it will be taken as proportional to the ultimate moment of the
                                    initial cross section.
        mmtLimit_hogging:       The smallest ultimate hogging moment the vessel is allowed to have before 
                                    repairs are undertaken to increase it.  If no limit is supplied
                                    it will be taken as proportional to the ultimate moment of the
                                    initial cross section.                
        interest_rate:          Estimated interest rate over the service life of the vessel (default=0.07)
        serice_extension:       The number of years passed the service life the vessel may
                                    have its operation extended.
        service_distribution:   The probability distribution representing the probability of
                                    extending the service life each year passed the initial
                                    estimate.
        mult:                   A multiplier for non-flat rate maintenance charges incured in 
                                    each year.  Can be used to estimate cost of full structure basd
                                    on number of transverse frames.  (default=1.)
        
        Returns
        -------
        total:                      The toal maintenance cost for the time period
        fatigue:                    The part of that cost due to fatigue damage
        corrosion:                  The part of that cost due to corrosion
        scheduled:                  The part of that cost due to schedule maintenance
        flat_rate:                  The part of that cost due to flat rate costs such as
                                        dry docking or underway repairs
        required_strength:          The part of that cost due to replacing panels to ensure
                                        that the ultimate strength of the vessel never falls
                                        below a required threshold
        yearly_fatigue:             An array of values representing the cost due to fatigue in
                                        each time step of the analysis
        yearly_corrosion:           An array of values representing the cost due to corrosion
                                        in each time step of the analysis
        yearly_scheduled:           An array of values representing the cost due to scheduled
                                        maintenance in each time step of the analysis
        yearly_flat_rates:          An array of values corresponding to the cost due to flat rate
                                        aspects of maintenance such as dry docking costs in each
                                        time step of the anaylsis
        average_cycles:             The average annual fatigue cycles
        average_fatigue_moment:     The average fatigue moment experienced by the section over
                                        its service life
        '''
        
        self.end = end
        self.service_distribution = service_distribution
        self.service_extension = service_extension
        
        #Some basic logistical input checks   
        if service_extension:
            if integrate.simps(service_distribution, np.linspace(0,service_extension,service_extension)) != 1:
                unity = False
                while not unity:
                    n = integrate.simps(service_distribution, np.linspace(0,service_extension,service_extension))
                    service_distribution /= n
                    
                    if integrate.simps(service_distribution, np.linspace(0,service_extension,service_extension)) != 1.0:
                        unity = True

        #ultimate_moment_sagging = self.get_ultimate_moment('s')
        #ultimate_moment_hogging = self.get_ultimate_moment('h')
        
        if type(end).__name__ != 'int':
            sys.exit("'end' argument must be integer value.  Got "+type(end).__name__+" value instead.")
        #if mmtLimit_sagging > ultimate_moment_sagging:
       #     print 'Warning! The supplied ultimate moment in sagging is larger than or equal to the \
       #             initial ultimate moment in sagging.  Using a proportional value instead!'
       #     mmtLimit_sagging = 0
            
       # if mmtLimit_hogging > ultimate_moment_hogging:
       #     print 'Warning! The supplied ultimate moment in hogging is larger than or equal to the \
       #             initial ultimate moment in sagging.  Using a proportional value instead!'
       #     mmtLimit_hogging = 0            
       
       ###Calculate time period
        
        time_period = (end + service_extension) - self.section_age        
        
        
        #If no minimum moment in sagging or hogging have been supplied take them as proportional
        #to the ultimate moment in the same direction
      #  if not mmtLimit_sagging:
      #      mmtLimit_sagging = ultimate_moment_sagging * 0.95
      #  if not mmtLimit_hogging:
      #      mmtLimit_hogging = ultimate_moment_hogging * 0.95
      #      
        if refresh_age:
            original_grillage = copy.deepcopy(self.grillages)
        
        #Initialize cost vectors
        yearly_fatigue = np.zeros(time_period)
        yearly_corrosion = np.zeros(time_period)
        yearly_scheduled = np.zeros(time_period)
        yearly_flat_rates = np.zeros(time_period)
        
        for year in range(int(self.section_age), end+service_extension):     
            if output:
                print ("Analyzing vessel in year "+str(year)+':')
                
            yearly_scheduled[year], yearly_flat_rates[year], replaced = self.performScheduledMaintenance()
            #Perform scheduled maintenance for this timestep
            if yearly_flat_rates[year]: 
                if output:
                    print ("Performing scheduled maintenance. Replaced panels: "+str(replaced))
            else:
                #If no scheduled maintenance is done then do any un-planned maintenance
                
                #Get the costs due to replacing individual panels that have corroded to an unacceptable degree
                yearly_corrosion[year], yearly_flat_rates[year], current_needs = self.replaceCorroded()
    
                if yearly_flat_rates[year]:
                    if output and yearly_flat_rates[year]==self.underway:
                        print ("Replaced panels "+str(current_needs)+". Performed repairs underway.")
                    if output and yearly_flat_rates[year]==self.emergency_dry_dock:
                        print ("Replaced panels "+str(current_needs)+". Needed to perform repairs at drydock.")
                
#                Get the costs due to replacing panels to ensure that the ultimate moment stays above a given threshold
#                regeneration_cost, repair_method = self.UM_regenerateSection(mmtLimit_sagging, mmtLimit_hogging, yearly_flat_rates[year]==self.emergency_dry_dock)
                regeneration_cost, repair_method = self.SM_regenerateSection(fracSM, fracUCS, yearly_flat_rates[year]==self.emergency_dry_dock)
                yearly_corrosion[year] += regeneration_cost
                if repair_method:
                    yearly_flat_rates[year] = repair_method
                    if output and yearly_flat_rates[year]==self.underway:
                        print ("Replaced panels to increase ultimate moment.  Repairs could be done underway")
                    if output and yearly_flat_rates[year]==self.emergency_dry_dock:
                        print ("Replaced panels to increase ultimate moment.  Repairs completed at dry dock")

            #Get the costs due to fatigue
            self.calculate_fatigue_loading(average_cycles, average_fatigue_moment)
            yearly_fatigue[year] = self.fatigueCost()
            if output:
                print ("Costs due to fatigue damage: "+str(yearly_fatigue[year]))
            
            #Apply multipler
            yearly_fatigue[year] *= (mult / ((1+discount_rate)**year)) 
            yearly_corrosion[year] *= (mult / ((1+discount_rate)**year)) 
            yearly_scheduled[year] *= (mult / ((1+discount_rate)**year)) 
            
            #Normalize costs by pdf of service life extension 
            if service_extension:
                if year > end:
                    yearly_fatigue[year] *= service_distribution[year-(end+1)]
                    yearly_corrosion[year] *= service_distribution[year-(end+1)]
                    yearly_scheduled[year] *= service_distribution[year-(end+1)]
                    yearly_flat_rates[year] *= service_distribution[year-(end+1)]
                    
                    
            #Age structure one year
            self.age_structure(1, None, self.fatigue_loaders)
            self.EIT, self._NAy, self.area, self.volume = self.section_data()
            self.EIT /= 10**9
        
        #Sum cost vectors
        fatigue = np.sum(yearly_fatigue)
        corrosion = np.sum(yearly_corrosion)
        scheduled = np.sum(yearly_scheduled)
        flat_rates = np.sum(yearly_flat_rates)
        total = fatigue + corrosion + scheduled + flat_rates

        if refresh_age:
            self.grillages = original_grillage
        
        
        return total, fatigue, corrosion, scheduled, flat_rates, yearly_fatigue, yearly_corrosion, yearly_scheduled, yearly_flat_rates """
