#this module optimizes the midship section
#author: Richard Thill
#date: 10/29/2024

import sys
sys.path.insert(0, 'C:/Users/rthill/Documents/MS-Thesis')
import midship_section_RTO
import TPanel_trans
import copy
import math
import Structures
import logging
from numpy import zeros
import midship_section
import Plate
import Cost
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
import Allowable_Permanent_Set
import nsga2_michigan_threaded as nsga2
import Box_Checker

class SD_Midship_Section_Test_Case(nsga2.Problem):

    def __init__(self, numObj, numConstraints, GeneNum, loBound, upBound):
        super(SD_Midship_Section_Test_Case, self).__init__(numObj, numConstraints, GeneNum, loBound, upBound)
        self.current_individual = None

    def Eval(self, individual_instance, metamodel=None):
        # Define the material properties for 5086-H116 Aluminum
        self.pmatl = Structures.EPMatl(207, 71000, 2660, 0.33) #MPa, MPa, kg/m^3, Poisson's ratio
        self.smatl = Structures.EPMatl(207, 71000, 2660, 0.33)
        self.tmatl = Structures.EPMatl(207, 71000, 2660, 0.33)

        self.Frame_Spacing = 2250/1000

        # Dimensions for bottom panel
        self.B_bot = 4433.08999/1000
        self.nstiff_bot = 8
        self.ntrans_bot = 1
        self.tp_bot = 9/1000
        self.tw_bot = 8/1000
        self.hw_bot = 100/1000
        self.tf_bot = 6/1000
        self.bf_bot = 100/1000
        self.twh_bot = 100/1000
        self.twt_bot = 8/1000
        self.tft_bot = 6/1000
        self.tfb_bot = 100/1000
        self.sloc_bot = [0,0,0]
        self.ornt_bot = -28
        self.qloc_bot = ['NA','NA','NA']
        self.eta_bot = 0.5

        # Dimensions for side shell panel
        self.B_side = 3615.29782/1000
        self.nstiff_side = 6
        self.ntrans_side = 1
        self.tp_side = 5.5/1000
        self.tw_side = 8/1000
        self.hw_side = 100/1000
        self.tf_side = 6/1000
        self.bf_side = 50/1000
        self.twh_side = 100/1000
        self.twt_side = 8/1000
        self.tft_side = 6/1000
        self.tfb_side = 50/1000
        self.sloc_side = [self.B_bot*math.cos(math.radians(-self.ornt_bot)),self.B_bot*math.sin(math.radians(-self.ornt_bot)),0]
        self.ornt_side = -71
        self.qloc_side = ['NA','NA','NA']
        self.eta_side = 0.5

        # Dimensions for sheer strake panel
        self.B_sheer = 2556.31049/1000
        self.nstiff_sheer = 4
        self.ntrans_sheer = 1
        self.tp_sheer = 5/1000
        self.tw_sheer = 8/1000
        self.hw_sheer = 100/1000
        self.tf_sheer = 6/1000
        self.bf_sheer = 50/1000
        self.twh_sheer = 100/1000
        self.twt_sheer = 8/1000
        self.tft_sheer = 6/1000
        self.tfb_sheer = 50/1000
        self.sloc_sheer = [self.sloc_side[0]+(self.B_side*math.cos(math.radians(-self.ornt_side))),self.sloc_side[1]+(self.B_side*math.sin(math.radians(-self.ornt_side))),0]
        self.ornt_sheer = -174
        self.qloc_sheer = ['NA','NA','NA']
        self.eta_sheer = 0.5

        # Dimensions for top panel
        self.B_top = 2566.89797/1000
        self.nstiff_top = 6
        self.ntrans_top = 1
        self.tp_top = 5/1000
        self.tw_top = 8/1000
        self.hw_top = 100/1000
        self.tf_top = 6/1000
        self.bf_top = 50/1000
        self.twh_top = 100/1000
        self.twt_top = 8/1000
        self.tft_top = 6/1000
        self.tfb_top = 50/1000
        self.sloc_top = [0,self.sloc_sheer[1]+(self.B_sheer*math.sin(math.radians((self.ornt_sheer+180)))),0]
        self.ornt_top = 180
        self.qloc_top = ['NA','NA','NA']
        self.eta_top = 0

        # Dimensions for 2nd deck plating
        self.B_deck = 4165.69171/1000
        self.nstiff_deck = 9
        self.ntrans_deck = 1
        self.tp_deck = 5/1000
        self.tw_deck = 8/1000
        self.hw_deck = 100/1000
        self.tf_deck = 6/1000
        self.bf_deck = 50/1000
        self.twh_deck = 100/1000
        self.twt_deck = 8/1000
        self.tft_deck = 6/1000
        self.tfb_deck = 50/1000
        self.sloc_deck = [0,2763.65281/1000,0]
        self.ornt_deck = 180
        self.qloc_deck = ['NA','NA','NA']
        self.eta_deck = 0

        # to be pulled from chromosome
        chromosome = individual_instance.chromosome
        self.nstiff_bot = chromosome[0]
        self.tp_bot = chromosome[1]
        self.tw_bot = chromosome[2]
        self.hw_bot = chromosome[3]
        self.tf_bot = chromosome[4]
        self.bf_bot = chromosome[5]

        self.nstiff_side = chromosome[6]
        self.tp_side = chromosome[7]
        self.tw_side = chromosome[8]
        self.hw_side = chromosome[9]
        self.tf_side = chromosome[10]
        self.bf_side = chromosome[11]

        self.nstiff_sheer = chromosome[12]
        self.tp_sheer = chromosome[13]
        self.tw_sheer = chromosome[14]
        self.hw_sheer = chromosome[15]
        self.tf_sheer = chromosome[16]
        self.bf_sheer = chromosome[17]

        self.nstiff_top = chromosome[18]
        self.tp_top = chromosome[19]
        self.tw_top = chromosome[20]
        self.hw_top = chromosome[21]
        self.tf_top = chromosome[22]
        self.bf_top = chromosome[23]

        self.nstiff_deck = chromosome[24]
        self.tp_deck = chromosome[25]
        self.tw_deck = chromosome[26]
        self.hw_deck = chromosome[27]
        self.tf_deck = chromosome[28]
        self.bf_deck = chromosome[29]

        # Creation of midship section
        self.bottom_panel = TPanel_trans.TPanel_trans(self.B_bot,self.Frame_Spacing,self.nstiff_bot,self.ntrans_bot,self.tp_bot,\
                                                self.tw_bot,self.hw_bot,self.tf_bot,self.bf_bot,self.twh_bot,self.twt_bot,\
                                                self.tft_bot,self.tfb_bot,self.sloc_bot,self.ornt_bot,self.qloc_bot,self.pmatl,self.smatl,self.tmatl,self.eta_bot)
        self.side_panel = TPanel_trans.TPanel_trans(self.B_side,self.Frame_Spacing,self.nstiff_side,self.ntrans_side,self.tp_side,\
                                                self.tw_side,self.hw_side,self.tf_side,self.bf_side,self.twh_side,self.twt_side,\
                                                self.tft_side,self.tfb_side,self.sloc_side,self.ornt_side,self.qloc_side,self.pmatl,self.smatl,self.tmatl,self.eta_side)
        self.sheer_panel = TPanel_trans.TPanel_trans(self.B_sheer,self.Frame_Spacing,self.nstiff_sheer,self.ntrans_sheer,self.tp_sheer,\
                                                self.tw_sheer,self.hw_sheer,self.tf_sheer,self.bf_sheer,self.twh_sheer,self.twt_sheer,\
                                                self.tft_sheer,self.tfb_sheer,self.sloc_sheer,self.ornt_sheer,self.qloc_sheer,self.pmatl,self.smatl,self.tmatl,self.eta_sheer)
        self.top_panel = TPanel_trans.TPanel_trans(self.B_top,self.Frame_Spacing,self.nstiff_top,self.ntrans_top,self.tp_top,\
                                                self.tw_top,self.hw_top,self.tf_top,self.bf_top,self.twh_top,self.twt_top,\
                                                self.tft_top,self.tfb_top,self.sloc_top,self.ornt_top,self.qloc_top,self.pmatl,self.smatl,self.tmatl,self.eta_top)
        self.deck_panel = TPanel_trans.TPanel_trans(self.B_deck,self.Frame_Spacing,self.nstiff_deck,self.ntrans_deck,self.tp_deck,\
                                                self.tw_deck,self.hw_deck,self.tf_deck,self.bf_deck,self.twh_deck,self.twt_deck,\
                                                self.tft_deck,self.tfb_deck,self.sloc_deck,self.ornt_deck,self.qloc_deck,self.pmatl,self.smatl,self.tmatl,self.eta_deck)

        #objective functions called here for box 
        self.structure = midship_section_RTO.Midship_Section([self.bottom_panel, self.side_panel, self.sheer_panel, self.top_panel, self.deck_panel],0)
        self.current_individual = individual_instance
        PC = self.structure.production_cost()
        data = self.structure.section_data()
        weight = data[3]

        #return the output of the objective functions
        obj = [PC, weight] 
        return (obj, None)
    
    #define constraints
    def constraint(self, individual_instance, metamodel=None):
        if individual_instance != self.current_individual:
            self.Eval(individual_instance)
        
        #determine if geometry is valid for individual panels
        valid_bot = self.bottom_panel.geoValid()
        valid_side = self.side_panel.geoValid()
        valid_sheer = self.sheer_panel.geoValid()
        valid_top = self.top_panel.geoValid()
        valid_deck = self.deck_panel.geoValid()

        #get section data for constraints
        section_data = self.structure.section_data()

        #evaluate attained SM against LR requirement
        SM = section_data[5]
        SM_R = 0.267 #m^3, required section modulus per LR calculations
        frac_SM = (SM_R-SM)/SM_R

        #evaluate beta against hull girder collapse
        HG_beta = self.structure.HG_reliability(section_data[6])[0]
        HG_beta_R = 5.5 #required beta against hull girder collapse
        HG_beta_frac = (HG_beta_R-HG_beta)/HG_beta_R

        #evaluate beta against individual panel collapse
        HG_stress = self.structure.HG_stress()
        HG_caps = self.structure.Hughes_Panel(2.4, 1025, 38.36, HG_stress)
        panel_beta = self.structure.Hughes_panel_reliability(HG_caps, HG_stress)[0]
        panel_beta_R = 4 #required beta against panel collapse
        panel_beta_frac = (panel_beta_R-panel_beta)/panel_beta_R

        #evaluate beta against bottom panel collapse
        stiff_spacing = self.B_bot/(self.nstiff_bot+1)
        s_t = stiff_spacing/self.tp_bot
        if s_t <= 80:
            aps = (1/100)*stiff_spacing*1000
        else:
            aps = (1/75)*stiff_spacing*1000
        p_allow = Allowable_Permanent_Set.Allowable_Permanent_Set(0, aps)._p_aps(self.bottom_panel)
        plating_beta = self.structure.plating_reliability(p_allow)[0]
        plating_beta_R = 3.5 #required beta against plating collapse
        plating_beta_frac = (plating_beta_R-plating_beta)/plating_beta_R

        #iterate through constraints
        constraints_empty = [0, 0, 0, 0, 0, 0, 0, 0, 0]
        constraints = [0, 0, 0, 0, 0, 0, 0, 0, 0]

        if valid_bot == bool(0):
            constraints[0] = 1

        elif valid_side == bool(0):
            constraints[1] = 1

        elif valid_sheer == bool(0):
            constraints[2] = 1

        elif valid_top == bool(0):
            constraints[3] = 1

        elif valid_deck == bool(0):
            constraints[4] = 1

        elif frac_SM > 0:
            constraints[5] = frac_SM

        elif HG_beta_frac > 0:
            constraints[6] = HG_beta_frac

        elif panel_beta_frac > 0:
            constraints[7] = panel_beta_frac

        elif plating_beta_frac > 0:
            constraints[8] = plating_beta_frac

        else:
            constraints = constraints_empty
        return (constraints, None)

# run the optimization
## lower and upper bounds for design variables
loBound = [1, 5/1000, 1/1000, 1/1000, 1/1000, 1/1000, 1, 5/1000, 1/1000, 1/1000, 1/1000, 1/1000, 1, 5/1000, 1/1000, 1/1000, 1/1000, 1/1000, 1, 5/1000, 1/1000, 1/1000, 1/1000, 1/1000, 1, 5/1000, 1/1000, 1/1000, 1/1000, 1/1000]
upBound = [15, 100/1000, 50/1000, 50/1000, 50/1000, 50/1000, 15, 100/1000, 50/1000, 50/1000, 50/1000, 50/1000, 15, 100/1000, 50/1000, 50/1000, 50/1000, 50/1000, 15, 100/1000, 50/1000, 50/1000, 50/1000, 50/1000, 15, 100/1000, 50/1000, 50/1000, 50/1000, 50/1000]
test_problem = SD_Midship_Section_Test_Case(2, 9, 30, loBound, upBound)
    # numObj, numConstraints, GeneNum, loBound, upBound
opt = nsga2.Optimizer(test_problem)
opt.run("SD_Midship_Section_Optimizer_Output_Case_1", "initial_test", 48109, 500, 500)
# num generations, pop size