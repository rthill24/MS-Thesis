#this module generates a stiffened box for testing
#author: Richard Thill
#date: 10/29/2024

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

class Box_Test_Case(nsga2.Problem):

    def __init__(self, numObj, numConstraints, GeneNum, loBound, upBound):
        super(Box_Test_Case, self).__init__(numObj, numConstraints, GeneNum, loBound, upBound)
        self.current_individual = None

    def Eval(self, individual_instance, metamodel=None):
        # Define the material properties
        self.pmatl = Structures.EPMatl(95e6, 70e9, 0.34)
        self.smatl = Structures.EPMatl(95e6, 70e9, 0.34)
        self.tmatl = Structures.EPMatl(95e6, 70e9, 0.34)



        # for bottom panel
        self.B_bot  =  1
        self.L_bot = 1
        self.ntrans_bot = 1
        self.tw_bot = 0.005
        self.hw_bot = 0.1133
        self.tf_bot = 0.0064
        self.bf_bot = 0.05
        self.twh_bot = 0.1133
        self.twt_bot = 0.005
        self.tft_bot = 0.0064
        self.tfb_bot = 0.05
        self.sloc_bot = [0,0,0]
        self.ornt_bot = 180
        self.qloc_bot = ['NA','NA','NA']
        self.eta_bot = 0.5

        # for side panels
        self.B_side = 1
        self.L_side = 1
        self.ntrans_side = 1
        self.tw_side = 0.005
        self.hw_side = 0.1133
        self.tf_side = 0.0064
        self.bf_side = 0.05
        self.twh_side = 0.1133
        self.twt_side = 0.005
        self.tft_side = 0.0064
        self.tfb_side = 0.05
        self.sloc_side = [self.B_bot,0,0]
        self.ornt_side = -90
        self.qloc_side = ['NA','NA','NA']
        self.eta_side = 0.5

        #for top panel
        self.B_top = 1
        self.L_top = 1
        self.ntrans_top = 1
        self.tw_top = 0.005
        self.hw_top = 0.1133
        self.tf_top = 0.0064
        self.bf_top = 0.05
        self.twh_top = 0.1133
        self.twt_top = 0.005
        self.tft_top = 0.0064
        self.tfb_top = 0.05
        self.sloc_top = [0,self.B_side,0]
        self.ornt_top = 180
        self.qloc_top = ['NA','NA','NA']
        self.eta_top = 0.5

        # to be pulled from chromosome
        chromosome = individual_instance.chromosome
        self.nstiff_top = chromosome[0]
        self.tp_top = chromosome[1]
        self.nstiff_bot = chromosome[2]
        self.tp_bot = chromosome[3]
        self.nstiff_side = chromosome[4]
        self.tp_side = chromosome[5]

        # Create the stiffened box
        self.test_panel_bot = TPanel_trans.TPanel_trans(self.B_bot,self.L_bot,self.nstiff_bot,self.ntrans_bot,self.tp_bot,\
                                        self.tw_bot,self.hw_bot,self.tf_bot,self.bf_bot,self.twh_bot,self.twt_bot,\
                                        self.tft_bot,self.tfb_bot,self.sloc_bot,self.ornt_bot,self.qloc_bot,self.pmatl,self.smatl,self.tmatl,self.eta_bot)
        
        self.test_panel_side = TPanel_trans.TPanel_trans(self.B_side,self.L_side,self.nstiff_side,self.ntrans_side,self.tp_side,\
                                        self.tw_side,self.hw_side,self.tf_side,self.bf_side,self.twh_side,self.twt_side,\
                                        self.tft_side,self.tfb_side,self.sloc_side,self.ornt_side,self.qloc_side,self.pmatl,self.smatl,self.tmatl,self.eta_side)
        
        self.test_panel_top = TPanel_trans.TPanel_trans(self.B_top,self.L_top,self.nstiff_top,self.ntrans_top,self.tp_top,\
                                        self.tw_top,self.hw_top,self.tf_top,self.bf_top,self.twh_top,self.twt_top,\
                                        self.tft_top,self.tfb_top,self.sloc_top,self.ornt_top,self.qloc_top,self.pmatl,self.smatl,self.tmatl,self.eta_top)

        #objective functions called here for box 
        self.structure = midship_section.Midship_Section([self.test_panel_bot,self.test_panel_side,self.test_panel_top],0)
        self.current_individual = individual_instance
        PC = self.structure.production_cost()
        data = self.structure.section_data()

        #return the output of the objective functions
        obj = [PC, data[3]] 
        return (obj, None)
    
    #define constraints
    def constraint(self, individual_instance, metamodel=None):
        if individual_instance != self.current_individual:
            self.Eval(individual_instance)
        SM_all = self.structure.section_modulii()
        SM = SM_all[0]
        SM_R = 0.0000303043931 #section modulus of reference box from "Box Test Case.py"
        frac_SM = (SM_R-SM)/SM_R
        constraints = []
        if frac_SM > 0:
            constraints.append(frac_SM)
        else:
            constraints.append(0)
        return (constraints, None)
    
test_problem = Box_Test_Case(2, 1, 6, [1,0.001,1,0.001,1,0.001], [8,0.012,8,0.012,8,0.012])
opt = nsga2.Optimizer(test_problem)
opt.run("test_sql", "initial_test", 48109, 100, 100)