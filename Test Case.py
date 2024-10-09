
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


pmatl = Structures.EPMatl(95e6, 70e9, 0.34)
smatl = Structures.EPMatl(95e6, 70e9, 0.34)
tmatl = Structures.EPMatl(95e6, 70e9, 0.34)
# Dimensions for TPanel_trans
_B = 1
_L = 1
_nstiff = 4
_ntrans = 1
_tp = 0.004
_tw = 0.005
_hw = 0.1133
_tf = 0.0064
_bf = 0.05
_twh = 0.1133
_twt = 0.005
_tft = 0.0064
_tfb = 0.05
_sloc = [0,0,0]
_ornt = 0
_qloc = ['NA','NA','NA']
_eta = 0
# Creation of TPanel_trans
test_panel_1 = TPanel_trans.TPanel_trans(_B,_L,\
    _nstiff,_ntrans,_tp,_tw,_hw,_tf,\
    _bf,_twh,_twt,_tft,_tfb,_sloc,\
    _ornt,_qloc,pmatl,smatl,tmatl,_eta)
     
#get the section modulii
structure = midship_section.Midship_Section([test_panel_1], 0)
SM = structure.section_modulii()
        
#get section data
structure = midship_section.Midship_Section([test_panel_1], 0)
data = structure.section_data()
        
#get production cost
structure = midship_section.Midship_Section([test_panel_1],0)
PC = structure.production_cost()

#produce plot
structure = midship_section.Midship_Section([test_panel_1],0)
plot = structure.plot_section()

print (SM)
print (data)
print (PC)
plot

#set up the smith collapse curve
#structure = Midship_Section([self.test_panel_3],0)
#collapse_data = structure.set_up_smith_collapse(True)
    
#get the curvature and moment curve    
#structure = Midship_Section([self.test_panel_3],0)
#curve, moment = structure.curvature_moment_curve(-1e-3,1e-3,nC=4)
    
#get the ultimate moment
#ddg = Hull.Hull('testCases/DTMB51_structure.txt', warnings=False)
#structure = Midship_Section(ddg.get_grillage_list(), 0)
#sag, hog = structure.get_ultimate_moment('s'), structure.get_ultimate_moment('h')

  