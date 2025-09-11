
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
import Box_Checker
import smith_collapse
import smith_v2

# Define the material properties
pmatl = Structures.EPMatl(207, 71000, 2660, 0.33) #MPa, MPa, kg/m^3, Poisson's ratio
smatl = Structures.EPMatl(207, 71000, 2660, 0.33)
tmatl = Structures.EPMatl(207, 71000, 2660, 0.33)

"""
        constructor
            B - Breadth of grillage
            L - Length of Grillage
            nstiff - number of longitudinal stiffeners
            ntrans - number of transverse members
            tp - thickness of the plating
            tw - thickness of longitudinal webs
            hw - height of longitudinal webs
            tf - thickness of longitudinal flanges
            bf - breadth of longitudinal flanges
            pamtl - basic elastic perfectly plastic class defining plating material
            smatl - basic elastic perfectly plastic class defining stiffener material
            tmatl - basic elactic perfectly plastic class defining transverse member material characteristics
            twh - transverse web height
            tht - transverse web thickness
            tft - transverse flange thickness
            tfb - transverse flange breath
            sloc - spacial location of the start of the panel
            ornt - orientation of the stiffener
            qloc - qualitative location used for corroision model
        """

# Dimensions for TPanel_trans_bot
B_bot = 1
L_bot = 1
nstiff_bot = 1
ntrans_bot = 1
tp_bot = 0.004
tw_bot = 0.005
hw_bot = 0.1133
tf_bot = 0.0064
bf_bot = 0.05
twh_bot = 0.1133
twt_bot = 0.005
tft_bot = 0.0064
tfb_bot = 0.05
sloc_bot = [0,0,0]
ornt_bot = 0
qloc_bot = ['NA','NA','NA']
eta_bot = 0

# Dimensions for TPanel_trans_side
B_side  =  1
L_side = 1
nstiff_side = 1
ntrans_side = 1
tp_side = 0.004
tw_side = 0.005
hw_side = 0.1133
tf_side = 0.0064
bf_side = 0.05
twh_side = 0.1133
twt_side = 0.005
tft_side = 0.0064
tfb_side = 0.05
sloc_side = [B_bot,0,0]
ornt_side = -90
qloc_side = ['NA','NA','NA']
eta_side = 0

# Dimensions for TPanel_trans_top
B_top  =  1
L_top = 1
nstiff_top = 1
ntrans_top = 1
tp_top = 0.004
tw_top = 0.005
hw_top = 0.1133
tf_top = 0.0064
bf_top = 0.05
twh_top = 0.1133
twt_top = 0.005
tft_top = 0.0064
tfb_top = 0.05
sloc_top = [0,B_side,0]
ornt_top = 180
qloc_top = ['NA','NA','NA']
eta_top = 0 

# Creation of TPanel_trans
test_panel_bot = TPanel_trans.TPanel_trans(B_bot,L_bot,nstiff_bot,ntrans_bot,tp_bot,\
                                        tw_bot,hw_bot,tf_bot,bf_bot,twh_bot,twt_bot,\
                                        tft_bot,tfb_bot,sloc_bot,ornt_bot,qloc_bot,pmatl,smatl,tmatl,eta_bot)
test_panel_side = TPanel_trans.TPanel_trans(B_side,L_side,nstiff_side,ntrans_side,tp_side,\
                                        tw_side,hw_side,tf_side,bf_side,twh_side,twt_side,\
                                        tft_side,tfb_side,sloc_side,ornt_side,qloc_side,pmatl,smatl,tmatl,eta_side)
test_panel_top = TPanel_trans.TPanel_trans(B_top,L_top,nstiff_top,ntrans_top,tp_top,\
                                        tw_top,hw_top,tf_top,bf_top,twh_top,twt_top,\
                                        tft_top,tfb_top,sloc_top,ornt_top,qloc_top,pmatl,smatl,tmatl,eta_top)

# Check the box for corner stiffener intersections
check_my_box = Box_Checker.boxchecker(tp_bot,hw_top,tf_top,B_side,nstiff_side,tw_side,bf_side, test_panel_bot, test_panel_side, test_panel_top)
check_my_box.check_corner_stiffs()

# Create the box girder
structure = midship_section.Midship_Section([test_panel_bot, test_panel_side, test_panel_top],0)

# Check if the geometry is valid for each individual panel
valid_bot = test_panel_bot.geoValid()
valid_side = test_panel_side.geoValid()
valid_top = test_panel_top.geoValid()
print(valid_bot, valid_side, valid_top)

# Produce plot of the box girder section
plot = structure.plot_section()
plt.show()

# Get section data for the box girder
'''
Returns
    EI:     Total EI of the section
    NAy:    Y location of the neutral axis
    area:   Total cross-sectional area of the section
    weight: weight of the structure
    I_NA:   Moment of inertia about the neutral axis
    SM_min: Minimum section modulus of the section
    My:     Yield moment of the section
'''
data = structure.section_data()

print ("here's EI: ", data[0])
print ("here's NAy: ", data[1])
print ("here's area: ", data[2])
print ("here's weight: ", data[3])
print ("here's I_NA: ", data[4])
print ("here's SM_min: ", data[5])
print ("here's My: ", data[6])
print ("here's Mult: ", data[7]) #this is the old ultimate moment, now replaced by Smith method
        
# Get production cost
PC = structure.production_cost()
print ("here's PC: ", PC)

# Run the Smith Method to get moment-curvature relationship and ultimate moment for the midship section
Smith = smith_v2.SmithMethod()
Smith.discretize([test_panel_bot, test_panel_side, test_panel_top])
Smith.getUltimateMoment()

# Get pressure for allowable permanent set, value comes from LR - Rules and Regulations for the Classification of Special Service Craft, July 2022 - Part 3 General Requirements and Constructional Arrangements - Chapter 1 General Regulations - Section 8 Building tolerances and associated repairs, Table 1.8.6
stiff_spacing = B_bot/(nstiff_bot+1)
s_t = stiff_spacing/tp_bot
if s_t <= 80:
    aps = (1/100)*stiff_spacing*1000
else:
    aps = (1/75)*stiff_spacing*1000
press_bot = Allowable_Permanent_Set.Allowable_Permanent_Set(0, aps)._p_aps(test_panel_bot)
print ("here's the pressure': ", press_bot)

# Determine the HG stress for each panel and calculate the ultimate axial compressive strength using the Hughes methods
HG_stress = structure.HG_stress()
Hughes_panel = structure.Hughes_Panel(2.4, 1025, 38.36, HG_stress)
print ("here's HG_stress for each panel: ", HG_stress)
print ("here's Hughes_panel for each panel: ", Hughes_panel)

# Determine the hull girder reliability and proabability of failure
beta_HG = structure.HG_reliability(data[6])[0]
P_F_HG = structure.HG_reliability(data[6])[1]
print ("here's beta_HG: ", beta_HG)
print ("here's P_F_HG: ", P_F_HG)

# Determine the bottom panel reliability and probability of failure against design pressure and allowable permanent set
beta_plate = structure.plating_reliability(press_bot)[0]
P_F_plate = structure.plating_reliability(press_bot)[1]
print ("here's beta_plate: ", beta_plate)
print ("here's P_F_plate: ", P_F_plate)

# Determine the panel reliability and probability of failure against axial compression using the Hughes method
beta_Hpanel = structure.Hansen_panel_reliability(Hughes_panel, HG_stress)[0]
P_F_Hpanel = structure.Hansen_panel_reliability(Hughes_panel, HG_stress)[1]
print ("here's beta_Hpanel: ", beta_Hpanel)
print ("here's P_F_Hpanel: ", P_F_Hpanel)



