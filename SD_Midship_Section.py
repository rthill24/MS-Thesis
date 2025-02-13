
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

# Define the material properties for 5086-H116 Aluminum
pmatl = Structures.EPMatl(207, 71000, 2660, 0.33) #MPa, MPa, kg/m^3, Poisson's ratio
smatl = Structures.EPMatl(207, 71000, 2660, 0.33)
tmatl = smatl

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
            twt - transverse web thickness
            tft - transverse flange thickness
            tfb - transverse flange breath
            sloc - spacial location of the start of the panel
            ornt - orientation of the stiffener
            qloc - qualitative location used for corroision model
        """

Frame_Spacing = 2250/1000

# Dimensions for bottom panel
B_bot = 4433.08999/1000 
nstiff_bot = 8 
ntrans_bot = 1 
tp_bot = 9/1000 
tw_bot = 8/1000 
hw_bot = 100/1000 
tf_bot = 6/1000 
bf_bot = 100/1000 
twh_bot = 100/1000 
twt_bot = 8/1000 
tft_bot = 6/1000 
tfb_bot = 100/1000 
sloc_bot = [0,0,0] 
ornt_bot = -28
qloc_bot = ['NA','NA','NA']
eta_bot = 0

# Dimensions for side shell panel
B_side = 3615.29782/1000
nstiff_side = 6 
ntrans_side = 1 
tp_side = 5.5/1000
tw_side = 8/1000 
hw_side = 100/1000 
tf_side = 6/1000 
bf_side = 50/1000 
twh_side = 100/1000 
twt_side = 8/1000 
tft_side = 6/1000 
tfb_side = 50/1000 
sloc_side = [B_bot*math.cos(math.radians(-ornt_bot)),B_bot*math.sin(math.radians(-ornt_bot)),0]
ornt_side = -71
qloc_side = ['NA','NA','NA']
eta_side = 0

# Dimensions for sheer strake panel
B_sheer = 2556.31049/1000
nstiff_sheer = 4
ntrans_sheer = 1
tp_sheer = 5/1000
tw_sheer = 8/1000
hw_sheer = 100/1000
tf_sheer = 6/1000
bf_sheer = 50/1000
twh_sheer = 100/1000
twt_sheer = 8/1000
tft_sheer = 6/1000
tfb_sheer = 50/1000
sloc_sheer = [sloc_side[0]+(B_side*math.cos(math.radians(-ornt_side))),sloc_side[1]+(B_side*math.sin(math.radians(-ornt_side))),0]
ornt_sheer = -174
qloc_sheer = ['NA','NA','NA']
eta_sheer = 0

# Dimensions for top panel
B_top = 2566.89797/1000
nstiff_top = 6
ntrans_top = 1
tp_top = 5/1000
tw_top = 8/1000
hw_top = 100/1000
tf_top = 6/1000
bf_top = 50/1000
twh_top = 100/1000
twt_top = 8/1000
tft_top = 6/1000
tfb_top = 50/1000
sloc_top = [0,sloc_sheer[1]+(B_sheer*math.sin(math.radians((ornt_sheer+180)))),0]
ornt_top = 180
qloc_top = ['NA','NA','NA']
eta_top = 0

# Dimensions for 2nd deck plating
B_deck = 4165.69171/1000
nstiff_deck = 9
ntrans_deck = 1
tp_deck = 5/1000
tw_deck = 8/1000
hw_deck = 100/1000
tf_deck = 6/1000
bf_deck = 50/1000
twh_deck = 100/1000
twt_deck = 8/1000
tft_deck = 6/1000
tfb_deck = 50/1000
sloc_deck = [0,2763.65281/1000,0]
ornt_deck = 180
qloc_deck = ['NA','NA','NA']
eta_deck = 0

# Creation of midship section
bottom_panel = TPanel_trans.TPanel_trans(B_bot,Frame_Spacing,nstiff_bot,ntrans_bot,tp_bot,\
                                        tw_bot,hw_bot,tf_bot,bf_bot,twh_bot,twt_bot,\
                                        tft_bot,tfb_bot,sloc_bot,ornt_bot,qloc_bot,pmatl,smatl,tmatl,eta_bot)
side_panel = TPanel_trans.TPanel_trans(B_side,Frame_Spacing,nstiff_side,ntrans_side,tp_side,\
                                        tw_side,hw_side,tf_side,bf_side,twh_side,twt_side,\
                                        tft_side,tfb_side,sloc_side,ornt_side,qloc_side,pmatl,smatl,tmatl,eta_side)
sheer_panel = TPanel_trans.TPanel_trans(B_sheer,Frame_Spacing,nstiff_sheer,ntrans_sheer,tp_sheer,\
                                        tw_sheer,hw_sheer,tf_sheer,bf_sheer,twh_sheer,twt_sheer,\
                                        tft_sheer,tfb_sheer,sloc_sheer,ornt_sheer,qloc_sheer,pmatl,smatl,tmatl,eta_sheer)
top_panel = TPanel_trans.TPanel_trans(B_top,Frame_Spacing,nstiff_top,ntrans_top,tp_top,\
                                        tw_top,hw_top,tf_top,bf_top,twh_top,twt_top,\
                                        tft_top,tfb_top,sloc_top,ornt_top,qloc_top,pmatl,smatl,tmatl,eta_top)
deck_panel = TPanel_trans.TPanel_trans(B_deck,Frame_Spacing,nstiff_deck,ntrans_deck,tp_deck,\
                                        tw_deck,hw_deck,tf_deck,bf_deck,twh_deck,twt_deck,\
                                        tft_deck,tfb_deck,sloc_deck,ornt_deck,qloc_deck,pmatl,smatl,tmatl,eta_deck)

# Create the midship section    
structure = midship_section.Midship_Section([bottom_panel, side_panel, sheer_panel, top_panel, deck_panel],0)
  
#get section data
'''
Returns
    EI:     Total EI of the section
    NAy:    Y location of the neutral axis
    area:   Total cross-sectional area of the section
    weight: weight of the structure
    I_NA:   Moment of inertia about the neutral axis
    SM_min: Minimum section modulus of the section
    My:     Yield moment of the section
    Mult:   Ultimate moment of the section based on specified yield strength
'''

data = structure.section_data()
HG_stress = structure.HG_stress()
Hughes_panel = structure.Hughes_Panel(2.4, 1025, 38.36, HG_stress)

beta_HG = structure.HG_reliability(data[6])[0]
P_F_HG = structure.HG_reliability(data[6])[1]
        
#get production cost
PC = structure.production_cost()

#produce plot
plot = structure.plot_section()
plt.show()

print ("here's the R: ", Hughes_panel)

#output data
print ("here's EI: ", data[0])
print ("here's NAy: ", data[1])
print ("here's area: ", data[2])
print ("here's weight: ", data[3])
print ("here's I_NA: ", data[4])
print ("here's SM_min: ", data[5])
print ("here's My: ", data[6])
print ("here's Mult: ", data[7])
print ("here's PC: ", PC)
print ("here's beta_HG: ", beta_HG)
print ("here's P_F_HG: ", P_F_HG)
print ("here's HG_stress for each panel: ", HG_stress)


#check if the geometry is valid for each individual panel
valid_bot = bottom_panel.geoValid()
valid_side = side_panel.geoValid()
valid_sheer = sheer_panel.geoValid()
valid_top = top_panel.geoValid()
valid_deck = deck_panel.geoValid()
print(valid_bot, valid_side, valid_sheer, valid_top, valid_deck)

#get pressure for allowable permanent set, value comes from LR - Rules and Regulations for the Classification of Special Service Craft, July 2022 - Part 3 General Requirements and Constructional Arrangements - Chapter 1 General Regulations - Section 8 Building tolerances and associated repairs, Table 1.8.6
stiff_spacing = B_bot/(nstiff_bot+1)
s_t = stiff_spacing/tp_bot
if s_t <= 80:
    aps = (1/100)*stiff_spacing*1000
else:
    aps = (1/75)*stiff_spacing*1000
press_bot = Allowable_Permanent_Set.Allowable_Permanent_Set(0, aps)._p_aps(bottom_panel)
print ("here's the pressure': ", press_bot)

#get plating reliability
beta_plate = structure.plating_reliability(press_bot)[0]
P_F_plate = structure.plating_reliability(press_bot)[1]
print ("here's beta_plate: ", beta_plate)
print ("here's P_F_plate: ", P_F_plate)
