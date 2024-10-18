
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

# Define the material properties
pmatl = Structures.EPMatl(95e6, 70e9, 0.34)
smatl = Structures.EPMatl(95e6, 70e9, 0.34)
tmatl = Structures.EPMatl(95e6, 70e9, 0.34)

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
nstiff_bot = 8
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
nstiff_side = 8
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
nstiff_top = 8
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

def check_corner_stiffs(tp_bot, hw_top, tf_top, B_side, nstiff_side, tw_side, bf_side):
    """
    Check that the corner stiffeners are not intersecting
    """
    stiff_spacing = (B_side)/(nstiff_side-1)

    h_bot_side_stiff = tp_bot + stiff_spacing

    h_bot_side_stiff_flange = h_bot_side_stiff - (tw_side/2) - (bf_side/2)

    h_bot_stiff = tp_bot + hw_top + tf_top

    print (stiff_spacing, h_bot_side_stiff, h_bot_side_stiff_flange, h_bot_stiff) 

    if h_bot_stiff >= h_bot_side_stiff_flange:
        print ("Corner stiffeners are intersecting") 
    else:
        print ("Corner stiffeners are not intersecting")

check_corner_stiffs(tp_bot, hw_top, tf_top, B_side, nstiff_side, tw_side, bf_side)


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

# Create the midship section    
structure = midship_section.Midship_Section([test_panel_bot,test_panel_side,test_panel_top],0)

#get the section modulii
SM = structure.section_modulii()
        
#get section data
'''
Returns
EI:     Total EI of the section
NAy:    Y location of the neutral axis
area:   Total cross-sectional area of the section
volume: Volume of the structure
'''
data = structure.section_data()
        
#get production cost
PC = structure.production_cost()

#produce plot
plot = structure.plot_section()

print (SM)
print (data)
print (PC)
plt.show()

#check if the geometry is valid for each individual panel
valid_bot = test_panel_bot.geoValid()
valid_side = test_panel_side.geoValid()
valid_top = test_panel_top.geoValid()
print(valid_bot, valid_side, valid_top)

#evaluate ABS buckling constraints
""" ABS_constraints_bot = test_panel_bot.constraints()
ABS_constraints_side = test_panel_side.constraints()
ABS_constraints_top = test_panel_top.constraints()
print(ABS_constraints_bot, ABS_constraints_side, ABS_constraints_top) """

#set up the smith collapse curve
""" collapse_data = structure.set_up_smith_collapse(True) """
    
#get the curvature and moment curve    
#curve, moment = structure.curvature_moment_curve(-1e-3,1e-3,nC=4)
    
#get the ultimate moment
#sag, hog = structure.get_ultimate_moment('s'), structure.get_ultimate_moment('h')

  