
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

# Define the material properties
pmatl = Structures.EPMatl(207, 71000, 0.33)
smatl = Structures.EPMatl(207, 71000, 0.33)
tmatl = Structures.EPMatl(207, 71000, 0.33)

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
sloc_bot = [0,100,0]
ornt_bot = -174
qloc_bot = ['NA','NA','NA']
eta_bot = 0

# Creation of TPanel_trans
test_panel_bot = TPanel_trans.TPanel_trans(B_bot,L_bot,nstiff_bot,ntrans_bot,tp_bot,\
                                        tw_bot,hw_bot,tf_bot,bf_bot,twh_bot,twt_bot,\
                                        tft_bot,tfb_bot,sloc_bot,ornt_bot,qloc_bot,pmatl,smatl,tmatl,eta_bot)

# Create the midship section    
structure = midship_section.Midship_Section([test_panel_bot],0)
        
#get section data
'''
Returns
EI:     Total EI of the section
NAy:    Y location of the neutral axis
area:   Total cross-sectional area of the section
weight: weight of the structure
I_NA:   Moment of inertia about the neutral axis
SM_min: Minimum section modulus of the section
'''
data = structure.section_data()

#produce plot
plot = structure.plot_section()
plt.show()