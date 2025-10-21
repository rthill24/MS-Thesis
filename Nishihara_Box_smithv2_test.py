
import TPanel_trans
import copy
import math
import Structures
import logging
from numpy import zeros
import midship_section_Nishihara
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
import smith_v2_Nishihara

# Define the material properties for the given box girder section
pmatl = Structures.EPMatl(245, 207000, 2660, 0.33) #MPa, MPa, kg/m^3, Poisson's ratio
smatl = Structures.EPMatl(245, 207000, 2660, 0.33)
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

Frame_Spacing = 900/1000

# Dimensions for bottom panel
B_bot = 720/1000 
nstiff_bot = 3
ntrans_bot = 1 
tp_bot = 3/1000 
tw_bot = 3/1000 
hw_bot = 50/1000 
tf_bot = 0/1000 
bf_bot = 0/1000 
twh_bot = 1/1000 
twt_bot = 1/1000 
tft_bot = 1/1000 
tfb_bot = 1/1000 
sloc_bot = [0,0,0] 
ornt_bot = 0
qloc_bot = ['NA','NA','NA']
eta_bot = 0

# Dimensions for right side panel
B_Rside = 720/1000
nstiff_Rside = 3
ntrans_Rside = 1 
tp_Rside = 3/1000
tw_Rside = 3/1000 
hw_Rside = 50/1000 
tf_Rside = 0/1000 
bf_Rside = 1/1000 
twh_Rside = 1/1000 
twt_Rside = 1/1000 
tft_Rside = 1/1000 
tfb_Rside = 1/1000 
sloc_Rside = [B_bot,0,0]
ornt_Rside = -90
qloc_Rside = ['NA','NA','NA']
eta_Rside = 0

# Dimensions for top panel
B_top = 720/1000
nstiff_top = 3
ntrans_top = 1
tp_top = 3/1000
tw_top = 3/1000
hw_top = 50/1000
tf_top = 0/1000
bf_top = 1/1000
twh_top = 1/1000
twt_top = 1/1000
tft_top = 1/1000
tfb_top = 1/1000
sloc_top = [B_bot,B_bot,0]
ornt_top = -180
qloc_top = ['NA','NA','NA']
eta_top = 0

# Dimensions for left side panel
B_Lside = 720/1000
nstiff_Lside = 3
ntrans_Lside = 1
tp_Lside = 3/1000
tw_Lside = 3/1000 
hw_Lside = 50/1000 
tf_Lside = 0/1000 
bf_Lside = 1/1000 
twh_Lside = 1/1000 
twt_Lside = 1/1000 
tft_Lside = 1/1000 
tfb_Lside = 1/1000 
sloc_Lside = [0,B_bot,0]
ornt_Lside = 90
qloc_Lside = ['NA','NA','NA']
eta_Lside = 0

# Creation of Nishihara box girder
bottom_panel = TPanel_trans.TPanel_trans(B_bot,Frame_Spacing,nstiff_bot,ntrans_bot,tp_bot,\
                                        tw_bot,hw_bot,tf_bot,bf_bot,twh_bot,twt_bot,\
                                        tft_bot,tfb_bot,sloc_bot,ornt_bot,qloc_bot,pmatl,smatl,tmatl,eta_bot)
Rside_panel = TPanel_trans.TPanel_trans(B_Rside,Frame_Spacing,nstiff_Rside,ntrans_Rside,tp_Rside,\
                                        tw_Rside,hw_Rside,tf_Rside,bf_Rside,twh_Rside,twt_Rside,\
                                        tft_Rside,tfb_Rside,sloc_Rside,ornt_Rside,qloc_Rside,pmatl,smatl,tmatl,eta_Rside)
top_panel = TPanel_trans.TPanel_trans(B_top,Frame_Spacing,nstiff_top,ntrans_top,tp_top,\
                                        tw_top,hw_top,tf_top,bf_top,twh_top,twt_top,\
                                        tft_top,tfb_top,sloc_top,ornt_top,qloc_top,pmatl,smatl,tmatl,eta_top)
Lside_panel = TPanel_trans.TPanel_trans(B_Lside,Frame_Spacing,nstiff_Lside,ntrans_Lside,tp_Lside,\
                                        tw_Lside,hw_Lside,tf_Lside,bf_Lside,twh_Lside,twt_Lside,\
                                        tft_Lside,tfb_Lside,sloc_Lside,ornt_Lside,qloc_Lside,pmatl,smatl,tmatl,eta_Lside)


# Create the box girder
structure = midship_section_Nishihara.Midship_Section([bottom_panel, Rside_panel, top_panel, Lside_panel],0)

# Produce plot of the midship section
plot = structure.plot_section()
plt.show()

# Get section data for the midship section
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

print ("here's NAy: ", data[1])
print ("here's area: ", data[2])
print ("here's I_NA: ", data[4])
print ("here's My: ", data[6])

# Run the Smith Method to get moment-curvature relationship and ultimate moment for the midship section
Smith = smith_v2_Nishihara.SmithMethod()
Smith.discretize([bottom_panel, Rside_panel, top_panel, Lside_panel])
data_smith = Smith.getOverallProperties()
print ("here's the area from smith method: ", data_smith[0])
print ("here's the NAy from smith method: ", data_smith[2])
print ("here's the Ixx from smith method: ", data_smith[3])

Smith.getCollapseCurve()
Smith.plotCollapseCurve()
Smith.getUltimateMoment()
print ("here's the ultimate moment from smith method: ", Smith.getUltimateMoment())
Myield, Myield_neg = Smith.getMomentatYldCrv()
print ("here's the positive moment at yield curvature: ", Myield)
print ("here's the negative moment at yield curvature: ", Myield_neg)




