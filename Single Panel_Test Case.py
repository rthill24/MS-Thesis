
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

# Dimensions for TPanel_trans
B = 1
L = 1
nstiff = 4
ntrans = 1
tp = 0.004
tw = 0.005
hw = 0.1133
tf = 0.0064
bf = 0.05
twh = 0.1133
twt = 0.005
tft = 0.0064
tfb = 0.05
sloc = [0,0,0]
ornt = 0
qloc = ['NA','NA','NA']
eta = 0

# Creation of TPanel_trans
Panel = TPanel_trans.TPanel_trans(B, L, nstiff, ntrans, tp, tw, hw, tf, bf, twh, twt, tft, tfb, sloc, ornt, qloc, pmatl, smatl, tmatl, eta)

# Create the midship section    
structure = midship_section.Midship_Section([Panel],0)
        
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
'''
data = structure.section_data()

print ("here's EI: ", data[0])
print ("here's NAy: ", data[1])
print ("here's area: ", data[2])
print ("here's weight: ", data[3])
print ("here's I_NA: ", data[4])
print ("here's SM_min: ", data[5])
print ("here's My: ", data[6])

#produce plot
plot = structure.plot_section()
plt.show()