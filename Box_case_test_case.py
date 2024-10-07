import TPanel_trans
import copy
import math
import Structures
import logging
from numpy import zeros
import midship_section_condensed
import Plate
import Cost
import sys
import shutil
import os
import Section
import scipy
import sys, os, fileinput, math, string, numpy, scipy, random
import matplotlib.pyplot as plt

"""
    class for defining a basic grillage of fixed length and breadth consisiting of
    longitudinally stiffened panels supported by transverse T_section members
    
    COORDINATE SYSTEM
    =================
    Origin: Center of bottom edge
    X:  Positive right
    Y:  Positive up
    Z:  Positive in
    Angle:  East positive clockwise
    """

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

pmatl = Structures.EPMatl(95e6, 70e9, 0.34)
smatl = Structures.EPMatl(95e6, 70e9, 0.34)
tmatl = Structures.EPMatl(95e6, 70e9, 0.34)


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
test_panel_1 = TPanel_trans.TPanel_trans(B,L,nstiff,ntrans,tp,tw,hw,tf,bf,twh,twt,tft,tfb,sloc,ornt,qloc,pmatl,smatl,tmatl,eta)

structure = midship_section_condensed.Midship_Section([test_panel_1,test_panel_1],0)
#data = structure.section_data()
SM = structure.section_modulii()
prod_cost = structure.production_cost()
