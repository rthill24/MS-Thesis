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
import Box_Post_Processing_NAME

plt.close()

# Define the material properties
pmatl = Structures.EPMatl(95e6, 70e9, 0.34)
smatl = Structures.EPMatl(95e6, 70e9, 0.34)
tmatl = Structures.EPMatl(95e6, 70e9, 0.34)

# for bottom panel
B_bot  =  1
L_bot = 1
ntrans_bot = 1
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
eta_bot = 0.5

# for side panels
B_side = 1
L_side = 1
ntrans_side = 1
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
eta_side = 0.5

#for top panel
B_top = 1
L_top = 1
ntrans_top = 1
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
eta_top = 0.5

nstiff_top = Box_Post_Processing_NAME.optvars[0][0]
tp_top = Box_Post_Processing_NAME.optvars[0][1]
nstiff_bot = Box_Post_Processing_NAME.optvars[0][2]
tp_bot = Box_Post_Processing_NAME.optvars[0][3]
nstiff_side = Box_Post_Processing_NAME.optvars[0][4]
tp_side = Box_Post_Processing_NAME.optvars[0][5]

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
plot = structure.plot_section()
plt.show()

