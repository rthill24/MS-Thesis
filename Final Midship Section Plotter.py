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
import SD_Midship_Section_Post_Processing

plt.close()

# Define the material properties
pmatl = Structures.EPMatl(207, 71000, 2660, 0.33) #MPa, MPa, kg/m^3, Poisson's ratio
smatl = Structures.EPMatl(207, 71000, 2660, 0.33)
tmatl = smatl

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

nstiff_bot = SD_Midship_Section_Post_Processing.optvars[0][0]
tp_bot = SD_Midship_Section_Post_Processing.optvars[0][1]
nstiff_side = SD_Midship_Section_Post_Processing.optvars[0][2]
tp_side = SD_Midship_Section_Post_Processing.optvars[0][3]
nstiff_sheer = SD_Midship_Section_Post_Processing.optvars[0][4]
tp_sheer = SD_Midship_Section_Post_Processing.optvars[0][5]
nstiff_top = SD_Midship_Section_Post_Processing.optvars[0][6]
tp_top = SD_Midship_Section_Post_Processing.optvars[0][7]
nstiff_deck = SD_Midship_Section_Post_Processing.optvars[0][8]
tp_deck = SD_Midship_Section_Post_Processing.optvars[0][9]

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
structure = midship_section.Midship_Section([bottom_panel,side_panel,sheer_panel,top_panel,deck_panel],0)
plot = structure.plot_section()
plt.show()

