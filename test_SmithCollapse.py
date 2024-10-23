import matplotlib.pyplot as plt

import HansenC as HC
import HansenCdata as HCdata
import Element as Element
import SmithCollapse as SC
import numpy as np
import Structures
import TPanel_trans as TPT


print ("Generating Cross Section...")
XSection = []

#Test Panel
pmatl = Structures.EPMatl(95e6, 70e9, 0.34)
smatl = Structures.EPMatl(95e6, 70e9, 0.34)
tmatl = Structures.EPMatl(95e6, 70e9, 0.34)


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


test_panel_bot = TPT.TPanel_trans(B_bot,L_bot,nstiff_bot,ntrans_bot,tp_bot,\
                                        tw_bot,hw_bot,tf_bot,bf_bot,twh_bot,twt_bot,\
                                        tft_bot,tfb_bot,sloc_bot,ornt_bot,qloc_bot,pmatl,smatl,tmatl,eta_bot)


H = HC.HansenC(test_panel_bot)
strn = H._strn
strs = H._strss
AE = H._AE

for i in range(0,200):
    XSection.append(Element.Element(test_panel_bot, strn, strs, AE, yloc = i)) #create an element object with the strain, stress, AE and zloc data


print ("Running Smith Progressive Collapse Method...")
UltMmt = SC.MUltimate(XSection) #Run Smith Method
print ("Curvature/moment data generated, ultimate moment on cross-section =',UltMmt,'N-m'")
print (UltMmt)
