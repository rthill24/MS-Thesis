################################################################################
# Panel data library for 24 panels from the following source:                  #
#                                                                              #
# S.E. Rutherford.  Stiffened compression panels: The analytical approach.     #
# Technical Report HSR 82/26/R2, Lloyd's Register of Shipping, 1984            #
################################################################################

import numpy as np
import Structures as ST

b = np.array([88.4, 147, 221, 236, 88.4, 147, 221, 236, 88.4, 147, 221, 236,
              88.4, 177, 265, 295, 88.4, 177, 265, 295, 88.4, 177, 265,
              295])/1000.
              
tp = np.array([3.07, 2.62, 2.54, 2.01, 3.07, 2.62, 2.54, 2.01, 3.07, 2.62, 2.54,
              2.01, 3.1, 3.05, 3.07, 2.57, 3.1, 3.05, 3.07, 2.57, 3.1, 3.05,
              3.07, 2.57])/1000.
              
hw = np.array([17.4, 30.4, 54.1, 43.6, 17.4, 30.4, 54.1, 43.6, 17.4, 30.4, 54.1,
               43.6, 26.4, 17.5, 34, 30.5, 26.4, 17.5, 34, 30.5, 26.4, 17.5, 34,
               30.5])/1000.

tw = np.array([4.88, 4.83, 4.9, 4.8, 4.88, 4.83, 4.9, 4.8, 4.88, 4.83, 4.9, 4.8,
               3.1, 4.85, 4.95, 4.9, 3.1, 4.85, 4.95, 4.9, 3.1, 4.85, 4.95,
               4.9])/1000.
               
bf = np.array([12.7, 12.7, 12.7, 12.7, 12.7, 12.7, 12.7, 12.7, 12.7, 12.7, 12.7,
               12.7, 0, 12.7, 12.7, 12.7, 0, 12.7, 12.7, 12.7, 0, 12.7, 12.7,
               12.7])/1000.

tf = np.array([6.17, 6.22, 6.1, 6.25, 6.17, 6.22, 6.1, 6.25, 6.17, 6.22, 6.1,
               6.25, 0, 6.15, 6.2, 6.12, 0, 6.15, 6.2, 6.12, 0, 6.15, 6.2,
               6.12])/1000.
               
L = np.array([244, 384, 638, 523, 488, 767, 1275, 1046, 732, 1151, 1913, 1570,
              262, 244, 422, 384, 523, 488, 843, 767, 785, 732, 1265,
              1151])/1000.

syp = np.array([250, 250, 256, 221, 225, 239, 270, 247, 230, 239, 239, 249, 253,
                242, 227, 244, 229, 229, 253, 261, 258, 242, 244, 239])*10.**6
                
sys = np.array([283, 262, 247, 250, 259, 259, 246, 259, 283, 258, 252, 266, 261,
                269, 267, 273, 256, 246, 266, 247, 262, 262, 262, 267])*10.**6

E = np.array([190, 190, 190, 190, 190, 190, 190, 190, 190, 190, 190, 190, 190,
              190, 190, 190, 190, 190, 190, 190, 190, 190, 190, 190])*10.**9
        
y0 = np.array([0, 0, 0, 0, 0.08, 0.76, 1.47, 2.73, 0, 0, 0, 0, 0, 0, 0, 0,
                0.10, 0.62, 2.0, 2.14, 0, 0, 0, 1.46])/1000.

w0 = np.array([0.05, 0.04, 0.21, 0.16, 1.20, 0.78, 1.08, 0.49, 1.78, 1.53, 3.37,
               1.13, 0.04, 0.07, 0.15, 0.07, 0.37, 0.19, 0.35, 0.24, 0.88, .40,
               0.81, 0.52])/1000.

e = np.array([1.9, -0.9, -2.6, -2.7, 1.9, -0.9, -2.6, -2.7, 1.9, -0.9, -2.6,
              -2.7, 0.11, -0.4, -1.5, -2.1, 0.11, -0.4, -1.5, -2.1, 0.11, -0.4,
              -1.5, -2.1])/1000.

phi = np.array([0.976, 0.733, 0.713, 0.567, 0.824, 0.75, 0.621, 0.515, 0.716,
                0.660, 0.494, 0.448, 0.988, 0.764, 0.569, 0.506, 0.822, 0.656,
                0.563, 0.455, 0.696, 0.515, 0.491, 0.384])

def getb(i):
    return float(b[i])

def gettp(i):
    return float(tp[i])

def gethw(i):
    return float(hw[i])

def gettw(i):
    return float(tw[i])
    
def getbf(i):
    return float(bf[i])
    
def gettf(i):
    return float(tf[i])
    
def getL(i):
    return float(L[i])
    
def getsyp(i):
    return float(syp[i])
    
def getsys(i):
    return float(sys[i])
    
def getE(i):
    return float(E[i])
    
def gety0(i):
    return float(y0[i])
    
def getw0(i):
    return float(w0[i])
    
def gete(i):
    return float(e[i])
    
def getphi(i):
    return float(phi[i])    

def HansenPanel(i):
    i = i - 1
    syp = getsyp(i)
    E = getE(i)
    sys = getsys(i)
    b = getb(i)
    tp = gettp(i)
    tw = gettw(i)
    hw = gethw(i)
    tf = gettf(i)
    bf = getbf(i)
    L = getL(i)
    w0 = getw0(i)
    e = gete(i)
    q = 0   
    pmatl = ST.EPMatl(syp, E, 0.3)
    smatl = ST.EPMatl(sys, E, 0.3)
    return ST.TPanel(b, tp, tw, hw, tf, bf, L, pmatl, smatl), w0, e, q
  
##validity check
#variables = [b, tp, hw, tw, bf, tf, L, syp, sys, E, y0, w0, e, phi]
#chars = ['b', 'tp', 'hw', 'tw', 'bf', 'tf', 'L', 'syp', 'sys', 'E', 'y0', 'w0',
#         'e', 'phi']
#count = 0
#for i in variables:    
#    char =  chars[count]    
#    print 'length of',char,'is',len(i),'(should be 24)'
#    print 'min/max of',char,'is',min(i),'/',max(i)
#    
#    print ''    
#    count += 1