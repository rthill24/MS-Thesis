# Smith ultimate strength method
# Author: Matthew Lankowski

################################################################################
#NOTES: 1)baseline (keel line) and centerline used as origin
#       2)suffix "i" = current iteration
################################################################################

import math
import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt

def MUltimate(XSection, output=1):
    '''Returns the curvature, moment, and neutral axis data for a cross-section'''    
    
    '''Define global variables, as required for certain functions'''
    global show_warnings
    show_warnings = output
    global Section    
    Section = XSection
    global NA0 
    
    '''Initialize effective area in each element (AE=A since strain=0)'''
    NAlocs = np.zeros(len(Section))
    for i in range(0,len(Section)):
        Section[i].AEi = max(Section[i]._AE)
        NAlocs[i] = Section[i].getYloc()         
    '''Establish neutral axis bounds for NA-location optimizer'''  
    NAmin = min(NAlocs)
    NAmax = max(NAlocs)
    
    '''Determine the curvature bounds for the ultimate strength optimizer'''
    NA0 = NAzeroStrain(Section)
    NAiTol = NA0*2/10000
    yloc_list = []
    for element in Section:
        yloc_list.append(element.getYloc())
    c1 = abs(NA0 - max(yloc_list))
    c2 = abs(NA0 - min(yloc_list))
    c = max(c1,c2)
    Ys = Section[0]._panel.getYsavg() #ATTN! currently just uses yield strength of first element
    E = Section[0]._panel._pmatl.getE()
    YldCrv = Ys/(c*E)
    
#    I=0
#    for el in Section:
#        el.MoI = el._panel.getINA() + el._panel.getArea()*(NA0 - el.getYloc())**2
#        I+=el.MoI
#    print 'My:', (Ys*I)/c

    crvMin = YldCrv/10
    crvMax = 2*YldCrv

#    crvStep = crvMax/50

    global NAi
    NAi = 0
    global getMomentCurrent
    global getMomentPrevious    
    getMomentCurrent = 0
    getMomentPrevious = 0
    
#    global showstats
#    showstats = 0
#################################################################################
#### THE FOLLOWING IS FOR DEBUGGING- PLOTS THE ENTIRE CURVE INSTEAD OF        ###
#### OPTIMIZING FOR MAX MOMENT                                                ###
################################################################################# 
#    crvtr_pos = [0]
#    posMom = [0]
#    ultcheck = 0
#    while ultcheck < 30:
#        maxMom = max(posMom)
#        crvtr_pos.append(crvtr_pos[-1]+crvStep)
#        posMom.append(-getMoment(crvtr_pos[-1],Section,'pos',NAmin,NAmax,NAiTol))
##        print crvtr_neg[-1], negMom[-1]
#        if posMom[-1] < maxMom:
#            ultcheck += 1    
#    
#    crvtr_neg = [0]
#    negMom = [0]
#    ultcheck = 0
#    while ultcheck < 30:
#        minMom = min(negMom)
#        crvtr_neg.append(crvtr_neg[-1]-crvStep)
#        negMom.append(getMoment(crvtr_neg[-1],Section,'neg',NAmin,NAmax,NAiTol))
##        print crvtr_neg[-1], negMom[-1]
#        if negMom[-1] > minMom:
#            ultcheck += 1
#    
#    posmoments = np.array(posMom)
#    negmoments = np.array(negMom)
#    print '\nSTEPPING:'
#    print 'max pos mom:',2*int(max(posMom))
#    print 'max neg mom:',2*int(min(negMom))
#    plt.plot(crvtr_neg, 2*negmoments, 'k')
#    plt.plot(crvtr_pos, 2*posmoments, 'k')
#    plt.xlabel('Curvature (1/m)')
#    plt.ylabel('Moment (N-m)')
#    plt.axhline(color = 'k')
#    plt.axvline(color = 'k')
#    plt.text(0.00025, 0.7e9, 'SAG')
#    plt.text(-0.00055, -1.8e9, 'HOG')
#    plt.grid()
#    plt.show()    
#    ult_mmt_neg = max(negMom)
#    print 'Crv at max pos mom:',crvtr_pos[posmoments.argmax()]
#    
#    showstats = 1
##    print 'HEY'
##    print posmoments
##    print 'argmax:',posmoments.argmax()
##    print 'crvtr_pos[argmax]',crvtr_pos[posmoments.argmax()]
#    print 'Max Pos Mom:',2*getMoment(crvtr_pos[posmoments.argmax()],Section,'pos',NAmin,NAmax,NAiTol)
##    print 'YO'
#    showstats = 0
#    
#################################################################################      
    
    '''Find the ultimate strength and NA for positive curvature'''
    find_pos = optimize.fminbound(getMoment, crvMin, crvMax, args=(Section,'pos',NAmin,NAmax,NAiTol), xtol=1e-09, disp=0)
    ult_mmt_pos = getMomentCurrent
#    ult_NA_pos = NAi
#    maxelmmts_pos = np.zeros(len(Section))    
#    for i in range(0, len(Section)):        
#        maxelmmts_pos[i] = Section[i].momenti #THESE ARE WRONG WHEN NOT USING OPTIMIZER (use last computed, not max)
    '''Find the ultimate strength and NA for negative curvature'''
    find_neg = optimize.fminbound(getMoment, -crvMax, -crvMin, args=(Section,'neg',NAmin,NAmax,NAiTol), xtol=1e-09, disp=0) 
    ult_mmt_neg = getMomentCurrent
#    print 'ult neg:', 2*ult_mmt_neg
#    ult_NA_neg = NAi
#    maxelmmts_neg = np.zeros(len(Section))    
#    for i in range(0, len(Section)):        
#        maxelmmts_neg[i] = Section[i].momenti   
    
    '''Store moment and NA data for ultimate curvature (either pos or neg)'''
    if abs(ult_mmt_pos) < abs(ult_mmt_neg):
        UltMmt = abs(ult_mmt_pos)
#        UltNA = ult_NA_pos
    else:
        UltMmt = abs(ult_mmt_neg)
#        UltNA = ult_NA_neg

#    '''store average pos/neg moments experienced by each element'''
#    nrml_max_el_mmts = (abs(maxelmmts_pos) + abs(maxelmmts_neg))/2

#    return UltMmt, UltNA, nrml_max_el_mmts #Use this and uncomment relevent lines above if moments on individual panels desired as output 
    return UltMmt


def getMoment(curvature, Section, crv_dir, NAmin, NAmax, NAiTol): 
    global showstats
    global NAi    
    global Curvature
    Curvature = curvature
    global getMomentCurrent

    if SumForces(NAmin)*SumForces(NAmax) >= 0:
        if show_warnings == 1:
            print('SumForces(a)*SumForces(b) >= 0; narrowing bracket...')
        NAmin = 0.25*(NAmin+NAmax)
        NAmax = 0.75*(NAmin+NAmax)
        if SumForces(NAmin)*SumForces(NAmax) >= 0:
            if show_warnings == 1:
                print('SumForces(a) * SumForces(b) still >= 0, using alternate NAi Finder...')     
            NAi = GoldenSectionSearch(NAmin, NAmax, NAiTol)
        else:
            NAi = optimize.brentq(SumForces, NAmin, NAmax, xtol=NAiTol)
    else:
        NAi = optimize.brentq(SumForces, NAmin, NAmax, xtol=NAiTol)

    getMomentCurrent = 0
    for i in range(0, len(Section)):
        Section[i].yi = Section[i].getYloc() - NAi
        if crv_dir == 'pos':
            Section[i].momenti = abs(Section[i].Fi * Section[i].yi) #get moment contribution from element
        else:
            Section[i].momenti = -abs(Section[i].Fi * Section[i].yi) #get moment contribution from element
        getMomentCurrent += Section[i].momenti 
    
    if crv_dir == 'pos':
        return -getMomentCurrent
    else:
        return getMomentCurrent

     
def SumForces(NAiGuess):
    '''Sums the forces above and below the neutral axis- used by BrentQ
    optimizer to find the instantaneous neutral axis'''
    curvaturei = Curvature
    SumF = 0
    for element in Section:
        straini = curvaturei * (element.getYloc() - NAiGuess)
        
#        global showstats
#        if showstats == 1:
#            stressel = element.getStress(straini)
#            print NAiGuess, element.getYloc(), straini, stressel 
#            
#            if (straini < 0 and stressel > 0) or (straini > 0 and stressel < 0):
#                print 'LKSJDFLSKDJF'

        #get effective area corresponding  to strain        
        element.AEi = element.getAE(straini)
#        if element.AEi == -1:
#            return -1
        #get stress corresponding to strain
        element.stressi = element.getStress(straini)
#        if element.stressi == -1:
#            return -1 
        #get force in element w.r.t. stress and location w.r.t. instant NA
        element.Fi = element.stressi * element.AEi 
        SumF += element.Fi
    return SumF
    
def NAzeroStrain(Section):
    '''returns the neutral axis of the cross-section in the zero-strain
    condition'''
    TotalMoment = 0
    TotalArea = 0
    for element in Section:
        TotalMoment += element.getYloc() * element.AEi #ttl mmnt about baseline
        TotalArea += element.AEi
    return TotalMoment / TotalArea

   
def GoldenSectionSearch(x1, x2, tolerance):
    
    tau = 0.3819660112501051
    x3 = tau*(x2 - x1) + x1
    x4 = tau*(x2 - x3) + x3
    f1 = abs(SumForces(x1))
    f2 = abs(SumForces(x2))
    f3 = abs(SumForces(x3))
    f4 = abs(SumForces(x4))
    
    x = np.array([x1,x2,x3,x4])
    y = np.array([f1,f2,f3,f4])
    x_closest = x[y.argmin()]    
    
    while abs(max(x) - min(x)) > tolerance:
        if f4 > f3:
            if x4 > x3:
                x2 = x4
                f2 = f4
                x4 = x3 - tau*(x3 - x1)
                f4 = abs(SumForces(x4))
            else:
                x1 = x4
                f1 = f4
                x4 = tau*(x2 - x3) + x3
                f4 = abs(SumForces(x4))
            
        if f4 < f3:
            if x4 > x3:
                x1 = x3
                f1 = f3
                x3 = x4
                f3 = f4
                x4 = tau*(x2 - x3) + x3
                f4 = abs(SumForces(x4))
            else:
                x2 = x3
                f2 = f3
                x3 = x4
                f3 = f4
                x4 = tau*(x3 - x1) + x1
                f4 = abs(SumForces(x4))
        
        x = np.array([x1,x2,x3,x4])
        y = np.array([f1,f2,f3,f4])
        x_closest = x[y.argmin()]
    
    return x_closest