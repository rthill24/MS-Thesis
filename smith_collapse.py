# -*- coding: utf-8 -*-
#Simplified implementation of the Smith collapse method for ultimate strength
# (c) 2025 Regents of the Univesity of Michigan

import TPanel_trans
import Element
import numpy as np
import HansenC as HC
import math
import matplotlib.pyplot as plt


class SmithCollapse(object):
    '''
    Class to implement the Smith collapse method for ultimate strength
    '''

    def __init__(self, t_panels):
        self.t_panels = t_panels
        
    def Hansen(self):
        '''get the Hansen data for each t_panel in the midship section'''
        strn_data = []
        strs_data = []
        AE_data = []
        HC_data = []

        for panel in self.t_panels:
            data_i = HC.HansenC(panel)
            HC_data.append(data_i)
            strn_data.append(-data_i._strn) #negative strain/stress for compression
            strs_data.append(-data_i._strss)
            AE_data.append(data_i._AE)
           
        return HC_data, strn_data, strs_data, AE_data
    
    def discretize(self):
        '''Discretize the cross section into elements (individual t_panels, then elements)'''
        self.HC_data = self.Hansen()[0]
        self.strn_data = self.Hansen()[1]
        self.strs_data = self.Hansen()[2]
        self.AE_data = self.Hansen()[3]
        elements = []
        XSection = []
        count = -1
        for panel in self.t_panels:
            elements = []
            count += 1
            stiff_spacing = panel.getb()
            numElements = panel.getnstiff()+2
            numStiff = panel.getnstiff()
            for j in range(numElements):
                if j < numStiff:
                    y_i = panel.get_bot()+((j+1)*(stiff_spacing)*math.sin(math.radians(-panel.getOrnt())))
                    #print ("here's the y_i: ", y_i)
                    elements.append(TPanel_trans.TPanel_trans(stiff_spacing,1,1,1,panel.gettp(),panel.gettw(),panel.gethw(),panel.gettf(),panel.getbf(),panel.gettwh() \
                                                            , panel.gettwt(), panel.gettft(),panel.gettfb(),[0,y_i,0],panel.getOrnt(),panel.getQloc(),panel.getpmatl(),panel.getsmatl(),panel.gettmatl(),panel.getEta()))
                elif j == numStiff:
                    y_0 = panel.get_bot()
                    y_last = panel.get_bot()+((numStiff)*(stiff_spacing)*math.sin(math.radians(-panel.getOrnt())))
                    span_0 = stiff_spacing
                    span_end = panel.getB()-(stiff_spacing*(numStiff))-span_0
                    #print ("y_0", y_0)
                    #print ("y_last", y_last)
                    elements.append(TPanel_trans.TPanel_trans(span_0,1,0,1,panel.gettp(),1e-6,1e-6,1e-6,1e-6,panel.gettwh() \
                                                            , panel.gettwt(), panel.gettft(),panel.gettfb(),[0,y_0,0],panel.getOrnt(),panel.getQloc(),panel.getpmatl(),panel.getsmatl(),panel.gettmatl(),panel.getEta()))
                    elements.append(TPanel_trans.TPanel_trans(span_end,1,0,1,panel.gettp(),1e-6,1e-6,1e-6,1e-6,panel.gettwh() \
                                                            , panel.gettwt(), panel.gettft(),panel.gettfb(),[0,y_last,0],panel.getOrnt(),panel.getQloc(),panel.getpmatl(),panel.getsmatl(),panel.gettmatl(),panel.getEta()))
                XSection.append(Element.Element(elements[j], self.strs_data[count], self.strn_data[count], self.AE_data[count], 0, y_i+elements[j].getoffset()))
        return XSection    
    
    def setup(self):
        '''Determine the curvature bounds for the ultimate strength calculator'''
        XSection = self.discretize()
        NA0 = self.NAzeroStrain()
        #print ("NA0=", NA0)
        yloc_list = []
        y_i = []
        for element in XSection:
            yloc_list.append(element.getYloc())
        c1 = abs(NA0 - max(yloc_list))
        c2 = abs(NA0 - min(yloc_list))
        c = max(c1,c2)
        #print ("c=", c)
        Ys = XSection[0].getPanel().getsmatl().getYld() #ATTN! currently just uses yield strength of first element
        E = XSection[0].getPanel().getsmatl().getE() #ATTN! currently just uses E of first element
        YldCrv = Ys/(c*E) #yield curvature in 1/m
        #print ("Yield curvature=", YldCrv, "1/m")
        
        num_crv_inc = 20

        #for sagging condition
        Crv_max_sag = -2*YldCrv
        Crv_min_sag = Crv_max_sag/10
        Crv_step_sag = (Crv_max_sag - Crv_min_sag)/num_crv_inc

        #for hogging condition
        Crv_max_hog = 2*YldCrv
        Crv_min_hog = Crv_max_hog/10
        Crv_step_hog = (Crv_max_hog - Crv_min_hog)/num_crv_inc

        return NA0, E, YldCrv, Crv_max_sag, Crv_min_sag, num_crv_inc, Crv_step_sag, Crv_max_hog, Crv_min_hog, Crv_step_hog
    
    def sumForce_and_moment(self):
        '''Sum the force and moment in the cross-section for a given NA location and curvature'''
        NA0 = self.setup()[0]
        crv_init = self.setup()[2]
        force = self.ApplyCurvature(NA0, crv_init)[0]
        moment = self.ApplyCurvature(NA0, crv_init)[1]
        print ("Ultimate force on cross-section=", force, "N")
        print ("Ultimate moment on cross-section=", moment, "kN-m")
        return force, moment

    def ApplyCurvature(self, NA_guess, curv_init, el_type = 'Norm', mirror = True):
        '''Takes a curvature and a proposed NA position, calcs resulting moment and force'''

        XSection = self.discretize()
        self.NA_guess =  NA_guess
        self.curv_init = curv_init
        self.el_type = el_type
        E = self.setup()[1]

        if mirror == True:
            factor = 2 
        else:
            factor = 1

        empty_force = []
        empty_moment = []
        empty_denom = []

        for element in XSection:
            element.strain = (-(element.getYloc() - self.NA_guess))/(1/self.curv_init)
            if self.el_type == 'PE':
                element.stress = element.strain*E
                element.force = element.stress*element.getPanel().getArea() * factor #in MN
            else:
                element.stress = element.getStress(element.strain)
                element.force = element.getForce(element.stress) * factor #in MN
            element.moment = element.force*(element.getYloc() - self.NA_guess) * factor #in MN*m
            empty_denom.append(E*element.getPanel().getArea()*factor*1e6) #in N
            empty_moment.append(element.moment)
            empty_force.append(element.force)
        total_force = sum(empty_force)*1e6 #in N
        total_moment = sum(empty_moment)*1e3 #in kN*m
        shift = total_force/(self.curv_init*sum(empty_denom)) #in m
        
        return total_force, total_moment, shift
    
    def plotCollapse(self):
        '''Iterates through curvature values and plots the collapse curve
        finding the equalibrium point on each'''
        force_tol = 10000 #in N
        NA0 = self.setup()[0]
        crv_array = []
        M_array = []
        num_crv_inc = self.setup()[5]

        #for sagging condition first
        Crv_min_sag = self.setup()[4]
        Crv_step_sag = self.setup()[6]
        empty_crv_sag = []
        empty_moment_sag = []
        NA0_sag = []

        for i in range(0,num_crv_inc+1):
            curv_applied_sag = Crv_min_sag + i*Crv_step_sag
            force, moment, shift = self.ApplyCurvature(NA0, curv_applied_sag)
            count = 0
            while abs(force) > force_tol:
                count += 1
                NA0 -= shift
                force, moment, shift = self.ApplyCurvature(NA0, curv_applied_sag)
                if count > 25:
                    #print ("NA0 not converging in sag")
                    break
            if abs(force) < force_tol:
                print ("Converged at NA0=", NA0)
                print ("Shift was =", shift)
            empty_crv_sag.append(curv_applied_sag)
            empty_moment_sag.append(moment)
            NA0_sag.append(NA0)
        print ("empty_crv_sag", empty_crv_sag)

        #for hogging condition
        NA0 = self.setup()[0]
        Crv_min_hog = self.setup()[8]
        Crv_step_hog = self.setup()[9]
        empty_crv_hog = []
        empty_moment_hog = []
        NA0_hog = []

        for i in range(0,num_crv_inc+1):
            curv_applied_hog = Crv_min_hog + i*Crv_step_hog
            force, moment, shift = self.ApplyCurvature(NA0, curv_applied_hog)
            count = 0
            while abs(force) > force_tol:
                count += 1
                NA0 += shift
                force, moment, shift = self.ApplyCurvature(NA0, curv_applied_hog)
                if count > 25:
                    print ("NA0 not converging in hog")
                    break
            if abs(force) < force_tol:
                print ("Converged at NA0=", NA0)
                print ("Shift was =", shift)
            empty_crv_hog.append(curv_applied_hog)
            empty_moment_hog.append(moment)
            NA0_hog.append(NA0)
        print ("empty_crv_hog", empty_crv_hog)
        
        crv_array = empty_crv_sag + empty_crv_hog
        M_array = empty_moment_sag + empty_moment_hog
        NA_array = NA0_sag + NA0_hog

        #print ("crv_array", crv_array)

        plt.plot(crv_array,NA_array, 'ro')
        plt.xlabel('Curvature [1/m]')
        plt.ylabel('VBM [kN*m]')
        plt.title('Collapse Curve for Section')
        plt.grid(True)
        plt.show()
        return crv_array, M_array, NA_array


        #M_array_sag.append(total_moment)
        #crv_array_sag.append(curv_applied)

        #max_m_sag = abs(M_array_sag[0])
        #for num in range(0,len(M_array_sag)):
        #    if abs(M_array_sag[num]) > max_m_sag:
        #        max_m_sag = abs(M_array_sag[num])
        #mult_sag = max_m_sag*1000 #in kN*m
        #print ("Ultimate moment on cross-section sagging=", mult_sag, "kN-m")

        '''Calculate strain in each element based on incremental curvature and guess NA location for sagging condition'''
        """ force_tol = 0.3
        empty_force_sag = []
        empty_moment_sag = []
        crv_array_sag = []
        M_array_sag = []
        NA_guess = NA0
        for i in range(0,num_crv_inc+1):
            curv_applied = Crv_min + i*Crv_step
            print ("curv_applied", curv_applied)
            for j in range(0,num_NA_inc+1):
                #print ("j", j)
                print ("NA_guess", NA_guess)
                empty_moment_sag = []
                empty_force_sag = []
                for element in XSection:
                    element.strain = (-(element.getYloc() - NA_guess))/(1/curv_applied)
                    element.stress = element.getStress(element.strain)
                    element.force = element.getForce(element.stress)
                    empty_force_sag.append(element.force)
                total_force = sum(empty_force_sag)*1e6 #in N
                #print ("total force:", total_force)
                if abs(total_force) <= force_tol:
                    print("converged at j =", j)
                    for element in XSection:
                        element.moment = element.force*(element.getYloc() - NA_guess)
                        empty_moment_sag.append(element.moment)
                    M_i = sum(empty_moment_sag)
                    M_array_sag.append(M_i)
                    crv_array_sag.append(curv_applied)
                    break
                else:
                    NA_guess = NAmin + j*NA_step

        max_m_sag = abs(M_array_sag[0])
        for num in range(0,len(M_array_sag)):
            if abs(M_array_sag[num]) > max_m_sag:
                max_m_sag = abs(M_array_sag[num])
        mult_sag = max_m_sag*1000 #in kN*m
        print ("Ultimate moment on cross-section sagging=", mult_sag, "kN-m") """

        '''Calculate strain in each element based on incremental curvature and guess NA location for hogging condition'''
        """ empty_force_hog = []
        empty_moment_hog = []
        crv_array_hog = []
        M_array_hog = []
        NA_guess = NAmin
        for i in range(0,num_crv_inc+1):
            curv_applied = Crv_min + i*Crv_step
            for j in range(0,num_NA_inc+1):
                empty_moment_hog = []
                empty_force_hog = []
                for element in XSection:
                    element.strain = ((element.getYloc() - NA_guess))/(1/curv_applied)
                    element.stress = element.getStress(element.strain)
                    element.force = element.getForce(element.stress)
                    empty_force_hog.append(element.force)
                total_force = sum(empty_force_hog)
                if abs(total_force) <= force_tol:
                    for element in XSection:
                        element.moment = element.force*(element.getYloc() - NA_guess)
                        empty_moment_hog.append(element.moment)
                    M_i = sum(empty_moment_hog)
                    M_array_hog.append(M_i)
                    crv_array_hog.append(curv_applied)
                    NA_guess = NAmin
                    break
                else:
                    NA_guess = NAmin + j*NA_step
                

        max_m_hog = abs(M_array_hog[0])
        for num in range(0,len(M_array_hog)):
            if abs(M_array_hog[num]) > max_m_hog:
                max_m_hog = abs(M_array_hog[num])
        mult_sag = max_m_hog*1000 #in kN*m
        print ("Ultimate moment on cross-section hogging=", mult_sag, "kN-m") """

       #return crv_array_sag, M_array_sag, #crv_array_hogM_array_hog




    def NAzeroStrain(self):
        '''returns the neutral axis of the cross-section in the zero-strain
        condition'''
        TotalMoment = 0
        TotalArea = 0
        XSection = self.discretize()
        for element in XSection:
            TotalMoment += element.getYloc() * element._panel.getArea() #ttl mmnt about baseline
            TotalArea += element._panel.getArea()
        NA_zero_strain = TotalMoment / TotalArea
        return NA_zero_strain
    
