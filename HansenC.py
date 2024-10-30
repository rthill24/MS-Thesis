#A.M. Hansen Panel Ultimate Strength Method for clamped-end panels
import math
import numpy as np
from scipy import optimize

import Structures

class HansenC:
    '''Class for calculating the stress-strain curve for clamped-end,
    conventional stiffened panels by the method presented by Anders Hansen, 
    as described in:
        "Reliability Methods for the Longitudinal Strength of Ships," Chpt. 2
        Technical University of Denmark, Jan. 1995'''

    def __init__(self, panel, w0=0, q=0, theta=0, output=0):
        self._output = output        
        
        self._tf = panel.gettf()
        self._bf = panel.getbf()
        self._tw = panel.gettw()
        self._hw = panel.gethw()
        self._tp = panel.gettp()
        self._b = panel.getb()
        self._A = panel.getArea()
        self._L = panel.geta()
        self._E = panel.getmatlP().getE()
        self._s0 = panel.getYsavg()
        self._NA = panel.getNA()
        self._matlP = panel.getmatlP()
        self._w0 = w0
        self._q = q
        self._theta = theta
            
        self._h = self._tp + self._hw + self._tf   #h = total height of panel       
        self._squashload = self._s0*self._A
        
        self.epanel = Structures.TPanel(self._b, self._tp, self._tw, self._hw, self._tf, self._bf, self._L, self._matlP)        
        
        self._LScurve()
        
    def _midbeamstress(self, P, z):
        '''returns the normal stress distribution over the cross-section of
        the beam-column at distance z from the neutral axis, given an in-plane
        load P, initial rotation theta, initial displacement w0, and transverse
        load q.'''
        IE = self.epanel.getINA()
        AE = self.epanel.getArea()        
        PE = ((math.pi/self._L)**2)*self._E*IE
        alpha = (P/(self._E*IE))**0.5
        return - P/AE - (self._w0*self._E*z*P*math.pi**2)/((self._L**2)*(P - PE)) - \
        (self._E*self._q*z)/P + \
        ((alpha*self._E*z)/(math.sin(alpha*self._L/2)))*((self._w0*P*math.pi)/(self._L*(P-PE)) + \
        (self._L*self._q)/(2*P) + self._theta)
        
    def _make_epanel(self, P):
        '''uses the Faulkner method to return an updated panel with effective
        breadth factor k, effective breadth b, effective moment of inertia I,
        and effective area A due to in-plane load P, initial rotation of theta,
        initial deformation w0, and transverse load q.'''
        b = self._b
        self.epanel.update(b=self._b)
        z = self._NA - 0.5*self._tp
        se = abs(self._midbeamstress(P, z))
        k1 = 1.0
        for i in range (1,1000):
            #calculate k2 using Faulkner's expression            
            beta = (b/self._tp)*((se/self._E)**0.5)
            if (beta > 1.0):
                k2 = (2.0/beta) - (1.0/(beta**2))
            else:
                k2 = 1.0
            #check for convergence
            if (abs(k1-k2) < 1.e-4):
                return
            else:               
                k1 = k2
                be = b*k1                
                self.epanel.update(b=be)
                z = self.epanel.getNA() - 0.5*self._tp               
                se = abs(self._midbeamstress(P, z))
            
        #If after i iterations, k does not converge, print error message   
        if self._output == 1:
            print('Error - Faulkner k method did not converge for panel at load',P)
            print('Setting k = 1.0 by default, but results may be inaccurate!!!')
        self.epanel.update(b=self._b)

    def _collapsecheck(self, P):
        '''checks for collapse of the panel based on an in-plane load P, initial
        angle of rotation theta, initial deformation w0, and transverse load q.
        Collapse defined as yielding first appearing in one of two z-locations:
        1)middle of the plate and 2)outer fibre of the stiffener.  Returns four
        solutions (positive and negative wxx)'''
        self._make_epanel(P)
        #check for yielding in middle of plate                    
        sy12 = self.epanel._pmatl.getYld()
        z0p12 = self.epanel.getNA() - 0.5*self._tp
        cc1 = sy12 - abs(self._midbeamstress(P, z0p12))
        cc2 = sy12 - abs(self._midbeamstress(P, -z0p12))
        #check for yielding in outer fibre of stiffener
        sy34 = self.epanel._smatl.getYld()
        z0p34 = ((self.epanel.gettp() + self.epanel.gethw() + self.epanel.gettf())
        - self.epanel.getNA())
        cc3 = sy34 - abs(self._midbeamstress(P, z0p34))
        cc4 = sy34 - abs(self._midbeamstress(P, -z0p34))

        if cc1 > 0 and cc2 > 0 and cc3 > 0 and cc4 > 0:
            return 0
        else:
            return 1
    
    def _Region2(self, P):
        dNA = (P*self._L)/(self._E*self._A)
        eNA = dNA / self._L
        return eNA

    def _Region3(self, P):
        L = self.epanel.geta()
        E = self.epanel._pmatl.getE()
        AE = self.epanel.getArea()
        dNA = (P*L)/(E*AE)
        eNA = dNA / L
        return eNA

    def _Region4(self, Pterm, Pc):
        '''Returns the neutral axis strain in the post-collapse region with the
        method from the following source:
            Rutherford & Caldwell, "Ultimate Longitudinal Strength of Ships: A
            Case Study."  SNAME, 1990'''
        global P
        P = Pterm
        global H
        H = self._h
        tol = H/1000
        c = 0
        if _Finda2(0, self.epanel)*_Finda2(H, self.epanel) > 0:
            if self._output == 1:
                print('_Finda2(a) * _Finda2(b) >= 0, using alternate c Finder...')     
            c = GoldenSectionSearch(0, H, tol)
        else:
            c = optimize.brentq(_Finda2, 0, H, args=(self.epanel), xtol=tol)
        b = self._R4breadth
        A = self._R4Area  
        Ap = b*self._tp
        Aw = self._hw*self._tw
        Af = self._bf*self._tf
        sm = P/A
        s1 = self._s0 - sm
        s2 = self._s0 + sm       
        
        #for plastic neutral axis in plate
        if c >= 0 and c <= self._tp:     
            Mp = s1*(Af*(self._h - self._tf/2) + Aw*(self._tp + self._hw/2) + \
            b*(self._tp**2 - c**2)/2) - 0.5*s2*b*c**2
#            print 'for P =',P,'plastic neutral axis in plate'
    
        #for plastic neutral axis in web
        if c > self._tp and c <= self._tp + self._hw:
            Mp = s1*(Af*(self._h - self._tf/2) + \
            self._tw*((self._hw + self._tp)**2 - c**2)/2) - \
            s2*(self._tw*(c**2 - self._tp**2) + Ap*self._tp)/2
#            print 'for P =',P,'plastic neutral axis in web'
            
        #for plastic neutral axis in flange
        if c > self._tp + self._hw and c <= self._h:       
            Mp = s1*self._bf*(self._h**2 - c**2)/2 - \
            s2*(self._bf*(self._tf - self._h + c)*(c - self._tf + self._h)/2 + \
            Aw*(self._tp + self._hw/2) + Ap*self._tp/2)
#            print 'for P =',P,'plastic neutral axis in flange'
            
        w = Mp/P
        z = self._NA - c
        eNA = Pc/(A*self._E) + 1 - (2*(self._L**2/4 - w**2)**0.5)/self._L + \
        4*z*w/self._L**2  
        
        return eNA
        
       
    def _LScurve(self):
        '''Returns the stress-strain-effective area relations for a stiffened
        panel'''
        
        Pincr = self._squashload / 50 #define the load-increase increment
        
        #Region 1   
        P = [-self._squashload]
        maxstrain = abs(self._Region2(P[0])*10) #how much strain to include in data
        self._strn = np.array([-maxstrain])
        self._AE = np.array([self._A])
        
        #Region 1/2 boundary (tension yield load)
        P.insert(1, -self._squashload)
        self._strn = np.append(self._strn, self._Region2(P[1]))
        self._AE = np.append(self._AE, self._A)
        ccheck = 0
        while ccheck == 0:
            Pterm = P[-1] + Pincr
            P.append(Pterm)
            
            #Region 2 (Elastic Tension Region)            
            if P[-1] <= 0:
                self._strn = np.append(self._strn, self._Region2(Pterm))
                self._AE = np.append(self._AE, self._A)
            
            #Region 3 (Elastic Compression Region)
            else:
                ccheck = self._collapsecheck(Pterm)
                self._strn = np.append(self._strn, self._Region3(Pterm))
                self._AE = np.append(self._AE, self.epanel.getArea())
        self._Pc = max(P)
        self._MaxR3strn = max(self._strn)
        self._R4breadth = self.epanel.getb()
        self._R4Area = self._AE[-1]
        
        #Region 4 (Post-collapse Compression Region)  
        Pterm = self._Pc - Pincr
        count = 0 #This counter is so there are at least 4 data points in Region 4 so it can be splined in post-processing (e.g. Element class object)
        newmax = maxstrain
        while self._strn[-1] < maxstrain and Pterm > 2*Pincr or count < 5:
            Pterm = Pterm - Pincr
            P.append(Pterm)
            PrevStrain = self._strn[-1]
            self._strn = np.append(self._strn, self._Region4(Pterm, self._Pc))
            self._AE = np.append(self._AE, self.epanel.getArea())
            count += 1
            
            #NOTE - THIS IS A PATCH TO AVOID CRASHING FOR UNREALISTIC PANEL GEOMETRIES
            #HERE IS AN EXAMPLE PANEL THAT FAILS FOR FUTURE DEBUGGING
            #smatl = Structures.EPMatl(475e6,200e9,0.26)
            #Panel = Structures.TPanel(.45, .0011816, 0.000392101080662, 0.34712179011, 0.00452926244859, 0.128898833475, 1.8288, smatl)     
            if self._strn[-1] <= PrevStrain:
                if PrevStrain >= maxstrain:                   
                    newmax += maxstrain*0.1
                    self._strn[-1] = newmax
                else:
                    self._strn[-1] = maxstrain
                self._AE[-1] = self._AE[len(self._AE)-2]
                
        #Add last data point to Region 4 so max negative strain = max positive strain
        if self._strn[-1] < maxstrain:
            P.append(P[-1])
            self._strn = np.append(self._strn, maxstrain)
            self._AE = np.append(self._AE, self._AE[-1])
            
        self._strss = np.array(P) / self._A
        self.P = np.array(P)
        
        
def _Finda1(a1guess, epanel):
    if a1guess <= epanel.gettp():
        A = a1guess*epanel.getb()
    elif a1guess > epanel.gettp() and a1guess <= epanel.gettp() + epanel.gethw():
        A = epanel.getb()*epanel.gettp() + epanel.gettw()*(a1guess - epanel.gettp())
    elif a1guess > epanel.gettp() + epanel.gethw() and a1guess <= H:
        A = epanel.getb()*epanel.gettp() + epanel.gettw()*epanel.gethw() + \
        epanel.getbf()*(a1guess - epanel.gettp() - epanel.gethw())  
    return 2*epanel.getYsavg()*A - epanel.getYsavg()*epanel.getArea() - P
   
   
def _Finda2(a2guess, epanel):
    if a2guess <= epanel.gettp():
        A = epanel.getArea() - a2guess*epanel.getb()
    elif a2guess > epanel.gettp() and a2guess <= epanel.gettp() + epanel.gethw():
        A = epanel.getArea() - epanel.gettp()*epanel.getb() - \
        epanel.gettw()*(a2guess - epanel.gettp())
    elif a2guess > epanel.gettp() + epanel.gethw() and a2guess <= H:
        A = epanel.getbf()*(H - a2guess)   
    return 2*epanel.getYsavg()*A - epanel.getYsavg()*epanel.getArea() - P
    

def GoldenSectionSearch(x1, x2, tolerance):
    
    tau = 0.3819660112501051
    x3 = tau*(x2 - x1) + x1
    x4 = tau*(x2 - x3) + x3
    f1 = abs(_Finda2(x1))
    f2 = abs(_Finda2(x2))
    f3 = abs(_Finda2(x3))
    f4 = abs(_Finda2(x4))
    
    x = np.array([x1,x2,x3,x4])
    y = np.array([f1,f2,f3,f4])
    x_closest = x[y.argmin()]    
    
    while abs(max(x) - min(x)) > tolerance:
        if f4 > f3:
            if x4 > x3:
                x2 = x4
                f2 = f4
                x4 = x3 - tau*(x3 - x1)
                f4 = abs(_Finda2(x4))
            else:
                x1 = x4
                f1 = f4
                x4 = tau*(x2 - x3) + x3
                f4 = abs(_Finda2(x4))
            
        if f4 < f3:
            if x4 > x3:
                x1 = x3
                f1 = f3
                x3 = x4
                f3 = f4
                x4 = tau*(x2 - x3) + x3
                f4 = abs(_Finda2(x4))
            else:
                x2 = x3
                f2 = f3
                x3 = x4
                f3 = f4
                x4 = tau*(x3 - x1) + x1
                f4 = abs(_Finda2(x4))
        
        x = np.array([x1,x2,x3,x4])
        y = np.array([f1,f2,f3,f4])
        x_closest = x[y.argmin()]
    
    return x_closest
