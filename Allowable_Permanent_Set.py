# Design by Allowable Permanent Set
# Author: Travis Zahradka
# Date: 03-13-10

class Allowable_Permanent_Set:
    """ 
    Class for backing out the compressive pressure on a conventional stiffened panel by method of
    Allowable Permanent Set as described in Hughe's : 
        Chapter 9 Plate Bending : 
            Section 9.3 : Plates Loaded Beyond the Elastic Limit                 pgs. 344-351
            Section 9.4 : Design of Plating Based On Allowable Permanent Set     pgs. 351-355 
    """
    
    # NOTE!!! - wpi is initialized in the _init_ method but should be a panel characteristic in the T_Panel class
    # NOTE!!! - _ar is a method in this class but should be part of the T_Panel Class
    # NOTE!!! - this class does not allow for initial permanent set with 'locked-in' stress
    #           we would ideally like to include this, how would we differentiate 'stress-free' from 'locked-in' permanent stress
        
    def __init__(self, wpi=0, wpt=2.0):
        """
        variable initialization
          wpi - initial 'stress-free' permanent set
          wpt - allowable permanent set
        """
        self._wpi = wpi
        self._wpt = wpt
    
    def _ar(self, a, b):
        """
        calculates aspect ratio of panel
          Panel - any class that implements the basic panel methods
        """
        return b/a
    
    def _wpo(self, beta, tp):       # MAYBE THIS SHOULD BE IN THE T_PANEL CLASS
        """
        calculates initial deflection at yield stress from eq. 9.3.18 in ref
          beta - panel slenderness ratio
          tp - panel thickness
          Note: Assumes a Poisson Ratio for steel of 0.30
        """
        return ((tp*0.07*(beta**2.0))/3.0)        
    
    def _Rw(self, wpi, wpt, wpo):
        # here you could import  alone wpi instead of using the calculated on in this class (wpi is a panel characteristic)
        """
        calculates deflection parameter Rw from eq. 9.4.4 in ref
          wpi - initial 'stress-free' permanent set
          wpt - allowable permament set
          wpo - initial deflection at yield stress
        """
        return ((wpt-wpi)/wpo)
        
    def _Qy(self, beta, v, ar):
        """
        calculates the non-dim load parameter to reach plate yeild stress as shown on pg. 251 in ref
          v - panel material poisson ratio
          beta - panel slenderness ratio
          ar - panel aspect ratio
        """        
        # calculate Qy
        return ((2.0/((1.0-v+v**2.0)**0.5))*(1.0+0.6*(ar**4.0))*(1.0/(beta**2.0)))
    
    def _dQ0(self, beta, v, ar):
        """
        calculates intercept of linear portion of Q/wpt curve from eq. 9.4.2 and shown in 
        figure 9.14 in ref
          ar - panel aspect ratio
          beta - panel slenderness ratio
          v - panel material poisson ratio
        """        
        return ((1.0+0.5*beta*ar*(1.0+ar*(3.3-(1.0/beta))))/(((1.0-v+(v**2.0))**0.5)))*(1/beta**2.0)
        
    def _dQ1(self, beta, ar):
        """
        calculates the further increment of load at the end of the transition zone from eq. 9.4.3 
        and shown in figure 9.14 in ref
          beta - panel slenderness ratio
          ar - panel aspect ratio        
        """        
        return (0.32*(ar/(beta**0.5))**1.5)
    
    def _TRw(self, Rw):   
        """
        calculates the parameter T(Rw) from 9.3.20 in ref
          Rw - deflection parameter calculated in this class
        """  
        if (Rw <= 1.0):                  
            Rw = ((1.0-(1.0-Rw)**3.0)**(1.0/3.0))
        else:
            Rw = 1.0
        return (Rw)
    
    def _Q(self, panel):
        """
        calculates the non dimensional load which causes the allowable permanent set 'wpa'
          b - panel width
          a - panel length
          beta - panel slenderness ratio
          tp - panel thickness
          v - panel material poisson ratio
          
          ar - panel aspect ratio
          wpo - deflection at panel yield stress
          Rw - panel deflection parameter
          
        """
        # Make local copies of needed variables
        self.b = panel.getb()
        self.a = panel.geta()
        self.beta = panel.getBeta()
        self.tp = panel.gettp()
        
        pmat = panel.getmatlP()
        self.v = pmat.getPoisson()
        
        # calculate aspect ratio [ar], yield stress deflection [wpo] and deflection parameter [Rw]
        self.ar = self._ar(self.a, self.b)
        self.wpo = self._wpo(self.beta, self.tp)
        
        Rw = self._Rw(self._wpi, self._wpt, self.wpo)
        
        # calculate Qy for this panel
        Qy = self._Qy(self.beta, self.v, self.ar)  
        # calculate dQ0 for this panel
        dQ0 = self._dQ0(self.beta, self.v, self.ar)
        # calculate dQ1 for this panel
        dQ1 = self._dQ1(self.beta, self.ar)
        # calculate T(Rw)
        TRw = self._TRw(Rw)
        
        # return Q value from eq. 9.4.1 in ref
        Q = Qy+TRw*(dQ0+dQ1*Rw)
        print(Q)
        return Q
        
    def _p_aps(self, panel):
        """
        calculates the dimensional out-of-plane pressure load on a panel from non-dimentional pressure load Q calculated above
          Yld - material yield stress
          E - panel material modulus of elasticity
        """
        Q = self._Q(panel)
        pmat = panel.getmatlP()
        E = pmat.getE()
        Yld = pmat.getYld()
        
        P = Yld**2*Q/E
        print(P)
        print(self.ar)
        print(self.beta)
        print(self.wpo)
        return P
    
    def APS_constraint(self, Panel, pressure):
        """
        calculates Allowable Permanent Set constraint violation
            if lateral pressure load specified is greater than the needed pressure to produce the permanent set
            
        ***Note: to evaluate this, you must set up an instance of this class specifying in the constructor the amount of allowable permanent set (default 2.0 mm)
                    (and initial permanent set can also be specified)
        """
        pressure_needed = self._p_aps(Panel)
        if pressure > pressure_needed:
            c1 = 1.0 - (pressure_needed/pressure)
        else:
            c1 = 0.0
            
        return c1
    
    
        
       
    
    
