# Method for calculating the cost of producing a stiffened panel with and without Transverse Web Frames
# Initial Author: Yi-Jen Wang

# Date Created: Fri Jul 23 16:04:00 2010
# Date Complete: __________

# Major Updates:
#    + 10/26/10 - Travis Zahradka - Added Transverse Web Calculation class [Cost_trans] inheriting [CostCal] class attributes
#        - *Note: a TPanel_trans can be passed to the original class [CostCal] - it will cost the panel only and inherit 
#                    ntrans from that class to calculate lstiff/trans intersection cost only - NO OTHER trans costs!
#
#    + 09/19/11 - Matthew Lankowski - Changed required input from millimeters to meters
#
#    + 02/08/13 - Doug Rigterink - Added simple weld pass calculator for cost of stiffeners

import math
from numpy import zeros

class CostCal:
    '''
    Class for calculating the cost of Tpanel is referenced from:
    SHIP STRUCTURES: IMPROVEMENT BY RATIONAL DESIGN OPTIMISATION
    M. K. Rahman and J. B. Caldwell 1992
    
    C[0] - Total Cost of TPanel WITHOUT transverse frames
        C[1] - cost of materials for hull plates
        C[2] - cost of material for longitudinal stiffeners
        C[3] - cost of material for transverse frames [excluded]
        C[4] - cost of welding for longitudinal stiffeners
        C[5] - cost of welding for transverse frames [excluded]
        C[6] - cost of intersections between lstiff and transverse frames [ntrans assumed = 1.0]
        C[7] - cost of preparation of brackets and joints [excluded]
        C[8] - cost of electricity and electrodes [total] [lstiff ONLY!]
        C[9] - fabrication cost of lstiff and transverse frames [excluded]
    '''
    def __init__(self,Pa=860,r=2660/907.185,Clm=1.05,Ps=27,Cls=1.2,Cis=0.6,Cbj=1.15,Cee=0.9,Cfb=1.5):
        '''
        set up the coefficient in order to be used later
        Pa = material price U.S$ / ton
        r = Specific weight of the material used ton/m^3
        Clm = material cost coefficient for longitudinal stiffeners
        Ps = labor rate US$/hr
        Cls = the labor hour required per meter welding of stiffeners to plate
        Cis = the labor hour required per intersection welding of stiffeners to 
                transverse frames
        Cbj = labor hour required per joint of stiffeners to transverse frames
        Cee = labor hour required per meter of stiffeners implementing electricity
        Cfb = labor hour required per meter of stiffeners for fabrication
        '''
        self.Pa = Pa
        self.r = r
        self.Clm = Clm
        self.Ps = Ps
        self.Cls = Cls
        self.Cis = Cis
        self.Cbj = Cbj 
        self.Cee = Cee
        self.Cfb = Cfb
                
    def Cplate(self,panel,length,B):
        '''
        cost of materials for hull plates
        '''
        #Total cross-sectional area of bottom plates
        B = B
        tp = panel.gettp()
        ATpb = B * tp
        #Total cross-sectional area of one-side plates
        ATps = 0
        #Total cross-sectional area of deck plates
        ATpd = 0
        #length of the cargo hold designed
        l = length
        #weight of plating
        Wp = (ATpb + 2*ATps + ATpd)*self.r*l
        return Wp * self.Pa
    
    def Cstiff(self,panel,nstiff,length):
        '''
        cost of materials for longitudinal stiffeners
        '''
        #Total cross-sectional area of longitudinal stiffeners in bottom
        ATsb = (panel.gettw() * panel.gethw() + panel.gettf() * panel.getbf())*nstiff
        #Total cross-sectional area of longitudinal stiffeners in one side
        ATss = 0
        #Total cross-sectional area of longitudinal stiffeners indeck
        ATsd = 0
        #Total cross-sectional area of longitudinal centre girder in bottom
        ATgb = 0
        #Total cross-sectional area of longitudinal centre girder deck
        ATgd = 0
        #weight of longitudinal stiffeners
        Wls = (ATsb + 2*ATss + ATsd + ATgb + ATgd)*self.r*length       
        return Wls*self.Clm*self.Pa
        
    def Cweld(self,panel,nstiff,length,Welding=0):                                     
        '''
        cost of welding for longitudinal stiffeners
        '''
        WeldPass=1.  #number of weld passes
        min_plate=0.
        self.Welding = Welding
        if self.Welding != 0:   #single weld pass is min plate thickness is less than 10.5mm
            min_plate=min(panel.gettw(),panel.gettp())
            WeldPass=1+(min_plate/5.-1)/7            
            #if min_plate>10.5:
             #   WeldPass=math.ceil(min_plate/10.5)
            

        return nstiff*length*self.Cls*self.Ps*WeldPass
        
    def Cintersect(self, nstiff, ntrans=1):
        '''
        cost of intersections between longitudinal stiffeners and transverse 
        frames and preparation of brackets and joints
        '''        
        return nstiff*ntrans*(self.Cbj + self.Cis)*self.Ps                          # IS THIS METHOD IMPLEMENTED CORRECTLY??
        
    def Celectrict_lstiff(self, nstiff, length):
        '''
        cost of electricity, electrodes and fabrication cost of longitudinal stiffeners  
        '''
        return nstiff*length*(self.Cee+self.Cfb)*self.Ps
    
    def TotalCost(self, panel, nstiff, length, B, ntrans=1.0):
        '''
        calculates and returns the total cost of the Tpanel without transverse members
        '''
        # initialize cost vector
        C = zeros((10))
        C[1] = self.Cplate(panel, length, B)
        C[2] = self.Cstiff(panel, nstiff, length)
        C[4] = self.Cweld(panel, nstiff, length)
        C[6] = self.Cintersect(nstiff, ntrans)
        C[8] = self.Celectrict_lstiff(nstiff, length)
        C[9] = 0.0
        C[0] = sum(C[1:10])
        self.C = C

        return C[0] 
                     
class Cost_trans(CostCal):
    """
    class for calculating the production cost of a cross stiffened panel of fixed length and breadth
    """
    
    def __init__(self, Panel, Welding=0, Cfm=1.4, Cwf=1.25):
        """
        constructor
            Panel - TPanel_trans member class - implementing tpanel functions with transverse web frames
                  - this is a '2D' panel with transverse and long. stiffeners
        """
        self.Panel = Panel
        self.Cfm = Cfm
        self.Cwf = Cwf
        self.Welding = Welding
        # run init method of derived class
        CostCal.__init__(self)
        
    def Cost_add(self):
        """
        this method implements the weight based cost function given in:
            - Ship Structures: Improvement by Rational Design and Optimizations - M.K. Rahman [1973] pgs. 88-91
        * This class assumes specific weight of material is constant for plating, longitudinal stiffeners and transverse members
        * To return total cost - Cost_add()[0] - i.e C[0] -or- Total_Cost_()
        
        C[0] - Total Cost of Panel WITH transverse frames
        C[1] - cost of materials for hull plates
        C[2] - cost of material for longitudinal stiffeners
        C[3] - cost of material for transverse frames
        C[4] - cost of welding for longitudinal stiffeners
        C[5] - cost of welding for transverse frames
        C[6] - cost of intersections between lstiff and transverse frames
        C[7] - cost of preparation of brackets and joints
        C[8] - cost of electricity and electrodes [total]
        C[9] - fabrication cost of lstiff and transverse frames
        """
        Atb = self.Panel.getta()                # area of tranverse frame
        B = self.Panel.getB()                   # length of transverse frame
        gamma = self.r                          # specific weight of material used
        Nw = self.Panel.getntrans()             # number of transverse frames
        
        self.Wtf = Atb*B*gamma*Nw               # total weight of transverse frames
        
        # get cost of TPanel without transverse members from derived class
        panel = self.Panel
        nstiff = self.Panel.getnstiff()
        length = self.Panel.getL()
        self.TotalCost(panel, nstiff, length, B, Nw)
        C = self.C
    
        Pa = self.Pa
        Ps = self.Ps
        
        WeldPass=1  #number of weld passes
        
        if self.Welding != 0:   #single weld pass is min plate thickness is less than 10.5mm
            min_plate=min(self.Panel.gettwt(),self.Panel.gettp())
            
            WeldPass=1+(min_plate/5.-1)/15
            #if min_plate>10.5:
             #   WeldPass=math.ceil(min_plate/10.5)
        
        C[3] = self.Wtf*self.Cfm*Pa        # cost of material for transverse frame
        C[5] = 2.0*B*self.Cwf*Ps*Nw*WeldPass        # cost of welding for transverse frames
        C[7] = 0.0
        C[8] = C[8]+2*B*Nw*self.Cee*Ps     # cost of electricity & electrodes for welding transverse members
        C[9] = 0.0                         # fabrication costs of lstiff & trans members          
        
        C[0] = sum(C[1:10])
        self.C = C
        
        return C
    
    def Total_Cost_(self):
        """
        returns the total cost of the TPanel_trans only
        """
        Total_Cost_ = self.Cost_add()[0]
        return Total_Cost_
        
    def get_all_cost(self):
        """
        returns plate cost
        """
        all_cost=self.Cost_add()[0:10]
        return all_cost