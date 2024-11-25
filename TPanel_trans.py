
# Basic Class for defining a grillage of fixed length and breadth with
# both stiffened panels with intermediate transverse T-shaped members 

# Author: Travis Zahradka

# Start Date: 09-20-10
# End Date:   ________

import copy
import math
import Structures
import logging
from numpy import zeros
import math
#LOG_FILENAME = 'example.log'
#logging.basicConfig(filename=LOG_FILENAME,level=logging.DEBUG)

class TPanel_trans(Structures.TPanel):
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
    
    def __init__(self, B, L, nstiff, ntrans, tp, tw, hw, tf, bf, \
                 twh, twt, tft, tfb, sloc=[0,0,0], ornt = 0, \
                 qloc=['NA','NA','NA'], pmatl = 0, smatl = 0, \
                 tmatl = 0, eta = 0):
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
        self.type = 'normal'
        self._B = B
        self._L = L
        if (nstiff != round(nstiff)):
#            logger.debug("number of stiffeners must be an integer: nstiff changed to int(nstiff)")
#            print 'number of stiffeners must be an integer: nstiff changed to int(nstiff)'
            self._nstiff = round(nstiff)
        else:
            self._nstiff = nstiff
        if (ntrans != round(ntrans)):
#            logger.debug("number of transverse members must be an integer: ntrans changed to int(ntrans)")
#            print 'number of transverse members must be an integer: ntrans changed to int(ntrans)'
            self._ntrans = round(ntrans)
        else:
            self._ntrans = ntrans
        self._b = B/(nstiff+1.0)    # breadth of long. stiff. panel
        self._a = L/(ntrans)    # length of long. stiff. panel
        self._tp = tp
        self._tw = tw
        self._hw = hw
        self._tf = tf
        self._bf = bf 
        self.sloc = sloc
        self.ornt = ornt
        self.qloc= qloc
        self._pmatl = pmatl
        if (smatl == 0):
            self._smatl = pmatl
        else:
            self._smatl = smatl
        self._eta = eta
        if (tmatl == 0):
            self._tmatl = pmatl
        else:
            self._tmatl = tmatl
        self._twh = twh
        self._twt = twt
        self._tft = tft
        self._tfb = tfb
        
        self._updateOverallProperties()
    
    def __str__(self):
        printPanel = 'Transverse T-Panel with properties:\n'
        try:
            printPanel += 'Start location: '+str(self.sloc)+'\n'
            printPanel += 'Orientation: '+str(self.ornt)+'\n'
            printPanel += 'Functional location: '+str(self.qloc)+'\n'
        except AttributeError:
            pass
        
        printPanel += 'Breadth: '+str(self._B)+'\n'
        printPanel += 'Number of stiffeners: '+str(self._nstiff)+'\n'
        printPanel += 'Plate thickness: '+str(self._tp)+'\n'
        printPanel += 'Web thickness: '+str(self._tw)+'\n'
        printPanel += 'Web height: '+str(self._hw)+'\n'
        printPanel += 'Flange thickness: '+str(self._tf)+'\n'
        printPanel += 'Flange breadth: '+str(self._bf)+'\n'
        printPanel += 'Transverse web height: '+str(self._twh)+'\n'
        printPanel += 'Transverse web thickness: '+str(self._twt)+'\n'
        printPanel += 'Transverse flange thickness: '+str(self._tft)+'\n'
        printPanel += 'Transverse flange breadth: '+str(self._tfb)+'\n'
        
        return printPanel
        
    def _updateOverallProperties(self):
        """
        updates the overall properties of the panel using inherited and self functions
        """
        self._upProp() 
        self.t_properties()   
        
    def update(self, B=-1, L=-1, nstiff=-1, ntrans=-1, tp=-1, tw=-1, hw=-1, tf=-1, bf=-1, \
                pmatl=-1, smatl=-1, eta=-1, tmatl=-1, twh=-1, twt=-1, tft=-1, tfb=-1, b=-1):
        """
        updates the TPanel_trans definition variables if specified
        ***NOTE: when updating the 'small' panel width [b], no other parameters are changed, i.e. 'B' and 'mstiff' remain the same
                    this property is used, for example, in Faulkner classes where effective panel properties need to be evaluated by replaceing 'b' with 'be'
                    in this way the total Fixed Length/Breadth Panel keeps its large scale properties [refer to author for more details]
        """
        if (B > 0):
            self._B = B
        if (L > 0):
            self._L = L
        if (nstiff >= 0):
            self._nstiff = nstiff
        if (ntrans > 0):
            self._ntrans = ntrans
        if (tp > 0):
            self._tp = tp
        if (tw > 0):
            self._tw = tw
        if (hw > 0):
            self._hw = hw
        if (tf > 0):
            self._tf = tf
        if (bf > 0):
            self._bf = bf
        if (tmatl > 0):
            self._tmatl = tmatl
        if (pmatl > 0):
            self._pmatl = pmatl
        if (smatl > 0):
            self._smatl = smatl
        if (eta > 0):
            self._eta = eta
        if (twh > 0):
            self._twh = twh
        if (twt > 0):
            self._twt = twt
        if (tft > 0):
            self._tft = tft
        if (tfb > 0):
            self._tfb = tfb
            
        self._b = self._B/(self._nstiff+1.0)
        self._a = self._L/(self._ntrans)
        
        # SPECIAL UPDATE - USE WITH UNDERSTANDING AND CAUTION
        if (b > 0):
            self._b = b
            
        self.t_properties()
        self._updateOverallProperties()
        
    def t_properties(self):
        """
        calculates additional properies of the grillage
        """
        # calculate area terms
        self.tpa = self._tp*self._a                  # transverse plate area
        self.twa = self._twt*self._twh               # transverse web area
        self.tfa = self._tft*self._tfb               # transverse flange area
        
        self.ta = self.tpa + self.twa + self.tfa     # transverse area
        self.tsa = self.twa + self.tfa               # transverse stiffener area
        
        # calculate first area moment about bottum of plate
        self.t_fmom = self.tpa*self._tp/2.0 + self.twa*(self._tp+self._twh/2.0) \
            + self.tfa*(self._tp+self._twh+self._tft/2.0)
        self.t_fmomStiff = self.twa*(self._tp+self._twh/2.0) + self.tfa*(self._tp+self._twh+self._tft/2.0)
        
        # calculate Neutral Axis of stiffener & stiffener/plate combination
        self._tNA = self.t_fmom/self.ta
        self._tNAStiff = self.t_fmomStiff/self.tsa
        
        # calculate distance to furthest fiber for section modulus calculation
        self.y_max = max(((self._tp+self._twh+self._tft)-self._tNA), self._tNA)
        
        #calculate spacing between stiffeners
        stiff_spacing = self._B/(self._nstiff+1)

        #calculate top of section considering orientation and stiffeners
        start_top = self.sloc[1]

        if self.ornt > 0 and self.ornt <=90:
            first_stiff_top = -((stiff_spacing)*math.sin(math.radians(self.ornt)))+((self._hw+self._tf+self._tp)*math.sin(math.radians(90-self.ornt)))+((self._bf/2)*math.sin(math.radians(self.ornt)))
            self.top_1 = start_top
            self.top_2 = start_top + first_stiff_top
            self.top = max(self.top_1, self.top_2)

        elif self.ornt > 90 and self.ornt <180:
            beta = self.ornt - 90
            first_stiff_top = -((stiff_spacing)*math.sin(math.radians(beta)))-((self._hw+self._tf+self._tp)*math.sin(math.radians(90-beta)))+((self._bf/2)*math.sin(math.radians(beta)))
            self.top_1 = start_top
            self.top_2 = start_top + first_stiff_top
            self.top = max(self.top_1, self.top_2)

        elif self.ornt == 180 or self.ornt == -180:
            self.top = start_top

        elif self.ornt < 0 and self.ornt >= -90:
            beta = abs(self.ornt)
            plate_only = self._B*math.sin(math.radians(beta))
            furthest_stiff_top = ((self._nstiff*stiff_spacing)*math.sin(math.radians(beta)))+((self._hw+self._tf+self._tp)*math.sin(math.radians(90-beta)))+((self._bf/2)*math.sin(math.radians(beta)))
            self.top_1 = start_top + plate_only 
            self.top_2 = start_top + furthest_stiff_top
            self.top = max(self.top_1, self.top_2) 

        elif self.ornt < -90 and self.ornt > -180:
            beta = 180 - abs(self.ornt)
            plate_only = self._B*math.sin(math.radians(beta))
            furthest_stiff_top = ((self._nstiff*stiff_spacing)*math.sin(math.radians(beta)))-((self._hw+self._tf+self._tp)*math.sin(math.radians(90-beta)))+((self._bf/2)*math.sin(math.radians(beta)))
            self.top_1 = start_top + plate_only
            self.top_2 = start_top + furthest_stiff_top
            self.top = max(self.top_1, self.top_2)
        
        elif self.ornt == 0:
            self.top = start_top + self._tp + self._hw + self._tf  

        #calculate bottom of section considering orientation and stiffeners
        start_bot = self.sloc[1]

        if self.ornt > 0 and self.ornt <=90:
            plate_only = -self._B*math.sin(math.radians(self.ornt))
            furthest_stiff_bot = -((self._nstiff*stiff_spacing)*math.sin(math.radians(self.ornt)))+((self._hw+self._tf+self._tp)*math.sin(math.radians(90-self.ornt)))-((self._bf/2)*math.sin(math.radians(self.ornt)))
            self.bot_1 = start_bot + plate_only
            self.bot_2 = start_bot + furthest_stiff_bot
            self.bot = min(self.bot_1, self.bot_2)

        elif self.ornt > 90 and self.ornt <180:
            beta = self.ornt - 90
            plate_only = -self._B*math.sin(math.radians(beta))
            furthest_stiff_bot = -((self._nstiff*stiff_spacing)*math.sin(math.radians(beta)))-((self._hw+self._tf+self._tp)*math.sin(math.radians(90-beta)))-((self._bf/2)*math.sin(math.radians(beta)))
            self.bot_1 = start_bot + plate_only
            self.bot_2 = start_bot + furthest_stiff_bot
            self.bot = min(self.bot_1, self.bot_2)

        elif self.ornt == 180 or self.ornt == -180:
            self.bot = start_bot - self._tp - self._hw - self._tf

        elif self.ornt < 0 and self.ornt >= -90:
            beta = abs(self.ornt)
            first_stiff_bot = ((stiff_spacing)*math.sin(math.radians(beta)))+((self._hw+self._tf+self._tp)*math.sin(math.radians(90-beta)))-((self._bf/2)*math.sin(math.radians(beta)))
            self.bot_1 = start_bot
            self.bot_2 = start_bot + first_stiff_bot
            self.bot = min(self.bot_1, self.bot_2)
        
        elif self.ornt < -90 and self.ornt > -180:
            beta = abs(self.ornt) - 90
            first_stiff_bot = ((stiff_spacing)*math.sin(math.radians(beta)))-((self._hw+self._tf+self._tp)*math.sin(math.radians(90-beta)))-((self._bf/2)*math.sin(math.radians(beta)))
            self.bot_1 = start_bot
            self.bot_2 = start_bot + first_stiff_bot
            self.bot = min(self.bot_1, self.bot_2)
        
        elif self.ornt == 0:
            self.bot = start_bot




        # calculate the moment of inertia of the cross section about the NA
        # this formulation uses the parallel axis theorem
        self.t_INA = 1.0/12.0*self._a*self._tp**3.0 + self.tpa*(self._tp/2.0 \
                            - self._tNA)**2.0 + \
                     1.0/12.0*self._twt*self._twh**3.0 + self.twa*(self._tp+self._twh/2.0 \
                            - self._tNA)**2.0 + \
                     1.0/12.0*self._tfb*self._tft**3.0 + self.tfa*(self._tp+self._twh \
                            + self._tft/2.0 - self._tNA)**2.0
        self.t_INAStiff = 1.0/12.0*self._twt*self._twh**3.0 + self.twa*(self._tp+self._twh/2.0 \
                            - self._tNAStiff)**2.0 + \
                     1.0/12.0*self._tfb*self._tft**3.0 + self.tfa*(self._tp+self._twh \
                            + self._tft/2.0 - self._tNAStiff)**2.0
                            
        self.t_rad_gyr = (self.t_INA/self.ta)**0.5
        
        # calculate actual section modulus of the transverse member
        self.t_asm = self.t_INA/self.y_max
        
        self.Volume_()
        
    def Volume_(self):
        """
        calculates the entire volume of the cross-stiffened panel
        """
        # calculate the volume of the transverse members
        v_trans = self.tsa*self._B*self._ntrans
        # calculate the volume of the longitudinal stiffeners
        stiff_area = self.getAreaStiff()
        v_stiff = stiff_area*self._nstiff*self._L
        # calculates volume of plating
        tp = self.gettp()
        v_plate = self._B*self._L*tp
        
        self.total_volume = v_trans + v_stiff + v_plate
        
    def constraints(self, C2=1.0, trans_stiff_ratio=2.0, w1=1.0, w2=1.0, w3=1.0, w4=1.0, w5=1.0, w6=1.0):                                                  
        """
        evaluates ABS stiffener buckling criteria: Part 3, Chapter 2, Section 4 [3-2-4]
            C2 - constant specified in ABS rules [3-2-4] - taken as 1.0 by default
            trans_stiff_ratio - allowable trans web height: stiff web height ratio (1.0 would mean they can be the same height)
            w# - respective constraint weights 
        """
    # EVALULATE STIFFENER w/ ABS
        c1 = zeros((2))
        hw = self.gethw()
        tw = self.gettw()
        smatl = self.getmatlS()
        E = smatl.getE()
        sys = smatl.getYld()
        # evaluate ABS stiff_web constraint
        if ((hw/tw) <= 1.5*((E/sys)**0.5)*C2):
            c1[0] = 0.0
        else:
            c1[0] = 1.0-(1.5*((E/sys)**0.5)*C2)/(hw/tw)
            
        # evaluate ABS stiff_flange constraint
        bf = self.getbf()/2.0
        tf = self.gettf()
        if ((bf/tf) <= (0.5*((E/sys)**0.5)*C2)):
            c1[1] = 0.0
        else:
            c1[1] = 1.0-((0.5*((E/sys)**0.5)*C2)/(bf/tf))
        
    # EVALUATE TRANSVERSE MEMBER w/ ABS
        c2 = zeros((2))
        twh = self.gettwh()     # transverse web height
        twt = self.gettwt()     # transverse web thickness
        tmatl = self.getmatlT() # transverse material
        Et = tmatl.getE()
        syst = tmatl.getYld()
        # evaluate ABS trans_web constraint
        if ((twh/twt) <= 1.5*((Et/syst)**0.5)*C2):
            c2[0] = 0.0
        else:
            c2[0] = 1.0-(1.5*((Et/syst)**0.5)*C2)/(twh/twt)
            
        # evaluate ABS trans_flange
        tfb = self.gettfb()
        tft = self.gettft()
        if ((tfb/tft) <= (0.5*((Et/syst)**0.5)*C2)):
            c2[1] = 0.0
        else:
            c2[1] = 1.0-((0.5*((Et/syst)**0.5)*C2)/(tfb/tft))
            
    # EVALUATE SUM(STIFF_FLANGEBREADTH*NSTIFF) < PANEL WIDTH
        B = self.getB()
        nstiff = self.getnstiff()
        flange_width_sum = bf*nstiff
        diff = B-flange_width_sum
        if (diff >= 0.0):
            c3 = 0.0
        else:
            c3 = 1.0-(B/flange_width_sum)
            
    # EVALUATE STIFF_HEIGHT TO TRANS_HEIGHT CONSTRAINT
        stiff_height = hw+tf
        if (twh/stiff_height > trans_stiff_ratio):
            c4 = 0.0
        else:
            c4 = 1.0 - (twh/stiff_height)/2.0
            
    # EVALUATE TRANSVERSE MEMBER FLANGE TO WEB ASPECT RATIO
        if (tfb > twh/4.0):
            c5 = 0.0
        else:
            c5 = 1.0 - (tfb/(twh/4.0))
    
    # EVALUATE STIFFENER MEMBER FLANGE TO WEB ASPECT RATIO
        if (bf > hw/4.0):
            c6 = 0.0
        else:
            c6 = 1.0 - (bf/(hw/4.0))
    
        
        panel_constraint_violation = w1*c1[0], w1*c1[1], w2*c2[0], w2*c2[1], w3*c3, w4*c4, w5*c5, w6*c6

        return panel_constraint_violation
    

    # easy function for additional parameter calls   
    def getB(self):
        """returns total breadth of grillage"""
        return self._B
    
    def getL(self):
        """returns total length of grillage"""
        return self._L
    
    def getnstiff(self):
        """returns the number of stiffeners"""
        return self._nstiff
    
    def getntrans(self):
        """returns number of transverse members"""
        return self._ntrans
    
    def getmatlT(self):
        """returns EPMat class information for the transvese member"""
        return self._tmatl
    
    def gettwh(self):
        """returns transverse web height"""
        return self._twh
    
    def gettwt(self):
        """returns transverse web thickness"""
        return self._twt
    
    def gettft(self):
        """returns transverse flange thickness"""
        return self._tft
    
    def gettfb(self):
        """returns transverse flange breadth"""
        return self._tfb
    
    def gettpa(self):
        """returns 'transverse' plate area (x-section of transverse memer)"""
        return self.tpa
    
    def gettwa(self):
        """returns transverse web area"""
        return self.twa
    
    def gettfa(self):
        """returns transverse flange area"""
        return self.tfa
    
    def getta(self):
        """returns total area of transverse member (x-section of transverse)"""
        return self.ta
    
    def gettsa(self):
        """returns transverse stiffener area (x-section of transverse)"""
        return self.tsa
    
    def getb(self):
        """returns stiffener spacing"""
        return self._b
    
    def geta(self):
        """returns length of panel (dist. between transverse members)"""
        return self._a
    
    def gett_fmom(self):
        """returns first area moment of inertia of transverse member about plate base"""
        return self.t_fmom
    
    def gett_fmomStiff(self):
        """return first area moment of inertia of trnasverse stiffener about plat base"""
        return self.t_fmomStiff
    
    def get_tNA(self):
        """returns Neutral Axis of plate/transverse member combination from base of plate"""
        return self._tNA
    
    def get_tNAStiff(self):
        """returns Neautral Axis of transverse member from base of plate"""
        return self._tNAStiff
    
    def gety_max(self):
        """returns distance to furthest fiber from Neutral Axis"""
        return self.y_max
    
    def gett_INA(self):
        """returns moment of inertial of plate/trans. member comb. about its neutral axis"""
        return self.t_INA
    
    def gett_INAStiff(self):
        """returns moment of inertia of trans. member about its neutral axis"""
        return self.t_INAStiff
    
    def gett_rad_gyr(self):
        """return radius of gyration of plate/trans. member comb."""
        return self.t_rad_gyr
    
    def gett_asm(self):
        """returns actual section modulus of transverse member"""
        return self.t_asm
    
    def getTotalVolume(self):
        """returns the total volume of the cross-stiffened panel"""
        # BE CAREFULE! - get_TotalVolume in derived class returns volume of individual panel
        return self.total_volume
        
    def getSinglePanel(self):
        '''returns a TPanel defined by the corresponding Transverse T-Panel Grillage'''
        return Structures.TPanel(self._b, self._tp, self._tw, self._hw, self._tf, self._bf, self._a, self._pmatl, self._smatl, self._eta)
        
    def geoValid(self):
        '''checks flange overlap and returns true for valid geometry'''
        if self._nstiff > 1:
            if self._b > self._bf:
                return bool(1)
            else:
                return bool(0)
        else:
            if 2*self._b > self._bf:
                return bool(1)
            else:
                return bool(0)
            
    def getXloc(self):
        ''' Returns the X location of the start of the T panel'''
        return self.sloc[0]
        
    def getYloc(self):
        ''' Returns the Y location of the start of the T panel'''
        return self.sloc[1]
        
    def getZloc(self):
        ''' Returns the Z location of the start of the T panel'''
        return self.sloc[2]     
    def getOrnt(self):
        ''' Retuns the Orientation of the panels, based on east as 0 and counterclockwise as positive'''
        return self.ornt
    def get_top(self):
        ''' Returns the top of the panel'''
        return self.top
    def get_bot(self):
        ''' Returns the bot of the panel'''
        return self.bot
    def getpmatl(self):
        ''' Returns the material properties of the plating'''
        return self._pmatl
    def getsmatl(self):
        ''' Returns the material properties of the stiffeners'''
        return self._smatl
    def gettmatl(self):
        ''' Returns the material properties of the transverse members'''
        return self._tmatl