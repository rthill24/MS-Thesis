#Tee-stiffened panel geometry plus simple elastic-plastic material
# (c) 2009 University of Michigan, Dr. Matthew Collette NA&ME Dept.

import math
import copy
import numpy.linalg as linalg
import logging


class EPMatl:
    """ Simple elastic-plastic material class with a yield or proof stress and
    and an elastic modulus"""
    def __init__(self, yld, E, Poisson=0.3):
        self._yld = yld
        self._E = E
        self._Poisson=Poisson
    
    def __str__(self):
        printMat = 'EPMatl instance with properties:\n'
        printMat += 'Yield strength: '+str(self._yld)+'\n'
        printMat += 'Modulus of elasticity: '+str(self._E)+'\n'
        printMat += "Poisson's ratio: "+str(self._Poisson)+'\n'
        
        return printMat

    def update(self, yld=-1, E=-1, Poisson=-1):
        '''
        Allows updates of material properties, two parameters yld for
        yield stress
        or proof stress and E for elastic modulus
        '''
        if yld > 0:
            self._yld = yld
        if E > 0:
            self._E=E
        if Poisson>0:
            self._Poisson=Poisson

    def getYld(self):
        '''returns the material's yield stress'''
        return self._yld

    def getE(self):
        '''retunrs the material's elastic modulus'''
        return self._E
    
    def getPoisson(self):
        '''returns the Poisson's ratio'''
        return self._Poisson

class TPanel:
    """
    A simple panel consisting of a span of plating and a single attached
    tee stiffener.  
    """
    def __init__(self, b, tp, tw, hw, tf, bf, a, pmatl, smatl=0, eta=0):
        #Class implicitly assumes one elastic modulus 
        if (smatl != 0 ):
            assert pmatl.getE() == smatl.getE() , \
               "Error stiffener and plate must have same modulus in TPanel"
        self._b = b
        self._tp = tp
        self._tw = tw
        self._hw = hw
        self._tf = tf
        self._bf = bf
        self._pmatl = pmatl;
        self._a = a;
        if (smatl == 0):
            self._smatl = pmatl
        else:
            self._smatl = smatl
        self._eta = eta
        self._upProp()
        
    def _upProp(self):
        # Calculate area in terms - we  will use these later for other
        #expressions
        Aplate = self._b*self._tp
        Aweb = self._tw*self._hw
        Aflange = self._tf*self._bf
        self._area =  Aplate + Aweb + Aflange
        self._areaStiff=Aweb+Aflange                            
        

        #Caclulate the area-averaged yield stress
        self._avgys = (Aplate*self._pmatl.getYld() + \
                       (Aweb+Aflange)*self._smatl.getYld())/self._area
        
        #Calculate first moment about bottom of plate
        fmom = Aplate*self._tp/2.0 + Aweb*(self._tp + self._hw/2.0) + \
               Aflange*(self._tp + self._hw + self._tf/2.0)
        fmonStiff=Aweb*(self._tp + self._hw/2.0) + Aflange*(self._tp + \
                  self._hw + self._tf/2.0)
        self._NA = fmom/self._area
        self._NAStiff=fmonStiff/self._areaStiff                         
        
        #Calculate the moment of intertia of the cross section about NA
        self._INA = 1.0/12.0*self._b*self._tp**3.0 + Aplate*(self._tp/2.0 \
                    - self._NA)**2.0 \
                    + 1.0/12.0*self._tw*self._hw**3.0 + Aweb*(self._tp + \
                      self._hw/2.0 - self._NA)**2.0 \
                    + 1.0/12.0*self._bf*self._tf**3.0 + \
                    Aflange*(self._tp + self._hw + self._tf/2.0 - self._NA)**2.0

        self._INAStiff=1.0/12.0*self._tw*self._hw**3.0 + Aweb*(self._tp + \
                        self._hw/2.0 - self._NAStiff)**2.0+ \
                        1.0/12.0*self._bf*self._tf**3.0 + \
                        Aflange*(self._tp + self._hw + self._tf/2.0 - \
                        self._NAStiff)**2.0

        #print self._INA
        
        rad_gyr = (self._INA/self._area)**0.5
        
        #Calculate the slenderness parameters beta and lambda
        #We have asserted Eplate=Estiffener
        matl_term_plate = (self._pmatl.getYld()/self._pmatl.getE())**0.5
        matl_term_section = (self._avgys/self._pmatl.getE())**0.5
        
        
        self._Beta = self._b/self._tp*matl_term_plate
        self._Lambda = self._a/(math.pi*rad_gyr)*matl_term_section
        #print self._a
        #print math.pi
        #print rad_gyr
        #print matl_term_section
        #print self._Lambda
    def update(self, b=-1, tp=-1, tw=-1, hw=-1, tf=-1, bf=-1, pmatl=0, \
               smatl=0, eta=-1, a=-1):
        """
        Update parameters of an existing panel.  Possible parameter are:
        b - plate width
        tp - plate thickness
        tw - web thickness
        hw -  height of web
        tf - flange thickness
        bf - flange breadth
        pmatl -  Plate material from EPMatl class
        smatl - Stiffener material from EPMatl class
        eta - residual stress tension block width
        a - panel length
        If parameters are not passed, original values are kept
        """
        if b > 0:
            self._b = b
        if tp > 0:
            self._tp = tp
        if tw > 0:
            self._tw = tw
        if hw > 0:
            self._hw = hw
        if tf > 0:
            self._tf = tf
        if bf > 0:
            self._bf = bf
        if eta > 0:
            self._eta = eta
        if a > 0:
            self._a = a

        #Material handling requires some complex logic
        if (pmatl != 0) and (smatl != 0):
            assert pmatl.getE() == smatl.getE(), \
            "Error stiffener and plate must have same modulus in TPanel"
            self._pmatl = pmatl
            self._smatl = smatl
        elif pmatl !=0 :
            assert self._smatl.getE() == pmatl.getE(), \
                   "Error stiffener and plate must have same modulus in TPanel"
            self._pmatl = pmatl
        elif smatl !=0:
            assert self._pmatl.getE() == smatl.getE(), \
                   "Error stiffener and plate must have same modulus in TPanel"
            self._smatl = smatl      
        #No else block - else is case of neither passed in        
        self._upProp()

    def getb(self):
        '''return plate width, b'''
        return self._b

    def gettp(self):
        '''returns plate thickness tp'''
        return self._tp

    def gettw(self):
        '''returns plate web thickness tw'''
        return self._tw

    def gethw(self):
        '''returns plate weh height hw'''
        return self._hw

    def gettf(self):
        '''returns plate flange thickness tf'''
        return self._tf

    def getbf(self):
        '''returns flange breadth bf'''
        return self._bf

    def getmatlS(self):
        '''returns current stiffener material'''
        return self._smatl

    def getmatlP(self):
        '''returns current plate material'''
        return self._pmatl

    def geta(self):
        '''returns overall length'''
        return self._a

    def getBeta(self):
        '''returns plate slenderness ratio beta'''
        return self._Beta

    def getLambda(self):
        '''returns column slenderness ratio lambda'''
        return self._Lambda

    def getArea(self):
        '''returns cross-sectional area'''
        return self._area
    
    def getAreaStiff(self):
        '''returns stiffner area'''
        return self._areaStiff

    def getNA(self):
        '''returns neutral axis'''
        return self._NA

    def getNAStiff(self):
        '''returns stiffner neutral axis'''
        return self._NAStiff

    def getINA(self):
        '''returns moment of intertia'''
        return self._INA

    def getINAStiff(self):
        '''returns moment of intertia of stiffners'''
        return self._INAStiff

    def getEta(self):
        '''return residual stress tension block width'''
        return self._eta

    def getYsavg(self):
        '''returns the area-average yield stress'''
        return self._avgys


#class ConstantThicknessPlate:
#    '''
#    A plate class for simple, constant thickness, metallic material.
#    
#    '''
#    def __init__(matl,tp, length, width, basept, widthv,lengthv):
#        '''
#        Defines the basic plate properties
#        matl - material that the plate is made out of
#        tp - plate thickness
#        length - plate length
#        width - plate width
#        basept- 3-d base point locating corner one in the structure
#        widthv- vector pointing along the width direction
#        lengthv - vector pointing along the length direction
#        widthv and lengthv must be in opposite directions
#        Corners are labeled as:
#              4 ------- 3
#              |         |          
#              L         |
#              E         |
#              N         |
#              G         |
#              T         |
#              H         |
#              |         |
#              1 -width- 2
#        '''
#        self._basemat = matl
#        self._tp = tp
#        self._length = length
#        self._width = width
#        self._basept = basept
#        #Normalize the width and span vectors
#        self._normwv = linalg.norm(widthv)
#        self._normlv = linalg.norm(lengthv)
#        
#        self._doCalcs()
#    
#    def _doCalcs(self):
#        '''
#        Function for updating overall plate properties - will print a 
#        warning if plate definition is off
#        '''
#        AngleTol = 0.001  #Tolerance in radians for width/length vector
#                          #angle
#        
#        #Check that angle between vectors is pretty much 90 deg
#        angle = math.acos(linalg.dot(self._normwv, slef._normlv))
#        if (math.abs(angle) - 1.57079633 > 0.001):
#            logging.ERROR("Error - not-rectangular plate detected")
#            logging.ERROR("Base point x" + str(self._basept[0]))
#            logging.ERROR("Base point y" + str(self._basept[1]))
#            logging.ERROR("Base point z" + str(self._basept[2]))
#            
#            
#        

    
class FixedWidthPanel(TPanel):
    """
    A tee-panel consisting of an overall width, and a variable number of 
    stiffeners, allowing the stiffener spacing to be computed automatically
    All of the original TPanel methods still work, and return properties
    associated with a single plate panel and stiffener. 
    """
    
    def __init__(self, B, nstiff, tp, tw, hw, tf, bf, a, pmatl, smatl=0, eta=0,sloc=[0,0,0], ornt = 0):
        # Class implicitly assumes one elastic modulus 
        if (smatl != 0 ):
            assert pmatl.getE() == smatl.getE() , \
               "Error stiffener and plate must have same modulus in TPanel"
        self._B = B
        self._nstiff = nstiff
        self._b = B/(nstiff+1.0)
        self._tp = tp
        self._tw = tw
        self._hw = hw
        self._tf = tf
        self._bf = bf
        self._pmatl = pmatl;
        self._a = a;
        self._L = a
        self.sloc = sloc
        self.ornt = ornt
        if (smatl == 0):
            self._smatl = pmatl
        else:
            self._smatl = smatl
        self._eta = eta
        self._updateOverallProp()
    
    def _updateOverallProp(self):
        '''
        Updates overall properties of the panel
        '''
        self._upProp()
        self._totalVolume = (self._B*self._tp+self._nstiff*(self._tw*self._hw+self._bf*self._tf)) * self._a
    
    def getBOverall(self):
        '''
        returns the overall width
        '''
        return self._B
    
    def getta(self):
        '''
        returns the area of one stiffener
        '''
        return self.getAreaStiff()
    def getNStiff(self):
        '''
        returns the number of stiffeners
        '''
        return self._nstiff
    def getnstiff(self):
        '''
        returns the number of stiffeners
        '''
        return self._nstiff
    
    def getTotalVolume(self):
        '''
        returns the total volume of the panel
        '''
        return self._totalVolume
    def getArea(self):
        '''
        returns total area of panel
        '''
        panelArea = self._B*self._tp
        singleStiffArea = self.getAreaStiff()
        return panelArea+self._nstiff*singleStiffArea
    
    def update(self, b=-1, B=-1, tp=-1, tw=-1, hw=-1, tf=-1, bf=-1, pmatl=0, \
               smatl=0, eta=-1, a=-1, nstiff=-1):
        """
        Update parameters of an existing panel.  Possible parameter are:
        B - total panel width
        tp - plate thickness
        tw - web thickness
        hw -  height of web
        tf - flange thickness
        bf - flange breadth
        pmatl -  Plate material from EPMatl class
        smatl - Stiffener material from EPMatl class
        eta - residual stress tension block width
        a - panel length
        nstiff - number of stiffeners
        If parameters are not passed, original values are kept
        plate spacing can not be directly set!
        """
        if b>0:
            self._b = b
            self.nstiff = math.ceil(B/b-1)
        if B > 0:
            self._B = B
            self.nstiff = math.ceil(B/b-1)
        if tp > 0:
            self._tp = tp
        if tw > 0:
            self._tw = tw
        if hw > 0:
            self._hw = hw
        if tf > 0:
            self._tf = tf
        if bf > 0:
            self._bf = bf
        if eta > 0:
            self._eta = eta
        if a > 0:
            self._a = a
        if nstiff > 0:
            self._nstiff = nstiff
            self._b = self._B/(nstiff + 1.0)


        #Material handling requires some complex logic
        if (pmatl != 0) and (smatl != 0):
            assert pmatl.getE() == smatl.getE(), \
            "Error stiffener and plate must have same modulus in TPanel"
            self._pmatl = pmatl
            self._smatl = smatl
        elif pmatl !=0 :
            assert self._smatl.getE() == pmatl.getE(), \
                   "Error stiffener and plate must have same modulus in TPanel"
            self._pmatl = pmatl
        elif smatl !=0:
            assert self._pmatl.getE() == smatl.getE(), \
                   "Error stiffener and plate must have same modulus in TPanel"
            self._smatl = smatl      
        #No else block - else is case of neither passed in        
        self._updateOverallProp()
    
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
    def getB(self):
        return self._B
    def getntrans(self):
        """returns number of transverse members"""
        return 0
    def getL(self):
        """returns total length of grillage"""
        return self._L

