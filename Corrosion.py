# -*- coding: utf-8 -*-
#Module for various corrosion models so that they can be interchangeably
#used with other code based on degredation
# (c) 2012 Regents of the Univesity of Michigan

#import util.errors as errors

class baseCorrosion:
    '''Basic corrosion class with common methods that should be overridden
    '''
    
    def __init__(self, idString):
        '''Base class initializer allows intelligent error messages
        
        Parameters
        ----------
        idString : string identifying the corrosion model for error messages
        '''
        self._modelIDString = idString
    
    def updateThickness(self, time, thickness):
        '''General thickness updating function
        
        Parameters
        ----------
        Time:   Time in years
        
        Thickness:  Original (NOT current) thickness to be reduced
        
        Returns
        -------
        Reduced thickness at the current time
        '''
        
        raise errors.SetupError(
         "Error update thickness called on corrosion model that does not"
         + "support it, corrosion model was " + self._modelIDString)
 
        return
    
    def updatePlateThickness(self, time, thickness):
        '''Plate thickness updating function
        
        Parameters
        ----------
        Time:   Time in years
        
        Thickness:  Original (NOT current) thickness to be reduced
        
        Returns
        -------
        Reduced thickness at the current time
        '''
        
        raise errors.SetupError(
         "Error update plate thickness called on corrosion model that does not"
         + "support it, corrosion model was " + self._modelIDString)
         
        return        

    def updateWebThickness(self, time, thickness):
        '''Web thickness updating function
        
        Parameters
        ----------
        Time:   Time in years
        
        Thickness:  Original (NOT current) thickness to be reduced
        
        Returns
        -------
        Reduced thickness at the current time
        '''
        
        raise errors.SetupError(
         "Error update web thickness called on corrosion model that does not"
         + "support it, corrosion model was " + self._modelIDString)
         
        return           

    def updateFlangeThickness(self, time, thickness):
        '''Web thickness updating function
        
        Parameters
        ----------
        Time:   Time in years
        
        Thickness:  Original (NOT current) thickness to be reduced
        
        Returns
        -------
        Reduced thickness at the current time
        '''
        
        raise errors.SetupError(
        "Error update flange thickness called on corrosion model that does not"
         + "support it, corrosion model was " + self._modelIDString)
         
        return


class paikCorrosion(baseCorrosion):
    '''Implements the Paik corrosion model based on M. Lankowskis' work
    
     Paik, J.K.; Lee, J.M.; Hwang, J.S.; Park, Y.I.  "A Time-Dependant 
     Corrosion Wastage Model for the Structures of Single- and Double-Hull 
     Tankers and FSOs and FPSOs."  Marine Technology, Vol. 40, No. 3, July 
     2003, pp. 201-217  
     '''
    def __init__(self, coatingLife, locationCodePlate=None,
                                    locationCodeWeb=None,
                                    locationCodeFlange=None,
                                    units='meters'):
        '''Sets up the object

        Sets up the location code and coating life.  Any location code can
        be omitted or passed as None which causes that component to have 
        zero corroion (all thickness returned as originals) 
        
        Parameters
        ----------
        coatingLife:            Life of the coating in years
        locationCodePlate:      A valid loaction code from Paik's paper:
        locationCodeWeb:        A valid loaction code from Paik's paper:
        locationCodeFlange:     A valid location code from Paik's paper:
        units:                  Units of the thickness function, Valid options
                                all as strings:
                                    'meters' :default
                                    'mm': millimeters
                                    'in': inches
        '''
        
        #Call base class constructor
        baseCorrosion.__init__(self, 'Paik Corrosion Model')        
        
        #Save coating life        
        self._coatingLife = coatingLife
        
        
        #Set up the unit conversion for different thickness values
        conv = {'meters':1., 'mm':1000.,  'in':39.3701}

        if not (units in conv):
            raise errors.SetupError("Error - invalid thickness unit passed"
                            + " to Paik corrosion model, unit was " + 
                            str(units) )                            
        
        #Process the location specific code
        
        #Catch any Nones
        if locationCodePlate == None:
            locationCodePlate = 'NA'
        if locationCodeWeb == None:
            locationCodeWeb = 'NA'
        if locationCodeFlange == None:
            locationCodeFlange = 'NA'
               
        #These are local dictionaries for corrosion rates from the paper
        plates = {'BSH': 0.0597e-3, 'ABH': 0.1084e-3, 'ABV': 0.0661e-3,
                       'BSV': 0.0622e-3, 'BLGB': 0.0619e-3, 'OBV': 0.1012e-3,
                       'BBH': 0.1408e-3, 'OSH': 0.0607e-3, 'AOH': 0.0581e-3,
                       'AOV': 0.0523e-3, 'OSV': 0.0423e-3, 'BLGC':0.0414e-3,
                       'OOV': 0.0577e-3, 'OOH':0.0405e-3, 'NA': 0.0}

        webs = {'BSLBW': 0.1367e-3, 'DLBW': 0.2403e-3, 'SSLBW': 0.1413e-3,
                'LBLBW': 0.1960e-3, 'BSLCW': 0.0466e-3, 'DLCW': 0.0716e-3,
                'SSLCW': 0.0420e-3, 'LBLCW': 0.0550e-3, 'BGLCW': 0.0377e-3,
                'DGLCW': 0.0477e-3, 'SSTLCW': 0.0261e-3, 'NA': 0.0}

        flanges = {'BSLBF': 0.1127e-3, 'SSLBF': 0.0882e-3, 'LBLBF': 0.1782e-3,
                   'BSLCF': 0.0437e-3, 'DLCF': 0.0588e-3, 'SSLCF': 0.0397e-3,
                   'LBLCF': 0.0508e-3, 'BGLCF': 0.0319e-3, 'DGLCF': 0.0449e-3,
                   'NA': 0.0}
                   
        #Check that we have a valid code for each
        if not (locationCodePlate in plates):
            raise errors.SetupError("Error invalid plate location code in "
                                     + " paik corrosion model" + 
                                     str(locationCodePlate))
        if not (locationCodeWeb in webs):
            raise errors.SetupError("Error invalid web location code in "
                                     + " paik corrosion model" + 
                                     str(locationCodeWeb))
        if not (locationCodeFlange in flanges):
            raise errors.SetupError("Error invalid flange location code in "
                                     + " paik corrosion model" + 
                                     str(locationCodeFlange))
        
        #Assign the variables
        self._pc = plates[locationCodePlate]*conv[units]
        self._wc = webs[locationCodeWeb]*conv[units]
        self._fc = flanges[locationCodeFlange]*conv[units]

        return        
        
    
    def _updater(self,time, thickness, rate):
        '''Internal utility function for doing updates
        '''
        #Start with the full thickness
        retVal = thickness
        #If beyond the coating life, start to decrease
        if time > self._coatingLife:
            retVal -= rate*(time - self._coatingLife)
#            print (time - self._coatingLife)
        #If over-corroded, set thickness to zero
        if retVal < 0:
            retVal = 0

        return retVal   
    
    def coating(self, time):
        '''Returns true/false if coating has failed at specified time
        
        Parameters
        ----------
        Time:   Time in years
        
        
        Returns
        -------
        True if coating is still OK, False if beyond coating life
        '''
        if time <= self._coatingLife:
            return True
        else:
            return False
        return

    def updatePlateThickness(self, time, thickness):
        '''Plate thickness updating function
        
        Parameters
        ----------
        Time:   Time in years
        
        Thickness:  Original (NOT current) thickness to be reduced
        
        Returns
        -------
        Reduced thickness at the current time
        '''
        
        return self._updater(time, thickness, self._pc)      


    def updateWebThickness(self, time, thickness):
        '''Web thickness updating function
        
        Parameters
        ----------
        Time:   Time in years
        
        Thickness:  Original (NOT current) thickness to be reduced
        
        Returns
        -------
        Reduced thickness at the current time
        '''
        
        return self._updater(time, thickness, self._wc)       

    def updateFlangeThickness(self, time, thickness):
        '''Web thickness updating function
        
        Parameters
        ----------
        Time:   Time in years
        
        Thickness:  Original (NOT current) thickness to be reduced
        
        Returns
        -------
        Reduced thickness at the current time
        '''
        
        return self._updater(time, thickness, self._fc)        
                        
        
     
    
    