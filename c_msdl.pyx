#Followed example from following web page:
#stackoverflow.com/questions/3046305/simple-wrapping-of-c-code-with-cython

import numpy as np
cimport numpy as np
import util.errors as errors
from cpython cimport bool

cdef extern from "C:/Users/Administrator/Documents/Michigan Documents/First Term/Master's Thesis/Code/Working/MS-Thesis/smith_collapse.h":
    
    ctypedef struct comp:
        pass
    
    ctypedef struct crossSection:
        pass
    
    crossSection* CreateEmptyCross()

    comp* CreateAndAddComponent( int id,
					double x_global, 
                    double y_global, 
                    int epp, 
                    double E,
                    double Ys,
                    double Area,
                    int numPts,
                    double stress[], 
                    double strain[],
                    crossSection *section)

    double momentPureBending(crossSection *section, 
                         double radcurvature,
                         double *NApos, double fTol,
                         double elementStrains[])

    void setupVerticalBending(crossSection *section)
    
    void cleanUp(crossSection* section)

    int findMaxMom(crossSection *section, double minCurv, double maxCurv,
				double forceTol, double curvTol, 
                double* NA, double *moment, double* actualCurv,
                double elementStrains[])
   


cdef class SmithCollapseC:
    '''Python/Cython front end to the underlying C library which implements
    a Smith-type progressive collapse analysis
    
    Major assumptions:
      
     Current only works for vertical bending, assumed to be about an axis where 
        yloc=constant. The underlying c 
        data structures can support arbitray axis bending as yeff 
        and yloc can be different, but the supporting functions are not 
        yet implemented for this
     * Sign convention is - distances positive up from keel (assumed origin)
     * tensile strains and stresses positive
     * hogging bending moments positive
     * Postive radius of curvature gives sagging (negative) response
    '''
    cdef crossSection* _section
    cdef double _forcetol
    cdef int _currID
    cdef int _numEls
    cdef bool _axisSet
    
    def __init__(self, name, forcetol=10.):
        self._section = CreateEmptyCross()
        self._forcetol = forcetol
        self._currID = 0
        self._numEls = 0 #Need to keep track of this for array sizing 
        self._axisSet = False 

    
    #Note - full defintion of stress/strain are required to be sure that we
    #get a nice, contigous block of memory in C-style (no FORTRAN style 
    #, discontinous slices, or the other nonsense that Numpy supports)
    def addElement(self, double xloc, double yloc, double area, epp, 
                    np.ndarray[double, ndim=1, mode="c"] stress not None ,
                    np.ndarray[double, ndim=1, mode="c"] strain not None,
                                E, ys):
        '''Adds an elastic perfectly plastic tension, given load-shortening
        curve in compression object to current cross section
                
        Parameters
        ----------
        xloc:           X location in global coordinates of the midship section
        
        yloc:           Y location in global coordinates of the midship section
        
        area:           Cross-sectional area of the component
        
        epp:            If the element is elastic-plastic in tension
                        True/False
        
        stress:         Numpy 1-D array of avearge effective stress over 
                        provided area of element
                        
        strain:         Numpy 1-D array of strain corrosponding to stress
                        values over the element
                        
        E:              Young's modulus of the element
        
        Ys:             Tensile yield strain of the element
       
        Returns
        -------
        id (interger) of added component to the section, will throw exception if 
        memory error or other problem
        '''
        self._currID += 1
        numpts = len(stress)
        
        #C uses a switch on epp
        if epp==True:
            eppflag = 1
        else:
            eppflag = 0
            
        result = CreateAndAddComponent(self._currID, xloc, yloc, eppflag, E, ys,
                                       area,
					     numpts, &stress[0], &strain[0], self._section)
        if result is NULL:
            print("Warning element addition failed in SmithCollapseC")
            raise errors.SetupError("SmithCollapseC failed componet add")
         
        self._numEls += 1 
        return self._currID
        
    def setupVerticalBending(self):
        '''
        Simple function to setup the C-data structures for bending about the 
        vertical direction, assumed to be about an axis where yloc=constant
        The c-data structures can support arbitray axis bending as yeff 
        and yloc can be different, but the supporting functions are not 
        yet implemented for this
        
        Parameters
        ----------
        None
        
        Returns
        -------
        None
        '''
        setupVerticalBending(self._section)
        self._axisSet = True 
        return
    
    def applyVerticalCurvature(self, double curvature, double NAguess):
        '''Applies a vertical curvature to the cross section
        
        Applied a vertical curvature to the cross section, will setup the 
        section for vertical bending
        
        Parameters
        ----------
        Curvature - a curvature (not radius of curvature as in C!) value.
                    see sign convention
        
        Returns
        -------
        Tuple where first value is the resisting moment, second is the found
        neutral axis position, and third is numpy array of all of the element
        strains at this curvature value
        '''
        setupVerticalBending(self._section)
        
        cdef np.ndarray[double, ndim=1, mode="c"] elements
        elements = np.zeros(self._numEls, dtype=np.double,order="c")
        
        cdef double result
        cdef double NApos 
        
        #Copy this over
        NApos = NAguess
        
        result = momentPureBending(self._section, 1./curvature, &NApos,
                                   self._forcetol, &elements[0])
                                   
        return (result, NApos, elements)
    

    def findMaxVerticalMom(self, double minCurv, double maxCurv, 
                           double curvTol, double NAguess):
        '''Find a maximum in the vertical curvature response 
        
        Finds a maximum in the vertical curvature response between two 
        curvature values. Uses a Brent minimization approach,  if there 
        are multiple extrema, then which one is returned is not defined. Will
        set up the section for vertical bending
        
        Parameters
        ----------
        minCurv - minimum curvature (not radius of curvature) to begin looking
                  at
                  
        maxCurv - maximum curvature (not radius of curvature) to stop looking
                  at
                  
        curvTol - tolerance for locating the peak in curvature domina
        
        NAguess - starting neutral axis guess can be for either curvature
        
        Returns
        -------
        Tuple where first value is the resisting moment, second is the found
        curvature, third is the corrosponding neutral axis position, and the 
        forth is a is numpy array of all of the element
        strains at this curvature value
        '''  
        
        cdef np.ndarray[double, ndim=1, mode="c"] elements
        elements = np.zeros(self._numEls, dtype=np.double,order="c")
        
        cdef double resultMom
        cdef int retCode 
        cdef double resultCurv
        cdef double NApos 
        
        #Copy this over
        NApos = NAguess
        
        setupVerticalBending(self._section)
        
        retCode = findMaxMom(self._section, minCurv,  maxCurv, self._forcetol,
                             curvTol, &NApos, &resultMom, &resultCurv, 
                             &elements[0])
        if retCode == 4:
            print ("Endpoints do not enclose a minimum")
            return [retCode]
        
        #Chek return code for a problem
        if retCode > 0:
            print("Warning peak finding routing in Smith method")
            raise errors.CalculationFail("SmithCollapseC failed peak find")
        
        return (resultMom, resultCurv, NApos, elements)
    
    def changeFTol(self, double newFtol):
        '''changes the force inbalance tolerance for curvature calcs.

        Parameters
        ----------
        newFtol - new value of Ftol                  
        
        Returns
        -------
        none
        '''
        self._forcetol = newFtol        

    def __dealloc__(self):
        '''Only calls the cleanup on the C library side - python handle
        by Cython we hope
        '''
        cleanUp(self._section)
        return

        

        
    
    
#myCrossSection = CreateEmptyCross()
#pan1 = CreateAndAddPanel(4.0, 5.0, 0, 210000., 235.,400, myCrossSection)
#pan2 = CreateAndAddPanel(5.0, 6.0, 0, 210000., 235.,400, myCrossSection)
#testPanels(myCrossSection)
