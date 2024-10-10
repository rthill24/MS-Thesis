#Fatigue Crack Occurence Probability for Time-Varying loading
#Based on original code by Mark Groden, extended by Matt Collette
# (c) 2012 Regents of the Univesity of Michigan

import math
from scipy import stats
import numpy as np

class CrackDetailVL:
    '''
    Class for calculting the probability of cracking at any time with time-
    varying loads.  Based on the derivation in the document
    "Lognormal Fatigue Models and Python Implementation"
    '''
    def lnPARAM(self, mean, std):
        
        zeta = math.sqrt(math.log(1+(std/mean)**2))
        lam = math.log(mean) - .5*zeta**2
        return (lam, zeta)

    def __init__(self, Amean=3*3.31*10**11, Astd=1.75*10**11, Kfmean=1., 
                 Kfstd=.1, Dmean=1., Dstd=.48, M=3, RefStress=100,
                 Ninc = 500, Nrepair=1):
        '''
        Sets up the crack growth object for subsequent loading
        
        Parameters
        ----------
        Amean:      Mean value of the S-N curve constant A (true mean)
        
        Astd:       Standard deviation of the S-N curve constant A

        Kfmean:     Means value of stress uncertainty factor 

        Kfstd:      Standard deviation of the stress uncertainty factor

        Dmean:      Mean value of the Palmgren-Miner damage summation

        Dstd:       Standard deviatin of the Palmgren-Miner damage summation

        M:          S-N curve slope paramater - deterministic
        
        RefStress:  Reference stress for converting variable amplitude stress
        
        Ninc:       Pre-allocated number of updating intervals (maximum) with
                    re-sizing all arrays (currently firm upper bound)
        
        Nrepair:    Number of repairs to consider, 0 or 1 are implemented
        '''
        #Copy over basic variables
        self.Nrepair = Nrepair
        if Nrepair != 0:
            print("BAD", Nrepair)
            self.RefStress = RefStress
        
        #Make numpy arrays for the crack probabilities and total equiv cycles
        self.RunningTotalCycles = np.empty(Ninc)
        self.NoRepairIncProb = np.empty(Ninc) #Incremental probability 
        #Initial cycles and cracking probability are zero
        self.RunningTotalCycles[0] = 0
        self.NoRepairIncProb[0] = 0
        self.Ninc = 0 #Keep track of number of applied integration pts
        self.NincMax = self.Ninc

        
        #Process fatigue parameters        
        self.M = M
        
        #Using lnParam function, creates arrays with respective lambda 
        #and zeta values 
        self.Kfparam = self.lnPARAM(Kfmean, Kfstd)           
        self.Aparam = self.lnPARAM(Amean, Astd)     
        self.Dcrparam = self.lnPARAM(Dmean, Dstd)
        
        
        #Make a lognormal distribution at the current reference stress
        self.lambdas = self.Dcrparam[0]+self.Aparam[0]-self.M*(self.Kfparam[0]+
                                            math.log(self.RefStress))
        self.zetas = (self.Dcrparam[1]**2 + self.Aparam[1]**2 + (self.M*
                                                self.Kfparam[1])**2)**(.5)
        self.lndist = stats.lognorm(self.zetas, loc=0, 
                                        scale=math.exp(self.lambdas))


    def CrackProbConstLoad(self, StressRange, Cycles):
        '''
        Adds a number of cycles at a constant-amplitude stress range 
        StressRange

        Parameters
        ----------
        StressRange:     Constant amplitude stress range
        
        Cycles:          Number of cycles at StressRange 
        
        Returns 
        -------
        Float of incremental probability of a crack occuring during this 
        load application.  If Nrepair has been set to 1.0, probability 
        includes probability of cracking, being fixed, and re-cracking once.
        '''
        #Increment the number of load application periods
        self.Ninc += 1
        #Convert applied stress to reference stress
        newCycles_Standard = Cycles*(StressRange/self.RefStress)**self.M

        #Add to the running stress total 
        self.RunningTotalCycles[self.Ninc] = (
            self.RunningTotalCycles[self.Ninc -1 ] + newCycles_Standard)
        
        #Calculate the current crack probability
        self.NoRepairIncProb[self.Ninc] = (
              self.lndist.cdf(self.RunningTotalCycles[self.Ninc]) -
              self.lndist.cdf(self.RunningTotalCycles[self.Ninc -1])
              )
        Prob = self.NoRepairIncProb[self.Ninc]

        #Calculate the previous repair probabilities
        if (self.Nrepair == 1) and (self.Ninc > 1):
            for i in range (1, self.Ninc):
                Prob += self.NoRepairIncProb[i]*self.lndist.cdf(
                   self.RunningTotalCycles[self.Ninc] -
                   self.RunningTotalCycles[i])
        return Prob
   
    def CrackProb(self, Moment, CycInc, E, EItot, Y, NA):        
        #Calculating stress range acting upon a crack, denoted delSig
        delSig = (E*Moment*abs(Y-NA))/(EItot)/(10**6)
            
        return self.CrackProbConstLoad(delSig, CycInc)

#Note: This is An ugly way to do this
def newInstance(Nrepair=1):
    return CrackDetailVL(Nrepair=Nrepair)
