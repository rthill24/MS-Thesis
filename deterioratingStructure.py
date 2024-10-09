'''
Module for representing deteriating structure

(c) 2012 The Regents of the University of Michigan
'''
import logging
#from methods.analysis import CrackDetailVL

class repairCost(object):
    """A basic class for returning repair cost for any type of structure.
    
    This class can store both a repair cost, and optionally be extended to 
    include dependencies (e.g. repair need to be performed in drydock) or
    different levels of repair
    """
    
    def __init__(self, cost):
        '''Simple init in the base class, only stores a cost associated 
        with repair
        
        Parameters
        ----------
        cost:   Cost to repair this object
        '''
        self._cost = cost
    
    def getRenewalCost(self):
        '''Returns the cost of renewing the object to new
        
        Parameters
        ----------
        None
        
        Returns
        -------
        Cost of repair
        '''
        return self._cost


class Loading(object):
    '''
    A base class for loading to apply to a dereriorating structure
    '''
    
    def __init__(self, load):
        '''
        Constructor takes a load value
        '''
        
        self._load = load
    
    def load(self):
        '''
        Returns the value of the loading object
        '''
        
        return self._load
    
    
class deterioratingStructure(object):
    """
    A base class for all member of structure that may deteriorate in service

    Methods
    -------
    The followings methods must be over-ridden in your derived problem class:
        
    __init__:    to copy data to base constructor
    
    repairNeeds:  returns a repair cost object without changing the structure    
    
    repair:     restores an object to as-new, returning a repairCost object
    
    age:    applies both time and loading to an object
    """

#Add repair cost to object     
    
    def __init__(self, RepairCostObject):
        '''Init - set up logger, stores repair cost
        
        Parameters
        ----------
        RepairCostObject:  Any object that satisfies the repairCost template
        '''
        self._logger = logging.getLogger('msdl')
        self._repaircost = RepairCostObject
        
        return
        
    def repairNeeds(self):
        '''Returns the repair needs of the object (reference) without changing
        the object's state
        
        Parameters
        ----------
        None
        
        Returns 
        -------
        An object that satifies the repair cost template
        '''    
        return self._repaircost
    
    def renew(self):
        '''Should be overriden to return object as new, returning the incurred
        repair cost
        
        Parameters
        ----------
        None
        
        Returns 
        -------
        An object that satifies the repair cost template        
        '''
        pass

    def age(self, time, loading):
        '''Applied an amount of aging to the object
        
        Applies aging to the object, in terms of both a time and a loading
        or loading components
        
        Parameters
        ----------
        time:   An amount of additional elasped time, measure in years

        loading: A loading class object        
        
        Returns
        -------
        None  - object internal state is updated       
        '''
        pass

class Fatigue_Loading(Loading):
    '''
    A class to represent fatigue loading for a TPanel
    '''
    
    def __init__(self, load, numCycles):
        '''
        Base constructor takes a load and a number of fatigue cycles
        '''
        self._numCycles = numCycles
        Loading.__init__(self, load)

    
    def get_fatigue_cycles(self):
        '''
        Returns the number of fatigue cycles of the loading object
        '''
        
        return self._numCycles
        
    
class TPanel_Repair(repairCost):
    '''
    A repair cost class designed for use with the TransTPanelDet class
    '''
    def __init__(self, cost, recoat_cost):
        '''
        Constructor takes a cost value, a cost to repair a fatigue crack value and
        a cost to recoat the panel value
        '''
        
        self.recoat_cost = recoat_cost
        repairCost.__init__(self, cost)
    
    def recoat(self):
        '''
        Returns the cost of recoating the panel.
        
        Parameters
        ----------
        None
        
        Returns
        -------
        self.recoat_cost - The cost to recoat the panel
        '''
        
        return self.recoat_cost
        
class TransTPanelDet(deterioratingStructure):
    '''
    Basic class that will take a Transverse T Panel and have it decay
    '''
    def __init__(self, RepairCost, TTransPanel, CorrosionModel, FatigueDetails = None,
                 limit_tp = -1., limit_tws = -1., limit_tfs = -1.,
                 limit_twf = -1., limit_tff = -1., age = -1.):
        '''
        Basic constructor taked a repair cost, TPanel_Trans, and corrosion 
        model
        
        Parameters
        ----------
        RepairCost:     A TPanel_Repair object for the repair cost of this panel 
        TTransPanel:    A TPanel_trans object with the basic, intact t-Panel
        CorrosionModel: A corrosion model which support plate, web, and flange
                            corrosion
        FatigueDetails: A list of fatigue detail objects corresponding to the fatigue details
                            of the Transverse T Panel (default=None) 
        limit_tp:       Limiting plate thickness before repair is needed
        limit_tws:      Limiting stiffener web thickness before repair needed
        limit_tfs:      Limiting stiffener flange thick. before repair needed
        limit_twf:      Limiting frame web thickness before repair needed
        limit_tff:      Limiting frame flange thick. before repair needed
        age:            Initial  age for both structure and coating, if not zero
        '''
       
        #Call base class constructor
        deterioratingStructure.__init__(self, RepairCost)

        #Store base models        
        self._CurrentTPanel = TTransPanel
        self._Corrosion = CorrosionModel    
        self._FatigueDetails = FatigueDetails
        
        #Store original dimensions for renewal
        self._origTp = TTransPanel.gettp()
        self._origTw_stiff = TTransPanel.gettw()
        self._origTf_stiff = TTransPanel.gettf()
        self._origTw_frame = TTransPanel.gettwt()
        self._origTf_frame = TTransPanel.gettft()
        self._origAge = age
        
        #Store current baseline values for models (allows re-coating)
        self._baseTp = TTransPanel.gettp()
        self._baseTw_stiff = TTransPanel.gettw()
        self._baseTf_stiff = TTransPanel.gettf()
        self._baseTw_frame = TTransPanel.gettwt()
        self._baseTf_frame = TTransPanel.gettft()
        self._baseNA = TTransPanel.getINA()
        self._baseArea = TTransPanel.getArea()
        self._baseSM = self._baseNA / TTransPanel.gety_max()
        
        #Store limiting corrosion values
        self._limit_tp = limit_tp
        self._limit_tws = limit_tws
        self._limit_tfs = limit_tfs
        self._limit_twf = limit_twf
        self._limit_tff = limit_tff
        
        #Store reference age and coating age
        self._age = age
        self._coatingAge = age
    
    def add_fatigue_detail(self, detail):
        '''
        Adds a fatigue detail object to the deteriorating T panel class.  These
        details will collect fatigue damage as the structure ages.
        
        Parameters
        ----------
        detail: Instance of the FatigueDetail class
        
        Return
        ------
        No return value
        '''
        if not self._FatigueDetails:
            self._FatigueDetails = [detail]
        else:
            self._FatigueDetails.append(detail)
        
    def recoat(self, newCorrosionModel ):
        '''Updates the corrosion model for a re-coating without steel renewal
        
        Parameters
        ----------
        newCorrosionModel  any corrosion model which support plate, web, and 
                            flange, corrosion
        
        Returns
        -------
        self._repaircost.recoat() - The cost to recoat the T-Panel
        '''
        #Reset model and age of coating        
        self._Corrosion = newCorrosionModel
        self._coatingAge = 0.
        
        #Current dimension are new baseline
        self._baseTp = self._CurrentTPanel.gettp()
        self._baseTw_stiff = self._CurrentTPanel.gettw()
        self._baseTf_stiff = self._CurrentTPanel.gettf()
        self._baseTw_frame = self._CurrentTPanel.gettwt()
        self._baseTf_frame = self._CurrentTPanel.gettft()
        
        return self._repaircost.recoat()
    
    
    def renew(self, newCorrosionModel, newFatigueModel=None, Nrepair=0):
        '''Returns the object to new, along with a new corrosion model
        
        Parameters
        ----------
        None
        
        Returns 
        -------
        An object that satifies the repair cost template with the renew cost        
        '''
        
        #Update the T-Panel to as-new
        self._CurrentTPanel.update(tp=self._origTp, tw=self._origTw_stiff,  
                                   tf=self._origTf_stiff, 
                                   twt=self._origTw_frame, 
                                   tft=self._origTf_frame)
        
        #Update the corrosion model
        self._Corrosion = newCorrosionModel
#        self._Fatigue = newFatigueModel
        
        #Reset the baseline values and coating age
        self._baseTp = self._origTp
        self._baseTw_stiff = self._origTw_stiff
        self._baseTf_stiff = self._origTf_stiff
        self._baseTw_frame = self._origTw_frame
        self._baseTf_frame = self._origTf_frame

        self._coatingAge = self._origAge
        self._age = self._origAge
        
        #Reset the fatigue details if present
        if self._FatigueDetails and newFatigueModel:
            for detail in self._FatigueDetails:
                detail.renew(newFatigueModel.newInstance(Nrepair=Nrepair))

                
        
        return self._repaircost

    def age(self, time, loading):
        '''Applied an amount of aging to the object
        
        Applies aging to the object, in terms of both a time and a loading
        or loading components
        
        Parameters
        ----------
        time:   An amount of _additional_ elasped time, measure in years

        loading: A loading class object        
        
        Returns
        -------
        None  - object internal state is updated       
        '''
        #Age the object
        self._age += time
        self._coatingAge += time
        
        #Figure out the corrosion for the object, and update T-panel properties
        newtp =  self._Corrosion.updatePlateThickness(self._coatingAge, 
                                                      self._baseTp)  
        newtw_stiff = self._Corrosion.updateWebThickness(self._coatingAge, 
                                                         self._baseTw_stiff) 
        newtf_stiff = self._Corrosion.updateFlangeThickness(self._coatingAge, 
                                                         self._baseTf_stiff) 
        newtw_frame = self._Corrosion.updateWebThickness(self._coatingAge, 
                                                         self._baseTw_frame) 
        newtf_frame = self._Corrosion.updateFlangeThickness(self._coatingAge, 
                                                         self._baseTf_frame) 
        self._CurrentTPanel.update(tp=newtp, 
                                   tw=newtw_stiff,  
                                   tf=newtf_stiff, 
                                   twt=newtw_frame, 
                                   tft=newtf_frame)
                                   
        #Age fatigue details if present
        if self._FatigueDetails and loading:
            for detail in self._FatigueDetails:
                detail.age(time, loading)
        return 
    
    def getTTPanRef(self):
        '''Returns a reference to the current Tpanel_trans with updated 
        thickness
        
        Parameters
        ----------
        None
        
        Returns
        -------
        Reference to the current TPanel_trans, changes made will be reflected
        in this call
        '''     
        return self._CurrentTPanel
        
    def needsRepair(self):
        '''Returns true or false if a repair is needed
        
        Parameters
        ----------
        None
        
        Returns
        -------
        True if one or more components below minimum thickness, false otherwise        
        '''
        retVal = False
        if self._CurrentTPanel.gettp() < self._limit_tp:
            retVal = True
        if self._CurrentTPanel.gettw() < self._limit_tws:
            retVal = True            
        if self._CurrentTPanel.gettf() < self._limit_tfs:
            retVal = True
        if self._CurrentTPanel.gettwt() < self._limit_twf:
            retVal = True
        if self._CurrentTPanel.gettft() < self._limit_tff:
            retVal = True

        return retVal
    
    def getCorrosionModel(self):
        '''
        Returns a reference to the current corrosion model.
        
        Parameters
        ----------
        No parameters
        
        Returns
        -------
        Reference to the CorrosionModel object assoaciated with this
        instance of TransTPanelDet.
        '''
        return self._Corrosion
    
    def getRepairObject(self):
        '''
        Returns a reference to this T Panels repairCost object. 
        
        Parameters
        ----------
        None
        
        Returns
        -------
        A reference to the repairCost object associated with this
        instance of TransTPanelDet.
        '''
        return self._repaircost
    
    def getFatigueDetails(self):
        '''
        Returns a reference to an associated list of fatigue details
        '''
        return self._FatigueDetails



class FatigueDetail(deterioratingStructure):
    '''
    Class that takes a single panel/stiffener and subjects it to fatigue
    '''
    
    def __init__(self, RepairCost, FatigueModel, age=0):
        '''
        Constructor that takes a repair cost object and a fatigue model
        that support constant stress range loading and returns a probability of fatigue damage in
        the timestep.
        
        Parameters
        -----------
        RepairCost:     A RepairCost object whose renewal cost is equal to the cost to
                            repair a fatigue crack
        FatigueModel:   An incemental fatigue costing model
        age:            The age of the detail (default=0)
        '''
        
        #Call base class constructor
        deterioratingStructure.__init__(self, RepairCost)
        
        #Store base fatigue model
        self._Fatigue = FatigueModel
        
        #Initialize probability of crack instantaneous/cumulative
        self._ProbCrack = 0.0
        self._cumProbCrack = 0.0
        self._age = age
        
    def age(self, time, loading):
        '''
        Applies an amount of age to the fatigue detail, in terms of a timestep and
        a fatigue loading object
        
        Parameters
        -----------
        time:       An incemental fatigue costing model
        loading:    A Fatigue_Loading object
        
        Returns
        -------
        No return value
        '''
        #Age the detail
        self._age += 1
        
        #Update the fatigue crack probability
        stress_range = loading.load()
        numCycles = loading.get_fatigue_cycles()
        self._ProbCrack = self._Fatigue.CrackProbConstLoad(stress_range, numCycles)
        self._cumProbCrack += self._ProbCrack
        
        return
    
    def renew(self, newFatigueModel):
        '''Returns the object to new, along with a new fatigue model
        
        Parameters
        ----------
        newFatigueModel:    A new fatigue model for the detail
        
        Returns 
        -------
        An object that satifies the repair cost template with the renew cost        
        '''  
        
        #Reset crack probabilities
        self._cumCrackProb = 0.0
        
        #Reset fatigue model
        self._Fatigue = newFatigueModel

        
        return self._repaircost
        
    def getCrackProb(self):
        '''
        Returns the probability that a crack occurs
        
        This method returns the probability that a crack will be found on the
        fatigue detail in the current timestep.
        
        Parameters
        -----------
        None
        
        Returns
        -------
        self._ProbCrack: The probability a crack occured in this time step
        '''
        return self._ProbCrack
        
    def getRepairObject(self):
        '''
        Returns a reference to this T Panels repairCost object. 
        
        Parameters
        ----------
        None
        
        Returns
        -------
        A reference to the repairCost object associated with this
        instance of TransTPanelDet.
        '''
        return self._repaircost        
        
        
    
    
        
        
       
#    
#class HardCorner(deterioratingStructure):
#    
#    def __init__(self, RepairCost):
#        '''Follows base class constructor in assigning repair cost
#        '''
#            deterioratingStructure.__init__(self, RepairCost)
        
