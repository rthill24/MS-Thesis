'''
NSGA-II Implementation with optional variability-fidelity metamodeling


(c) 2011 The Regents of the University of Michigan
'''


import math
import copy
import os
import sqlite3
from random import Random
from operator import itemgetter, attrgetter
import datetime
import numpy as np
import kriging
import logging
import sys
import time
import mpi4py

try:
    from mpi4py import MPI
    entry = {'mpi_imported':True}
    globals().update(entry)
except ImportError:
    entry = {'mpi_imported':False}
    globals().update(entry)


def pretty_print(string, cNum):
    if cNum == 0:
        print(string)
        
        
__author__ = "Yi-Jen Wang, Dylan Temple, Matt Collette"
__copyright__ = "Copyright 2011, The Regents of the University of Michigan"
__license__ = "See license.dat in source distribution"
__status__ = "Development"


class Problem:
    """
    A base class explaining the attributes and  methods required to be present
    to use with NSGA-II optimizer 


    Attributes
    ----------
    metamodel:  Stores the metamodel currently in use for the problem. Does
                not need to be set by the user/derived class, will be written
                by NSGA-II
    
    GeneNum:    Integer, Number of genes in the chromosome for problem
    
    numObj:     Integer, number of objective functions
    
    numConstraints: Integer, number of constraints
    
    upBound:    List of floats, as long as GeneNum, listing the upper bounds 
                of each gene in the chromosome

    loBound:    List of floats, as long as GeneNum, listing lower bounds of 
                of each gene in the chromosome
    

    Methods
    -------
    The followings methods must be over-ridden in your derived problem class:
    __init__:    to copy data to base constructor
    
    Eval:        to handle function evalautions
    
    The following method is only required if you have constraints - otherwise
    base class will pass a "None/None" tuple back successfully     
    constraint:  to handle constraint evaluations
    
    The following methods are only required if make a VFO run:
    lowfidelity: handles low-fidelity calls in VFO
    
    highfidelity: handles high-fidelity calls in VFO


    """
    def __init__(self, numObj, numConstraints, GeneNum, loBound, upBound):
        self.numObj = numObj
        self.numConstraints = numConstraints
        
        if mpi_imported:
            self.comm = MPI.COMM_WORLD
            if self.comm.rank == 0:
                self.rank_flag = 0
            else:
                self.rank_flag = 1
        else:
            self.rank_flag = 0
        
        if self.rank_flag == 0:
            logger = logging.getLogger('msdl')         
            #Do some very basic error checking here
            if (len(loBound) != GeneNum):
                logger.error('NSGA-II problem class- number of genes and lower' \
                             'bound vector length mis-match')
                raise Exception('NSGA-II Problem length mis-match genes/loBound')
            if (len(upBound) != GeneNum):
                logger.error('NSGA-II problem class- number of genes and upper' \
                             'bound vector length mis-match')
                raise Exception('NSGA-II Problem length mis-match genes/loBound')
        
        self.GeneNum = GeneNum
        self.LoBound = loBound
        self.UpBound = upBound
        #Compatibility with old code - remove when fixed
        #TODO - remove old code
        self.LoBond = loBound
        self.UpBond = upBound
        
        #We will have no metamodel until NSGA-II sets one if requested, but
        #field must be populated
        self.metamodel = None 
        
    def Eval(self, individual_instance, metamodel=None):
        '''
        This function evaluates all the objective functions in the problem
        for an passed instantce of the individual class. 
        
        Parameters
        ----------
        individual_instance: An instance of the individual class, which defines
                             a particular location in design space via the 
                             genes in its chromosome
        
        metamodel:  A metamodel, currently only Kriging is supported. May be
                    left blank or passed as None if making a non-VFO run
       
        Returns
        -------
        Tuple of two lists, one is the objective function values for each
        objective function, second is the metamodel scaling values for each
        objective function (Note this is to support future revisions with
        mulitple meta-models, current code can only support one meta-model).
        Second list may be safely returned as None if not a VFO problem
        
        '''
        pass
    
    def constraint(self, individual_instance, metamodel=None):
        '''
        This function evaluates all the constraints in the problem
        for an passed instantce of the individual class. 

        Parameters
        ----------
        individual_instance: An instance of the individual class, which defines
                             a particular location in design space via the 
                             genes in its chromosome
        
        metamodel:  A metamodel, currently only Kriging is supported. May be
                    left blank or passed as None if making a non-VFO run
       
        Returns
        -------
        Tuple of two lists, one is the contraint violation values for each
        constraint function, second is the metamodel scaling values for each
        constraint function.  Contraints that are *not* violated should be 
        returned as zero, violations should be > 0
        (Note list of scaling factors is to support future revisions with
        mulitple meta-models, current code can only support one meta-model).
        Second list may be safely returned as None if not a VFO problem               
        '''
        return (None, None)

    def lowfidelity(self, individual_instance):
        '''
        This function evaluates the low-fidelilty version of the metamodel
        constraint or objective function. Only required for VFO runs,
        may be left as-is for non-VFO runs. 

        Parameters
        ----------
        individual_instance: An instance of the individual class, which defines
                             a particular location in design space via the 
                             genes in its chromosome
           
        Returns
        -------
        Low fidelity version of the constraint or meta-model evaluation.           
        '''        
        pass
    
    def highfidelity(self, individual_instance):
        '''
        This function evaluates the high-fidelilty version of the metamodel
        constraint or objective function. Only required for VFO runs,
        may be left as-is for non-VFO runs. 

        Parameters
        ----------
        individual_instance: An instance of the individual class, which defines
                             a particular location in design space via the 
                             genes in its chromosome
           
        Returns
        -------
        High fidelity version of the constraint or meta-model evaluation.           
        '''        
        pass

class Individual:
    def __init__(self,chromosome,prob,idnum,Generation,addupdist,rank):
        '''
        Store identities of each chromosome here
        '''
        self.prob=prob
        self.chromosome=chromosome
        self.id=idnum
        self.Generation=Generation
        self.addupdist = addupdist
        self.current=False
        self.current_const=False
        self.obj_correction = []
        self.const_correction = []
        self.rank = None #Rank of the individual
        self.crowding_distance = None #Crowding distance 
        self.violation_sum = 0 #Sum of all constraint violations
        
    def getobj(self):
        '''Updates/set the obejective function vector and any VFO correction
        vector for the objective functions if not already set
        '''
        if self.current==False:
            (self.__obj,self.obj_correction)=  \
                      self.prob.Eval(self, self.prob.metamodel)
            self.current=True
        return self.__obj
    def violation(self):
        '''Updates/sets the constraint function vector and any VFO correction
        terms for the constraints if not already set
        '''
        if (self.current_const==False):
            (self.constraint, self.const_correction) = \
                   self.prob.constraint(self, self.prob.metamodel)
            #It is possible for problems to have no constraints in which
            #case a sum of 0 and an empty list is required 
            if (self.constraint == None):
                self.violation_sum = 0
                self.constraint = []
            #Otherwise, the sum is just the sum of the individual violations
            else:
                self.violation_sum = math.fsum(self.constraint)
            self.current_const=True
        return self.violation_sum
    def getid(self):
        return self.id
    def getGeneration(self):
        return self.Generation
    def getchromosome(self):
        return self.chromosome
    def getVFOCorrections(self):
        '''Returns a tuple of two lists - objective function VFO scale factors
        and constraint function VFO scale factors
        '''
        #Make sure we are up-to-date
        self.violation()
        self.getobj()
        
        
        return (self.obj_correction, self.const_correction)
        
    def dom(self,i2):
        '''
        return 1 if this individual dominate i2
        or 0 otherwise
        '''
        if self.violation()<i2.violation():
            return 1
        elif self.violation()>i2.violation():
            return 0
        else:
            flag=1
            o1=self.getobj()
            o2=i2.getobj()
#            print "01", o1
#            print "o2", o2
            for i in range(0,len(o1)):
                if o1[i]>o2[i]:
                    flag=0
                elif o1[i]<o2[i]:
                    flag*=2
                '''
                Comparing two chromosome, if one of the chromosome has at least one
                object value less than the other chromosome, then it's nondominated
                '''
            if flag>1:
                return 1
            else:
                return 0
    def checkdom(self,i2):
        '''
        return 1 if self dominate i2 (self is better)
        return 0 if self and i2 are nondominated 
        return -1 if i2 dominate self (i2 is better)
        '''
        if self.dom(i2)==1:
            return 1
        elif i2.dom(self)==1:
            return -1
        else:
            return 0        
    def __repr__(self):
        return repr(self.chromosome)
      
class Optimizer:
    #Todo add database documentation
    def __init__(self, prob,  pc=0.8, sbx_ex = 4.0,
                 p_mut = 0.001, mut_ex = 2.0, use_metamodel=False, offset=-1, 
                 density=-1,spacing=-1, max_correction = 10.e12,
                 krig_theta_guess = 10, tol = 0.00, thread=False):
       
        '''Initializes an NSGA-II optimizer, either with or without a VFO model
        for either a constraint or objective function.  Currently real-coded
        only
        
        Based on the NSGA-II algorithm proposed by Deb in 
        K. Deb, A. Pratap, S. Agarwal, and T. Meyarivan, 
        "A fast and elitist multiobjective genetic algorithm: NSGA-II,"
        IEEE Transactions on Evolutionary Computation, vol. 6, no. 2, pp. 
        182-97, Apr. 2002. Crossover and mutation are based on Deb's SBX
        framework with user-specified exponents presented in:
        K. Deb, Multi-objective optimization using evolutionary algorithms. 
        Chichester: John Wiley & Sons, 2001.


        
        Parameters
        ----------
        prob:  An implementation of the problem class - see problem class
               documentation for more details
        
        pc: Probability of crossover per gene in chromosome during 
            reproduction.  Defaults to 0.8
        
        sbx_ex: Exponent on the SBX crossover method, default to 4.0
    
        p_mut: Probability of mutation per crossed-over gene. Defaults to 
                0.001

        mut_ex: Exponent on the mutation operator, default to 2.0
        
        use_metamodel: True/false if metamodel should be used. Default to false.
        
        offset: Only for VFO runs, must be >1, generation number when metamodel
                will first be constructed
        
        density: Must be > 1, number of individuals to choose from population
                 to update metamodel
        
        spacing: spacing of metamodel updates, in generations
        
        max_correction: maximum multiplictive correction factor for kriging
                        model Useful if low-fidelty model either fails or 
                        returns zeros.  Points above this factor will be dropped
                        before Krigning model build. Default 1.0e12
                        
        krig_theta_guess:  Initial guess for the Kriging theta value when 
                           solving the Kriging model. supplied vector is the
                           unit vector time this factor.
        
        tol:    Independent variable space normalized tolerance for VFO 
                kringing models.  Will not allow points to be added if within
                a normalized tol radius of existing points. Defaults to 0.00
                which prohibts only strict duplicate points. Normalized by
                (high bound - low bound) for each dimension from problem class
        
        thread: A flag indicating whether or not the optimizer will be run across
                multiplt processors in parallel
        
        Returns
        -------
        No return value
        
        '''
        
        self.thread = thread
        self.prob=prob
        self.use_metamodel = use_metamodel
        self.pc = pc
        self.sbx_ex = sbx_ex
        self.p_mut = p_mut
        self.mut_ex = mut_ex
        self.metamodel = None   #This is the instatiated metamodel -set later
        self.offset = offset
        self.density = density
        self.spacing = spacing
        self.metamodel_count = 0
        self.current_generation=0
        self.max_correction = max_correction
        self.krig_theta_guess =  krig_theta_guess
        self.TOL = tol
        
        self._krigingpts = []  #List of all individuals in Kriging metamodel
        self._krigingx = []   #List of list of independent variable for kmm
        self._krigingy = []   #List of dependent variable for kmm
    
    def __createGATables(self, runtitle):
        '''
        creates basic GA tables for further writing later
        numgenes is the number of genes in the chromosome
        '''
        
        if self.prob.rank_flag == 0:
            #Create initial informational table
            self.cur.execute('''CREATE TABLE description(title, value)''')
            #Create individual tables        
            self.cur.execute( '''CREATE TABLE individuals (IndID,born,died)''' )
            #create chromosome tables
            self.cur.execute( '''CREATE TABLE chromosomes (IndID,ChromID,Value)''' )
            #Create objective function value table
            self.cur.execute( 
              '''CREATE TABLE objfun (IndID,ObjID,Generation,Value,Correction)''' )
            #Creat constraint table
            self.cur.execute(
              '''CREATE TABLE constraints (IndID,ConstID,Generation,Value,Correction)''')
            #Create generation table
            self.cur.execute( 
               '''CREATE TABLE generation (GenNum, IndID, Front, CrowdDist)''' )  
            #Create a list of metamodels
            self.cur.execute('''CREATE TABLE metamodels(modelID, Generation)''')
            #Create a list of metamodel underlying data  - must retrieve chromosome
            #from chromosomes table 
            self.cur.execute('''CREATE TABLE mmData(modelID, IndID, Value)''')
            #Create a list of metamodel parameters           
            self.cur.execute(
                          '''CREATE TABLE mmParameters(modelID, ParamID, Value)''')
        
            #store key values - starting with run title
            self.cur.execute('insert into description values (?,?)', ('run title',
                        runtitle))
            #Time of execution
            self.cur.execute('insert into description values (?,?)', ('run time',
                        datetime.datetime.now().isoformat()))
            #
            if self.metamodel == None:
                self.cur.execute('insert into description values (?,?)', ('Metamodel',
                        'None'))
            else:
               self.cur.execute('insert into description values (?,?)', ('Metamodel',
                        self.metamodel.getType())) 
                        
            self.conn.commit()
        
    
    def __resetPop(self, pop):
        '''Simple function that resets both constraints and objective funtion
        to unevaluated upon meta-model update
        '''
        for indiv in pop:
            indiv.current=False
            indiv.current_const=False
    
    def __getVFOSamples(self, sample_pts):
        '''Goes to the problem class and gets both the high and low fidelity
        variants needed for the metamodel.  Returns a tuple of two numpy 
        arrays, one is the chromosomes, one for the ratio of high/low fidelity,
        and a modified sample_pts, with any pts above max_correction thrown
        out
        '''
        #Will be a Numpy array for the chromosome
        x_val_list = None
        #Will be a Numpy array for high/low fidelity
        y_val_list = []
        return_sample_pts = []
        
        
        kount = 0
        for indiv in sample_pts:
            #Get high/low from problems
            low_fi = self.prob.lowfidelity(indiv)
            high_fi = self.prob.highfidelity(indiv)
            #This is Metamodel term - 
            y = high_fi/low_fi
            if y < self.max_correction:
                if (x_val_list == None):
                    x_val_list = np.array(indiv.chromosome)
                    y_val_list.append(y)
                else:
                    temp_x = np.array(indiv.chromosome)              
                    x_val_list = np.vstack((x_val_list, temp_x))
                    y_val_list.append(y)
                kount += 1
                return_sample_pts.append(indiv)
        
        y_val_list = np.array(y_val_list)
        y_val_list = np.transpose(y_val_list)  
  
        return (x_val_list, y_val_list, return_sample_pts)
    
    def __formModel(self, xpts, ypts, sample_pts):
        '''Forms a new metamodel.
        '''
        logger = logging.getLogger('msdl')
        logger.info("In Form Model")
        self.metamodel_count += 1

        #Thought needs to be put in for how to make this generic later
        num_sita = self.prob.GeneNum
        sita0 = np.ones((num_sita))*self.krig_theta_guess 

        #Setup and solve the model
        
        ypts = np.transpose(np.matrix(ypts))
        self.prob.metamodel=kriging.Kriging(xpts,ypts, sita0)
        logger.info("About to solve Kriging model")
        theta0=self.prob.metamodel.solve()
        logger.info("Kriging model solve, theta vector " + str(theta0))
        
        #Now do the database updates
        metamodels = (self.metamodel_count, self.current_generation)
        self.cur.execute('insert into metamodels values (?,?)', metamodels)
        for kount in range(0, len(sample_pts)):
            data_tuple = (self.metamodel_count, sample_pts[kount].getid(),
                          ypts[kount,0])
            self.cur.execute('insert into mmData values (?,?,?)', data_tuple)
        kount = 0
        for value in theta0:
            kount += 1
            theta_tuple = (self.metamodel_count, kount, value)
            self.cur.execute('insert into mmParameters values (?,?,?)', 
                             theta_tuple)
        self.conn.commit()
        return

    def __writeChromsomestoDB(self,pop):        
        '''Writes a list of individuals/chromosome to the DB.  Does not have 
        to be sorted, only tables individuals and chromosomes will be touched
        death generation set to none initially'''
        if self.prob.rank_flag == 0:
            for indiv in pop:
                dbwrite = (indiv.id, self.current_generation, None)
                self.cur.execute('insert into individuals values (?,?,?)', dbwrite)
                kount = 0
                for gene in indiv.chromosome:
                    kount += 1
                    dbwrite = (indiv.id, kount, gene)
                    self.cur.execute('insert into chromosomes values (?,?,?)', 
                                     dbwrite)
            self.conn.commit()
        return
    
    def __writeGenObjFunctiontoDB(self,pop):
        '''This function writes just the objective functions/constraints to 
           the objfun and constraints table in the database - population does 
           *NOT* need to be sorted at this time
        '''
        if self.prob.rank_flag == 0:
            for i in pop:
                #Get obj values
                i.violation()
                if i.violation_sum == 0:
                    obj_vals = i.getobj()
                else:
                    obj_vals = []
                    for j in range(self.prob.numObj):
                        obj_vals.append('Infeasible')
                
                #Get violations individually - first call sum to make sure is 
                #current, then pull straight from class 
                constraints = i.constraint
                if self.use_metamodel:
                    (obj_corr, const_corr) = i.getVFOCorrections()
                else:
                    obj_corr = None
                    const_corr = None

                for k in range(0,self.prob.numObj):
                    #Make up the database write statement- Fields are
                    #IndID,ObjID,Generation,Value,Correction
                    if (obj_corr == None):
                        mm_corr = 1.0
                    else:
                        mm_corr = obj_corr[k]
                    dbwrite = (i.id, k+1, self.current_generation, obj_vals[k],
                                           mm_corr)
                    self.cur.execute('insert into objfun values (?,?,?,?,?)', 
                                         dbwrite)
                for j in range(0,self.prob.numConstraints):
                    #Make up the database write statement- Fields are
                    #IndID,ConstID,Generation,Value,Correction  
                    if (const_corr == None):
                        mm_corr = 1.0
                    else:
                        mm_corr = const_corr[k]                
                    dbwrite = (i.id, j+1, self.current_generation, constraints[j],
                                           mm_corr)
                    self.cur.execute('insert into constraints values (?,?,?,?,?)', 
                                         dbwrite)
            self.conn.commit()
        return
        
    def __writeSortedPoptoDB(self,pop):
        '''This function write a sorted population to the generations
        table in the database.  Popoulation should be sorted and have
        crowding distance metric done
        '''
        if self.prob.rank_flag == 0:
            for i in pop:       
                    dbwrite = (self.current_generation, i.id, i.rank, 
                                i.crowding_distance)
                    self.cur.execute('insert into generation values (?,?,?,?)', 
                                         dbwrite)
            self.conn.commit()
        return 
        
    
    def __okToAdd(self, indiv, TOL, NormalVector):
        '''
        Checks if a point is within a normalized tolerance TOL of an existing
        point in the Kriging model point list.  The normalization is 
        performed by dividing the distance between each gene in the 
        chromosome by the corrosponding entry in the list NormalVector. 
        Returns true if OK to add, or false if point is within TOL
        '''
        for existingPt in self._krigingpts:
            total = 0
            #Add up the distance normalized
            for k in range(0, len(NormalVector)):
                    total += ((indiv.chromosome[k] - existingPt.chromosome[k])/
                                    NormalVector[k])**2.
            total = total**0.5
            #As soon as any point fails, abort and return with False
            if total < TOL:
                return False
        
        #If we make it through the full list, point passes
        return True
        
    def __selKrigingPts(self, f, first=False):
        '''
        Returns a list of points to sample for the Kriging model.  May
        be the first Kriging model (first=True) in which case the list will be
        randomly shuffled before selection (to avoid any correlation from 
        crossover) or an update (first=False)
        in which case the points will be selected by predicted error.  If 
        the selection runs out of potential points before hitting the target
        density, a warning will be printed to the logger.info channel, but 
        no other errors will be generated. f is the list of sorted fronts
        from the main NSGA loop
        
        NOTE! This function updates self._krigingpts as it goes so that the 
        proper distance check can be made.
        '''
        sample_pts = []  #List of individuals to sample at
        curr_front = 0
        indiv_still_required = self.density
        #Loop over fronts filling as we go
        while ((indiv_still_required > 0) and (curr_front < len(f))):           
            #Length of current front and counter 
            lcurrFront = len(f[curr_front])
            i = 0
            #Make either a random list or sorted list 
            if (first):
                sel_list= self.rnd.sample(f[curr_front], lcurrFront)
            else:
                #Compute COV at every pt
                krig_cov=[]
                for indiv in f[curr_front]:
                    xpts = np.matrix(indiv.chromosome)
                    (yp,mse)=self.prob.metamodel.predictor(xpts)
                    krig_cov.append([((mse[0,0]**0.5)/yp[0,0]),indiv])
                #Sorted list of errors
                #reverse = False means sorting cov from low to high 
                sorted_error = sorted(krig_cov, key = itemgetter(0),
                                              reverse = True)
                #Go through and make the selection list - this is ugly
                sel_list = []
                for entry in sorted_error:
                    sel_list.append(entry[1])
            #Go through the list, if points are OK, then add them
            while ((i < lcurrFront) and (indiv_still_required > 0)):
                indiv = sel_list[i]
                if (self.__okToAdd(indiv, self.TOL, self.NormVector)):
                    indiv_still_required -= 1
                    sample_pts.append(indiv)
                    self._krigingpts.append(indiv)
                i += 1
            curr_front += 1
            
        #If we didn't get enough points - issue a warning
        if (len(sample_pts) < self.density):
            print('Warning - VFO unable to completely fill update')
            print(('Number of points left ', indiv_still_required))
        
        #Return the list of sample points 
        return sample_pts
                
    
    
    def __removeExistingIndividuals(self, sorted_fronts):
        '''This function goes through all the fronts in a sorted_fronts 
        structures (the one returend from non-dominated sorting function and
        checks if any individual is already in the class Kringing pts master
        list, if so it removes the point so that duplicate points don't get
        added to Kriging metamodel '''
        
        no_krig_duplicates = []   #This will be filled in as return structure        
        
        for front in sorted_fronts:
            new_front = []
            for individual in front:
                if self._krigingpts.count(individual) == 0: #not in list yet
                    new_front.append(individual)
            if len(new_front) > 0:
                no_krig_duplicates.append(new_front)
        return no_krig_duplicates
 
   
    def run(self, dbname, runtitle, seed, num_generations, pop_size):
 
        '''Start a run of the NSGA-II algorithm, either with or without 
        metamodel updates. All runs start from scratch. 
        
        Parameters
        ----------
        dbname:  Name of sqlite database file to write with the results
        
        runtitle: String title of run, will be copied into database
        
        seed:  Random number generator seed, used for tournmanent selection,
               crossover, etc. Runs with identical settings and seed should
               always be identical
               
        num_generations: Number of generations to run
        
        pop_size: Desired population size, should be even for 2-pass tournament
        
        Returns
        -------
        No return value
        
        '''       
        
        #@TODO - work in pop_size in place of problem - store on class
     
        '''
        Perform initial setup
        '''
        if self.prob.rank_flag == 0:
            logger = logging.getLogger('msdl')
        self.rnd=Random(seed)
        self.num_generations = num_generations
        if pop_size%2 != 0:
            if self.prob.rank_flag == 0:
                logger.error("NSGA II passed non-even pop size")
                raise Exception('NSGA II pop error')
        self.pop_size = pop_size        
        #Set up the database, overwriting if required and store key values
        if self.prob.rank_flag == 0:
            if (os.path.exists(dbname)):
                os.remove(dbname)
            self.conn = sqlite3.connect(dbname)
            self.cur = self.conn.cursor()
        self.__createGATables(runtitle)
        
        '''
        Create the first randomly-distributed population
        '''
        pretty_print('GENERATION '+str(self.current_generation), self.prob.rank_flag)
        self.id=0
        self.addupdist = 0
        self.current_generation=0
        pop=[]          
        #First population        
        for x in range(self.pop_size):
            self.id = self.id+1
            self.rank = 0
            self.addupdist = 0
            #Initialize distance
            chromosome = []
#            UpBound = [0.75,5.0]
#            LoBound = [0.39,0]
            for y in range(self.prob.GeneNum):
                chromosome.append(self.rnd.uniform(self.prob.UpBound[y],
                                                   self.prob.LoBound[y]))
            ind = Individual(chromosome,self.prob,self.id,
                             self.current_generation, 
                             self.addupdist, self.rank)
            pop.append(ind)
        #sort and calculate crowding distance    
        pop = self.evaluate_pop_parallel(pop)


             
        if self.prob.rank_flag == 0:
            f = self.fast_nondominated_sort(pop)
            self.__writeGenObjFunctiontoDB(pop)
            self.crowding_distance(pop, f)


        '''
        Store the first generation of Individual and it's indentities into 
        individuals, chromosome, and objective function, and generation database
        '''
        self.__writeChromsomestoDB(pop)
        self.__writeGenObjFunctiontoDB(pop)
        self.__writeSortedPoptoDB(pop)
        
        if self.prob.rank_flag == 0:
            logger.info('Generation:'+ str(self.current_generation))     
            k=0
            for front in f:
                logger.debug('front number: ' + str(k))
                k+=1
                for i in front:
                    logger.debug('id ' + str(i.id) + 'chromosome ' +
                                    str(i.chromosome))
        
        '''
        Run through the required number of generations
        '''
        for self.current_generation in range(1,self.num_generations+1):
            if self.current_generation > 1:
                pretty_print('TIME TO FOR GENERATION '+str(self.current_generation)+' '+str(fTime-sTime), self.prob.rank_flag)
            sTime = time.time()
            pretty_print('GENERATION '+str(self.current_generation), self.prob.rank_flag)
            if self.prob.rank_flag == 0:
                logger.info('Generation:'+ str(self.current_generation))   
            
            #See if we have to do a metalmodel update
            if (self.use_metamodel) and \
                  (self.current_generation == self.offset):
                logger.info('Setting up initial metamodel')          
                
                #Need a normalized distance vector for checking duplicates/
                #tightly-spaced points
                self.NormVector = []
                for i in range(0, len(self.prob.LoBound)):
                    self.NormVector.append(self.prob.UpBound[i] -
                                           self.prob.LoBound[i])
                              
                #First flag all members of the population as stale
                self.__resetPop(pop)
                
                sample_pts = self.__selKrigingPts(f, True)                
                
                xpts, ypts, sample_pts = self.__getVFOSamples(sample_pts)
                #And form a model
                self.__formModel(xpts, ypts, sample_pts)
                #Record the individuals used so we can append to them later
                #selKrigingPts takes care of list on self (krigingpts)
                self._krigingx = xpts
                self._krigingy = ypts
                
                #And re-sort, re-scale the population with new metamodel
                #NOTE: These ranks are not written to the DB - only used
                #for making children - 2-pass tournament needs rank and 
                #crowding distance metrics! Objective function values
                #will be stored in DB, however, for all members 
                f = self.fast_nondominated_sort(pop)
                self.crowding_distance(pop, f)      
            elif ((self.use_metamodel) and 
                    (self.current_generation > self.offset) and 
                    ((self.current_generation-self.offset)%self.spacing==0)):
                
                logger.info('Refining metamodel')                   
                #First flag all members of the population as stale
                self.__resetPop(pop)
                
                sample_pts = self.__selKrigingPts(f, False)
                
                #Get just the new high-fidelity runs
                new_xpts, new_ypts, sample_pts = \
                                    self.__getVFOSamples(sample_pts) 
                
                #update our list for the kriging model
                self._krigingx = np.append(self._krigingx, new_xpts, axis=0)
                self._krigingy = np.append(self._krigingy, new_ypts, axis=0)      

                #And form a model
                self.__formModel(self._krigingx, self._krigingy, 
                                             self._krigingpts)

                #And re-sort, re-scale the population with new metamodel
                #NOTE: These ranks are not written to the DB - only used
                #for making children - 2-pass tournament needs rank and 
                #crowding distance metrics! Objective function values
                #will be stored in DB, however, for all members 
                f = self.fast_nondominated_sort(pop)
                self.crowding_distance(pop, f)                                     
            else:
                if self.prob.rank_flag == 0:
                    logger.info('No model or no model update')
            
            #Write the parents objective functions to the DB
            self.__writeGenObjFunctiontoDB(pop)
            
            #Make the children
            child = []
            if self.prob.rank_flag == 0:
                child = self.two_pass_tournament_selection(pop)               
            child = self.evaluate_pop_parallel(child)
	    
            #store the new create child
            self.__writeChromsomestoDB(child)
            self.__writeGenObjFunctiontoDB(child)    
            #combine the current parents and children
            if self.prob.rank_flag == 0:
                pop = pop+child
            else:
                pop = None
            
            #Run fast_nondominated_sort and crowding distance to delete worse 
            #individual
            if self.prob.rank_flag == 0:
                f = self.fast_nondominated_sort(pop)
                self.crowding_distance(pop, f)

            '''
            fill the next generation by sort
            '''
            #Make a sorted vector by rank, then crowding distance, reject all
            #entities that are over the rank
            Identities = []
            Finalpop = []
            if self.prob.rank_flag == 0:
                for i in pop:
                    if i.crowding_distance != 'Infeasible':
                        Ind = [i.rank, -i.crowding_distance, i.id, i]
                    elif i.crowding_distance == 'Infeasible':
                        Ind = [i.rank, i.crowding_distance, i.id, i]
                    Identities.append(Ind)
                newpop = sorted(Identities, key = itemgetter(0,1),reverse = False)
                
                #Log the sorted list to the log file for detailed checking
                logger.debug('Sorted population list for generation' + 
                                    str(self.current_generation))
                logger.debug('Data is: rank, -crodwing distance, id')
                for i in newpop:
                    logger.debug( str(i[0]) + ',' + str(i[1]) + ',' + str(i[2]))
            #population is now twice as long as needed
            for kount in range(0, 2*self.pop_size):
                if kount >= self.pop_size:
                    if self.prob.rank_flag == 0:
                        ID = newpop[kount][2]
                        logger.debug('Removing individual ' + str(ID))
                        #store when the indivdual die 
                        self.cur.execute(
                         "UPDATE individuals SET died = ? WHERE IndID = ?", 
                         (self.current_generation , ID))
                else:
                    if self.prob.rank_flag == 0:
                        Finalpop.append(newpop[kount][3])
                    else:
                        Finalpop = []
            if self.prob.rank_flag == 0:
                self.conn.commit()            
            #Run fast_nondominated_sort  and crowding distance
            #@TODO - replace with a more efficient scheme - individual ranks
            #are still known and correct, just re-form f without doing the 
            #full sort 
            pop = Finalpop
            if self.prob.rank_flag == 0:
                f = self.fast_nondominated_sort(pop)
                self.crowding_distance(pop, f)
            
            self.__writeSortedPoptoDB(pop)
            fTime = time.time()
    
    def evaluate_pop_parallel(self,pop):
        '''
        This method retrieves the objective values and constraint values for each
        individual in the population.
        
        This method uses the mpi4py module to scatter a popluation across multiple
        processors and evaluate the objective function value and constraint violation
        in parallel.  The fully updated individuals are then recombined into a single
        population on the root machine.
        
        Parameters
        ----------
        pop - A list of instances of the 'individual' class to be updated
        
        Returns
        -------
        pop_indis (root processor) - A list of individuals with fully updated 
                                    objective values and constraint violations
                                    
        None (remote processors) - A place holder empty object
        '''
        if self.thread:
	    
            #Set up local array sizes for each processor
            num_2_update = len(pop)
            local_indis = math.ceil(num_2_update/self.prob.comm.size)
            local_indis = int(local_indis)
            local_size = (self.prob.comm.size, local_indis)
            
            #Initialize lists on each processor
            if self.prob.comm.rank == 0:
                com_indis = np.empty(local_size, dtype=object)
                k = 0
                for i in range(local_indis):
                    for j in range(self.prob.comm.size):
                        try:
                            #Fill empty spaces with objects
                            com_indis[j,i] = pop[k]
                            k += 1
                        except IndexError:
                            continue
            else:
                #Give all non-root processors a place holder
                com_indis = None

            #Scatter the data across all the processors
            com_indis = self.prob.comm.scatter(com_indis, root=0)
            self.prob.comm.Barrier()
            
            #Evaluate the objective function and constraint for each 
            #individual on all processors
            for i in com_indis:
                if type(i).__name__ == 'instance':
                    i.violation()
                    if i.violation_sum == 0:
                        i.getobj()
                else:
                    print(type(i).__name__)
                    #Used for unused array indices on processors
                    continue
                
            #Gather all the individuals on the root processor
            com_indis = self.prob.comm.gather(com_indis)
            now_updated = []
            self.prob.comm.Barrier()
            
            #Place all newly updated individuals in list with no empty places
            if self.prob.comm.rank == 0:
                for i in com_indis:
                    try:
                        for j in i:
                            now_updated.append(j)
                    except TypeError:
                        #Do not keep empty array entries
                        continue
                
                pop = now_updated
                #Ensure new population is only comprised of individuals on root
                #machine
                pop_indis = []
                for i in pop:
                    if not i:
                        continue
                    else:
                        pop_indis.append(i)
                #return fully updated population on root processor
                return pop_indis
            else:
                #Return place holder on remote processors
                return None
        else:
            return pop
    
    def slow_nondominated_sort(self,pop):
        self.pop=pop
        for i in range(0,len(pop)):
            print('\n"',i,'"Front')
            nonID=[]
            ID=[]
            con=0
            for i in pop:
                print('\n',i.getid(),'\'th Chromosome is:' ,i.getchromosome())
                con = con+1
                pop2=copy.copy(pop)
                pop2.remove(i)
                #remove chromosome itself when comparing dominating 
                counta=countb=countc=count=0
                for j in pop2:
                    if i.checkdom(j) == 1:
                        counta=0
                    elif i.checkdom(j)==0:
                        countb=0
                    else:
                        countc=1
                    count=counta+countb+countc
                if count >=1:
                    print('dominated')
                else:
                    print('nondominated')
                    nonID.append(i.id)
                    ID.append(con)
                    #In order not to ruin the oringianl ID creat another  ID list
            print('\nNondominated individual are:',nonID)
            ID.sort()
            ID.reverse()
            for i in ID:
                j = i-1
                print('\nNondominated Chromosome are:',self.pop.pop(j))
    #Careful! pop out all the individual from the population
    def fast_nondominated_sort(self, P):
        '''
        NSGA2 non-dominated sorting, based on Michal Fita 2007 (py-nsga2)
        Apache License, Version 2.0
        Notations match Deb IEEE TRANSACTIONS ON EVOLUTIONARY, VOL.6,NO.2,APRIL.
        2002
        Note - you should *NOT* write other code using single-letter variable 
        names - always use descriptive names.  The names used here match the 
        reference for consistency 
        '''
        F = []
        F.append([])
        S = dict()
	
        for p in P:
            S[hash(p)] = []
            p.n = 0
            for q in P:
                if q.id != p.id:
                    domflag = p.checkdom(q)
                    if domflag == 1:  #p dominates q
                        S[hash(p)].append(q)
                    elif domflag == -1: #q doominates p 
                        p.n += 1
            if p.n == 0:
                p.rank = 1
                F[0].append(p)
        i = 0
        while len(F[i]) != 0:
            Q = []
            for p in F[i]:
                for q in S[hash(p)]:
                    q.n -= 1
                    if q.n == 0:
                        q.rank = i + 1
                        Q.append(q)
            i += 1
            F.append([])
            F[i] = Q
        for i in range(0,len(F)-1):
            for j in F[i]:
                j.rank = i
        #store back the ranks
        return F
        
        
    def sbx(self,x1,x2):
        '''
        only do the chromosome crossover
        '''
        Pc = self.pc
        eta_c = self.sbx_ex
        X1=[]
        X2=[]
        for i in range(0,len(x1)):
            pc = self.rnd.random()
            if pc >= Pc:
                a = x1[i]
                b = x2[i]
            else:
                u=self.rnd.random()
                if u<0.5:
                    Bq = (2*u)**(1/(eta_c+1))
                else:
                    Bq = (1/(2*(1-u)))**(1/(1+eta_c))
                a=0.5*((1.0+Bq)*x1[i]+(1.0-Bq)*x2[i])
                b=0.5*((1.0-Bq)*x1[i]+(1.0+Bq)*x2[i])
                #control the children not to exceed the boundary
            X1.append(a)
            X2.append(b)
        for i in range(self.prob.GeneNum):
            if X1[i] > self.prob.UpBound[i]:
                X1[i] = self.prob.UpBound[i]
            elif X1[i] < self.prob.LoBound[i]:
                X1[i] = self.prob.LoBound[i]
            if X2[i] > self.prob.UpBound[i]:
                X2[i] = self.prob.UpBound[i]
            elif X2[i] < self.prob.LoBound[i]:
                X2[i] = self.prob.LoBound[i]
        return [X1,X2]
        
    def mutation(self,mut):
        '''
        do the mutation of chromosome after the crossover is done
        '''
        n = self.mut_ex
        for i in mut:
        #when passing in mutation, two individuals form a pair in each list
        #entry, so unpack them into j 
            for j in i:
                for k in range(0,len(j)):
                #each individual has multi chromosome
                    Max_mut = (self.prob.UpBound[k]-self.prob.LoBound[k])
                    pm = self.rnd.random()
                    u = self.rnd.random()
                    if pm >= 1.00-self.p_mut:
                    #1% of chance to do the mutation
                        if u < 0.5:
                            delta = (2*u)**(1/(1+n))-1
                        else:
                            delta = 1-(2*(1-u))**(1/(n+1))
                    else:
                        delta = 0
                    j[k] = j[k] + delta*Max_mut
                    #control the children not to exceed the boundary
                for i in range(self.prob.GeneNum):
                    if j[i] > self.prob.UpBound[i]:
                        j[i] = self.prob.UpBound[i]
                    elif j[i] < self.prob.LoBound[i]:
                        j[i] = self.prob.LoBound[i]
        return mut
    
    
    def crowding_distance(self, pop, f):
        '''Perfoms the Deb (2002) crowding distance calculation on the provided
        population, storing results on individual.crowding_distance
        
        
        Parameters
        ----------
        pop:  Complete population
        
        f: Population sorted into fronts 
       
        Returns
        -------
        No return value
        
        '''        

        logger = logging.getLogger('msdl')        
        logger.info('Starting crowding distance calc')
        
        #num_objs = len(pop[0].getobj())
        num_objs = self.prob.numObj

        #for front in f:
        #Does not work as f seems to have a trailing front with length 0 at the 
        #end. Other code may rely on this so, sort not fixed
        # TODO - fix trailing zero-length front and find all code that relies
        # on non-standard front indexing and update 
        for front_kount in range(0,len(f)-1):
            front = f[front_kount]
            logger.debug("Working on front" + str(front[0].rank))
            #This will be a list-of-lists, column 0 is the individual class
            #column 2...n are the objective function values
            front_list = []
            
            #First make up the front list
            for indiv in front:
                    if indiv.violation_sum == 0:
                        #Reset crowding distance
                        indiv.crowding_distance = 0
                        front_member = [indiv]
                        #Add objective function values
                        front_member.extend(indiv.getobj())
                        #Add to the big list
                        front_list.append(front_member)
                        logger.debug("Appended to front list" + str(front_member))
            
            #Figure out how long front_list is
            num_pts = len(front_list)
            #Now loop over each objective function, sorting the list from
            #min to max
            if front_list:
                for current_obj in range(1, num_objs+1):
                    #Sort from min to max
                    sorted_list = sorted(front_list, key = itemgetter(current_obj))
                    #We will normalize by the spread of values
                    min_val = sorted_list[0][current_obj]
                    max_val = sorted_list[-1][current_obj]
                    normalize_factor = max_val - min_val
                    #Check for a spread of zero - single pt front
                    if (normalize_factor == 0):
                        normalize_factor = 1
                    
                    #Log some values
                    logger.debug('Working on objective number: ' + str(current_obj))
                    logger.debug('min val this generation: ' + str(min_val))
                    logger.debug('max val this generation: ' + str(max_val))
                    logger.debug('Crowding distance denominator: ' + 
                                                 str(normalize_factor))  
                    logger.debug('Print sorted list' + str(sorted_list))
                    
                    #Go through and do the crowding distance calculation
                    #Ends get infinite values
                    sorted_list[0][0].crowding_distance += float("inf")
                    sorted_list[-1][0].crowding_distance += float("inf")
                    
                    #Everyone else gets one dimensions of the cubeoid 
                    for kount in range(1, num_pts-1):
                        sorted_list[kount][0].crowding_distance += math.fabs(
                            sorted_list[kount - 1][current_obj] -
                            sorted_list[kount + 1][current_obj])/math.fabs(
                            normalize_factor)
            
            for indiv in front:
                if indiv.violation_sum != 0:
                    indiv.crowding_distance = 'Infeasible'
                    
        return
    
    def crowding_distance_older(self, pop, f):
        #Start up the logging facilities
        logger = logging.getLogger('msdl')        
        logger.info('Starting crowding distance calc')
        Inf = 1e30000 #create an infinity number
        obj = []
        Indi2 = []
        g = len(f)
        #Total fronts in this pop is g-1
        d = len(pop[0].getobj())
        #number of objective function
        for j in range(0,g-1):
        #Calculate distance front by front
            logger.debug('Crowding distance for front: ' + str(j))
            h = j #copy the front number in order to be used later
            distance = [] #distance stores boundary distance in each front 
            k = len(f[j])
            #Number of objectives in each front
            for e in range(0,d):
                obj1 = []
                for i in f[j]:
                    obj1.append(i.getobj()[e])
                obj1.sort()
                min_val = obj1[0]
                max_val = obj1[len(f[j])-1]
                diff = max_val - min_val
                '''
                If same individual appear more than twice and happen to have 
                zero boundary distance, we need to set up a small number in 
                order to prevent the case 'divided by zero' happen
                '''
                if diff == 0:
                    distance.append(1)
                else:
                    distance.append(diff)
                logger.debug('Working on objective number: ' + str(e))
                logger.debug('min val this generation: ' + str(min_val))
                logger.debug('max val this generation: ' + str(max_val))
                logger.debug('Crowding distance denominator: ' + 
                                             str(distance[-1]))
            obj = []
            ID = []
            #Get IDS and objective functions values as a list
            for i in f[j]:
                obj.append(i.getobj())
                ID.append(i.id)
            for e in range(0,d):
            #calculate distance objective by objective
                obj2 = sorted(obj,key = itemgetter(e))
                Indi = []
                ID2 = []
                for i in f[h]: #h equals j at this point
                    Ind2 = [i.getobj()[e], i.id]
                    Indi.append(Ind2)  
                #Sort by current objective functions     
                Ind3 = sorted(Indi, key = itemgetter(0))
                for i in Ind3:
                    ID2.append(i[1])
                #ID after sortng is ID2                
                #sort obj by the different objective value
                for i in range(0,k): #k is number of individuals in front
                #calculate distance in each objective function
                    addupdist = 0
                    if i == 0:
                        id_val = ID2[0] 
                        addupdist = Inf
                    elif i == k-1:
                        id_val = ID2[k-1] 
                        addupdist = Inf
                        #The distances are infinity at the boundary points
                    else:
                        id_val = ID2[i]
                        left = obj2[i-1]
                        right = obj2[i+1]
                        distance2 = [abs(m - n) for m,n in zip(right,left)]
                        #distance between the neigbors of each chromosome is distance2
                        normalize = [m/n for m,n in zip(distance2,distance)]
                        #normalize the distance
                        for j in normalize:
                            addupdist = addupdist + j
                    Ind4 = [id_val, addupdist]
                    Indi2.append(Ind4)
                    Indi2.sort()
        Indi3 = []
        #add up different distance by different using objective value
        for j in range(1,len(pop)+1):
            j = d*j
            p = j-1
            aa = Indi2[p]
            distance3 = aa[1]
            for i in range(2,d+1):
                ji = j-i
                bb = Indi2[ji]
                distance3 = distance3 + bb[1]
            Ind5 = [aa[0], distance3]
            Indi3.append(Ind5[1])
        #The bigger the crowded distance means less crowded and is preferred        
        #store crowding distance back to pop
        x=[]
        Indi3.reverse()
        for i in pop:
            i.crowding_distance = Indi3.pop()
            x.append(i.crowding_distance)
        return x
        
    def Individual_comparison(self,x,y):
        '''
        input individual directly with all those identities
        comparing two individuals and return the best one base 
        on the ranking first then crowding distance
        '''
        if x.rank < y.rank or (x.rank == y.rank and x.crowding_distance  > y.crowding_distance):
            return x
        elif x.rank == y.rank and x.crowding_distance == y.crowding_distance:
        #If both chromosomes have same criterion, depends on random number
            p = self.rnd.random()
            if p >= 0.5:
                return x
            else:
                return y
        else:
            return y
    def two_pass_tournament_selection(self, pop):
        '''
        Using the comparism operator to do the selection
        combing with sbx and mutation, create a list of child
        '''
        #Selected members of the population for crossover
        chosen = []
        #List of children before mutation
        child = []
        #List of children after mutation
        mut = []
        for i in range(0,2):        
            #Do the tournament selection twice
            #shuffle the population first
            pop2 = copy.copy(pop)
            self.rnd.shuffle(pop2)
            for j in range(0,(len(pop2)//2)):
                x = pop2.pop()
                y = pop2.pop()
                #each time take two chromosome out to compare
                a = self.Individual_comparison(x,y)
                #chosen chromosome after doing the comparism
                chosen.append(a)
        #Now go through and make children
        for i in range(0,(len(chosen)//2)):
                x = chosen.pop()
                y = chosen.pop()
                b = self.sbx(x.chromosome, y.chromosome)
                #return the children chromosome - this is a list of two 
                #chromosomes
                mut.append(b)
        aftmut = self.mutation(mut) 
        for i in aftmut:
            for j in i:
                self.id = self.id + 1
                #since each bracket has two individual 
                Indi = Individual(j,self.prob,self.id,self.current_generation, 
                                                           0, 0)
                child.append(Indi)
        return child
