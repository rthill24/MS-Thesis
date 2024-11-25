### Unit test code for NSGA-II optimizer with variable fidelity option

### imports - unittest framework, nsga2_michigan_vfo module, random module, copy module, logging module,
### sqlite3 module, and math module
import unittest
from random import Random
import copy
import math
import sqlite3
import nsga2_michigan_threaded as nsga2
import os
import kriging
import numpy as np
import logging

class testNSGA2(unittest.TestCase):
    """
    Unit test framework for NSGA-II optimizer
    
    This class tests various methods of the NSGA-II algorithm with a variable fidelity option.  This class does not
    test the methods associated with the variable fidelity.  The data and calculations used to generate the assertions
    are taken from:
        Course notes for Dr. Matthew Collette's NAVARCH 570 class at the University of Michigan
        Example calculations in "Multi-Objective Optimization using Evolutionary Algorithms" by K. Deb
    
    Hand calculations supporting the assertions can be seen in various '_handCalcs' files.
    """


    def setUp(self):
        """
        Sets up the data for the tests, called before each test case

        This sets up the data for a unit test on the NSGA-II framework with a variable fidelity option.  It currently
        creates a test population that can be sorted into fronts or bred in order to test the various methods within
        the optimizer.  It also sets up a second set of test data that is used to test the crowding distance operator.

        Parameters
        ----------
        No parameters

        Returns
        -------
        No return value

        
        """

        ### Test population data [number, function 1 value, function 2 value]
        ### This test data is from Dr. Matthew Collete's NAVARCH 570 (Advanced
        ### Marine Design) class.
        self.function_table = [[1,5,7],
                          [2,7,4],
                          [3,2,7],
                          [4,3,8],
                          [5,10,2.5],
                          [6,1,9],
                          [7,7,10],
                          [8,3,6],
                          [9,8,4]]
                

        ### Second test population data (from hand calcs in Deb's "Multi Objective Optimization using 
        ### Evolutionary Algorithms).  Used for the crowding ditance unit test.
        self.population_data = [[1,0.31,0.89,0.31,6.10],
                                [2,0.43,1.92,0.43,6.79],
                                [3,0.22,0.56,0.22,7.09],
                                [4,0.59,3.63,0.59,7.85],
                                [5,0.66,1.41,0.66,3.65],
                                [6,0.83,2.51,0.83,4.23],
                                ['a',0.21,0.24,0.21,5.90],
                                ['b',0.79,2.14,0.79,3.97],
                                ['c',0.51,2.32,0.51,6.51],
                                ['d',0.27,0.87,0.27,6.93],
                                ['e',0.58,1.62,0.58,4.52],
                                ['f',0.24,1.05,0.24,8.54]]                



        ### Set up data in terms of 'Individual' class for use with NSGA-II algorithm
        self.testPop = []
        

        for i in range(1,10):
            ind = nsga2.Individual(self.function_table[i-1][0],0,
                                   self.function_table[i-1][0],0,0,0)
            
            self.testPop.append(ind)
            self.testPop[i-1].current = True
            self.testPop[i-1].current_const = True
            self.testPop[i-1].constraint = None
            self.testPop[i-1]._Individual__obj = [self.function_table[i-1][1],
                         self.function_table[i-1][2]]
        
        self.testPop2 = []
        
        for i in range(1,13):
            ind = nsga2.Individual(self.population_data[i-1][0],0,
                                   self.population_data[i-1][0],0,0,0)
            
            self.testPop2.append(ind)
            self.testPop2[i-1].current = True
            self.testPop2[i-1].current_const = True
            self.testPop2[i-1].constraint = None
            self.testPop2[i-1]._Individual__obj = [self.population_data[i-1][3],
                          self.population_data[i-1][4]]

    def test_fast_nondominated_sort(self):
        '''
        Test the fast_nondominated_sort method using the above test population

        This test uses the population data defined in in the setUp portion of the test and sorts it into fronts
        (minimizing both objective functions).  It test to make sure that the 'fast_nondominated_sort' method of the
        nsga2_michigan_vfo optimizer script properly sorts the data.  

        Parameters
        ----------
        No parameters

        Returns
        -------
        No return value
        '''
        
        ### Initiate instance of optimizer and run nondominated sorting method
        self.Opt = nsga2.Optimizer(0)
        self.fronts = self.Opt.fast_nondominated_sort(self.testPop)
        
        ### Assert that the proper individuals are in each front
        self.assertTrue(self.testPop[0] in self.fronts[1], \
                        "Individual 1 is in front 2")
        self.assertTrue(self.testPop[1] in self.fronts[0], \
                        "Individual 2 is in front 1")
        self.assertTrue(self.testPop[2] in self.fronts[0], \
                        "Individual 3 is in front 1")
        self.assertTrue(self.testPop[3] in self.fronts[1], \
                        "Individual 4 is in front 2")
        self.assertTrue(self.testPop[4] in self.fronts[0], \
                        "Individual 5 is in front 1")
        self.assertTrue(self.testPop[5] in self.fronts[0], \
                        "Individual 6 is in front 1")
        self.assertTrue(self.testPop[6] in self.fronts[2], \
                        "Individual 7 is in front 3")
        self.assertTrue(self.testPop[7] in self.fronts[0], \
                        "Individual 8 is in front 1")
        self.assertTrue(self.testPop[8] in self.fronts[1], \
                        "Individual 9 is in front 2")
        


    def test_minMax(self):
        '''tests the min-max internal function for adding to the VFO model'''
        
        #reference population        
        ref1 = nsga2.Individual([0.,0.],None,0,0,0,0)
        ref2 = nsga2.Individual([5.,5.],None,0,0,0,0)
        
        #simple initial population
        init1 = nsga2.Individual([1.,1.],None,0,0,0,0)
        init2 = nsga2.Individual([1.,2.],None,0,0,0,0) 
        init3 = nsga2.Individual([3.,3.],None,0,0,0,0)  
        
        reflist = [ref1, ref2]
        
        initlist = [init1, init2, init3]
       
        normalvect = [1., 1.]
        
        #optimizer
        opt = nsga2.Optimizer(None) 
        result = opt._Optimizer__minMax(reflist, initlist, normalvect)
        
        self.assertEqual(result, init3, "minMax simple case")

        #Another population for testing scaling 
        pop21 = nsga2.Individual([0.1,10.],None,0,0,0,0)        
        pop22 = nsga2.Individual([0.2,0.1],None,0,0,0,0)  
        pop2 = [pop21, pop22]
        result = opt._Optimizer__minMax(reflist, pop2, normalvect) 
        self.assertEqual(result, pop21, "minMax  case 2 ")

        shrinkvec = [1., 1000.]
        result = opt._Optimizer__minMax(reflist, pop2, shrinkvec)
        self.assertEqual(result, pop22, "minMax  case 3 - shrink vector")       

        return
        
    def test_newFrontSelection(self):
        '''test the new first-front selection for the Kriging model'''
        front1_1 = nsga2.Individual([0.,0.],None,0,0,0,0)
        front1_2 = nsga2.Individual([5.,5.],None,0,0,0,0)
        front1 = [front1_1, front1_2]
        
        front2_1 = nsga2.Individual([0.,0.],None,0,0,0,0)
        front2_2 = nsga2.Individual([1.,20.],None,0,0,0,0)
        front2_3 = nsga2.Individual([1.,1.],None,0,0,0,0)
        front2_4 = nsga2.Individual([5.,1.],None,0,0,0,0)
        front2_5 = nsga2.Individual([0.3,0.2],None,0,0,0,0)
        front2 = [front2_1, front2_2, front2_3, front2_4, front2_5]

        f = [front1, front2]
        
        f2 =[front2]

        #Two test problems with different scaling, manually doing densities
        noScale = nsga2.Problem(2,0,2,[0.,0.],[1.,1.])
        bigScale = nsga2.Problem(2,0,2,[0.,0.],[1.,1000.])
        optNo = nsga2.Optimizer(noScale)
        optNo.density = 5
        optNo.NormVector = [1.,1.]
        optBig = nsga2.Optimizer(bigScale)
        optBig.density = 2
        optBig.NormVector = [1.,2000.]
        
        
        #Do the first front selection
        result1 = optNo._Optimizer__selKrigingPts(f, True)
        result2 = optBig._Optimizer__selKrigingPts(f2, True)
        
       
        #Go ahead and check 
        self.assertEqual(len(result1), 5, "length first front selection 1")        
        self.assertEqual(result1[0], front1_1, "First front selection 1-1/5")
        self.assertEqual(result1[1], front1_2, "First front selection 1-2/5")
        self.assertEqual(result1[2], front2_1, "First front selection 1-3/5") 
        self.assertEqual(result1[3], front2_2, "First front selection 1-4/5")
        self.assertEqual(result1[4], front2_4, "First front selection 1-5/5")
        
        #Now make one to test the objective function selection methods. 
        front2_1.current=True
        front2_1._Individual__obj = [40,100]
        front2_2.current=True
        front2_2._Individual__obj = [40,80]
        front2_3.current=True
        front2_3._Individual__obj = [40,20]
        front2_4.current=True
        front2_4._Individual__obj = [40,50]
        front2_5.current=True
        front2_5._Individual__obj = [40,60]

        #Make an optimizer just to check this
        optMethod2 = nsga2.Optimizer(noScale, krig_select=1,
                 krig_obj = 2)  
        optMethod2.density = 5
        optMethod2.NormVector = [1.,1.]
        result2 = optMethod2._Optimizer__selKrigingPts(f, True)
        self.assertEqual(len(result1), 5, 
                         "krig select = 1 length first front selection 1")        
        self.assertEqual(result2[0], front1_1, 
                         "krig select = 1  First front selection 1-1/5")
        self.assertEqual(result2[1], front1_2, 
                         "krig select = 1  First front selection 1-2/5")
        self.assertEqual(result2[2], front2_3, 
                         "krig select = 1  First front selection 1-3/5") 
        self.assertEqual(result2[3], front2_5, 
                         "krig select = 1  First front selection 1-4/5")
        self.assertEqual(result2[4], front2_1, 
                         "krig select = 1  First front selection 1-5/5")        
        return
        
    def test_sbx(self):
        '''
        Test the sbx method using the test population

        This test uses the first two individuals in the test population ([5,7] and [7,4]) and uses the simulated binary crossover
        method in the NSGA-II optimizer to birth two children.  An sbx exponent and a random seed of 2 are both used.  The test
        ensures that the sbx function creates the proper children based on the information provided.  The calculations supporting the
        assertions in this test can be seen in 'SBX_handCalcs.py'.

        Parameters
        ----------
        No parameters

        Returns
        -------
        No return value
        '''
        
        ### Pick two data points from first test population and identify second two values as genes
        parent_1 = self.function_table[0][1:3]
        parent_2 = self.function_table[1][1:3]
        
        ### Initiate optimizer and problem instance, then run crossover method
        Prob = nsga2.Problem(2,0,2,[-10,-10],[10,10])
        Opt = nsga2.Optimizer(Prob,pc=1.0,sbx_ex=2)
        Opt.rnd = Random(2)
        children = Opt.sbx(parent_1,parent_2)

        ### Assert that crossover method births proper children
        self.assertAlmostEqual(children[0][0],3.87589450248,7, \
                          "returned child 1,1")
        self.assertAlmostEqual(children[0][1],6.3305314095,7, \
                          "returned child 1,2")
        self.assertAlmostEqual(children[1][0],8.12410549752,7, \
                          "returned child 2,1")
        self.assertAlmostEqual(children[1][1],4.6694685905,7, \
                          "retruend child 2,2")

    def test_mutation(self):
        '''
        Test the mutation method using the test population
        
        This test also uses the first two individuals in the first test population ([5,7] and [7,4]) and uses the polynomial mutation 
        operator as programmed in the mutation method within the NSGA-II optimizer.  A mutation exponent and random seed of 2 are used.
        The test ensures that the mutation function properly mutates the individuals put into it.  The calculations supporting the assertions
        in this test can be seen in 'MUT_handCalcs.py'.
        
        Parameters
        ----------
        No parameters
        
        Returns
        -------
        No return value
        '''
        
        ### Indentify two individuals two mutate
        individual_1 = self.function_table[0][1:3]
        individual_2 = self.function_table[1][1:3]
        mut = [[individual_1,individual_2]]

        ### Initiate optimizer and problem isntance, then run mutaton method
        Prob = nsga2.Problem(2,0,2,[0,0],[10,10])
        Opt = nsga2.Optimizer(Prob,p_mut=1.00,mut_ex=2)
        Opt.rnd = Random(2)
        mutation = Opt.mutation(mut)

        ### Assert that genes are correctly mutated
        self.assertEqual(mutation[0][0][0],10, "returned first mutated gene")
        self.assertAlmostEqual(mutation[0][0][1],2.536876063310209,7, \
                               "returned second mutated gene")
        self.assertAlmostEqual(mutation[0][1][0],8.9172137042557083,7, \
                               "returned third mutated gene")
        self.assertAlmostEqual(mutation[0][1][1],2.5098981126016193,7, \
                               "returned fourth mutated gene")
    
    def test_crowding_distance(self):
        '''
        Test the crowding distance method using the second test population
        
        This test uses the second set of test population data defined in the set
        up method to test the crowding distance tournament
        method in the NSGA-II algorithm.  The data is taken from " Multi-
        Objective Optimization using Evolutionary Algorithms" by K. Deb.
        The process has been somewhat althered from the example in the book and 
        the hand calcs script (crowding_handCalcs.py) can be consulted to
        find the data used for the assertions.
        
        Parameters
        ----------
        No Parameters
        
        Returns
        -------
        No return value
        '''
        
        ### Initiate optimizer instance, sort by fronts, and calculate crowding distances for
        ### second test population
        Opt = nsga2.Optimizer(0)
        fronts = Opt.fast_nondominated_sort(self.testPop2)
        Opt.crowding_distance(self.testPop2,fronts)
        
        logger = logging.getLogger('msdl')
        logname = 'NSGA_II_RegTest.log'
        if (os.path.exists(logname)):
                    os.remove(logname)
        hdlr = logging.FileHandler(logname)
        formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
        hdlr.setFormatter(formatter)
        logger.addHandler(hdlr) 
        logger.setLevel(logging.DEBUG)
        logger.debug("test message")        
        
        ### Assert that each individual is assigned correct crowding distance
        self.assertAlmostEqual(self.testPop2[0].crowding_distance,
                               3.72199730094/2.0,
                               7,"returned crowding distance for element '1'" +
                               "expected 3.7219973004 got " + 
                               str(self.testPop2[0].crowding_distance))
        self.assertAlmostEqual(self.testPop2[1].crowding_distance,
                               1.85724959692/2.0,
                               7,"returned crowding distance for element '2'" +
                               "expected 1.85724959692 got " + 
                               str(self.testPop2[1].crowding_distance))
        self.assertTrue(math.isinf(self.testPop2[2].crowding_distance), 
                                "returned crowding distance for element '3'")
        self.assertTrue(math.isinf(self.testPop2[3].crowding_distance), \
                                "returned crowding distance for element '4'")      
        self.assertTrue(math.isinf(self.testPop2[4].crowding_distance), \
                                "returned crowding distance for element '5'")
        self.assertTrue(math.isinf(self.testPop2[5].crowding_distance), \
                                "returned crowding distance for element '6'")        
        self.assertTrue(math.isinf(self.testPop2[6].crowding_distance), \
                                "returned crowding distance for element 'a'")        
        self.assertTrue(math.isinf(self.testPop2[7].crowding_distance), \
                                "returned crowding distance for element 'b'")    
        self.assertAlmostEqual(self.testPop2[8].crowding_distance,
                               2.54386723819/2.0,
                               7, "returned crowding distance for element 'c'")         
        self.assertAlmostEqual(self.testPop2[9].crowding_distance,
                               0.9504048583/2.0,
                               7, "returned crowding distance for element 'd'")            
        self.assertEqual(self.testPop2[10].crowding_distance,4.0/2.0,  
                                "returned crowding distance for element 'e'")     
        self.assertTrue(math.isinf(self.testPop2[11].crowding_distance), 
                                "returned crowding distance for element 'f'")   
            
    
    def test_Individual_comparison(self):
        '''
        Test the individual comparison method using the first set of test population data
        
        This test uses the first set of test data defined in the set up method to test the individual comparison function in the 
        NSGA-II algorithm.  The test data is shuffled twice and the adjacent pairs are compared on each shuffle to test the method.  
        
        Parameters
        ----------
        No Parameters
        
        Returns
        -------
        No return value
        '''
        
        ### Initiate optimizer instance, calculate fronts, and crowding distance
        rnd = Random(2)
        Opt = nsga2.Optimizer(0)
        fronts = Opt.fast_nondominated_sort(self.testPop)
        Opt.crowding_distance(self.testPop,fronts)
        
        ### Copy population to be shuffled and compared
        Pop2  = copy.copy(self.testPop)
        winners = []
        
        ### shuffle the population twice, each time comparing different individuals and
        ### assigning victor to 'winners' list
        for i in range(0,2):
            j = 0
            rnd.shuffle(Pop2)
            for k in range(0,4):
                x = Pop2[j]
                y = Pop2[j+1]
                z = Opt.Individual_comparison(x,y)
                winners.append(z)
                j += 2
        
        ### Assert that the correct individuals were victorious in each comparison
        self.assertEqual(winners[0],self.testPop[5], \
                            "6 victorious over 2")
        self.assertEqual(winners[1],self.testPop[2], \
                            "3 victorious over 4")
        self.assertEqual(winners[2],self.testPop[4], \
                            "5 victorious over 7")
        self.assertEqual(winners[3],self.testPop[7], \
                            "8 victorious over 1")
        self.assertEqual(winners[4],self.testPop[2], \
                            "3 victorious over 9")
        self.assertEqual(winners[5],self.testPop[5], \
                            "6 victorious over 1")
        self.assertEqual(winners[6],self.testPop[1], \
                            "2 victorious over 4")
        self.assertEqual(winners[7],self.testPop[4], \
                            "5 victorious over 8") 
    
    
    def test_run_and_tournament(self):
        '''
        A simple test case problem to test the two_pass_tournament_selection and run methods
        
        This method performs a unit test meant to test both the run and two_pass_tournament methods of the NSGA-II
        algorithm.  It assumes that the other methods tested above are working properly.  This test uses the objective 
        functions found in the 'Min-Ex' example problem in "Multi-Objective Optimization using Evolutionary Algorithms" 
        by K. Deb.  The method ensures that the proper individuals, after cross over (and assuming no mutation in the population)
        get carried over into the second generation of the optimizer.  It also makes sure that the individuals in the second genration
        have the proper objective function values.  The hand calcs supporting this can be seen in the attached .PDF file entitled
        'run_and_tournament_handCalcs'.
        
        Parameters
        ----------
        No parameters
        
        Returns
        -------
        No return value
        '''
        
        ### Set up a problem class to be ran through optimizer with 'Min-Ex' objective functions
        class testProblem(nsga2.Problem):
            def __init__(self):
                nsga2.Problem.__init__(self,2,0,2,[0.1,0], [1,5])
            
            def Eval(self, Individual, metamodel=None):
                f = []
                f.append(Individual.chromosome[0])
                f.append((1+Individual.chromosome[1])/Individual.chromosome[0])
                
                return (f,None)
        
        ### Initiate problem and optimizer, then run the optimizer using run method
        Problem = testProblem()
        
        Opt = nsga2.Optimizer(Problem, pc=1.0,p_mut=0.0)
        Opt.run('unittest','Unit Test Problem',2,2,4)
        
        ### Read database produced by the run method and extract data for assertions
        newFront = []
        
        conn = sqlite3.connect('unittest')
        cur = conn.cursor()
        
        cur.execute('select IndID from generation where GenNum=1')
        for i in cur:
            newFront.append(i[0])
        
        f1_new = []
        f2_new = []
        
        cur.execute('select * from objfun where Generation=1')
        for i in cur:
            if i[0] >= 5:
                if i[1] == 1:
                    f1_new.append(i[3])
                elif i[1] == 2:
                    f2_new.append(i[3])
        
        ### Assert that proper individuals got moved through to second generation
        self.assertTrue(1 in newFront, "Individual 1 in second generation")
        self.assertFalse(2 in newFront, "Individual 2 not in second generation")
        self.assertFalse(3 in newFront, "Individual 3 not in second generation")
        self.assertFalse(4 in newFront, "Individual 4 not in second generation")
        self.assertTrue(5 in newFront, "Individual 5 in second generation")
        self.assertTrue(6 in newFront, "Individual 6 in second generation")
        self.assertTrue(7 in newFront, "Individual 7 in second generation")
        self.assertFalse(8 in newFront, "Individual 8 not in second generation")
                
        ### Assert that the new individuals created within two pass tournament are correct
        self.assertAlmostEqual(f1_new[0], 0.2795697782777084, 7, \
                            "returned individual 5; objective function 1")
        self.assertAlmostEqual(f1_new[1], 0.10805038670546248, 7, \
                            "returned individual 6, objective function 1")
        self.assertAlmostEqual(f1_new[2], 0.30528623524565457, 7, \
                            "returned individual 7, objective function 1")
        self.assertAlmostEqual(f1_new[3], 0.78338668909989306, 7, \
                            "returned individual 8, objective function 1")
        self.assertAlmostEqual(f2_new[0], 8.2552363669662423, 7, \
                            "returned individual 5; objective function 2")
        self.assertAlmostEqual(f2_new[1], 11.78244759167945, 7, \
                            "returned individual 6; objective function 2")
        self.assertAlmostEqual(f2_new[2], 4.2322449818998278, 7, \
                            "returned individual 7; objective function 2")                    
        self.assertAlmostEqual(f2_new[3], 7.0775474346107181, 7, \
                            "returned individual 8; objective function 2")

    
    def getVFOSamples(self):
        '''
        Simple utility function to get some individuals for testing the 
        Kriging model functions
        
        Parameters
        ----------
        No parameters
        
        Returns
        -------
        List of individual for testing Kriging VFO model and problem class 
        instance 
        '''
        I1_chromo = [1,1]
        I2_chromo = [1,2]
        I3_chromo = [2,3]
        class simpleKriging(nsga2.Problem):
            def __init__(self):
                nsga2.Problem.__init__(self,1,0,2,[0.0,0.0], [4.0,4.0])
            
            def Eval(self, Individual, metamodel=None):
                return (Individual.chromosome[0]**2.0 + 
                             Individual.chromosome[1]**2.0 ,None)
            
            def lowfidelity(self, individual_instance):
                temp  =  2.*individual_instance.chromosome[0] + \
                          1.*individual_instance.chromosome[1]
                return temp
            
            def highfidelity(self, individual_instance):
                return individual_instance.chromosome[0]**2.0 +  \
                             individual_instance.chromosome[1]**2.0
        
        testProb = simpleKriging()
        Ind1 = nsga2.Individual(I1_chromo, testProb, 1, 0, 0, 0)
        Ind2 = nsga2.Individual(I2_chromo, testProb, 2, 0, 0, 0)
        Ind3 = nsga2.Individual(I3_chromo, testProb, 3, 0, 0, 0)
        
        return([Ind1, Ind2, Ind3], testProb)
    
    
    def testgetVFOSamples(self):
        '''
        Tests the process of getting VFO samples
        
        Parameters
        ----------
        No parameters
        
        Returns
        -------
        No return value
        '''   
        temp  = self.getVFOSamples()
        indivs = temp[0]
        testProb = temp[1]
        #Data is
        #x1   x2    low   high  dif
        #1    1     3     2      2/3 
        #1    2     4     5      5/4 
        #2    3     7     13    13/7.
        Prob = nsga2.Optimizer(testProb)
        (xvals, yvals, samples) =  Prob._Optimizer__getVFOSamples(indivs)
        self.assertEqual(xvals[0][1], 1, "x value check 1")
        self.assertEqual(xvals[1][0], 1, "x value check 2")
        self.assertEqual(xvals[1][1], 2, "x value check 3")        
        self.assertEqual(xvals[2][1], 3, "x value check 4")  
        self.assertAlmostEqual(yvals[0], 2./3., 7, "yval check 1")
        self.assertAlmostEqual(yvals[1], 5./4., 7, "yval check 2")        
        self.assertAlmostEqual(yvals[2], 13./7., 7, "yval check 3") 

    
    
    
    def testKrigingSelect(self):
        '''
        Test the selection process for Krining model updates
        '''
       
        #Get some parameters we will need for the testing
        temp  = self.getVFOSamples()
        indivs = temp[0]
        testProb = temp[1]                
        opt = nsga2.Optimizer(testProb)
        
        #Make up some individuals
        i1 = nsga2.Individual([1.,1.], testProb, 1, 0, 0, 0)
        i2 = nsga2.Individual([2.,3.], testProb, 2, 0, 0, 0)
        i3 = nsga2.Individual([1.5,2.0], testProb, 3, 0, 0, 0)
        i4 = nsga2.Individual([100.,50.4], testProb, 4, 0, 0 ,0)
        
        #Make into a false list of fronts
        f1 = [i1, i2]
        f2 = [i3, i4]
        f = [f1, f2]        
        
        #Make sample metamodel 
        #Pts are clustered low - should cause i3 to be prefered over i4
        #Class assumes a DB has been set up, so make a simple one
        dbname = 'unittest_testformModel.db'        
        if (os.path.exists(dbname)):
            os.remove(dbname)
        opt.conn = sqlite3.connect(dbname)
        opt.cur = opt.conn.cursor()
        opt._Optimizer__createGATables('test case') 
        opt.rnd=Random(1001)         
        
        (xvals, yvals, samples) =  opt._Optimizer__getVFOSamples(indivs)  
        opt._Optimizer__formModel(xvals, yvals, indivs)
        #Compute the errors:
       # print('estimated COVs')
       # for indiv in  [i1,i2,i3,i4]:
        #    xpts = np.matrix(indiv.chromosome)
         #   (yp,mse)=testProb.metamodel.predictor(xpts)
        #    print((mse[0,0]**0.5)/yp[0,0])
        
        #Try some selects
        opt.density = 3
        opt.TOL = 0.1
        opt.NormVector = [1.,1.]

        first_list = opt._Optimizer__selKrigingPts(f, True)
        self.assertEqual(len(first_list), 3, "Returns 3 points on first call")
        self.assertEqual(i1, first_list[0], "First call, item 1")
        self.assertEqual(i2, first_list[1], "First call, item 2")
        self.assertEqual(i3, first_list[2], "First call, item 3")
        
        opt._krigingpts=[]
        second_list = opt._Optimizer__selKrigingPts(f, False)
        self.assertEqual(len(second_list), 3, "Returns 3 points on 2nd call")
        self.assertEqual(i1, second_list[0], "Second call, item 1")
        self.assertEqual(i2, second_list[1], "Second call, item 2")        
        self.assertEqual(i4, second_list[2], "Second call, item 3")

        #Now make one pt fail distance by including it in list
        opt._krigingpts=[i1]
        third_list = opt._Optimizer__selKrigingPts(f, True)
        self.assertEqual(len(third_list), 3, "Returns 3 points on 3rd call")
        self.assertEqual(i2, third_list[0], "third call, item 1")
        self.assertEqual(i3, third_list[1], "third call, item 2")
        self.assertEqual(i4, third_list[2], "third call, item 3")  

        #now try to over-fill list
        opt._krigingpts=[]
        opt.density = 5
        fourth_list = opt._Optimizer__selKrigingPts(f, False)                      
        self.assertEqual(len(fourth_list), 4, "Returns 4 points on 4th call")
        self.assertEqual(i1, fourth_list[0], "fourth call, item 1")
        self.assertEqual(i2, fourth_list[1], "fourth call, item 2")
        self.assertEqual(i4, fourth_list[2], "fourth call, item 3")
        self.assertEqual(i3, fourth_list[3], "fourth call item 4")        
    
    def testokToAdd(self):
        '''
        Test the distance measuring function to see if points should be added
        to the Kriging model
        '''

        temp  = self.getVFOSamples()
        indivs = temp[0]
        testProb = temp[1]                
        opt = nsga2.Optimizer(testProb)
        
        #Make up some individuals
        i1 = nsga2.Individual([1.,1.,1.], testProb, 1, 0, 0, 0)
        i2 = nsga2.Individual([2.,3.,4.], testProb, 1, 0, 0, 0)
        opt._krigingpts = [i1, i2]
        
        #Made up normalization vector and TOL
        normVect = [2.,1.,4.]
        TOL = 0.05
                
        #Test adding an identical point - should fail
        self.assertFalse(opt._Optimizer__okToAdd(i1, TOL, normVect,
                                                opt._krigingpts),
                         'Add identical point 1/2 returns false')
        self.assertFalse(opt._Optimizer__okToAdd(i2, TOL, normVect,
                                                 opt._krigingpts),
                         'Add identical point 2/2 returns false')
        
        #Addition points
        #Just fails on first pt
        i3 = nsga2.Individual([1.05,1.04,1.06], testProb, 1, 0, 0, 0)
        #Just passes on second pt
        i4 = nsga2.Individual([1.97,2.95,4.05], testProb, 1, 0, 0, 0)
    
        self.assertFalse(opt._Optimizer__okToAdd(i3, TOL, normVect,
                                                 opt._krigingpts),
                         "Point that just fails")
        self.assertTrue(opt._Optimizer__okToAdd(i4, TOL, normVect,
                                                opt._krigingpts),
                         "Point that just passes")
    
    def testformModel(self):
        '''
        Tests the process of forming and storing a model
        
        Parameters
        ----------
        No parameters
        
        Returns
        -------
        No return value
        '''  
        temp  = self.getVFOSamples()
        indivs = temp[0]
        testProb = temp[1]
        Prob = nsga2.Optimizer(testProb)
    
        #Class assumes a DB has been set up, so make a simple one
        dbname = 'unittest_testformModel.db'        
        if (os.path.exists(dbname)):
            os.remove(dbname)
        Prob.conn = sqlite3.connect(dbname)
        Prob.cur = Prob.conn.cursor()
        Prob._Optimizer__createGATables('test case')  

        (xvals, yvals, pts) =  Prob._Optimizer__getVFOSamples(indivs)
        Prob._Optimizer__formModel(xvals, yvals, indivs)

        #temporary testing code to get the theta values out
#        sita0 = np.ones((2))*10.0
#        
#        #Setup and solve the model
#        testKrig=kriging.Kriging(xvals,np.transpose(np.matrix(yvals)), sita0)
#        theta0=testKrig.solve()                                 
#        print('-----')
#        print(theta0)
#       [ 12.58792896   1.91314248]

        #Test that the values are coming back ok from the database
        Prob.cur.execute('select * from metamodels')
        numrow = 0
        for row in Prob.cur:
            numrow += 1
            output = row
        self.assertEqual(numrow, 1, "length of returned cursor")
        self.assertEqual(output[0], 1, "first model id")
        self.assertEqual(output[1], 0, "generation")
        
        #mmData(modelID, IndID, Value)
        Prob.cur.execute('select * from mmData')
        #Correct answers
        ids = [1,2,3]
        vals = [2./3., 5./4., 13./7.]
        kount = 0
        for row in Prob.cur:
            self.assertEqual(row[0], 1, "model id in data")
            self.assertEqual(row[1], ids[kount], "ind id in data")
            self.assertAlmostEqual(row[2], vals[kount], 7, "values in data")
            kount += 1
        
        Prob.cur.execute('select * from mmParameters')
        #modelID, ParamID, Value
        vals = [ 12.58792896 ,   1.91314248]
        kount = 0
        for row in Prob.cur:
            self.assertEqual(row[0], 1, "model id")
            self.assertEqual(row[1], ids[kount], "sita/theta")
            self.assertAlmostEqual(row[2], vals[kount], 7, 
                "sita/theta value was " +str(row[2]) + " should be " +
                str(vals[kount]))
            kount += 1
        Prob.cur.close()
            
    def testwriteChromsomestoDB(self):
        '''
        Tests the process of writing a DB model
        
        Parameters
        ----------
        No parameters
        
        Returns
        -------
        No return value
        '''           
        temp  = self.getVFOSamples()
        indivs = temp[0]
        testProb = temp[1]
        Prob = nsga2.Optimizer(testProb)
    
        #Class assumes a DB has been set up, so make a simple one
        dbname = 'unittest_testformModel.db'        
        if (os.path.exists(dbname)):
            os.remove(dbname)
        Prob.conn = sqlite3.connect(dbname)
        Prob.cur = Prob.conn.cursor()
        Prob._Optimizer__createGATables('test case')  
        Prob._Optimizer__writeChromsomestoDB(indivs)
        
        #IndID,born,died
        Prob.cur.execute('select * from individuals')
        ids = [1,2,3]
        born = [0,0,0]
        died =[None, None, None]
        for row, i, b, d in zip(Prob.cur, ids, born, died):  
            self.assertEqual(row[0], i, "indiv id")
            self.assertEqual(row[1], b, "born gen")
            self.assertEqual(row[2], d, "died gen")
            
        Prob.cur.execute('select * from chromosomes')
        #IndID,ChromID,Value
        ids = [1,1,2,2,3,3]
        chromid = [1,2,1,2,1,2]
        val = [1,1,1,2,2,3]
        for row, i,c, v in zip(Prob.cur, ids, chromid, val):  
            self.assertEqual(row[0], i, "indiv id")
            self.assertEqual(row[1], c, "chrom id")
            self.assertEqual(row[2], v, "value")        


    def oldtestVFORun(self):
        '''
        This no longer works SVN version > 180 as NSGA-II changed to allow
        minimum distance filtering in setting up Kriging Models
        
        Tests a simple VFO run, making a database and comparing the values
        
        Parameters
        ----------
        No parameters
        
        Returns
        -------
        No return value
        '''                                 
        class SCH_vfo(nsga2.Problem):
            def __init__(self):
                nsga2.Problem.__init__(self,2,0,1,[-1], [6])
                self.low_fi_cov = 0.1
                self.high_fi_cov = 0.05
            
            def Eval(self, indiv, metamodel):
                obj = []
                corr = []
                if (metamodel == None):
                    obj.append((indiv.chromosome[0]-0.5)**2.0)
                    obj.append((indiv.chromosome[0]-2)**2.0)
                    #Tuple, none in place of VFO list and 
                    return (obj, None)
                else:
                    (yp,mse)=metamodel.predictor(np.matrix(indiv.chromosome))
                    cov_k=(mse[0,0]**0.5)/yp[0,0]
                    if (cov_k**2.0+self.high_fi_cov**2.0)**0.5<=self.low_fi_cov:
                        obj.append(yp[0,0]*(indiv.chromosome[0]-0.5)**2.0)
                        corr.append(yp[0,0])
                    else:
                        obj.append((indiv.chromosome[0]-0.5)**2.0)
                        corr.append(1.)
                    obj.append((indiv.chromosome[0]-2)**2.0)
                    corr.append(1.)
                    return (obj, corr)    
    
            def lowfidelity(self, indiv):
                return (indiv.chromosome[0]-0.5)**2.0
            
            def highfidelity(self, indiv):
                return (indiv.chromosome[0])**2.0
 
        #Make the run
        myProb = SCH_vfo()
        myOpt = nsga2.Optimizer(myProb, use_metamodel=True, offset=1, density=3, 
                              spacing=2)
        myOpt.run('NSGA_regtest_sch_vfo.db', 
                         'SCH Test Problem VFO', 1024, 3, 6)
        
        conn = sqlite3.connect('NSGA_regtest_sch_vfo.db')
        cur = conn.cursor()
        #Correct data from carefully-followed and hand-calced run
        #see file nsgaII_vfo_regression_test.ods/pdf in repo
        cd= [ [1,1,0,0.00629816834717053,1], 
            [1,2,0,2.49438116884182,1], 
            [2,1,0,4.47537371111352,1], 
            [2,2,0,0.378849985480589,1], 
            [3,1,0,7.69510544210702,1], 
            [3,2,0,1.62308952173273,1], 
            [4,1,0,0.0935559787172721,1], 
            [4,2,0,3.26116363368504,1], 
            [5,1,0,0.351689223812039,1], 
            [5,2,0,4.38079096747606,1], 
            [6,1,0,23.010770800631,1], 
            [6,2,0,10.8699078242255,1], 
            [1,1,1,0.176947536470543,28.0950788732151], 
            [1,2,1,2.49438116884182,1], 
            [2,1,1,6.8334970266069,1.52691092805], 
            [2,2,1,0.378849985480589,1], 
            [3,1,1,10.7168801624676,1.39268789012632], 
            [3,2,1,1.62308952173273,1], 
            [4,1,1,2.69024454520605,28.7554529607992], 
            [4,2,1,3.26116363368504,1], 
            [5,1,1,0.351689223812039,1], 
            [5,2,1,4.38079096747606,1], 
            [6,1,1,23.010770800631,1], 
            [6,2,1,10.8699078242255,1], 
            [7,1,1,5.12516093185667,1], 
            [7,2,1,0.583515447387534,1], 
            [8,1,1,1.48423937611872,28.6183566926124], 
            [8,2,1,2.98506795301668,1], 
            [9,1,1,0.351689223812039,1], 
            [9,2,1,4.38079096747606,1], 
            [10,1,1,0.176947536470543,28.0950788732151], 
            [10,2,1,2.49438116884182,1], 
            [11,1,1,10.5048297082298,1], 
            [11,2,1,3.03148319717361,1], 
            [12,1,1,2.04912586748409,1], 
            [12,2,1,8.59355630787509,1], 
            [1,1,2,0.176947536470543,28.0950788732151], 
            [1,2,2,2.49438116884182,1], 
            [2,1,2,6.8334970266069,1.52691092805], 
            [2,2,2,0.378849985480589,1], 
            [10,1,2,0.176947536470543,28.0950788732151], 
            [10,2,2,2.49438116884182,1], 
            [7,1,2,5.12516093185667,1], 
            [7,2,2,0.583515447387534,1], 
            [3,1,2,10.7168801624676,1.39268789012632], 
            [3,2,2,1.62308952173273,1], 
            [5,1,2,0.351689223812039,1], 
            [5,2,2,4.38079096747606,1], 
            [13,1,2,0.243381105233774,1], 
            [13,2,2,3.97339122246405,1], 
            [14,1,2,6.39828735411599,1], 
            [14,2,2,1.05983651174744,1], 
            [15,1,2,5.12516093185667,1], 
            [15,2,2,0.583515447387534,1], 
            [16,1,2,5.12516093185667,1], 
            [16,2,2,0.583515447387534,1], 
            [17,1,2,0.520336511770751,28.3342426830189], 
            [17,2,2,2.67490827057957,1], 
            [18,1,2,4.71611405720991,1], 
            [18,2,2,0.451129290160243,1], 
            [1,1,3,0.176939762407044,28.0938445360125], 
            [1,2,3,2.49438116884182,1], 
            [2,1,3,6.84215490824614,1.52884548864719], 
            [2,2,3,0.378849985480589,1], 
            [10,1,3,0.176939762407044,28.0938445360125], 
            [10,2,3,2.49438116884182,1], 
            [18,1,3,7.13922464980715,1.51379389115767], 
            [18,2,3,0.451129290160243,1], 
            [7,1,3,7.6407314847675,1.4908276220703], 
            [7,2,3,0.583515447387534,1], 
            [13,1,3,0.243381105233774,1], 
            [13,2,3,3.97339122246405,1], 
            [19,1,3,0.176939762407044,28.0938445360125], 
            [19,2,3,2.49438116884182,1], 
            [20,1,3,0.176939762407044,28.0938445360125], 
            [20,2,3,2.49438116884182,1], 
            [21,1,3,6.84215490824614,1.52884548864719], 
            [21,2,3,0.378849985480589,1], 
            [22,1,3,0.176939762407044,28.0938445360125], 
            [22,2,3,2.49438116884182,1], 
            [23,1,3,6.7348830669166,1.53518161731945], 
            [23,2,3,0.353457686222586,1], 
            [24,1,3,0.223116290784778,1], 
            [24,2,3,3.89017190615878,1]]
    
        cur.execute('select * from objfun')
        numrow = 0
        for row in cur:
            self.assertEqual(row[0],cd[numrow][0], "ind id row " + 
                            str(numrow))
            self.assertEqual(row[1],cd[numrow][1], "obj id row " + 
                            str(numrow))
            self.assertEqual(row[2],cd[numrow][2], "generation row " + 
                            str(numrow))
            self.assertAlmostEqual(row[3],cd[numrow][3], 7,
                                   "value row " + str(numrow))
            self.assertAlmostEqual(row[4],cd[numrow][4], 7,
                          "correction row " + str(numrow))                                
            numrow += 1

def testKriging():
    '''Simple test function for the vfo test case above that
    build the Kriging model and shows the estimated and true errors
    '''
    
    #First Kriging model data
    xvals = np.matrix([0.420638999835117, 2.61550790854431, 3.27400530679143])
    yvals = np.matrix([28.0934326345495, 1.52856097864412, 1.39297776093415])
    sita = np.array([1.])
    km = kriging.Kriging(np.transpose(xvals),np.transpose(yvals), sita)
    actual_sita = km.solve()
    
    print(("Site for KM 1 " + str(actual_sita)))
    
    
    #Generation 0 gets to use the KM before the selection pool, so check
    #them
    pts_gen0 = [0.420638999835117, 2.61550790854431, 3.27400530679143,
                0.19413078167741, -0.093033914554673, 5.2969543254685]
    print("Generation 0")
    for pt in pts_gen0:
        (yp, mse) = km.predictor(np.matrix([pt]))
        print(('Point ', pt, ' y pred ', yp[0,0] , ' mse ', mse[0,0], ','))
    
    
    #Generation 1 
    
    pts_gen1 = [2.76388182815638, 0.272265080223046, -0.093033914554673
                , 0.420638999835117, 3.74111550368538 ,-0.931476813463666]
    print("Generation 1 Children")                
    for pt in pts_gen1:
        (yp, mse) = km.predictor(np.matrix([pt]))
        print(('Point ', pt, ' y pred ', yp[0,0] , ' mse ', mse[0,0], ','))
    
    #Generation 2
    pts_gen2 = [0.00666329425657597, 3.02948361412285, 2.76388182815638,
                2.76388182815638, 0.364485319362871, 2.67166158901655]
                
    print("Generation 2 Children")                  
    for pt in pts_gen2:
        (yp, mse) = km.predictor(np.matrix([pt]))
        print(('Point ', pt, ' y pred ', yp[0,0] , ' mse ', mse[0,0], ','))

    #Generation 3 before KM update
    pts_gen3_parents=[0.420638999835117, 2.61550790854431, 0.420638999835117,
              2.67166158901655,2.76388182815638,0.00666329425657597]
    print("Generation 3 Parents for Kriging model update ")                  
    
    for pt in pts_gen3_parents:
        (yp, mse) = km.predictor(np.matrix([pt]))
        print(('Point ', pt, ' y pred ', yp[0,0] , ' mse ', mse[0,0], ',')) 

    #Generation three kringing model
    xvals = np.matrix([0.420638999835117, 2.61550790854431, 3.27400530679143,
                       0.420638999835117 ,2.67166158901655, 2.76388182815638])
    yvals = np.matrix([28.0934326345496,1.52856097864412,1.39297776093415,
                      28.0934326345496,1.51348664592078,1.49049812514779])
    sita = np.array([1.])                      
    km2 = kriging.Kriging(np.transpose(xvals),np.transpose(yvals), sita)
    actual_sita = km2.solve()
    
    print(("Site for KM 2 " + str(actual_sita)))
    
    #Gen 3 parents re-sort after KMM
    pts_gen3_parents_post =[0.420638999835117, 2.61550790854431, 
                            0.420638999835117,
              2.67166158901655,2.76388182815638,0.00666329425657597]
    print("Generation 3 Parents after Kriging model update ")   
    
    for pt in pts_gen3_parents_post:
        (yp, mse) = km2.predictor(np.matrix([pt]))
        print(('Point ', pt, ' y pred ', yp[0,0] , ' mse ', mse[0,0], ',')) 
        
    #Gen 3 children
    pts_gen3_children = [0.420638999835117,0.420638999835117, 2.61550790854431,
                        0.420638999835117,2.59452307459222,0.0276481282086652]
                        
    print("Generation 3 children after Kriging model update ")                         
    for pt in pts_gen3_children:
        (yp, mse) = km2.predictor(np.matrix([pt]))
        print(('Point ', pt, ' y pred ', yp[0,0] , ' mse ', mse[0,0], ','))         
        
        
#This is currently used to load and run the tests - eventually
#we can replace this with a better test procedure
suite = unittest.TestLoader().loadTestsFromTestCase(testNSGA2)
unittest.TextTestRunner(verbosity=2).run(suite)


        

        
