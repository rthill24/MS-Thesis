#General post-processing module for optimizers.  Contains several classes
#designed to post-process different databases formats
# (c) 2011 University of Michigan, Dr. Matthew Collette NA&ME Dept

import unittest
import sqlite3
import post_process
import os
import numpy as np

class testPP(unittest.TestCase):
    """
    Unit test framework for post-processor
    """

    def setUp(self):
        """
        Setups up the database used in the test
        """
        
        dbname = 'test_post_process.sql3'
        if (os.path.exists(dbname)):
            os.remove(dbname)
        conn = sqlite3.connect(dbname)        
        cur = conn.cursor()       
        #Create individual tables        
        cur.execute( '''CREATE TABLE individuals (IndID,born,died)''' )
        #create chromosome tables
        cur.execute( '''CREATE TABLE chromosomes (IndID,ChromID,Value)''' )
        #Create objective function value table
        cur.execute( 
          '''CREATE TABLE objfun (IndID,ObjID,Generation,Value,Correction)''' )
        #Create generation table
        cur.execute( 
           '''CREATE TABLE generation (GenNum, IndID, Front, CrowdDist)''' ) 
        conn.commit()
        
        ##Add some data  - totally fictitous, not real fronts
        cur.execute('insert into individuals values (?,?,?)', (1,1,1))
        cur.execute('insert into individuals values (?,?,?)', (2,1,None))
        cur.execute('insert into individuals values (?,?,?)', (3,2,2))
        cur.execute('insert into individuals values (?,?,?)', (4,2,None))
        conn.commit()
        
        cur.execute('insert into chromosomes values (?,?,?)', (1,1,1))
        cur.execute('insert into chromosomes values (?,?,?)', (2,1,2))
        cur.execute('insert into chromosomes values (?,?,?)', (3,1,3))
        cur.execute('insert into chromosomes values (?,?,?)', (4,1,4))
        cur.execute('insert into chromosomes values (?,?,?)', (1,2,10))
        cur.execute('insert into chromosomes values (?,?,?)', (2,2,20))
        cur.execute('insert into chromosomes values (?,?,?)', (3,2,30))
        cur.execute('insert into chromosomes values (?,?,?)', (4,2,40))        
        
        conn.commit()
        
        cur.execute('insert into objfun values (?,?,?,?,?)', (1,1,1,5,1.2))
        cur.execute('insert into objfun values (?,?,?,?,?)', (1,2,1,4,1.0))

        cur.execute('insert into objfun values (?,?,?,?,?)', (2,1,1,1,1.0))
        cur.execute('insert into objfun values (?,?,?,?,?)', (2,1,2,8,1.0))
        cur.execute('insert into objfun values (?,?,?,?,?)', (2,2,1,3,1.0))
        cur.execute('insert into objfun values (?,?,?,?,?)', (2,2,2,2,1.4))
        
        cur.execute('insert into objfun values (?,?,?,?,?)', (3,1,2,11,1.1))
        cur.execute('insert into objfun values (?,?,?,?,?)', (3,2,2,14,3.2))  
        
        cur.execute('insert into objfun values (?,?,?,?,?)', (4,1,2,21,1.0))
        cur.execute('insert into objfun values (?,?,?,?,?)', (4,2,2,24,1.1)) 
        conn.commit()

        cur.execute('insert into generation values (?,?,?,?)', (1,1,1,1.1))        
        cur.execute('insert into generation values (?,?,?,?)', (1,2,0,1.2)) 
        cur.execute('insert into generation values (?,?,?,?)', (2,2,0,1.4)) 
        cur.execute('insert into generation values (?,?,?,?)', (2,3,1,1.5)) 
        cur.execute('insert into generation values (?,?,?,?)', (2,4,0,1.6)) 
                                                               

        #Some partial data to test infinite crowding distances
        cur.execute('insert into generation values (?,?,?,?)', (3,1,0,1.0))        
        cur.execute('insert into generation values (?,?,?,?)', (3,2,0,1.2)) 
        cur.execute('insert into generation values (?,?,?,?)', (3,2,0,
                    float('Inf'))) 
        cur.execute('insert into generation values (?,?,?,?)', (3,3,0,
                    float('Inf')))          
        
                                                               
        conn.commit()  
        
        conn.close()
        self.dbname = dbname
    
    
    def test_ParetoStats(self):
        testobj = post_process.NSGA_PostProcess(self.dbname) 
        numpts, avg_c, span = testobj.ParetoStats([1,2], [5.5, 7.8])
        self.assertEqual(numpts[0],1, "numpt gen 1")
        self.assertEqual(numpts[1], 2, "numpt gen 2")
        self.assertEqual(avg_c[0],1.2, "avg crowd gen 1")        
        self.assertAlmostEqual(avg_c[1],1.5,7, "avg crowd gen 2")
        self.assertAlmostEqual(span[0],0.,7, "span gen 1")
        self.assertAlmostEqual(span[1],6.66666666666667,7, "span gen 2")  
        numpts, avg_c, span = testobj.ParetoStats([3], None)
        self.assertEqual(numpts[0],4, "numpt gen 3")
        self.assertAlmostEqual(avg_c[0],1.1,7,"avg crowd with infs")        
        
    
    
    def test_KrigStats(self):
        testobj = post_process.NSGA_PostProcess(self.dbname) 
        (mx, mn, avg, std, per) = testobj.KrigingStats([1,2], 1)
        self.assertEqual(mx[0],1.2, "maximum gen 1 obj 1")
        self.assertEqual(mn[0],1.2, "min gen 1 obj 1")
        self.assertEqual(avg[0],1.2, "avg gen 1 obj 1")  
        self.assertEqual(std[0],0., "st dev gen 1 obj 1")   
        self.assertEqual(per[0],0.5, "percent updated gen 1 obj 1")           
        self.assertEqual(mx[1],1.1, "maximum gen 2 obj 1")
        self.assertEqual(mn[1],1.1, "min gen 2 obj 1")
        self.assertEqual(avg[1],1.1, "avg gen 2 obj 1")  
        self.assertEqual(std[1],0., "st dev gen 2 obj 1")   
        self.assertEqual(per[1],1./3., "percent updated gen 2 obj 1")       
        (mx, mn, avg, std, per) = testobj.KrigingStats([1,2], 2)
        self.assertEqual(mx[0],0., "maximum gen 1 obj 2")
        self.assertEqual(mn[0],0., "min gen 1 obj 2")
        self.assertEqual(avg[0],0., "avg gen 1 obj 2")  
        self.assertEqual(std[0],0., "st dev gen 1 obj 2")   
        self.assertEqual(per[0],0., "percent updated gen 1 obj 2")           
        self.assertEqual(mx[1],3.2, "maximum gen 2 obj 2")
        self.assertEqual(mn[1],1.1, "min gen 2 obj 2")
        self.assertAlmostEqual(avg[1],1.90, 7, "avg gen 2 obj 2")  
        self.assertAlmostEqual(std[1],1.13578166916005, 7, "st dev gen 2 obj 2")   
        self.assertEqual(per[1],1., "percent updated gen 2 obj 2")           
    
    
    def test_getIndVariables(self):
        #Test the independent variable retrivial
        testobj = post_process.NSGA_PostProcess(self.dbname)
        retind = testobj.getFront(2,0,[2,])
        retchromosome = testobj.getIndVariables(retind)
        self.assertEqual(len(retchromosome), 2, "Length of chromosome list")
        self.assertEqual(retchromosome[0][0],2,"Returned gene 1-1")
        self.assertEqual(retchromosome[0][1],20,"Returned gene 1-2")
        self.assertEqual(retchromosome[1][0],4,"Returned gene 2-1")
        self.assertEqual(retchromosome[1][1],40,"Returned gene 2-2")        
        
    def testDom(self):
        #Test the domination sorter
        i1 = [1,2,3]
        i2 = [1,1,3]
        i3 = [0,0,6]
        i4 = [0,0,0]
        self.assertEqual(post_process.isDom(i1, i1), 0)
        self.assertEqual(post_process.isDom(i1, i2),1)
        self.assertEqual(post_process.isDom(i2, i1),-1)        
        self.assertEqual(post_process.isDom(i1, i3),0)
        self.assertEqual(post_process.isDom(i1, i4),1)
        return

    def testcombinePareto(self):
        #Test the domination sorter
        i1 = [1,1,2,3]
        i2 = [2,0,4,3]
        i3 = [3,4,0,2]
        f1 = [i1, i2, i3]
        
        i4 = [4,1,2,2]
        i5 = [5,1,4,3]
        i6 = [6,3,-1,2]        
        f2 = [i4, i5, i6]
        
        #final list should be (i4, i2, i6)
        result = post_process.combinePareto([f1, f2])
        self.assertEqual(result[0][0][0], 0)
        self.assertEqual(result[0][0][1], 1)
        self.assertEqual(result[0][0][2], 3)
        self.assertEqual(result[1][0][0], 0)
        self.assertEqual(result[1][0][1], 2) 
        self.assertEqual(result[1][2][0], 1)
        self.assertEqual(result[1][2][1], 6)         
        return
        
        
        
    def test_getFront(self):
        #Test function for get front method
        testobj = post_process.NSGA_PostProcess(self.dbname)
        retlist = testobj.getFront(1,1, [1,2])
        self.assertEqual(len(retlist),1, "Length of returned list")
        self.assertEqual(retlist[0][0], 1, "returned ID")
        self.assertEqual(retlist[0][1], 5, "First obj function")
        self.assertEqual(retlist[0][2], 4, "Second obj function")
        retlist = testobj.getFront(2,0,[2,])
        self.assertEqual(len(retlist),2, "Length of 2nd returned list")
        self.assertEqual(retlist[0][0], 2, "returned ID - 1/2nd try")
        self.assertEqual(retlist[0][1], 2, " obj function 1-2/ 2nd try")
        self.assertEqual(retlist[1][0], 4, "returned ID - 2/2nd try")
        self.assertEqual(retlist[1][1], 24, " obj function 2-2/ 2nd try") 
        #Try with VFO on and off
        retlist = testobj.getFront(2,0,[2,], (2, True))
        self.assertEqual(len(retlist),2, "Length of 2nd returned list VFO")
        self.assertEqual(retlist[0][0], 2, "returned ID - 1/2nd try VFO")
        self.assertEqual(retlist[0][1], 2, " obj function 1-2/ 2nd try VFO")
        self.assertEqual(retlist[1][0], 4, "returned ID - 2/2nd try VFO")
        self.assertEqual(retlist[1][1], 24, " obj function 2-2/ 2nd try VFO")
        retlist = testobj.getFront(2,0,[2,], (1, False))
        self.assertEqual(len(retlist),2, "Length of 2nd returned list VFO")
        self.assertEqual(retlist[0][0], 2, "returned ID - 1/2nd try VFO")
        self.assertEqual(retlist[0][1], 2, " obj function 1-2/ 2nd try VFO")
        self.assertEqual(retlist[1][0], 4, "returned ID - 2/2nd try VFO")
        self.assertEqual(retlist[1][1], 24, " obj function 2-2/ 2nd try VFO")
        retlist = testobj.getFront(2,0,[2,], (2, False))
        self.assertEqual(len(retlist),0, "Length of 2nd returned list  No VFO")
    
    def test_frontPlotFormat(self):
        testobj = post_process.NSGA_PostProcess(self.dbname)
        retlist = testobj.frontPlotFormat(2,0,[2,])        
        self.assertEqual(retlist[0][0], 2, " obj function 1-2")
        self.assertEqual(retlist[1][0], 24, " obj function 2-2")
        return 
    
    def test_reduceFront(self):
        #Test function for reduced front method
        testobj = post_process.NSGA_PostProcess(self.dbname)
        retlist = testobj.getFront(1,1, [1,2])
        shortlist = testobj.reduceFront(retlist)
        self.assertEqual(len(shortlist),1, "Length of returned list")
        self.assertEqual(shortlist[0][0], 5, "First obj function")
        self.assertEqual(shortlist[0][1], 4, "Second obj function")     
        retlist = testobj.getFront(2,0,[2,])
        shortlist = testobj.reduceFront(retlist)
        self.assertEqual(len(shortlist),2, "Length of 2nd returned list")
        self.assertEqual(shortlist[0][0], 2, " obj function 1-2/ 2nd try")
        self.assertEqual(shortlist[1][0], 24, " obj function 2-2/ 2nd try") 

    def test_windowFront(self):
        #Assemble a list of 4 points        
        pt1 = [4,8,23]
        pt2 = [6,8,8]
        pt3 = [-10,7,8]
        pt4 = [7,10,11]
        pts_list = [pt1, pt2, pt3, pt4]
        window_list = [[5,15], [5,10], [7,12]]
        
        #Make up and test the object
        testobj = post_process.NSGA_PostProcess(self.dbname)
        retlist = testobj.windowFront(pts_list, window_list)
        self.assertEqual(len(retlist), 2, "Length of returned list")
        self.assertEqual(retlist[0], pt2, "First point in returned list")
        self.assertEqual(retlist[1], pt4, "Second point in returned list")
            
    def test_UtilAppendFront(self):
        test = np.array( [[3,4,5],[6,7,8]])
        tuple_1 = (50,1,2,3)
        tuple_2 = (50,8,9,7)
        test_list = []
        test_list.append(tuple_1)
        test_list.append(tuple_2)
        testobj = post_process.NSGA_PostProcess(self.dbname)
        ret_test = testobj._UtilAppendFront(None,test_list, [1,1,1])
        self.assertEqual(ret_test[0, 0], 1, "Build with no list 1")
        self.assertEqual(ret_test[1, 2], 7, "Build with no list 2")
        ret_test = testobj._UtilAppendFront(test,test_list, [-0.5, -0.5, 1])
        self.assertEqual(ret_test[0, 0], 3, "Build with  list 3")
        self.assertEqual(ret_test[1, 2], 8, "Build with  list 4")        
        self.assertEqual(ret_test[2, 0], -0.5, "Build with  list 5")
        self.assertEqual(ret_test[3, 2], 7, "Build with  list 6")        
        
      
        
suite = unittest.TestLoader().loadTestsFromTestCase(testPP)
unittest.TextTestRunner(verbosity=2).run(suite)