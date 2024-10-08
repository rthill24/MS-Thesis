# Test Class for tee-stiffened panel with transverse members (Length & Breadth are Set)
# Unit test code for TPanel_trans
# utilizing unittest framework

# Author: Travis Zahradka
# Date Created: 10-05-10
# Date Completed: _______

import unittest
import TPanel_trans as Panel
import Structures

class test_TPanel_trans(unittest.TestCase):
    """
    Unit test framework for tee-stiffened panel with transverse members and fixed length and breadth
    """
    
    def setUp(self):
        """
        sets up the panels used in each test, called before each test
        """
        self.matl = Structures.EPMatl(262.0, 70000.0, 0.3)
        self.smatl = Structures.EPMatl(200.0, 70000.0,0.3)
        self.Panel = Panel.TPanel_trans(3000.0, 2000.0, 9.0, 3.0, 5.0, 4.0, 80.0, 10.0, 30.0, \
                200.0, 10.0, 20.0, 75.0, pmatl = self.matl, smatl = self.smatl, tmatl = self.matl, eta = 4.5)
        
    def test_transFunctions(self):
        """
        tests all basic get functions in TPanel_trans class
        """
        self.assertEqual(self.Panel.getB(), 3000.0)
        self.assertEqual(self.Panel.getL(), 2000.0)
        self.assertEqual(self.Panel.getnstiff(), 9.0)
        self.assertEqual(self.Panel.getntrans(), 3.0)
        self.assertEqual(self.Panel.getb(), 300.0)
        self.assertAlmostEqual(self.Panel.geta(), 666.66666666666663, 7)
        
        self.assertEqual(self.Panel.getmatlT(), self.matl)
        self.assertEqual(self.Panel.gettwh(), 200.0)
        self.assertEqual(self.Panel.gettwt(), 10.0)
        self.assertEqual(self.Panel.gettft(), 20.0)
        self.assertEqual(self.Panel.gettfb(), 75.0)
        self.assertAlmostEqual(self.Panel.gettpa(), 3333.333333333333, 7)
        self.assertEqual(self.Panel.gettwa(), 2000.0)
        self.assertEqual(self.Panel.gettfa(), 1500.0)
        self.assertAlmostEqual(self.Panel.getta(), 6833.333333333333, 7)
        self.assertEqual(self.Panel.gettsa(), 3500.0)   
        
        self.assertAlmostEqual(self.Panel.gett_fmom(), 540833.3333333333, 7) 
        self.assertAlmostEqual(self.Panel.gett_fmomStiff(), 532500.0, 7)
        self.assertAlmostEqual(self.Panel.get_tNA(), 79.146341463414643, 7)
        self.assertAlmostEqual(self.Panel.get_tNAStiff(), 152.1428571, 7)  
        self.assertAlmostEqual(self.Panel.gety_max(), 145.85365853658536, 7)   
        self.assertAlmostEqual(self.Panel.gett_INA(), 55326964.769647688, 7)
        self.assertAlmostEqual(self.Panel.gett_INAStiff(), 17088095.238095239, 7)
        self.assertAlmostEqual(self.Panel.gett_rad_gyr(), 89.981270221530821, 7)
        self.assertAlmostEqual(self.Panel.gett_asm(), 379332.0327015979, 7)
        self.assertAlmostEqual(self.Panel.getTotalVolume(), 72660000.0, 7)
        
        
#       Test update function of TPanel_trans class
        self.Panel.update(B=2000.0, L=3000.0, nstiff=7.0, ntrans=4.0, tp=4.0, tw=-1, hw=-1, tf=-1, bf=-1, \
                pmatl=-1, smatl=-1, eta=-1, tmatl=-1, twh=250.0, twt=5.0, tft=2.0, tfb=25.0)
        
        self.assertEqual(self.Panel.getB(), 2000.0)
        self.assertEqual(self.Panel.getL(), 3000.0)
        self.assertEqual(self.Panel.getnstiff(), 7.0)
        self.assertEqual(self.Panel.getntrans(), 4.0)
        self.assertEqual(self.Panel.getb(), 250.0)
        self.assertAlmostEqual(self.Panel.geta(), 750.0, 7)
        self.assertEqual(self.Panel.gettp(), 4.0)        
        
        self.assertEqual(self.Panel.getmatlT(), self.matl)
        self.assertEqual(self.Panel.gettwh(), 250.0)
        self.assertEqual(self.Panel.gettwt(), 5.0)
        self.assertEqual(self.Panel.gettft(), 2.0)
        self.assertEqual(self.Panel.gettfb(), 25.0)
        self.assertEqual(self.Panel.gettpa(), 3000.0)
        self.assertEqual(self.Panel.gettwa(), 1250.0)
        self.assertEqual(self.Panel.gettfa(), 50.0)
        self.assertEqual(self.Panel.getta(), 4300.0)
        self.assertEqual(self.Panel.gettsa(), 1300.0)  
        
    def test_constraints(self):
        """
        test the constaint violation method in TPanel_trans class
        """
        self.assertAlmostEqual(self.Panel.constraints()[0], 0.0, 7)
        self.assertAlmostEqual(self.Panel.constraints()[1], 0.0, 7)
        # test with panel parameters that violate constraints
        self.Panel.update(B=2000.0, L=3000.0, nstiff=7.0, ntrans=4.0, tp=4.0, tw=1.2, hw=35.0, tf=-1, bf=-1, \
                pmatl=-1, smatl=-1, eta=-1, tmatl=-1, twh=250.0, twt=5.0, tft=2.0, tfb=25.0)
        self.assertAlmostEqual(self.Panel.constraints()[0], 0.037859529115272239, 7)
        
        self.Panel.update(B=2000.0, L=3000.0, nstiff=7.0, ntrans=4.0, tp=4.0, tw=4.0, hw=80.0, tf=0.5, bf=40.0, \
                pmatl=-1, smatl=-1, eta=-1, tmatl=-1, twh=250.0, twt=5.0, tft=2.0, tfb=25.0)
        self.assertAlmostEqual(self.Panel.constraints()[1], 0.76614641332662869, 7)
        
    
    def test_EPmatl(self):
        """
        test the derived elastic-plastic material class
        """
        self.matl.update()
        self.failUnlessEqual(self.matl.getYld(), 262.0, "Returned yield")
        self.assertEqual(self.matl.getE(), 70000.0, \
                         "Returned elastic modulus")
        self.matl.update(yld=245.0, E=71000.0)
        self.assertEqual(self.matl.getYld(), 245.0, \
                         "Updated yield")
        self.assertEqual(self.matl.getE(), 71000.0, \
                         "Updated elastic modulus")
        
    def test_TPanel(self):
        """
        Test derived TPanel class methods
        """
        self.assertEqual(self.Panel.gettp(), 5.0)
        self.assertEqual(self.Panel.gettw(), 4.0)
        self.assertEqual(self.Panel.gethw(), 80.0)
        self.assertEqual(self.Panel.gettf(), 10.0)
        self.assertEqual(self.Panel.getbf(), 30.0)
        self.assertAlmostEqual(self.Panel.geta(), 666.66666666666663, 7)
        self.assertEqual(self.Panel.getmatlP(), self.matl)
        self.assertEqual(self.Panel.getmatlS(), self.smatl)
        self.assertEqual(self.Panel.getEta(), 4.5)
        self.assertAlmostEqual(self.Panel.getYsavg(), 243.867924528, 7)
        self.assertAlmostEqual(self.Panel.getBeta(),3.67073368,7)
        self.assertAlmostEqual(self.Panel.getLambda(),0.38009580746767235,7)
        self.assertAlmostEqual(self.Panel.getArea(), 2120.0, 7)
        self.assertAlmostEqual(self.Panel.getNA(), 21.2971698, 7)
        self.assertAlmostEqual(self.Panel.getINA(), 2302099.44968553, 7)

suite = unittest.TestLoader().loadTestsFromTestCase(test_TPanel_trans)
unittest.TextTestRunner(verbosity=2).run(suite)