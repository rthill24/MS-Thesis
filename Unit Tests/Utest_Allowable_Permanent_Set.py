# Test panel for the Allowable Permanent Set Approach from Hughe's
# Using Python unit test

#Import unit test framework and Faulkner code, plus structures
import unittest
import Allowable_Permanent_Set as APS
import Structures as Panel_
import TPanel_trans

class testAllowable_Permanent_Set(unittest.TestCase):
    """
    unit test framework for Allowable Permanent Set
    """
          
    def setUp(self):
        """
        sets up a panel and material for future use
        """
        self.pmat = Panel_.EPMatl(260.0, 70000.0, 0.3)
        self.smat = Panel_.EPMatl(300.0, 70000.0, 0.3)
        self.tP = Panel_.TPanel(300.0, 6.0, 5.0, 80.0, 8.0, 30.0, 300.0, \
                                    self.pmat, self.smat, 3.5)
        self.PanelT = TPanel_trans.TPanel_trans(3000.0, 3000.0, 9.0, 11.0, 6.0, 5.0, 80.0, 8.0, 30.0, \
                200.0, 10.0, 20.0, 75.0, self.pmat, self.smat, self.smat, 4.5)
        
    
    def test_ar(self):
        testobj = APS.Allowable_Permanent_Set()
        self.assertAlmostEqual(testobj._ar(10., 50., \
                                           ), 5., 7)
    def test_wpo(self):
        testobj = APS.Allowable_Permanent_Set()
        self.assertAlmostEqual(testobj._wpo(10., 3., \
                                           ),7 , 7)
    
    def test_Rw(self):
        testobj = APS.Allowable_Permanent_Set()
        self.assertAlmostEqual(testobj._Rw(1.0, 2., 2., \
                                           ), 0.5000000, 7)
        
    def test_Qy(self):
        testobj = APS.Allowable_Permanent_Set()
        self.assertAlmostEqual(testobj._Qy(10., 0.3, 2., \
                                           ), 0.238518635, 7)
                                           
    def test_dQ0(self):
        testobj = APS.Allowable_Permanent_Set()
        self.assertAlmostEqual(testobj._dQ0(10., 0.3, 2., \
                                           ), 0.843815926, 7)
                                         
    def test_dQ1(self):
        testobj = APS.Allowable_Permanent_Set()
        self.assertAlmostEqual(testobj._dQ1(10., 2., \
                                           ), 0.1609514789, 7)
        
    def test_TRw(self):
        testobj = APS.Allowable_Permanent_Set()
        # test case Rw <= 1
        self.assertAlmostEqual(testobj._TRw(0.5, \
                                           ), 0.9564656, 7)
        # test case Rw > 1
        self.assertAlmostEqual(testobj._TRw(1.5, \
                                            ), 1.0000000, 7)
    
    def test_Q(self):
        testobj = APS.Allowable_Permanent_Set(0, 12.)
        self.assertAlmostEqual(testobj._Q(self.tP, \
                                           ), 2.522844178 , 7)
    
    def test_p_aps(self):
        testobj = APS.Allowable_Permanent_Set(0, 12.)
        self.assertAlmostEqual(testobj._p_aps(self.tP, \
                                           ), 2.4363466633384592, 7)
    
    def test_graphs(self):
        """
        this method tests a bunch of beta ratio/aspect ratio pairs and compares them to the 
        graphs (figure 9.15) in Hughe's Ch. 9 - Plate benfding (pg. 352)
        ***Note: aspect ratio is defined as b/a - graphs plot a/b
        """
        
        # test for a beta ratio of ~3.04
        # aspect ratio (ar = 1) : wpt/t = 2
        self.pmat = Panel_.EPMatl(260.0, 70000.0, 0.3)
        self.smat = self.pmat
        self.tP = Panel_.TPanel(300.0, 6.0, 5.0, 80.0, 8.0, 30.0, 300.0, \
                                    self.pmat, self.smat, 3.5)
        testobj = APS.Allowable_Permanent_Set(0, 12.)
        self.assertAlmostEqual(testobj._Q(self.tP, \
                                         ), 2.5, 1)
        
        # aspect ratio (ar = 1) : wpt/t = 3
        self.pmat = Panel_.EPMatl(260.0, 70000.0, 0.3)
        self.smat = self.pmat
        self.tP = Panel_.TPanel(300.0, 6.0, 5.0, 80.0, 8.0, 30.0, 300.0, \
                                    self.pmat, self.smat, 3.5)
        testobj = APS.Allowable_Permanent_Set(0, 18.)
        self.assertAlmostEqual(testobj._Q(self.tP, \
                                         ), 3.2, 1)
        
        # aspect ratio (ar = 1) : wpt/t = 3
        self.pmat = Panel_.EPMatl(260.0, 70000.0, 0.3)
        self.smat = self.pmat
        self.tP = Panel_.TPanel(300.0, 6.0, 5.0, 80.0, 8.0, 30.0, 300.0, \
                                    self.pmat, self.smat, 3.5)
        testobj = APS.Allowable_Permanent_Set(0, 24.)
        self.assertAlmostEqual(testobj._Q(self.tP, \
                                         ), 3.8, 1)
        
        # aspect ratio (ar = 0.3333) : wpt/t = 2
        self.pmat = Panel_.EPMatl(260.0, 70000.0, 0.3)
        self.smat = self.pmat
        self.tP = Panel_.TPanel(300.0, 6.0, 5.0, 80.0, 8.0, 30.0, 900.0, \
                                    self.pmat, self.smat, 3.5)
        testobj = APS.Allowable_Permanent_Set(0, 12.)
        self.assertAlmostEqual(testobj._Q(self.tP, \
                                         ), 0.7, 1)
        
        # aspect ratio (ar = 0.3333) : wpt/t = 5
        self.pmat = Panel_.EPMatl(260.0, 70000.0, 0.3)
        self.smat = self.pmat
        self.tP = Panel_.TPanel(300.0, 6.0, 5.0, 80.0, 8.0, 30.0, 900.0, \
                                    self.pmat, self.smat, 3.5)
        testobj = APS.Allowable_Permanent_Set(0, 30.)
        self.assertAlmostEqual(testobj._Q(self.tP, \
                                         ), 1.1, 1)
        
        # aspect ratio (ar = 0.3333) : wpt/t = 10
        self.pmat = Panel_.EPMatl(260.0, 70000.0, 0.3)
        self.smat = self.pmat
        self.tP = Panel_.TPanel(300.0, 6.0, 5.0, 80.0, 8.0, 30.0, 900.0, \
                                    self.pmat, self.smat, 3.5)
        testobj = APS.Allowable_Permanent_Set(0, 60.)
        self.assertAlmostEqual(testobj._Q(self.tP, \
                                         ), 1.7, 1)
        
        # test for a beta ratio of ~2.28
        # aspect ratio (ar = 1) : wpt/t = 0.5
        self.pmat = Panel_.EPMatl(260.0, 125036.0, 0.3)
        self.smat = self.pmat
        self.tP = Panel_.TPanel(300.0, 6.0, 5.0, 80.0, 8.0, 30.0, 300.0, \
                                    self.pmat, self.smat, 3.5)
        testobj = APS.Allowable_Permanent_Set(0, 3.)
        self.assertAlmostEqual(testobj._Q(self.tP, \
                                         ), 2.6, 1)
        
        # aspect ratio (ar = 1) : wpt/t = 1.0
        self.pmat = Panel_.EPMatl(260.0, 125036.0, 0.3)
        self.smat = self.pmat
        self.tP = Panel_.TPanel(300.0, 6.0, 5.0, 80.0, 8.0, 30.0, 300.0, \
                                    self.pmat, self.smat, 3.5)
        testobj = APS.Allowable_Permanent_Set(0, 6.)
        self.assertAlmostEqual(testobj._Q(self.tP, \
                                         ), 3.3, 1)
        
        # aspect ratio (ar = 0.5) : wpt/t = 0.5
        self.pmat = Panel_.EPMatl(260.0, 125036.0, 0.3)
        self.smat = self.pmat
        self.tP = Panel_.TPanel(300.0, 6.0, 5.0, 80.0, 8.0, 30.0, 600.0, \
                                    self.pmat, self.smat, 3.5)
        testobj = APS.Allowable_Permanent_Set(0, 3.)
        self.assertAlmostEqual(testobj._Q(self.tP, \
                                         ), 1.2, 1)
        
        # aspect ratio (ar = 0.5) : wpt/t = 1.5
        self.pmat = Panel_.EPMatl(260.0, 125036.0, 0.3)
        self.smat = self.pmat
        self.tP = Panel_.TPanel(300.0, 6.0, 5.0, 80.0, 8.0, 30.0, 600.0, \
                                    self.pmat, self.smat, 3.5)
        testobj = APS.Allowable_Permanent_Set(0, 9.)
        self.assertAlmostEqual(testobj._Q(self.tP, \
                                         ), 1.7, 1)
        
        # aspect ratio (ar = 0.2) : wpt/t = 1.5
        self.pmat = Panel_.EPMatl(260.0, 125036.0, 0.3)
        self.smat = self.pmat
        self.tP = Panel_.TPanel(300.0, 6.0, 5.0, 80.0, 8.0, 30.0, 1500.0, \
                                    self.pmat, self.smat, 3.5)
        testobj = APS.Allowable_Permanent_Set(0, 9.)
        self.assertAlmostEqual(testobj._Q(self.tP, \
                                         ), 0.9, 1)
        
    def test_TPanel_trans(self):
        """
        test using TPanel with transverse members (TPanel_trans class)
        """
        testobj = APS.Allowable_Permanent_Set(0, 12.)
        self.assertAlmostEqual(testobj._p_aps(self.PanelT, \
                                           ), 2.4363466633384592, 7)
        
        
    # Run test_Allowable_Permanent_Set
    
suite = unittest.TestLoader().loadTestsFromTestCase(testAllowable_Permanent_Set)
unittest.TextTestRunner(verbosity=2).run(suite)
       
        
                                           
        