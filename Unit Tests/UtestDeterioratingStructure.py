import unittest
from methods.analysis import Corrosion
import deterioratingStructure
import Structures
import TPanel_trans as Panel
from methods.analysis  import CrackDetailVL as CrackDetail

class test_detStructures(unittest.TestCase):
    """
    Unit for TransTPanelDet - deteriorating transverse panel
    """
    
    def setUp(self):
        
        #Set up the transverse T-Panel
        self.matl = Structures.EPMatl(262.0, 70000.0, 0.3)
        self.smatl = Structures.EPMatl(200.0, 70000.0,0.3)
        self.Panel = Panel.TPanel_trans(3000.0, 2000.0, 9.0, 3.0, 5.0, 4.0, 
                                        80.0, 10.0, 30.0, \
                200.0, 10.0, 20.0, 75.0, pmatl=self.matl, smatl=self.smatl, tmatl=self.matl, eta=4.5)
        
        #Set up the corrosion model
        corr = Corrosion.paikCorrosion(0, 'OSH', 'BSLCW', 'DLCF')
        
        #Set up the cost model
        costModel = deterioratingStructure.TPanel_Repair(50000., 10000.)
        fatigue_cost = deterioratingStructure.repairCost(100000.)
        
        #Set up the test classes
        self.testDetStruct = deterioratingStructure.TransTPanelDet(costModel, self.Panel, corr,age=0)
        self.testFatigue = deterioratingStructure.FatigueDetail(fatigue_cost, CrackDetail.CrackDetailVL())
        
    def test_add_detail(self):
        
        self.testDetStruct.add_fatigue_detail(self.testFatigue)
        self.assertTrue(self.testFatigue in self.testDetStruct.getFatigueDetails())
        
    def test_recoat(self):
        
        recoat_cost = self.testDetStruct.recoat(Corrosion.paikCorrosion(5., 'OSH', 'BSLCW', 'DLCF'))
        self.assertEqual(recoat_cost, 10000.)
    
    def test_renew(self):
        
        renew_cost = self.testDetStruct.renew(Corrosion.paikCorrosion(5., 'OSH', 'BSLCW', 'DLCF'))
        self.assertEqual(renew_cost.getRenewalCost(), 50000.)
    
    def test_age(self):
        
        self.testDetStruct.age(2, None)
        panel = self.testDetStruct.getTTPanRef()
        self.assertEqual(self.testDetStruct._age, 2)
        self.assertAlmostEqual(panel.gettp(), 4.9998786, 7)
        self.assertAlmostEqual(panel.gettf(), 9.9998824, 7)
        self.assertAlmostEqual(panel.gettw(), 3.9999068, 7)
    
    def test_needsRepair(self):
        
        panel = self.testDetStruct.getTTPanRef()
        self.testDetStruct._limit_tp = panel.gettp()
        self.testDetStruct.age(2, None)
        
        repair = self.testDetStruct.needsRepair()
        self.assertTrue(repair)
    
    def test_fatigueAge(self):
        loader = deterioratingStructure.Fatigue_Loading(1000, 1000)
        self.testFatigue.age(1, loader)
        self.assertAlmostEqual(self.testFatigue.getCrackProb(), 0.576942023142, 7)
    
    def test_fatigueRenew(self):
        
        crack_cost = self.testFatigue.renew(None)
        self.assertEqual(crack_cost.getRenewalCost(), 100000.)
        
    
    
    
    
    
        
        
#This is currently used to load and run the tests - eventually
#we can replace this with a better test procedure
suite = unittest.TestLoader().loadTestsFromTestCase(test_detStructures)
unittest.TextTestRunner(verbosity=2).run(suite)  