
### Unit test for midship_section class

from __future__ import division
import unittest
import numpy as np
import matplotlib.pyplot as plt
import sys
import TPanel_trans as Panel
from midship_section import *
import Structures

class test_midship_section(unittest.TestCase):
    '''
    Unit test framework for both the midship_section.py class
    '''
    
    def setUp(self):
        
        
        self.matl = Structures.EPMatl(262.0, 70000.0, 0.3)
        self.smatl = Structures.EPMatl(200.0, 70000.0,0.3)
        self.test_panel_1 = Panel.TPanel_trans(3000.0, 2000.0, 9.0, 3.0, 5.0, 4.0, 80.0, 10.0, 30.0, \
                200.0, 10.0, 20.0, 75.0, pmatl = self.matl, smatl = self.smatl, tmatl = self.matl, eta = 4.5)

        # Materials for TPanel_trans
        self.pmatl = Structures.EPMatl(95e6, 70e9, 0.34)
        self.smatl = Structures.EPMatl(95e6, 70e9, 0.34)
        self.tmatl = Structures.EPMatl(95e6, 70e9, 0.34)
        # Dimensions for TPanel_trans
        self._B = 1
        self._L = 1
        self._nstiff = 4
        self._ntrans = 1
        self._tp = 0.004
        self._tw = 0.005
        self._hw = 0.1133
        self._tf = 0.0064
        self._bf = 0.05
        self._twh = 0.1133
        self._twt = 0.005
        self._tft = 0.0064
        self._tfb = 0.05
        self._sloc = [0,0,0]
        self._ornt = 0
        self._qloc = ['NA','NA','NA']
        self._eta = 0
        # Creation of TPanel_trans
        self.test_panel_2 = TPanel_trans.TPanel_trans(self._B,self._L,\
            self._nstiff,self._ntrans,self._tp,self._tw,self._hw,self._tf,\
            self._bf,self._twh,self._twt,self._tft,self._tfb,self._sloc,\
            self._ornt,self._qloc,self.pmatl,self.smatl,self.tmatl,self._eta)

        matl = Structures.EPMatl(262e6, 70e9, 0.3)
        smatl = Structures.EPMatl(200e6, 70e9, 0.3)
        self.test_panel_3 = Panel.TPanel_trans(1.0, 1.0, 4, 3, 0.005, 0.004, 0.06, 0.008, 0.03, \
                0.2, 0.01, 0.02, 0.075, pmatl=matl, smatl=smatl, tmatl=matl, eta=4.5)
                
    #def test_section_modulii(self):
        
        #structure = Midship_Section([self.test_panel_1,self.test_panel_1], 0)
        #SM = structure.section_modulii()
        #for i in range(len(SM)):
            #self.assertAlmostEqual(SM[i], 15783.62498948, 7)
    
    def test_section_data(self):
        
        structure = Midship_Section([self.test_panel_2], 0)
        data = structure.section_data()
        data_need = (2.303563408195345, 0.0362301749271137, 0.015092000000000001, 44.8609, 3.244455504500486e-05, 0.0003709228298784142, 0.07678102578483174)
        for i in range(len(data)):
            self.assertAlmostEqual(data[i], data_need[i], 7)
        
    def test_production_cost(self):
        
        structure = Midship_Section([self.test_panel_3],0)
        self.assertAlmostEqual(structure.production_cost(), 1497.961716, 7)
    
    def test_fatigue_details(self):
        
        structure = Midship_Section([self.test_panel_3],0)
        structure.generate_fatigue_details(10000.)
        details = structure.grillages[0].getFatigueDetails()
        self.assertEqual(len(details), self.test_panel_3.getnstiff()+2)
        for detail in details:
            self.assertEqual(detail.getRepairObject().getRenewalCost(),10000.)
        
        structure.calculate_fatigue_loading(1000,1000)
        k = 0
        for l in structure.fatigue_loaders:
            for lo in l:
                self.assertEqual(lo.get_fatigue_cycles(), 1000)
                if k != 0 and k != len(structure.fatigue_loaders[0])-1:
                    self.assertAlmostEqual(lo.load(),14.5245386982, 7)
                else:
                    self.assertAlmostEqual(lo.load(), 13.3442381279, 7)
                
                k += 1
        structure.calculate_fatigue_loading(10000,10000)
        for i in range(10):
            structure.age_structure(1, None, structure.fatigue_loaders)
        fatigue_cost = structure.fatigueCost()
        self.assertAlmostEqual(fatigue_cost, 509.17371172201035, 7)
    
    def test_age(self):
        
        structure = Midship_Section([self.test_panel_3],0)
        structure.generate_fatigue_details(10000.)
        details = structure.grillages[0].getFatigueDetails()
        structure.calculate_fatigue_loading(1000,1000)

        structure.age_structure(1, None, structure.fatigue_loaders)
        details = structure.grillages[0].getFatigueDetails()
        self.assertEqual(structure.section_age, 1)
        self.assertEqual(structure.grillages[0]._age, 1)
        
        k = 0
        for detail in details:
            if k == 0 or k == len(details):
                self.assertAlmostEqual(detail.getCrackProb(), 8.83785656515e-112, 7)
            else:
                self.assertAlmostEqual(detail.getCrackProb(), 1.74352432076e-107, 7)
                
            k += 1
    
    def test_schedule_maintenance(self):
        
        structure = Midship_Section([self.test_panel_3],0)
        schedule = [1,2,3,4,5]
        for i in range(len(structure.schedule)):
            self.assertEqual(structure.schedule[i], schedule[i])
    
    def test_repair_needs(self):
                
        structure = Midship_Section([self.test_panel_3],0)
        structure.grillages[0]._limit_tp = self.test_panel_3.gettp() + 1
        needs = structure.repair_needs()
        self.assertEqual(needs[0], 0)
    
    def test_renew_grillages(self):
        
        structure = Midship_Section([self.test_panel_3],0)
        cost = structure.renew_grillages([0])
        self.assertAlmostEqual(cost, 1497.961716, 7)
    
    def test_set_up_smith_collapse(self):
        
        structure = Midship_Section([self.test_panel_3],0)
        collapse_data = structure.set_up_smith_collapse(True)
        correct = [0.018554054054054053, 0.018554054054054053, 0.00010983999999999999, 0.0059199999999999999]
        for i in range(len(collapse_data)):
            if i != 0:
                self.assertAlmostEqual(collapse_data[i], correct[i-1])
    
    def test_curvature_moment_curve(self):
        
        structure = Midship_Section([self.test_panel_3],0)
        curve, moment = structure.curvature_moment_curve(-1e-3,1e-3,nC=4)
        correct_curve = [-0.001, -0.00033333,0.00033333,0.001]
        correct_moment = [15.43416421,5.14472126,-5.14472116,-15.43416348]
        for i in range(len(curve)):
            self.assertAlmostEqual(curve[i], correct_curve[i], 7)
        for i in range(len(moment)):
            self.assertAlmostEqual(moment[i], correct_moment[i], 7)
    
    def test_get_ultimate_moment(self):
        
        ddg = Hull.Hull('testCases/DTMB51_structure.txt', warnings=False)
        structure = Midship_Section(ddg.get_grillage_list(), 0)
        sag, hog = structure.get_ultimate_moment('s'), structure.get_ultimate_moment('h')
        self.assertAlmostEqual(sag, 1555110080.8531084, 7)
        self.assertAlmostEqual(hog, 1486806990.0672166, 7)
    
    
    def test_predict_repairs(self):
        
        structure = Midship_Section([self.test_panel_3],0)
        structure.grillages[0]._limit_tp = self.test_panel_3.gettp() + 1
        predicted = structure.predictRepair(2)
        self.assertEqual(predicted[0], 0)


    def test_perform_scheduled(self):        
        
        structure = Midship_Section([self.test_panel_3],0)
        structure.grillages[0]._limit_tp = self.test_panel_3.gettp() + 1
        structure.schedule_mainteance(1)        
        structure.age_structure(1, None, None)
        structure.end = 30
        structure.service_distribution = None
        structure.service_extension = 0
        maint_cost, repair_cost, repaired = structure.performScheduledMaintenance()
        self.assertAlmostEqual(maint_cost, 1497.961716)
        self.assertEqual(repair_cost, 2.5e5)
        self.assertEqual(repaired[0], 0)
    
    def test_replace_corroded(self):

        structure = Midship_Section([self.test_panel_3],0,coatLife=1)
        structure.grillages[0]._limit_tp = self.test_panel_3.gettp() + 1
        structure.age_structure(2, None, None)
        cost, rtype, repaired = structure.replaceCorroded()    
        self.assertAlmostEqual(cost, 1497.961716)
        self.assertEqual(rtype, 50000.)
        self.assertEqual(repaired[0], 0)
    
    def test_regenerate_section(self):

        ddg = Hull.Hull('testCases/DTMB51.txt', warnings=False)
        structure = Midship_Section(ddg.get_grillage_list(), 0, coatLife=1)
        slimit, hlimit = structure.get_ultimate_moment('s'), structure.get_ultimate_moment('h')
        structure.age_structure(2, None, None)
        cost, repair = structure.regenerateSection(slimit, hlimit)
        self.assertAlmostEqual(cost, 86184.806647186328, 7)
        self.assertEqual(repair, 500000.)
    
    def test_maintenance_cost(self):
        ddg = Hull.Hull('testCases/DTMB51.txt', warnings=False)
        structure = Midship_Section(ddg.get_grillage_list(), 0, coatLife=1)
        structure.generate_fatigue_details(ddg.get_crack_repair_cost())
        total, fatigue, corrosion, scheduled, flat =  structure.maintenance_cost(10, output=False,average_cycles=ddg.get_cycle_increment(),average_fatigue_moment=ddg.get_fatigue_moment())[:5]     
        self.assertAlmostEqual(total, 1063260.4314784333, 7)
        self.assertAlmostEqual(fatigue, 0.231563653243369, 7)
        self.assertAlmostEqual(corrosion, 13260.199914780049, 7)
        self.assertEqual(scheduled, 0.)
        self.assertEqual(flat, 1050000.0)

suite = unittest.TestLoader().loadTestsFromTestCase(test_midship_section)
unittest.TextTestRunner(verbosity=2).run(suite)