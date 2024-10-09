"""
Unittest file for class Section
designed to check subfunctions within the class
Created by Tom Devine
September 2012
University of Michigan MSDL
"""

import unittest
import sys, os, fileinput, math, string, numpy, scipy, random
import Structures, TPanel_trans, Plate
import Section

class testSection(unittest.TestCase):
    
    def setUp(self):
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
        self._qloc = [0,0,0]
        self._eta = 0
        # Creation of TPanel_trans
        self._grillage = TPanel_trans.TPanel_trans(self._B,self._L,\
            self._nstiff,self._ntrans,self._tp,self._tw,self._hw,self._tf,\
            self._bf,self._twh,self._twt,self._tft,self._tfb,self._sloc,\
            self._ornt,self._qloc,self.pmatl,self.smatl,self.tmatl,self._eta)
        # Creation of section
        self.section = Section.section()
        self.section.Append_Panels(self._grillage)
        self.section.Explode()
        self.section._upCalcs()
            
    def test_SectionFunctions(self):
        self.assertAlmostEqual(self.section.getSectionArea(),0.007546,7)
        self.assertAlmostEqual(self.section.getYCentroid(),0.0362301749271,7)
        self.assertAlmostEqual(self.section.getSectionYMOI(),1.62222775225e-05, 7)
     
#This is currently used to load and run the tests - eventually
#we can replace this with a better test procedure
suite = unittest.TestLoader().loadTestsFromTestCase(testSection)
unittest.TextTestRunner(verbosity=2).run(suite)
