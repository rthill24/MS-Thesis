# -*- coding: utf-8 -*-
import unittest
import Corrosion
import util.errors as errors

class testPaikCorrosion(unittest.TestCase):
    """
    Unit for Paik Corrosion model
    """

    def testObjectSetupErrors(self):
        self.assertRaises(errors.SetupError, Corrosion.paikCorrosion,
                              10., '#$%^', 'BSLBW', 'BSLBF', 'mm')
        self.assertRaises(errors.SetupError, Corrosion.paikCorrosion,
                              10., 'BSH', '#$%^', 'BSLBF', 'mm')
        self.assertRaises(errors.SetupError, Corrosion.paikCorrosion,
                              10., 'BSH', 'BSLBW', '#$%^', 'mm')
        self.assertRaises(errors.SetupError, Corrosion.paikCorrosion,
                              10., 'BSH', 'BSLBW', 'BSLBF', '#$%^')
        return

    def testBasicCorrosion(self):
        #Boring object that should never have corrosion
        testObj = Corrosion.paikCorrosion(10.)

        self.assertEqual(10., testObj.updatePlateThickness(4, 10.),
                        "Plate thickness before coating life - no corrosion")        
        self.assertEqual(10., testObj.updateWebThickness(4, 10.),
                        "Web thickness before coating life - no corrosion")
        self.assertEqual(10., testObj.updateFlangeThickness(4, 10.),
                        "Flange thickness before coating life - no corrosion")        
        self.assertEqual(10., testObj.updatePlateThickness(12, 10.),
                        "Plate thickness after coating life - no corrosion")        
        self.assertEqual(10., testObj.updateWebThickness(12, 10.),
                        "Web thickness after coating life - no corrosion")
        self.assertEqual(10., testObj.updateFlangeThickness(12, 10.),
                        "Flange thickness after coating life - no corrosion")
        
        testObj2 = Corrosion.paikCorrosion(5., 'OSH', 'BSLCW', 'DLCF')
 
        self.assertEqual(10., testObj2.updatePlateThickness(4, 10.),
                        "Plate thickness before coating life -  corrosion")        
        self.assertEqual(10., testObj2.updateWebThickness(4, 10.),
                        "Web thickness before coating life -  corrosion")
        self.assertEqual(10., testObj2.updateFlangeThickness(4, 10.),
                        "Flange thickness before coating life -  corrosion")
        self.assertAlmostEqual(10 - 7.*0.0607e-3, 
                               testObj2.updatePlateThickness(12, 10.),7, 
                        "Plate thickness after coating life - corrosion")          
        self.assertAlmostEqual(10 - 7.*0.0466e-3, 
                               testObj2.updateWebThickness(12, 10.),7, 
                        "Web thickness after coating life - corrosion")         
        self.assertAlmostEqual(10 - 7.*0.0588e-3 , 
                               testObj2.updateFlangeThickness(12, 10.),7, 
                        "Flange thickness after coating life - corrosion")         
        
        return
        
    def testUnitConversion(self):
        testObj2 = Corrosion.paikCorrosion(5., 'OSH', 'BSLCW', 'DLCF', 
                                           'meters')
        self.assertAlmostEqual(10 - 7.*0.0607e-3, 
                               testObj2.updatePlateThickness(12, 10.),7, 
            "Plate thickness after coating life - corrosion +meters")          
        self.assertAlmostEqual(10 - 7.*0.0466e-3, 
                               testObj2.updateWebThickness(12, 10.),7, 
            "Web thickness after coating life - corrosion+meters")         
        self.assertAlmostEqual(10 - 7.*0.0588e-3 , 
                               testObj2.updateFlangeThickness(12, 10.),7, 
            "Flange thickness after coating life - corrosion+meters")         
        
        testObj3 = Corrosion.paikCorrosion(5., 'OSH', 'BSLCW', 'DLCF', 
                                           'mm')
        self.assertAlmostEqual(10 - 7.*0.0607e-3*1000., 
                               testObj3.updatePlateThickness(12, 10.),7, 
            "Plate thickness after coating life - corrosion +mm")          
        self.assertAlmostEqual(10 - 7.*0.0466e-3*1000., 
                               testObj3.updateWebThickness(12, 10.),7, 
            "Web thickness after coating life - corrosion+mm")         
        self.assertAlmostEqual(10 - 7.*0.0588e-3*1000., 
                               testObj3.updateFlangeThickness(12, 10.),7, 
            "Flange thickness after coating life - corrosion+mm")         

        testObj4 = Corrosion.paikCorrosion(5., 'OSH', 'BSLCW', 'DLCF', 
                                           'in')
        self.assertAlmostEqual(10 - 7.*0.0607e-3*39.3701, 
                               testObj4.updatePlateThickness(12, 10.),7, 
            "Plate thickness after coating life - corrosion +in")          
        self.assertAlmostEqual(10 - 7.*0.0466e-3*39.3701, 
                               testObj4.updateWebThickness(12, 10.),7, 
            "Web thickness after coating life - corrosion+in")         
        self.assertAlmostEqual(10 - 7.*0.0588e-3*39.3701, 
                               testObj4.updateFlangeThickness(12, 10.),7, 
            "Flange thickness after coating life - corrosion+in")  
        
        return

suite = unittest.TestLoader().loadTestsFromTestCase(testPaikCorrosion)
unittest.TextTestRunner(verbosity=2).run(suite)


