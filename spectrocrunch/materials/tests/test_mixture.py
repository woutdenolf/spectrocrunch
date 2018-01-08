# -*- coding: utf-8 -*-
#
#   Copyright (C) 2016 European Synchrotron Radiation Facility, Grenoble, France
#
#   Principal author:   Wout De Nolf (wout.de_nolf@esrf.eu)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

import unittest

from ..compoundfromformula import CompoundFromFormula as compound
from ..mixture import Mixture as mixture
from ..types import fractionType
from ...geometries import xrf as xrfgeometries
from ...detectors import xrf as xrfdetectors
from ...common import instance

import numpy as np

class test_mixture(unittest.TestCase):

    def test_molefractions(self):
        compound1 = compound("La2O3",0)
        compound2 = compound("SrO",0)
        compound3 = compound("Co2O3",0)
        compound4 = compound("Fe2O3",0)

        mole = np.array([1.3,2,0.2,3],dtype=float)
        m = mixture([compound1,compound2,compound3,compound4],mole,fractionType.mole)

        # Test compound mole fractions
        nfrac1 = m.molefractions(total=True)
        nfrac2 = mole
        labels = ["La2O3","SrO","Co2O3","Fe2O3"]
        for i in range(len(labels)):
            self.assertAlmostEqual(nfrac1[labels[i]],nfrac2[i])

        # Test elemental mole fractions
        nfrac1 = m.elemental_molefractions(total=True)
        nLa = 2*mole[0]
        nSr = 1*mole[1]
        nCo = 2*mole[2]
        nFe = 2*mole[3]
        nO = 3*mole[0]+1*mole[1]+3*mole[2]+3*mole[3]
        nfrac2 = np.array([nSr,nCo,nFe,nO,nLa])
        labels = ["Sr","Co","Fe","O","La"]
        for i in range(len(labels)):
            self.assertAlmostEqual(nfrac1[labels[i]],nfrac2[i])
        
    def test_addcompound(self):
        c1 = compound("Co2O3",1.5)
        c2 = compound("Fe2O3",1.6)
        m1 = mixture([c1,c2],[2,3],fractionType.mole)
        m2 = mixture([c1],[2],fractionType.mole)
        m2.addcompound(c2,3,fractionType.mole)

        n1 = m1.molefractions(total=True)
        n2 = m2.molefractions(total=True)
        self.assertEqual(n1,n2)

    def test_tocompound(self):
        c1 = compound("Co2O3",1.5)
        c2 = compound("Fe2O3",1.6)
        c3 = mixture([c1,c2],[2,3],fractionType.mole).tocompound("mix")
        c4 = compound("Co4O15Fe6",c3.density)
        for e in c3.elements:
            self.assertEqual(c3.elements[e],c4.elements[e])

    def _spectrum(self,lstlines,thickness):
        geometry = xrfgeometries.factory("sdd120",detectorposition=-15.)
        detector = xrfdetectors.factory("leia",geometry=geometry)
    
        out = {}
        for spectrum,d in zip(lstlines,thickness):
            
            for line,prob in spectrum.probabilities:
                intensity = prob*detector.solidangle*d
                energy = line.energy(scatteringangle=detector.scatteringangle)
                energy = instance.asarray(energy)
                intensity = instance.asarray(intensity)
                for en,inten in zip(energy,intensity):
                    if en in out:
                        out[en] += intensity
                    else:
                        out[en] = intensity
        return out
    
    def _spectrum_equal(self,spectrum1,spectrum2):
        np.testing.assert_allclose(sorted(spectrum1.keys()),sorted(spectrum2.keys()))
        
        for k in spectrum1:
            np.testing.assert_allclose(spectrum1[k],spectrum2[k])
    
    def test_cross_sections(self):
        c = [compound("Co2O3",1.5),compound("Fe2O3",1.6)]
        rhoc = np.array([x.density for x in c])
        dc = np.array([0.1,0.2])
        
        rhom = sum(rhoc*dc)/sum(dc)
        wc = rhoc*dc
        wc /= sum(wc)

        c3 = mixture(c,wc,fractionType.weight)
        c4 = c3.tocompound("mix")

        np.testing.assert_allclose(c3.density,rhom)
        np.testing.assert_allclose(c4.density,rhom)
        
        energy = 10.
        muc = [x.mass_att_coeff(energy) for x in c]
        
        muL = sum(muc*rhoc*dc)
        np.testing.assert_allclose(c3.density*c3.mass_att_coeff(energy)*sum(dc),muL)
        np.testing.assert_allclose(c4.density*c4.mass_att_coeff(energy)*sum(dc),muL)
        
        spectrum = self._spectrum([x.xrayspectrum(energy) for x in c],dc)
        spectrum3 = self._spectrum([c3.xrayspectrum(energy)],[sum(dc)])
        spectrum4 = self._spectrum([c3.xrayspectrum(energy)],[sum(dc)])

        self._spectrum_equal(spectrum3,spectrum4)
        self._spectrum_equal(spectrum,spectrum3)
        
def test_suite_all():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_mixture("test_molefractions"))
    testSuite.addTest(test_mixture("test_addcompound"))
    testSuite.addTest(test_mixture("test_tocompound"))
    testSuite.addTest(test_mixture("test_cross_sections"))
    return testSuite
    
if __name__ == '__main__':
    import sys

    mysuite = test_suite_all()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)
