# -*- coding: utf-8 -*-

import os, sys
sys.path.insert(1,os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import matplotlib.pyplot as plt
import spectrocrunch.materials.element
import spectrocrunch.materials.compoundfromformula
import spectrocrunch.materials.mixture
import spectrocrunch.materials.types
import spectrocrunch.materials.multilayer
import spectrocrunch.detectors.xrf
import spectrocrunch.geometries.xrf
import spectrocrunch.geometries.source


element = spectrocrunch.materials.element.Element("Ca")

compound1 = spectrocrunch.materials.compoundfromformula.CompoundFromFormula("PbSO4",density=6.29)

compound2 = spectrocrunch.materials.compoundfromformula.CompoundFromFormula("CaSO4",density=2.32)

mixture = spectrocrunch.materials.mixture.Mixture([compound1,compound2],[0.5,0.5],spectrocrunch.materials.types.fractionType.weight)

source = spectrocrunch.geometries.source.factory("synchrotron")
detector = spectrocrunch.detectors.xrf.factory("leia")
geometry = spectrocrunch.geometries.xrf.factory("sxm120",detectorposition=-15.,detector=detector,source=source)


sample = spectrocrunch.materials.multilayer.Multilayer(material=[element,compound1,mixture],\
                                            thickness=[1e-4,1e-4,1e-4],\
                                            geometry = geometry)

energy = [4.9,5]
weights = [5,1]



fig,axs = plt.subplots(2,2,figsize=(10,10))

plt.axes(axs[0][0])
spectrum = element.xrayspectrum(energy,weights=weights,emin=2,emax=5.1)
spectrum.plot()

plt.axes(axs[0][1])
spectrum = compound1.xrayspectrum(energy,weights=weights,emin=2,emax=5.1)
spectrum.plot()

plt.axes(axs[1][0])
spectrum = mixture.xrayspectrum(energy,weights=weights,emin=2,emax=5.1)
spectrum.plot()

plt.axes(axs[1][1])
spectrum = sample.xrayspectrum(energy,weights=weights,emin=2,emax=5.1)
spectrum.geometry = sample.geometry
spectrum.plot()

plt.tight_layout()
plt.show()

