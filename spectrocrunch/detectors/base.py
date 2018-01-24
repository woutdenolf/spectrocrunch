# -*- coding: utf-8 -*-
#
#   Copyright (C) 2017 European Synchrotron Radiation Facility, Grenoble, France
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

from ..common import instance

import numpy as np
import fisx


class Material(object):

    DETMATERIALLABEL = "Detector"

    def __init__(self,attenuators=None):
        #TODO: separate attenuators, beamfilters and detector material
        
        if attenuators is None:
            attenuators = {}
        self.attenuators = attenuators

    def addattenuator(self,name,material,thickness):
        self.attenuators[name] = {"material":material,"thickness":thickness}

    def __str__(self):
        s1 = "\n ".join(["{} = {} cm".format(attinfo["material"],attinfo["thickness"]) for attinfo in self.attbefore()])
        s2 = "\n ".join(["{} = {} cm".format(attinfo["material"],attinfo["thickness"]) for attinfo in self.attafter()])

        if s1:
            s1 = "Beam filters:\n {}".format(s1)
        else:
            s1 = "Beam filters: None"
            
        if s2:
            s2 = "Attenuators:\n {}".format(s2)
        else:
            s2 = "Attenuators: None"
        
        if self.hasmaterial:
            s3 = " Material = {} ({} cm)".format(self.material,self.thickness)
        else:
            s3 = " Material = None" 

        return "\n".join([s3,s1,s2])
        
    @property
    def hasmaterial(self):
        return self.DETMATERIALLABEL in self.attenuators
        
    @property
    def material(self):
        return self.attenuators[self.DETMATERIALLABEL]["material"]

    @property
    def thickness(self):
        return self.attenuators[self.DETMATERIALLABEL]["thickness"]

    @thickness.setter
    def thickness(self,value):
        self.attenuators[self.DETMATERIALLABEL]["thickness"] = value

    def addtopymca(self,setup,cfg): 
        for attlabel,attinfo in self.attenuators.items():
            self.addtopymca_attenuator(setup,cfg,attlabel,attinfo["material"],attinfo["thickness"])

    def addtopymca_attenuator(self,setup,cfg,attlabel,material,thickness):
        matname = setup.addtopymca_material(cfg,material)
        cfg["attenuators"][attlabel] = [1,matname,material.density,thickness,1.0]

    def loadfrompymca(self,config):
        self.attenuators = {}
        for attlabel,attinfo in self.attenuators["attenuators"].items():
            self.loadfrompymca_attenuator(config,attlabel,attinfo)

    def loadfrompymca_attenuator(self,config,attlabel,attinfo):
        matname,density,thickness = attinfo[1],attinfo[2],attinfo[3]
        material = self.loadfrompymca_material(config,matname,density)
        self.attenuators[attlabel] = {"material":material,"thickness":thickness}
        
    def loadfrompymca_material(self,config,matname,density):
        if matname in config["materials"]:
            material = mixture.frompymca(config["materials"][matname])
        else:
            material = compoundfromformula.CompoundFromFormula(matname,density)
        return material

    def isbefore(self,name):
        return "BeamFilter" in name

    def isafter(self,name):
        return "BeamFilter" not in name and name != self.DETMATERIALLABEL

    def attbefore(self):
        return [self.attenuators[att] for att in self.attenuators if self.isbefore(att)]
        
    def attafter(self):
        return [self.attenuators[att] for att in self.attenuators if self.isafter(att)]

    def addtofisx(self,setup,cfg):
        for attlabel,attinfo in self.attenuators.items():
            cfg.addtofisx_material(attinfo["material"])

        atts = [[attinfo["material"].pymcaname,attinfo["material"].density,attinfo["thickness"],1.0] for attinfo in self.attbefore()]
        if atts:
            setup.setBeamFilters(atts)

        atts = [[attinfo["material"].pymcaname,attinfo["material"].density,attinfo["thickness"],1.0] for attinfo in self.attafter()]
        if atts:
            setup.setAttenuators(atts)

    def filter_absorbance(self,energy,source=False):
        """Absorbance by filters (before or after sample)
        """
        if source:
            atts = self.attbefore()
        else:
            atts = self.attafter()
        if atts:
            return sum([att["material"].mass_att_coeff(energy)*(att["thickness"]*att["material"].density) for att in atts])
        else:
            return np.zeros_like(energy,dtype=float)

    def filter_transmission(self,energy,source=False):
        """Transmission through filters (before or after sample)
        """
        return np.exp(-self.filter_absorbance(energy,source=source))
        
    def filter_attenuation(self,energy,source=False):
        """Attenuation by filters (before or after sample)
        """
        return 1-self.filter_transmission(energy,source=source)

    def absorbance(self,energy):
        """Absorbance by detector material
        """
        if self.hasmaterial:
            return self.material.mass_att_coeff(energy)*(self.thickness*self.material.density)
        else:
            return np.full_like(3,np.inf,dtype=float)

    def transmission(self,energy):
        """Transmission through detector material
        """
        return np.exp(-self.absorbance(energy))

    def attenuation(self,energy):
        """Attenuation by detector material
        """
        return 1-self.transmission(energy)


class SolidState(Material):

    def __init__(self,ehole=None,**kwargs):
        self.ehole = ehole
        super(SolidState,self).__init__(**kwargs)

    def __str__(self):
        return " Ionization energy = {} eV\n{}".format(self.ehole,super(SolidState,self).__str__())
        

class SolidAngle(SolidState):

    def __init__(self,solidangle=None,**kwargs):
        self._solidangle = solidangle # srad
        super(SolidAngle,self).__init__(**kwargs)

    @property
    def solidangle(self):
        return self._solidangle
    
    def __str__(self):
        return " Solid angle = 4*pi*{} srad\n{}".format(self.solidangle/(4*np.pi),super(SolidAngle,self).__str__())

    def efficiency(self,energysource,energydet):
        """Detector efficiency = S/cos(ain)*T(energysource)*T(energydet)*A(energydet)
            S: solid angle detector
            ain: angle of beam with surface normal
            T: transmission by filters (before sample and detector)
            A: attenuation by the detector crystal
            
        Args:
            energysource: n0
            energydet: n1

        Returns:
            array: n0 x n1
        """
        energysource = instance.asarray(energysource)
        energydet = instance.asarray(energydet)
        
        g = self.solidangle/self.geometry.cosnormin
        T0 = super(SolidAngle,self).filter_transmission(energysource,source=True)
        T1 = super(SolidAngle,self).filter_transmission(energydet,source=False)
        A = self.attenuation(energydet)
        
        # the cosine term is put here for convenience (comes from integration over sample thickness)
        
        return (g*T0)[:,np.newaxis]*(T1*A)[np.newaxis,:]
        
        
class PointSourceCentric(SolidAngle):

    def __init__(self,activearea=None,**kwargs):
        self.activearea = activearea # cm^2
        super(PointSourceCentric,self).__init__(**kwargs)

    def __str__(self):
        return " Active area = {} cm^2\n{}".format(self.activearea,super(PointSourceCentric,self).__str__())

    @classmethod
    def solidangle_calc(cls,activearea=None,distance=None):
        r2 = activearea/np.pi # squared disk radius
        return 2.*np.pi*(1.-(distance/np.sqrt(r2+distance**2.)))

    @property
    def solidangle(self):
        return self.solidangle_calc(activearea=self.activearea,distance=self.geometry.distance)

    @solidangle.setter
    def solidangle(self,value):
        solidanglefrac = value/(4*np.pi)
        if solidanglefrac>=0.5:
            raise ValueError("Solid angle must be < 2.pi")
        r2 = self.activearea/np.pi # squared disk radius
        self.geometry.distance = np.sqrt(r2)*(0.5-solidanglefrac)/np.sqrt((1-solidanglefrac)*solidanglefrac)

    def addtofisx(self,setup,cfg):
        super(PointSourceCentric,self).addtofisx(setup,cfg)
        
        if self.hasmaterial:
            detector = fisx.Detector(self.material.name, self.material.density, self.thickness)
            detector.setActiveArea(self.activearea)
            detector.setDistance(self.geometry.distance)
            #detector.setMaximumNumberOfEscapePeaks(0)
            setup.setDetector(detector)

    def addtopymca(self,setup,cfg): 
        super(PointSourceCentric,self).addtopymca(setup,cfg)
        cfg["concentrations"]["area"] = self.activearea
        
