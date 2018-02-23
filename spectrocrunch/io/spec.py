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

from PyMca5.PyMcaCore import SpecFileDataSource
import re
import numpy as np
from collections import OrderedDict
from ..common.instance import isarray
from .. import ureg
import logging

logger = logging.getLogger(__name__)


def zapline_values(start,end,npixels):
    """Values of the pixel centers
    """
    inc = (end-start)/np.float(npixels)
    return start + inc/2 + inc*np.arange(npixels)

def zapline_range(start,end,npixels):
    """First and last pixel center
    """
    inc = (end-start)/np.float(npixels)
    return [start + inc/2, end - inc/2]
    
def zapline_pixelsize(start,end,npixels):
    """Pixel size
    """
    return (end-start)/np.float(npixels-1)

def zapline_scansize(start,end,npixels):
    """Distance between last and first pixel center
    """
    inc = (end-start)/np.float(npixels)
    return end - start - inc
    
def ascan_values(start,end,nsteps):
    """Values of the pixel centers
    """
    inc = (end-start)/np.float(nsteps)
    return start + inc*np.arange(nsteps+1)

def ascan_range(start,end,nsteps):
    """First and last pixel center
    """
    return [start,end]

def ascan_pixelsize(start,end,nsteps):
    """Pixel size
    """
    return (end-start)/np.float(nsteps)

def ascan_scansize(start,end,npixels):
    """Distance between last and first pixel center
    """
    return end-start
    
def zapimage_submap(header,cmdlabel,scanrange,currentpos,microntounits):
    """
    Args:
        header(dict): edf header
        cmdlabel(str): scan command header key
        scanrange(num): in micron
        currentpos(dict): current motor positions
        microntounits(dict): factors to convert microns to motor units 
    """

    # Old map motor values
    o = cmd_parser()
    result = o.parsezapimage(header[cmdlabel])
    if result["name"]!='zapimage':
        raise RuntimeError("Cannot extract zapimage information from edf header")
    fastvalues = zapline_values(result["startfast"],result["endfast"],result["npixelsfast"])
    slowvalues = ascan_values(result["startslow"],result["endslow"],result["nstepsslow"])
    
    # New map motor values
    pfasta = currentpos[result["motfast"]]-scanrange*microntounits[result["motfast"]]/2.
    pfastb = currentpos[result["motfast"]]+scanrange*microntounits[result["motfast"]]/2.
    pslowa = currentpos[result["motslow"]]-scanrange*microntounits[result["motslow"]]/2.
    pslowb = currentpos[result["motslow"]]+scanrange*microntounits[result["motslow"]]/2.
    
    ifasta = (np.abs(fastvalues-pfasta)).argmin()
    ifastb = (np.abs(fastvalues-pfastb)).argmin()
    islowa = (np.abs(slowvalues-pslowa)).argmin()
    islowb = (np.abs(slowvalues-pslowb)).argmin()
    
    result["startfast"] = fastvalues[ifasta]
    result["endfast"] = fastvalues[ifastb]
    result["npixelsfast"] = abs(ifastb-ifasta)+1
    
    result["startslow"] = slowvalues[islowa]
    result["endslow"] = slowvalues[islowb]
    result["nstepsslow"] = abs(islowb-islowa)
    
    startpositions = {mot:header[mot] for mot in currentpos}
    startpositions[result["motfast"]] = result["startfast"]
    startpositions[result["motslow"]] = result["startslow"]

    d = (result["endfast"]-result["startfast"])/(result["npixelsfast"]-1.)

    scancmd = "zapimage {} {} {} {} {} {} {} {} {} 0".format(\
              result["motfast"],result["startfast"]+d/2.,result["endfast"]+d/2.,result["npixelsfast"],\
              result["motslow"],result["startslow"],result["endslow"],result["nstepsslow"],\
              int(result["time"].to("ms").magitude))
    
    mvcmd = "mv "+" ".join("{} {}".format(mot,pos) for mot,pos in startpositions.items())
    
    if ifasta>ifastb:
        ifast = -1
    else:
        ifast = 1
    if islowa>islowb:
        islow = -1
    else:
        islow = 1
    
    for k,v in currentpos.items():
        if k!=result["motfast"] and k!=result["motslow"]:
            if v!=startpositions[k]:
                logger.warning("Current position of {} ({}) is ignored and set to {}".format(k,v,startpositions[k]))
    
    return scancmd,mvcmd,[[ifasta,ifastb+1,ifast],[islowa,islowb+1,islow]]
        
class cmd_parser(object):
    def __init__(self):
        self.fnumber = "(?:[+-]?[0-9]*\.?[0-9]+)"
        self.inumber = "\d+"
        self.blanks = "\s+"
        self.motor = "[a-zA-Z]+"
        self.motornum = "[a-zA-Z0-9]+"
        
    def parse(self,cmd):
        scanname = cmd.split(' ')[0]
        if scanname=="zapimage":
            return self.parsezapimage(cmd)
        elif scanname=="ascan":
            return self.parseascan(cmd)
        elif scanname=="zapline":
            return self.parsezapline(cmd)
        elif scanname=="zapenergy":
            return self.parsezapenergy(cmd)
        elif scanname=="mesh":
            return self.parsemesh(cmd)
        elif scanname=="puzzle":
            return self.parsepuzzle(cmd)
        else:
            return {'name':'unknown'}

    def match(self,cmd,patterns):
        for pattern in patterns:
            m = re.match(pattern,cmd)
            if m:
                return m
        return m
        
    def patternzapimage(self,name="zapimage"):
        pat1 = "(?P<name>" + name + ")" + self.blanks +\
                   "(?P<motfast>"+ self.motor +")" + self.blanks +\
                   "(?P<startfast>"+ self.fnumber +")" + self.blanks +\
                   "(?P<endfast>"+ self.fnumber +")" + self.blanks +\
                   "(?P<npixelsfast>"+ self.inumber +")" + self.blanks +\
                   "(?P<time>"+ self.inumber +")" + self.blanks +\
                   "(?P<motslow>"+ self.motor +")" + self.blanks +\
                   "(?P<startslow>"+ self.fnumber +")" + self.blanks +\
                   "(?P<endslow>"+ self.fnumber +")" + self.blanks +\
                   "(?P<nstepsslow>"+ self.inumber +")"
                   
        pat2 = "(?P<name>" + name + ")" + self.blanks +\
               "(?P<motfast>"+ self.motor +")" + self.blanks +\
               "(?P<startfast>"+ self.fnumber +")" + self.blanks +\
               "(?P<endfast>"+ self.fnumber +")" + self.blanks +\
               "(?P<npixelsfast>"+ self.inumber +")" + self.blanks +\
               "(?P<motslow>"+ self.motor +")" + self.blanks +\
               "(?P<startslow>"+ self.fnumber +")" + self.blanks +\
               "(?P<endslow>"+ self.fnumber +")" + self.blanks +\
               "(?P<nstepsslow>"+ self.inumber +")" + self.blanks +\
               "(?P<time>"+ self.inumber +")"

        return [pat1,pat2]
    
    def castzapimage(self,m):
        if m:
            result = m.groupdict()
            result['name'] = str(result['name'])
            
            result['motfast'] = str(result['motfast'])
            result['startfast'] = np.float(result['startfast'])
            result['endfast'] = np.float(result['endfast'])
            result['npixelsfast'] = np.int(result['npixelsfast'])
            
            result['motslow'] = str(result['motslow'])
            result['startslow'] = np.float(result['startslow'])
            result['endslow'] = np.float(result['endslow'])
            result['nstepsslow'] = np.int(result['nstepsslow'])
            
            result['time'] = ureg.Quantity(np.float(result['time']),"ms")
        else:
            result = {'name':'unknown'}
        return result

    def patternpuzzle(self,name="puzzle"):
        return ["(?P<name>" + name + ")" + self.blanks +\
                   "(?P<motfast>"+ self.motornum +")" + self.blanks +\
                   "(?P<startfast>"+ self.fnumber +")" + self.blanks +\
                   "(?P<endfast>"+ self.fnumber +")" + self.blanks +\
                   "(?P<npixelsfast>"+ self.inumber +")" + self.blanks +\
                   "(?P<motslow>"+ self.motornum +")" + self.blanks +\
                   "(?P<startslow>"+ self.fnumber +")" + self.blanks +\
                   "(?P<endslow>"+ self.fnumber +")" + self.blanks +\
                   "(?P<nstepsslow>"+ self.inumber +")" + self.blanks +\
                   "(?P<time>"+ self.inumber +")"]
                   
    def castpuzzle(self,m):
        return self.castzapimage(m)
        
    def patternmesh(self,name="mesh"):
        return ["(?P<name>" + name + ")" + self.blanks +\
               "(?P<motfast>"+ self.motor +")" + self.blanks +\
               "(?P<startfast>"+ self.fnumber +")" + self.blanks +\
               "(?P<endfast>"+ self.fnumber +")" + self.blanks +\
               "(?P<nstepsfast>"+ self.inumber +")" + self.blanks +\
               "(?P<motslow>"+ self.motor +")" + self.blanks +\
               "(?P<startslow>"+ self.fnumber +")" + self.blanks +\
               "(?P<endslow>"+ self.fnumber +")" + self.blanks +\
               "(?P<nstepsslow>"+ self.inumber +")" + self.blanks +\
               "(?P<time>"+ self.fnumber +")"]

    def castmesh(self,m):
        if m:
            result = m.groupdict()
            result['name'] = str(result['name'])
            
            result['motfast'] = str(result['motfast'])
            result['startfast'] = np.float(result['startfast'])
            result['endfast'] = np.float(result['endfast'])
            result['nstepsfast'] = np.int(result['nstepsfast'])

            result['motslow'] = str(result['motslow'])
            result['startslow'] = np.float(result['startslow'])
            result['endslow'] = np.float(result['endslow'])
            result['nstepsslow'] = np.int(result['nstepsslow'])
            
            result['time'] = ureg.Quantity(np.float(result['time']),"s")
        else:
            result = {'name':'unknown'}
        return result
        
    def patternzapenergy(self,name="zapenergy"):
        pat1 = "(?P<name>" + name + ")" + self.blanks + \
               "SUM" + self.blanks +\
               "(?P<repeats>"+ self.inumber +")" + self.blanks +\
               "(?P<time>"+ self.fnumber +")"
        pat2 = "(?P<name>" + name + ")" + self.blanks + \
               "SUM2" + self.blanks +\
               "(?P<repeats>"+ self.inumber +")" + self.blanks +\
               "(?P<time>"+ self.fnumber +")"
        return [pat1,pat2]
    
    def castzapenergy(self,m):
        if m:
            result = m.groupdict()
            result['name'] = str(result['name'])
            
            result['repeats'] = np.int(result['repeats'])
            result['time'] = ureg.Quantity(np.float(result['time']),"ms")
        else:
            result = {'name':'unknown'}
        return result
        
    def patternzapline(self,name="zapline"):
        return ["(?P<name>" + name + ")" + self.blanks + \
               "(?P<motfast>"+ self.motor +")" + self.blanks +\
               "(?P<startfast>"+ self.fnumber +")" + self.blanks +\
               "(?P<endfast>"+ self.fnumber +")" + self.blanks +\
               "(?P<npixelsfast>"+ self.inumber +")" + self.blanks +\
               "(?P<time>"+ self.fnumber +")"]
  
    def castzapline(self,m):
        if m:
            result = m.groupdict()
            result['name'] = str(result['name'])
            
            result['motfast'] = str(result['motfast'])
            result['startfast'] = np.float(result['startfast'])
            result['endfast'] = np.float(result['endfast'])
            result['npixelsfast'] = np.int(result['npixelsfast'])
            
            result['time'] = ureg.Quantity(np.float(result['time']),"ms")
        else:
            result = {'name':'unknown'}
        return result
        
    def patternascan(self,name="ascan"):
        return ["(?P<name>" + name + ")" + self.blanks + \
               "(?P<motfast>"+ self.motor +")" + self.blanks +\
               "(?P<startfast>"+ self.fnumber +")" + self.blanks +\
               "(?P<endfast>"+ self.fnumber +")" + self.blanks +\
               "(?P<nstepsfast>"+ self.inumber +")" + self.blanks +\
               "(?P<time>"+ self.fnumber +")"]
               
    def castascan(self,m):
        if m:
            result = m.groupdict()
            result['name'] = str(result['name'])
            
            result['motfast'] = str(result['motfast'])
            result['startfast'] = np.float(result['startfast'])
            result['endfast'] = np.float(result['endfast'])
            result['nstepsfast'] = np.int(result['nstepsfast'])
            
            result['time'] = ureg.Quantity(np.float(result['time']),"s")
        else:
            result = {'name':'unknown'}
        return result
        
    def parsezapimage(self,cmd,name="zapimage"):
        patterns = self.patternzapimage(name=name)
        m = self.match(cmd,patterns)
        return self.castzapimage(m)
        
    def parsepuzzle(self,cmd,name="puzzle"):
        patterns = self.patternpuzzle(name=name)
        m = self.match(cmd,patterns)
        return self.castpuzzle(m)
    
    def parsemesh(self,cmd,name="mesh"):
        patterns = self.patternmesh(name=name)
        m = self.match(cmd,patterns)
        return self.castmesh(m)
        
    def parsezapenergy(self,cmd,name="zapenergy"):
        # Only SUM is called zapenergy, otherwise its called zapline
        patterns = self.patternzapenergy(name=name)
        m = self.match(cmd,patterns)
        return self.castzapenergy(m)
        
    def parsezapline(self,cmd,name="zapline"):
        patterns = self.patternzapline(name=name)
        m = self.match(cmd,patterns)
        return self.castzapline(m)
        
    def parseascan(self,cmd,name="ascan"):
        patterns = self.patternascan(name=name)
        m = self.match(cmd,patterns)
        return self.castascan(m)
        

class edfheader_parser(object):

    def __init__(self,fastlabel=None,slowlabel=None,\
                    timelabel=None,timeunit=None,\
                    energylabel=None,energyunit=None,\
                    speclabel=None):
        self.fastlabel = fastlabel
        self.slowlabel = slowlabel
        self.timelabel = timelabel
        if timeunit is None:
            self.timeunit = "s"
        else:
            self.timeunit = timeunit
        self.energylabel = energylabel
        if energyunit is None:
            self.energyunit = "keV"
        else:
            self.energyunit = energyunit
        self.speclabel = speclabel
        self.specparser = cmd_parser()

    def parse(self,header):
        try:
            r = self.specparser.parse(str(header[self.speclabel]))
        except:
            r = {}
            
        try:
            r['motfast'] = str(header[self.fastlabel+"_mot"])
            r['startfast'] = np.float(header[self.fastlabel+"_start"])
            r['endfast'] = np.float(header[self.fastlabel+"_end"])
            r['npixelsfast'] = np.int(header[self.fastlabel+"_nbp"])
        except:
            r['name'] = 'unknown'
            
        try:
            r['motslow'] = str(header[self.slowlabel+"_mot"])
            r['startslow'] = np.float(header[self.slowlabel+"_start"])
            r['endslow'] = np.float(header[self.slowlabel+"_end"])
            r['nstepsslow'] = np.int(header[self.slowlabel+"_nbp"])
        except:
            r['name'] = 'unknown'

        try:
            r['time'] = ureg.Quantity(np.float(header[self.timelabel]),self.timeunit)
        except:
            pass

        try:
            r['energy'] = ureg.Quantity(np.float(header[self.energylabel]),self.energyunit)
        except:
            pass
                     
        if 'name' not in r:
            r['name'] = 'zapimage'
        
        return r
        
        
class spec(SpecFileDataSource.SpecFileDataSource):
    """An interface to a spec file
    """

    def __init__(self, filename):
        """
        Args:
            filename(str): file name

        Raises:
            ValueError: file cannot be loaded
        """
        SpecFileDataSource.SpecFileDataSource.__init__(self,filename)
        self.parser = cmd_parser()

    def getdata(self, scannumber, labelnames):
        """Get counters + info on saved data

        Args:
            scannumber(num): spec scan number
            labelnames(list(str)): list of counter names

        Returns:
            (np.array, dict): counters (nCounter x npts), info on collected data

        Raises:
            KeyError: scan number doesn't exist
            TypeError: unknown scan type
            ValueError: no data corresponding to the labelnames
        """

        # Get data object
        scan = self.getDataObject("{:d}.1".format(scannumber))

        # Extract xia data names
        info = {"DIRECTORY":"", "RADIX":"", "ZAP SCAN NUMBER":"", "ZAP IMAGE NUMBER":""}
        for s in scan.info["Header"]:
            if s.startswith("#C "):
                tmp = s[2:].split(":")
                if len(tmp)==2:
                    tmp = [s.strip() for s in tmp]
                    if tmp[0] in info:
                        info[tmp[0]] = tmp[1]

        # Extract data
        ind = []
        labels = scan.info["LabelNames"]
        for i in range(len(labelnames)):
            try:
                j = labels.index(labelnames[i])
                ind.append(j)
            except:
                pass
        if len(ind)>0:
            data = scan.data[:,ind]
        else:
            data = None

        return data,info
        
    def getdata2(self, scannumber, labelnames):
        """Get counters

        Args:
            scannumber(num): spec scan number
            labelnames(list(str)): list of labels

        Returns:
            np.array: counters (nCounter x npts)

        Raises:
            KeyError: scan number doesn't exist
            TypeError: unknown scan type
            ValueError: no data corresponding to the labelnames
        """

        # Get data object
        try:
            scan = self.getDataObject("{:d}.1".format(scannumber))
        except TypeError:
            raise TypeError("Error loading scan number {}".format(scannumber))
        
        # Extract data
        ind = []
        labels = scan.info["LabelNames"]
        for i in range(len(labelnames)):
            try:
                j = labels.index(labelnames[i])
                ind.append(j)
            except:
                raise RuntimeError("Label not in list: {}".format(labels))
        if len(ind)>0:
            data = scan.data[:,ind]
        else:
            data = None

        return data

    def haslabel(self,scannumber,label):
        try:
            scan = self.getDataObject("{:d}.1".format(scannumber))
        except TypeError:
            raise TypeError("Error loading scan number {}".format(scannumber))
        return label in scan.info["LabelNames"]

    def getmotorvalues(self,scannumber,motors):
        """Get start positions for the specified motors
        """
        try:
            info = self.getKeyInfo("{:d}.1".format(scannumber))
        except:
            msg = "Failed to retrieve scan number {} from {}".format(scannumber,self.sourceName)
            raise KeyError(msg)
            
        names = info["MotorNames"]
        values = info["MotorValues"]
        return [values[names.index(mot)] for mot in motors]

    def getxialocation(self,scannumber):
        """Get info on saved data
        """
        info = self.getKeyInfo("{:d}.1".format(scannumber))

        ret = {"DIRECTORY":"", "RADIX":"", "ZAP SCAN NUMBER":"", "ZAP IMAGE NUMBER":""}
        for s in info["Header"]:
            if s.startswith("#C "):
                tmp = s[2:].split(":")
                if len(tmp)==2:
                    tmp = [s.strip() for s in tmp]
                    if tmp[0] in ret:
                        ret[tmp[0]] = tmp[1]
        return ret

    def scancommand(self,scannumber):
        info = self.getKeyInfo("{:d}.1".format(scannumber))
        return self.parser.parse(info["Command"])

    def getdimensions(self,scannumber,motors):
        """Get scan dimensions for the specified motors
        """
        info = self.getKeyInfo("{:d}.1".format(scannumber))
        names = info["MotorNames"]
        values = info["MotorValues"]
        cmd = info["Command"]

        # Parse motor positions
        ret = {mot:values[names.index(mot)] for mot in motors}

        # Parse command
        result = self.parser.parse(cmd)
        if result['name']=="zapimage":
            if result['motfast'] in motors:
                ret[result['motfast']] = np.array(zapline_range(result['startfast'],result['endfast'],result['npixelsfast']))
            if result['motslow'] in motors:
                ret[result['motslow']] = np.array(ascan_range(result['startslow'],result['endslow'],result['nstepsslow']))
        elif result['name']=="ascan":
            if result['motfast'] in motors:
                ret[result['motfast']] = np.array(ascan_range(result['startfast'],result['endsfast'],result['nstepsfast']))
        elif result['name']=="mesh":
            if result['motfast'] in motors:
               ret[result['motfast']] = np.array(ascan_range(result['startfast'],result['endfast'],result['npixelsfast']))
            if result['motslow'] in motors:
                ret[result['motslow']] = np.array(ascan_range(result['startslow'],result['endslow'],result['nstepsslow']))

        return ret

    @staticmethod
    def addunit(value,mot,units):
        u = units.get(mot,None)
        if units is None:
            return value
        else:
            return ureg.Quantity(value,u)

    def getimages(self,motors=None,units={}):
        """Get list of all images
        """
        ret = OrderedDict()
        for k in self.getSourceInfo()["KeyList"]:

            info = self.getKeyInfo(k)
            cmd = info["Command"]
            if cmd.startswith("zapimage") or cmd.startswith("puzzle"):# or cmd.startswith("mesh"): NO LINK TO DATA FOR MESH
                result = self.parser.parse(cmd)
                if result["name"] == "unknown":
                    continue

                add = {}
                add["specnumber"] = k.split('.')[0]
                add["scansize"] = "{} x {}".format(result['npixelsfast'],result['nstepsslow']+1)

                ffast = lambda op: self.addunit(op(result['startfast'],result['endfast'],result['npixelsfast']),result['motfast'],units)
                fslow = lambda op: self.addunit(op(result['startslow'],result['endslow'],result['nstepsslow']),result['motslow'],units)
                
                ufast = units.get(result['motfast'],None)
                uslow = units.get(result['motslow'],None)

                add["scansize_units"] = "{} x {}".format(ffast(zapline_scansize),ffast(ascan_scansize))
                
                add["pixelsize"] = "{} x {}".format(ffast(zapline_pixelsize),ffast(ascan_pixelsize))

                add["puzzle"] = result['name']=="puzzle"
                add["mesh"] = result['name']=="mesh"
                
                if motors is None:
                    add["motors"] = []
                else:
                    names = info["MotorNames"]
                    values = info["MotorValues"]
                    add["motors"] = [values[names.index(mot)] if mot in names else np.nan for mot in motors]

                if result["name"] != "mesh":
                    h = info["Header"]
                    h = [i.split(":")[1].strip() for i in h if "#C " in i]
                    add["scandir"] = h[0]
                    add["scanname"] = h[1]
                    add["scannumber"] = h[2]

                    ret['{}_{}'.format(h[1],int(h[2]))] = add

        return ret

    def getregexscans(self,pattern):
        ret = OrderedDict()
        fmt = re.compile(pattern)
        for k in self.getSourceInfo()["KeyList"]:
            info = self.getKeyInfo(k)
            m = fmt.match(info["Command"])
            if m:
                ret[int(k.split('.')[0])] = m.groupdict()
        return ret
        
    def extractxanesinfo(self,skip=None,nrmin=None,nrmax=None):
        """Get list of all ID21 XANES
        """
        ret = []

        lst = self.getSourceInfo()["KeyList"]
        if skip is None:
            skip = []

        for k in lst:
            scannumber = int(k.split('.')[0])
            if nrmin is not None:
                if scannumber < nrmin:
                    continue
            if nrmax is not None:
                if scannumber > nrmax:
                    continue
            if scannumber in skip:
                continue

            info = self.getKeyInfo(k)
            if info["Command"].startswith("zapline mono "):
                if info["Lines"]==0:
                    continue
                ret += [{"scannumber":scannumber,"repeats":1}]
                
            elif info["Command"].startswith("zapenergy SUM "):
                if info["Lines"]==0:
                    continue
                ret += [{"scannumber":scannumber,"repeats":int(info["Command"].split(' ')[2])}]
            elif info["Command"].startswith("zapenergy SUM2 "):
                # This used to be the sum fo the repeats, energy interpolated
                if info["Lines"]==0:
                    continue
                ret[-1] = {"scannumber":scannumber,"repeats":int(info["Command"].split(' ')[2])}

        return ret

    def extractxanesginfo(self,keepsum=False,sumingroups=False,keepindividual=False,skip=None,nrmin=None,nrmax=None):
        """Get list of all ID21 XANES, grouping repeats
        """
        data = self.extractxanesinfo(skip=skip,nrmin=nrmin,nrmax=nrmax)

        ret = []

        bproc = [True]*len(data)

        # Groups: [rep1,rep2,...,(sum)]
        for i in range(len(data)):
            if not bproc[i]:
                continue

            n = data[i]["repeats"]
            if n>1:
                # [rep1,rep2,....]
                i0 = i
                i1 = i
                while (data[i0-1]["repeats"] if i0>0 else 0) ==1:
                    i0 -= 1
                rng = range(max(i0,i-n),i1)
                add = [data[k]["scannumber"] for k in rng if bproc[k]]
                for l in rng:
                    bproc[l] = keepindividual

                # [rep1,rep2,...,(sum)]
                if keepsum:
                    if sumingroups:
                        add += [data[i1]["scannumber"]]
                        bproc[i1] = keepindividual
                else:
                    bproc[i1] = False

                ret += [add]

        # Add each scan as an individual group
        for i in range(len(data)):
            if bproc[i]:
                ret += [[data[i]["scannumber"]]]

        return ret

    def extractxanes(self,scannumbers,labelnames):
        """Get list of specific ID21 XANES with counters
        """
        ret = {}

        for scannumber in scannumbers:
            k = "{:d}.1".format(scannumber)

            info = self.getKeyInfo(k)
            if info["Command"].startswith("zapline mono"):
                result = self.parser.parse(info["Command"])
                if result['name']!="zapline" or info["Lines"]==0:
                    continue
                ret[scannumber] = {"repeats":1, \
                                    "data":self.getdata2(scannumber, labelnames),\
                                    "time":result["time"],\
                                    "labels":["{}.{}".format(s,scannumber) for s in labelnames]}
            elif info["Command"].startswith("zapenergy SUM"):
                result = self.parser.parse(info["Command"])
                if result['name']!="zapenergy" or info["Lines"]==0:
                    continue
                ret[scannumber] = {"repeats":int(info["Command"].split(' ')[2]),\
                                    "data":self.getdata2(scannumber, labelnames),\
                                    "time":result["time"],\
                                    "labels":["{}.{}".format(s,scannumber) for s in labelnames]}

        return ret

