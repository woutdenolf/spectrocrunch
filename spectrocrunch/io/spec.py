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

def zapline_values(start,end,npixels):
    inc = (end-start)/np.float(npixels)
    return start + inc/2 + inc*np.arange(npixels)

def zapline_range(start,end,npixels):
    inc = (end-start)/npixels
    return [start + inc/2, end - inc/2]

def ascan_values(start,end,nsteps):
    inc = (end-start)/np.float(nsteps)
    return start + inc*np.arange(nsteps+1)

def ascan_range(start,end,nsteps):
    return [start,end]

class cmd_parser(object):
    def __init__(self):
        self.fnumber = "(?:[+-]?[0-9]*\.?[0-9]+)"
        self.inumber = "\d+"
        self.blanks = "\s+"
        self.motor = "[a-zA-Z]+"

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
        else:
            return {'name':'unknown'}

    def parsezapimage(self,cmd,name="zapimage"):
        #scanname = cmd.split(' ')[0]
        expr = name + self.blanks +\
                   "("+ self.motor +")" + self.blanks +\
                   "("+ self.fnumber +")" + self.blanks +\
                   "("+ self.fnumber +")" + self.blanks +\
                   "("+ self.inumber +")" + self.blanks +\
                   "("+ self.inumber +")" + self.blanks +\
                   "("+ self.motor +")" + self.blanks +\
                   "("+ self.fnumber +")" + self.blanks +\
                   "("+ self.fnumber +")" + self.blanks +\
                   "("+ self.inumber +")"
        result = re.findall(expr,cmd)
        if len(result)==1:
            return {'name':name,\
                    'motfast':str(result[0][0]),\
                    'startfast':np.float(result[0][1]),\
                    'endfast':np.float(result[0][2]),\
                    'npixelsfast':np.int(result[0][3]),\
                    'time':np.float(result[0][4])/1000,\
                    'motslow':str(result[0][5]),\
                    'startslow':np.float(result[0][6]),\
                    'endslow':np.float(result[0][7]),\
                    'nstepsslow':np.int(result[0][8])}
        else:
            return {'name':'unknown'}

    def parsemesh(self,cmd,name="mesh"):
        #scanname = cmd.split(' ')[0]
        expr = name + self.blanks +\
                   "("+ self.motor +")" + self.blanks +\
                   "("+ self.fnumber +")" + self.blanks +\
                   "("+ self.fnumber +")" + self.blanks +\
                   "("+ self.inumber +")" + self.blanks +\
                   "("+ self.motor +")" + self.blanks +\
                   "("+ self.fnumber +")" + self.blanks +\
                   "("+ self.fnumber +")" + self.blanks +\
                   "("+ self.inumber +")" + self.blanks +\
                   "("+ self.inumber +")"
        result = re.findall(expr,cmd)
        if len(result)==1:
            return {'name':name,\
                    'motfast':str(result[0][0]),\
                    'startfast':np.float(result[0][1]),\
                    'endfast':np.float(result[0][2]),\
                    'nstepsfast':np.int(result[0][3]),\
                    'motslow':str(result[0][4]),\
                    'startslow':np.float(result[0][5]),\
                    'endslow':np.float(result[0][6]),\
                    'nstepsslow':np.int(result[0][7]),\
                    'time':np.float(result[0][8])}
        else:
            return {'name':'unknown'}

    def parsezapenergy(self,cmd,name="zapenergy"):
        # Only SUM is called zapenergy, otherwise its called zapline
        exprSUM = name + self.blanks + \
                   "SUM" + self.blanks +\
                   "("+ self.inumber +")" + self.blanks +\
                   "("+ self.fnumber +")"
        result = re.findall(exprSUM,cmd)

        if len(result)==0:
            exprSUM = name + self.blanks + \
                       "SUM2" + self.blanks +\
                       "("+ self.inumber +")" + self.blanks +\
                       "("+ self.fnumber +")"
            result = re.findall(exprSUM,cmd)

        if len(result)==1:
            return {'name':name,\
                    'repeats':np.int(result[0][0]),\
                    'time':np.float(result[0][1])}
        else:
            return {'name':'unknown'}

    def parsezapline(self,cmd):
        ret = self.parseascan(cmd,name="zapline")
        if "time" in ret:
            ret["time"] /= 1000
        return ret

    def parseascan(self,cmd,name="ascan"):
        expr = name + self.blanks + \
                   "("+ self.motor +")" + self.blanks +\
                   "("+ self.fnumber +")" + self.blanks +\
                   "("+ self.fnumber +")" + self.blanks +\
                   "("+ self.inumber +")" + self.blanks +\
                   "("+ self.fnumber +")"
        result = re.findall(expr,cmd)
        if len(result)==1:
            return {'name':name,\
                    'motfast':str(result[0][0]),\
                    'startfast':np.float(result[0][1]),\
                    'endfast':np.float(result[0][2]),\
                    'nstepsfast':np.int(result[0][3]),\
                    'time':np.float(result[0][4])}
        else:
            return {'name':'unknown'}

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
                pass
        if len(ind)>0:
            data = scan.data[:,ind]
        else:
            data = None

        return data

    def getmotorvalues(self,scannumber,motors):
        """Get start positions for the specified motors
        """
        info = self.getKeyInfo("{:d}.1".format(scannumber))
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

    def getzapimages(self):
        """Get list of all zapimages
        """
        ret = {}
        for k in self.getSourceInfo()["KeyList"]:

            info = self.getKeyInfo(k)
            if info["Command"].startswith("zapimage"):
                
                h = info["Header"]
                h = [i.split(":")[1].strip() for i in h if "#C " in i]
                add = {"specnumber":k.split('.')[0],"scanname":h[1],"scannumber":h[2],"scandir":h[0],"scansize":""}

                result = self.parser.parse(info["Command"])
                if result['name']=="zapimage":
                    add["scansize"] = "{} x {}".format(result['npixelsfast'],result['nstepsslow']+1)

                ret['{}_{}'.format(h[1],int(h[2]))] = add

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

