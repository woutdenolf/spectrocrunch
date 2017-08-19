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
        else:
            return {'name':'unknown'}

    def parsezapimage(self,cmd,name="zapimage"):
        scanname = cmd.split(' ')[0]
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
                    'time':np.float(result[0][4]),\
                    'motslow':str(result[0][5]),\
                    'startslow':np.float(result[0][6]),\
                    'endslow':np.float(result[0][7]),\
                    'nstepsslow':np.int(result[0][8])}
        else:
            return {'name':'unknown'}

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
                    'nstepsfast':np.int(result[0][2]),\
                    'time':np.float(result[0][2])}
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
        """
        Args:
            scannumber(num): spec scan number
            labelnames(list(str)): list of labels

        Returns:
            (np.array, dict): first dimension are the labelnames, information on the real data

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
        """
        Args:
            scannumber(num): spec scan number
            labelnames(list(str)): list of labels

        Returns:
            np.array: first dimension are the labelnames

        Raises:
            KeyError: scan number doesn't exist
            TypeError: unknown scan type
            ValueError: no data corresponding to the labelnames
        """

        # Get data object
        scan = self.getDataObject("{:d}.1".format(scannumber))

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

        return ret

    def getzapimages(self):
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

    def extractxanesinfo(self):
        ret = {}
        for k in self.getSourceInfo()["KeyList"]:
            scannumber = int(k.split('.')[0])
            info = self.getKeyInfo(k)
            if info["Command"].startswith("zapline mono"):
                ret[scannumber] = {"scannumber":scannumber,"repeats":1}
            elif info["Command"].startswith("zapenergy SUM"):
                ret[scannumber] = {"scannumber":scannumber,"repeats":int(info["Command"].split(' ')[2])}
        return ret

    def extractxanes(self,scannumbers,labelnames):
        ret = {}
        for scannumber in scannumbers:
            k = "{:d}.1".format(scannumber)

            info = self.getKeyInfo(k)
            if info["Command"].startswith("zapline mono"):
                ret[scannumber] = {"repeats":1, \
                                    "data":self.getdata2(scannumber, labelnames),\
                                    "labels":["{}.{}".format(s,scannumber) for s in labelnames]}
            elif info["Command"].startswith("zapenergy SUM"):
                ret[scannumber] = {"repeats":int(info["Command"].split(' ')[2]),\
                                    "data":self.getdata2(scannumber, labelnames),\
                                    "labels":["{}.{}".format(s,scannumber) for s in labelnames]}

        return ret

