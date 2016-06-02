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

    def getdata(self, scannumber, labelnames):
        """
        Args:
            filename(str): file name
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
        
    def getmotorvalues(self,scannumber,motors):
        info = self.getKeyInfo("{:d}.1".format(scannumber))
        names = info["MotorNames"]
        values = info["MotorValues"]
        return [values[names.index(mot)] for mot in motors]

