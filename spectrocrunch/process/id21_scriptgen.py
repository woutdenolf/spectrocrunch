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

import fnmatch
import os
from glob import glob
from PyMca5.PyMcaIO import EdfFile
from spectrocrunch.io.spec import spec

def specinfo(specfolder):
    ret = {}
    files = glob(os.path.join(specfolder,"*.dat"))
    for filename in files:
        try:
            f = spec(filename)
        except:
            continue
        ret[os.path.basename(filename)] = f.getzapimages()
    return ret

def findmapinspec(info,scanname,scannumber):
    name = '{}_{}'.format(scanname,scannumber)
    for specfile in info:
        if name in info[specfile]:
            return specfile,info[specfile][name]["specnumber"],info[specfile][name]["scansize"]
    return "not found",0,""

def extractscaninfo(filename):
    basename = os.path.basename(filename)
    basename,_ = os.path.splitext(basename)

    tmp = basename.split("_")
    scanname = "_".join(tmp[:tmp.index("arr")])
    scannumber = int(tmp[-2])

    f = EdfFile.EdfFile(filename)
    try:
        header = f.GetHeader(0)
        energy = float(header["DCM_Energy"])
    except:
        energy = 0
    return energy,scanname,scannumber

def fluoXASscript(rootdir,name,specfolder="",exclude=["test"],minnum=0,fltr='*fluoXAS*',ctr="*iodet*"):

    if specfolder=="":
        specfolder = os.path.join(rootdir,"spec")
    info = specinfo(specfolder)

    scans = []
    for root, dirnames, filenames in os.walk(rootdir):
        for fluoXASdir in fnmatch.filter(dirnames, fltr):
            path = os.path.join(root, fluoXASdir)

            files = sorted(glob('/'.join((path,ctr))))
            n = len(files)
            if n>=2:
                for i in range(n):
                    energya,scanname,a = extractscaninfo(files[i])
                    if energya != 0: 
                        break
                for i in range(n-1,-1,-1):
                    energyb,scanname,b = extractscaninfo(files[i])
                    if energyb != 0: 
                        break
                
                if not any([e in scanname for e in exclude]) and (b-a+1 > minnum):
                    specfile,specnumbera,scansize = findmapinspec(info,scanname,a)
                    specfile,specnumberb,scansize = findmapinspec(info,scanname,b)
                    scans.append({'sourcepath':path,'scanname':scanname,\
                                    'a':a,'b':b,'ea':energya,'eb':energyb,\
                                    'specfile':specfile,'specnumbera':specnumbera,\
                                    'specnumberb':specnumberb,'scansize':scansize})

    sortlist = [s["scanname"] for s in scans]
    scans = [s for (l,s) in sorted(zip(sortlist,scans))]
    
    with open("{}.py".format(name),"w") as f:
        f.write("import align\n\n")

        for s in scans:
            f.write("def {}_{}():\n".format(name,s["scanname"]))
            f.write("    print(\"Processing {} ...\")\n".format(s["scanname"]))
            f.write("    sourcepath =  [\"{}\"]\n".format(s["sourcepath"]))
            f.write("    scanname = [\"{}\"]\n".format(s["scanname"]))
            f.write("    scannumbers = [range({},{})]\n".format(s["a"],s["b"]+1))
            f.write("    cfgfile = None\n")
            f.write("    alignmethod = None\n")
            f.write("    alignreference = \"arr_absorp2\"\n")
            f.write("    refimageindex = 100\n")
            f.write("    default = \"arr_absorp2\"\n")
            f.write("    try:\n")
            f.write("        align.do(sourcepath,scanname,scannumbers,cfgfile,alignmethod,alignreference,refimageindex,default,\"{}\")\n".format(name))
            f.write("    except:\n")
            f.write("        print(\"{} failed!\")\n".format(s["scanname"]))
            f.write("\n")

        f.write("if __name__ == '__main__':\n")
        for s in scans:
            f.write("    # Number of maps: {} or size {}\n".format(s["b"]-s["a"]+1,s["scansize"]))
            f.write("    # Scan numbers: {} - {}\n".format(s["a"],s["b"]))
            f.write("    # Energy range: {} - {} keV\n".format(s["ea"],s["eb"]))
            f.write("    # Spec: {} - {} in {}\n".format(s["specnumbera"],s["specnumberb"],s["specfile"]))
            f.write("    {}_{}()\n".format(name,s["scanname"]))
            f.write("\n")

    with open("{}.show.py".format(name),"w") as f:
        for s in scans:
            f.write("    #lst += [shape_hdf5(\"{{}}{}/{}/{}.norm.h5\".format(dirname),\\\n".format(name,s["scanname"],s["scanname"]))
            f.write("    #        [\"/counters/{{}}\".format(ctr)],range({}),offset=offset,noimage=noimage,ROIs=ROIs,name=\"{}\")]\n".format(s["b"]-s["a"]+1,s["scanname"]))
            f.write("    #inc += 1\n")
            f.write("\n")

