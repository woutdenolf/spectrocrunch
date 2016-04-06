# -*- coding: utf-8 -*-
#
#   Copyright (C) 2015 European Synchrotron Radiation Facility, Grenoble, France
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

import os
import json

from spectrocrunch.xrfxas.create_hdf5_imagestacks import create_hdf5_imagestacks as makestack

def createconfig():
    path = os.path.dirname(os.path.abspath(__file__))
    scanname = "fXANES5"

    config = {
            # Input
            "sourcepath": os.path.join(path,"testdata","xrfxanes","id21",scanname),
            "scanname": scanname,
            "scannumbers": [80],
            "counters": ["arr_iodet","arr_idet","arr_absorp1","arr_absorp2","arr_absorp3","xmap_x1","xmap_x1c","xmap_x2","xmap_x2c","xmap_x3","xmap_x3c"],

            # Meta data
            "metacounter": "arr_iodet",
            "stacklabel": "DCM_Energy",
            "fastlabel": "fast",
            "slowlabel": "slow",

            # Deadtime correction
            "dtcor": True,

            # Configuration for fitting
            "detectorcfg": [os.path.join(path,"testdata","xrfxanes","id21",scanname+".cfg")],
            "fit": True,
            "fastfitting": True,
            "addbeforefitting": True,
            "addafterfitting": False,
            "exclude_detectors":[],

            # Output directories
            "outdatapath": os.path.join(path,"testresults",scanname+"_data"),
            "outfitpath": os.path.join(path,"testresults",scanname+"_fit"),
            "hdf5output": os.path.join(path,"testresults",scanname+".h5")
    }

    # Create configuration file
    jsonpath = os.path.join(path,"testresults")
    if not os.path.exists(jsonpath):
        os.makedirs(jsonpath)
    jsonfile = os.path.join(jsonpath,"fXANES5.json")
    with open(jsonfile,'w') as f:
        json.dump(config,f,indent=2)

    return jsonfile, config

if __name__ == '__main__':
    makestack(createconfig())

        
    


    

