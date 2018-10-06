# -*- coding: utf-8 -*-
#
#   Copyright (C) 2018 European Synchrotron Radiation Facility, Grenoble, France
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
import errno
import shutil
import tempfile
import random
import string

def randomstring(size=6, chars=string.ascii_letters + string.digits):
    # Number of combinations: n^size  (default: 62^6)
    return ''.join(random.choice(chars) for _ in range(size))
    
class Copy(object):
    
    def __init__(self,filename,copyname):
        self.filename = filename
        self.copyname = copyname

    def __enter__(self):
        shutil.copy2(self.filename, self.copyname)
        return self.copyname

    def __exit__(self,exc_type, exc_val, exc_tb):
        if exc_type:
            if os.path.isfile(self.copyname):
                os.remove(self.copyname)
        
class TemporaryCopy(object):
    
    def __init__(self,filename,ext=".tmp"):
        self.filename = filename
        self.tmpfilename = None
        self.ext = ext

    def __enter__(self):
        temp_dir = tempfile.gettempdir()
        temp_name = next(tempfile._get_candidate_names())
        self.tmpfilename = os.path.join(temp_dir,temp_name+self.ext)
        shutil.copy2(self.filename, self.tmpfilename)
        return self.tmpfilename

    def __exit__(self,exc_type, exc_val, exc_tb):
        if os.path.isfile(self.tmpfilename):
            os.remove(self.tmpfilename)
        self.tmpfilename = None

class TemporaryFilename(object):
    
    def __init__(self,path,suffix='.tmp',prefix=''):
        self.tmpfilename = os.path.join(path,prefix+randomstring()+suffix)
        
    def __enter__(self):
        return self.tmpfilename

    def __exit__(self,exc_type, exc_val, exc_tb):
        if os.path.exists(self.tmpfilename):
            os.remove(self.tmpfilename)

def mkdir(path):
    try:
        os.makedirs(path)
    except OSError as e:
        if e.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise e

