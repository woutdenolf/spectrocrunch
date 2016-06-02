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

"""Initialize cctbx library. To use it, add to following line to the top of yout script:
    execfile("initcctbx.py")
"""

import os, sys

LIBTBX_BUILD = "/usr/local/cctbx-dev-643/build"
os.environ["LIBTBX_BUILD"] = LIBTBX_BUILD

FONTCONFIG_PATH="/usr/local/cctbx-dev-643/base/etc/fonts"
os.environ["FONTCONFIG_PATH"] = FONTCONFIG_PATH

if "LD_LIBRARY_PATH" not in os.environ:
    os.environ["LD_LIBRARY_PATH"] = ":".join([os.path.abspath(os.path.join(LIBTBX_BUILD,'lib')), \
                                os.path.abspath(os.path.join(LIBTBX_BUILD,'bin')), \
                                os.path.abspath(os.path.join(LIBTBX_BUILD,'..','base','lib'))])
    try:
        os.execl(sys.executable, 'python', __file__, *sys.argv[1:])
        sys.exit()
    except Exception, exc:
        print 'Failed re-exec:', exc
        sys.exit(1)

sys.path.append(os.path.abspath(os.path.join(LIBTBX_BUILD,'lib')))
sys.path.append(os.path.abspath(os.path.join(LIBTBX_BUILD,'..','modules','cctbx_project','libtbx','pythonpath')))
sys.path.append(os.path.abspath(os.path.join(LIBTBX_BUILD,'..','modules','cctbx_project','boost_adaptbx')))
sys.path.append(os.path.abspath(os.path.join(LIBTBX_BUILD,'..','modules')))
sys.path.append(os.path.abspath(os.path.join(LIBTBX_BUILD,'..','modules','cctbx_project')))
