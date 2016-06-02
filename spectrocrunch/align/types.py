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

from spectrocrunch.common.Enum import Enum

dataType = Enum(['h5','h5ext','singlefile','nparray'])
alignType = Enum(['full','usetransfo','calctransfo'])
transformationType = Enum(['translation','orthogonal','rigid','linear similarity','similarity','linear','affine','homography'])

# Affine transformation:
#   Coordinate transformation: X' = R.X + T
#   Change of reference frame: X = R^-1.X' - R^-1.T
#
# rotation: R = [[cos,-sin],[sin,cos]]
# reflection: R = [[cos,sin],[sin,-cos]]
# isotropic scaling: R = [[s,0],[0,s]]
# scaling: R = [[sx,0],[0,sy]]
# shear parallel to x-axis: R = [[1,sx],[0,1]]
# shear parallel to y-axis: R = [[1,0],[sy,1]]
#
# orthogonal:           rotation + reflection
# rigid:                rotation + reflection + translation
# linear similarity:    rotation + reflection + isotropic scaling
# similarity:           rotation + reflection + isotropic scaling + translation
# linear:               rotation + reflection + scaling + shear
# affine:               rotation + reflection + scaling + shear + translation
#
# For an orthogonal basis:
# orthogonal and rigid: R^T.R = id
# rotation + reflection + isotropic scaling: R = s.[[cos,-r.sin],[sin,r.cos]] (with r=1 or -1)

