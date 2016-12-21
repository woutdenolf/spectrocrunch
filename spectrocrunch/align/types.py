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
transformationType = Enum(['translation','rigid','similarity','affine','homography']) # B-spline, moving least-squares

# Affine transformation:
#   Coordinate transformation: X' = R.X + T
#   Change of reference frame: X = R^-1.X' - R^-1.T
#
# Homogeneous coordinates: M = [[R,T],[P,1]]
#   [[Y'],[1]] = M . [[X],[1]]
#   X' = Y'/Z
#
# rotation: R = [[cos,-sin],[sin,cos]]
# reflection: R = [[-1,0],[0,-1]]
# scaling: R = [[s,0],[0,s]]
# aspect: R = [[a,0],[1/a,0]]
# shear: R = [[1,s],[s,1]]
#
# rigid(eucledian):        rotation + translation          dof = 3      R = [[a,sqrt(1-a^2)],[sqrt(1-a^2),a]]  T = [[tx],[ty]]
# similarity:              rigid + scaling                 dof = 4      R = [[a,-b],[b,a]]  T = [[tx],[ty]]
# affine:                  similarity + aspect + shear     dof = 6      R = [[a,b],[c,d]] T = [[tx],[ty]]
# homography(projection):  affine + projection             dof = 8      R = [[a,b],[c,d]] T = [[tx],[ty]] P = [px,py]
#



