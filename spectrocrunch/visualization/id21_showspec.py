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

import h5py
import numpy as np
import fabio
import os
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.animation as animation
import warnings

from ..math.utils import logscale
from ..io.spec import spec
from ..utils import instance


class shape_object(object):

    def __init__(self, nframes, noborder=False, noROIborder=False, notitle=False, noimage=False, notext=False, static=True, ROIs=[],
                 cmap="jet", rellim=[[0., 1.]], abslim=[], log=[], offset=[0, 0], name="", figindex=0, xanesnormalize={}):
        self.noborder = noborder
        self.noROIborder = noROIborder
        self.notitle = notitle
        self.noimage = noimage
        self.notext = notext
        self.offset = offset
        self.cmap = cmap
        self.rellim = rellim
        self.abslim = abslim
        self.log = log
        self.img = None
        self.im = None
        self.markers = []
        self.markerstext = []
        self.ROIplts = []
        self.nframes = nframes
        self.static = static
        self.name = name
        self.figindex = figindex
        self.xanesnormalize = xanesnormalize

        self.calcROIs(ROIs)

    def unloadimage(self):
        self.img = None

    def setimage(self, images, shape):
        if shape is None:
            return
        images = [np.zeros(shape) if img is None else img for img in images]
        self.img = np.dstack(images)
        if len(self.img) == 2:
            self.img = self.img[..., np.newaxis]

    def get2ddims(self):
        dim1 = np.array(self.dim1)
        dim2 = np.array(self.dim2)
        shape = 2
        if len(dim1) == 1:
            dim1 = np.array([dim1, dim1])
            shape -= 1
        if len(dim2) == 1:
            dim2 = np.array([dim2, dim2])
            shape -= 1
        return dim1.tolist(), dim2.tolist(), shape

    def extendname(self, name):
        if name == "":
            return self.name
        if self.name == "":
            return name
        return '_'.join((self.name, name))

    def calcROIs(self, ROIs):
        self.nROI = len(ROIs)
        self.nmarkers = self.nROI
        self.ROIs = ROIs
        if self.nROI == 0:
            self.loadimage(0)
            return

        mask = [None]*self.nROI
        for j in range(self.nROI):
            if "minmax" in ROIs[j]:
                self.loadimage(ROIs[j]["frame"])
                n1, n2, nchan = self.img.shape
                minmax = ROIs[j]["minmax"]
                if nchan > 1:
                    tmp = self.img[:, :, 0]
                else:
                    tmp = self.img

                mi = np.nanmin(tmp)
                ma = np.nanmax(tmp)
                mi2 = mi + minmax[0]*(ma-mi)
                ma2 = mi + minmax[1]*(ma-mi)
                mask[j] = (tmp >= mi2) & (tmp <= ma2) & (~np.isnan(tmp))

        for i in range(self.nframes):
            self.loadimage(i)
            n1, n2, nchan = self.img.shape

            for j in range(self.nROI):
                if "range" in ROIs[j]:
                    ROI = ROIs[j]["range"]
                    if nchan > 1:
                        avg = np.nanmean(np.nanmean(
                            self.img[ROI[0]:ROI[1], ROI[2]:ROI[3], 0], axis=0), axis=0)
                        navg = 1
                    else:
                        avg = np.nanmean(np.nanmean(
                            self.img[ROI[0]:ROI[1], ROI[2]:ROI[3]], axis=0), axis=0)
                        navg = len(avg)
                elif "minmax" in ROIs[j]:
                    minmax = ROIs[j]["minmax"]
                    if nchan > 1:
                        tmp = self.img[:, :, 0]
                    else:
                        tmp = self.img
                    avg = np.nanmean(tmp[mask[j]])
                    navg = 1
                else:
                    raise Exception("No known ROI type")

                if i == 0 and j == 0:
                    self.nROI *= navg
                    self.ROIx = np.full(self.nframes, np.nan)
                    self.ROIy = np.full((self.nframes, self.nROI), np.nan)
                    self.ROInames = [""]*self.nROI

                self.ROIx[i] = self.energy
                self.ROIy[i, j*navg:j*navg+navg] = avg
                self.ROInames[j*navg:j*navg +
                              navg] = [self.extendname(ROIs[j]["name"])]*navg

        if self.nROI != 0 and len(self.xanesnormalize) != 0:
            preedge = np.average(
                self.ROIy[self.xanesnormalize['pre'][0]:self.xanesnormalize['pre'][1], :], axis=0)
            postedge = np.average(
                self.ROIy[self.xanesnormalize['post'][0]:self.xanesnormalize['post'][1], :], axis=0)
            self.ROIy -= preedge
            self.ROIy /= postedge-preedge

    def applytransform(self, extent, transform):
        for c in transform:
            if c == 't':  # Transpose
                extent = (extent[2], extent[3], extent[0], extent[1])
            elif c == 'v':  # Flip vertical
                extent = (extent[0], extent[1], extent[3], extent[2])
            elif c == 'h':  # Flip horizontal
                extent = (extent[1], extent[0], extent[2], extent[3])
        return extent

    def rangeind_to_coord(self, i, j, n):
        # arr[i:j]        -> arr[i], arr[j-1]
        # arr[None:None]  -> arr[0], arr[n-1]
        # i < 0           -> n+i
        # j < 0           -> n-1+j
        if i is None:
            i = 0
        elif i < 0:
            i += n

        if j is None:
            j = n
        elif j < 0:
            j += n
        j -= 1

        return i, j

    def getimage(self, transform, origin):
        # Image
        img = self.img
        if img is None:
            nchan = 0
            n1 = 0
            n2 = 0
        else:
            n1, n2, nchan = self.img.shape
            # Transform
            for c in transform:
                if c == 't':  # Transpose
                    img = img.transpose((1, 0, 2))
                elif c == 'v':  # Flip vertical
                    img = img[::-1, :, :]
                elif c == 'h':  # Flip horizontal
                    img = img[:, ::-1, :]

        # Plot range
        dim1, dim2, shape = self.get2ddims()
        if n1 <= 1:
            d1 = 0.
        else:
            d1 = (dim1[1]-dim1[0])/(n1-1.)
        if n2 <= 1:
            d2 = 0.
        else:
            d2 = (dim2[1]-dim2[0])/(n2-1.)

        halfpixel = np.array([-d2, d2, -d1, d1])/2

        # extent = [x0,x1,y0,y1]
        extent = np.array([dim2[0]-origin[1]+self.offset[1],
                           dim2[1]-origin[1]+self.offset[1],
                           dim1[0]-origin[0]+self.offset[0],
                           dim1[1]-origin[0]+self.offset[0]])

        markers = []
        for i in range(self.nROI):
            if "range" not in self.ROIs[i]:
                continue
            ROI = self.ROIs[i]["range"][:]
            ROI[0], ROI[1] = self.rangeind_to_coord(ROI[0], ROI[1], n1)
            ROI[2], ROI[3] = self.rangeind_to_coord(ROI[2], ROI[3], n2)

            ROI = np.array([extent[0]+ROI[2]*d2, extent[0]+ROI[3]
                            * d2, extent[2]+ROI[0]*d1, extent[2]+ROI[1]*d1])
            # ROI = [x0,x1,y0,y1]

            ROI += halfpixel
            ROI = self.applytransform(ROI, transform)
            markers.append(
                {"extent": ROI, "name": self.extendname(self.ROIs[i]["name"])})

        extent += halfpixel
        extent = self.applytransform(extent, transform)

        if self.nROI == 0 and not self.noborder:
            markers.append({"extent": extent, "name": self.extendname("")})

        return img, extent, markers, nchan, shape

    def plot(self, ax, ax2, origin, frame, transform=""):
        bfirst = self.im is None

        # Load the image
        frame = min(frame, self.nframes-1)
        self.loadimage(frame)

        # Get image and border
        img, extent, markers, nchan, shape = self.getimage(transform, origin)

        # Clear previous plot when not static
        if not self.static:
            if len(self.markers) != 0:
                for o in self.markers:
                    if o is not None:
                        o.remove()
                del self.markers
                self.markers = []
            if len(self.markerstext) != 0:
                for o in self.markerstext:
                    if o is not None:
                        o.remove()
                del self.markerstext
                self.markerstext = []
            if self.im is not None:
                o = self.im
                o.remove()
                del o
                self.im = None

        # Plot image (if any)
        if not self.noimage and ax is not None:
            # Images scaling
            for i in range(nchan):
                if len(self.log) != 0:
                    if self.log[i]:
                        img[..., i] = logscale(img[..., i])

                if len(self.abslim) == 0:
                    j = min(len(self.rellim), i)
                    mi = np.nanmin(img[..., i])
                    ma = np.nanmax(img[..., i])
                    d = ma-mi
                    ma = mi + d*self.rellim[j][1]
                    mi += d*self.rellim[j][0]
                else:
                    j = min(len(self.abslim), i)
                    mi = self.abslim[j][0]
                    ma = self.abslim[j][1]
                img[..., i] -= mi
                img[..., i] /= ma-mi
                img[..., i] = np.clip(img[..., i], 0, 1)

            # Plot image (if any)
            if self.im is None:
                if nchan == 1:
                    self.im = ax.imshow(np.squeeze(img), extent=extent, origin='lower',
                                        interpolation='nearest', aspect=1, cmap=cm.get_cmap(self.cmap), vmin=0, vmax=1)
                elif nchan == 2:
                    self.im = ax.imshow(np.concatenate((img, np.zeros(
                        img.shape[0:2])[..., np.newaxis]), axis=2), extent=extent, origin='lower', interpolation='nearest', aspect=1, vmin=0, vmax=1)
                elif nchan == 3:
                    self.im = ax.imshow(
                        img, extent=extent, origin='lower', interpolation='nearest', aspect=1, vmin=0, vmax=1)
            else:
                if nchan == 1:
                    self.im.set_array(np.squeeze(img))
                elif nchan == 2:
                    self.im.set_array(np.concatenate(
                        (img, np.zeros(img.shape[0:2])[..., np.newaxis]), axis=2))
                elif nchan == 3:
                    self.im.set_array(img)

        # Plot ROI's
        bredocolor = True
        if ax2 is not None and self.nROI != 0:
            if len(self.ROIplts) == 0:
                self.ROIplts = [None]*self.nROI
                self.color = [None]*self.nROI
                bredocolor = False
                for i in range(self.nROI):
                    self.ROIplts[i] = ax2.plot(
                        self.ROIx[:frame+1], self.ROIy[:frame+1, i], color=self.color[i], label=self.ROInames[i])[0]
                    self.color[i] = self.ROIplts[i].get_color()
            else:
                for i in range(self.nROI):
                    self.ROIplts[i].set_data(
                        self.ROIx[:frame+1], self.ROIy[:frame+1, i])

        # Plot markers (rectangles, lines, points)
        nmarkers = len(markers)
        if nmarkers != 0 and (not self.static or bfirst) and not self.noROIborder and ax is not None:
            if len(self.markers) == 0:
                self.markers = [None]*nmarkers
                self.markerstext = [None]*nmarkers
                if bredocolor:
                    self.color = [None]*nmarkers
                for i in range(nmarkers):
                    marker = markers[i]["extent"]
                    name = markers[i]["name"]
                    j = i*nchan
                    color = self.color[j]
                    if shape == 0:
                        self.markers[i] = ax.scatter(
                            marker[0], marker[2], marker='o', color=color)[0]
                    elif shape == 1:
                        if marker[0] == marker[1]:
                            self.markers[i] = ax.plot([marker[0], marker[0]], [
                                                      marker[2], marker[3]], color=color, linewidth=2)[0]
                        else:
                            self.markers[i] = ax.plot([marker[0], marker[1]], [
                                                      marker[2], marker[2]], color=color, linewidth=2)[0]
                    else:
                        self.markers[i] = ax.plot([marker[0], marker[0], marker[1], marker[1], marker[0]],
                                                  [marker[2], marker[3], marker[3], marker[2], marker[2]], color=color, linewidth=2)[0]
                    self.color[j] = self.markers[i].get_color()
                    if not self.notext:
                        pos = [np.max(marker[:2]), np.min(marker[2:])]
                        self.markerstext[i] = ax.annotate(
                            name, xy=pos, xytext=pos, color=self.color[j])
            else:
                for i in range(nmarkers):
                    marker = markers[i]["extent"]
                    if shape == 0:
                        self.markers[i].set_data(marker[0], marker[2])
                    elif shape == 1:
                        if marker[0] == marker[1]:
                            self.markers[i].set_data([marker[0], marker[0]], [
                                                     marker[2], marker[3]])
                        else:
                            self.markers[i].set_data([marker[0], marker[1]], [
                                                     marker[2], marker[2]])
                    else:
                        self.markers[i].set_data([marker[0], marker[0], marker[1], marker[1], marker[0]],
                                                 [marker[2], marker[3], marker[3], marker[2], marker[2]])
                    if not self.notext:
                        pos = [np.max(marker[:2]), np.min(marker[2:])]
                        self.markerstext[i].set_position(pos)

        # Free memory image
        # self.unloadimage()


class shape_hdf5(shape_object):

    def __init__(self, filename, subpaths, subindices, coordinatesindex=None, ignoremotorpositions=False, **kwargs):
        self.filename = filename
        self.subpaths = subpaths
        self.subindices = subindices
        self.coordinatesindex = coordinatesindex
        self.ignoremotorpositions = ignoremotorpositions

        shape_object.__init__(self, len(subindices), **kwargs)

    def loadimage(self, frame):
        oh5 = h5py.File(self.filename)

        # Prepare global coordinates
        dim1off = 0.
        dim1name = "samz"
        dim1mult = 1
        dim2off = 0.
        dim2name = "samy"
        dim2mult = 1

        if self.ignoremotorpositions:
            frame2 = 0
        else:
            if self.coordinatesindex is not None:
                frame2 = self.coordinatesindex
            else:
                frame2 = frame

        try:
            ocoord = oh5["stackinfo"]
        except KeyError:
            warnings.warn(
                "\"coordinates\" is deprecated and should be replaced by \"stackinfo\"", DeprecationWarning)
            ocoord = oh5["coordinates"]
        for f in ocoord:
            if f == "samz":
                v = instance.asarray(ocoord[f].value)
                if len(v) == 1:
                    dim1off = v[0]*1000
                else:
                    dim1off = v[self.subindices[frame2]]*1000

                dim1name = "sampz"
                dim1mult = 1
            if f == "sampz":
                v = instance.asarray(ocoord[f].value)
                if len(v) == 1:
                    dim1off = v[0]
                else:
                    dim1off = v[self.subindices[frame2]]

                dim1name = "samz"
                dim1mult = 1000
            if f == "samy":
                v = instance.asarray(ocoord[f].value)
                if len(v) == 1:
                    dim2off = v[0]*1000
                else:
                    dim2off = v[self.subindices[frame2]]*1000

                dim2name = "sampy"
                dim2mult = 1
            if f == "sampy":
                v = instance.asarray(ocoord[f].value)
                if len(v) == 1:
                    dim2off = v[0]
                else:
                    dim2off = v[self.subindices[frame2]]

                dim2name = "samy"
                dim2mult = 1000

        # Get image with axes in micron
        n = len(self.subpaths)
        s = None
        images = [None]*n
        self.energy = 0
        for i in range(n):
            if self.subpaths[i] in oh5:
                ogrp = oh5[self.subpaths[i]]
                odset = ogrp[ogrp.attrs["signal"]]
                self.dim1 = dim1off + ogrp[dim1name].value[[0, -1]]*dim1mult
                self.dim2 = dim2off + ogrp[dim2name].value[[0, -1]]*dim2mult

                tmp = ogrp.attrs["axes"]
                if instance.isstring(tmp):
                    tmp = tmp.split(":")
                else:
                    tmp = [t for t in tmp]

                idim1 = tmp.index(dim1name)
                idim2 = tmp.index(dim2name)
                if idim2 != 0 and idim1 != 0:
                    img = odset[self.subindices[frame], ...]
                    self.energy = ogrp[tmp[0]][self.subindices[frame]]
                elif idim2 != 1 and idim1 != 1:
                    img = odset[:, self.subindices[frame], :]
                    self.energy = ogrp[tmp[1]][self.subindices[frame]]
                else:
                    img = odset[..., self.subindices[frame]]
                    self.energy = ogrp[tmp[2]][self.subindices[frame]]
                if idim1 > idim2:
                    img = img.T
                s = img.shape
                images[i] = img
        self.setimage(images, s)

        print("Energy = {} keV".format(self.energy))
        print(" range = [{}, {}]".format(
            np.nanmin(self.img), np.nanmax(self.img)))

        oh5.close()


class shape_spec(shape_object):

    def __init__(self, specfile, scannumbers, labels=[], olddir="", newdir="", **kwargs):
        self.specfile = specfile
        self.scannumbers = scannumbers
        self.labels = labels
        self.olddir = olddir
        self.newdir = newdir

        shape_object.__init__(self, len(scannumbers), **kwargs)

    def loadimage(self, frame):
        motors = ["samz", "sampz", "samy", "sampy", "Energy MONO"]

        f = spec(self.specfile)

        p = f.getdimensions(self.scannumbers[frame], motors)
        self.dim1 = p["samz"]*1000 + p["sampz"]
        self.dim2 = p["samy"]*1000 + p["sampy"]
        self.energy = p["Energy MONO"]

        # Get images
        n = len(self.labels)
        if n != 0:
            xia = f.getxialocation(self.scannumbers[frame])
            xia["DIRECTORY"] = xia["DIRECTORY"].replace(
                self.olddir, self.newdir)
            xia["ZAP SCAN NUMBER"] = int(xia["ZAP SCAN NUMBER"])
            xia["ZAP IMAGE NUMBER"] = int(xia["ZAP IMAGE NUMBER"])

            images = [None]*n
            s = None
            for i in range(n):
                filename = "%s/%s_%s_%04d_%04d.edf" % (
                    xia["DIRECTORY"], xia["RADIX"],
                    self.labels[i], xia["ZAP SCAN NUMBER"],
                    xia["ZAP IMAGE NUMBER"])
                if os.path.isfile(os.path.join(filename)):
                    images[i] = fabio.open(filename).data
                    s = images[i].shape
                    print("  {}".format(filename))
                else:
                    print("  File not found: {}".format(filename))
            self.setimage(images, s)


class lstPlot(object):

    def __init__(self, lst, limits=[], figsize=None, nframes=None, legendloc=0, transform=""):
        self.lst = lst
        tmp = [l.figindex for l in self.lst if l.figindex is not None]
        if len(tmp) == 0:
            self.naxes = 0
        else:
            self.naxes = max(tmp)+1
        self.nframes = max([l.nframes for l in self.lst])
        if nframes is not None:
            self.nframes = min(self.nframes, nframes)
        self.legendloc = legendloc
        self.transform = transform

        self.prepare_axes(figsize, limits)

    def prepare_axes(self, figsize, limits):
        nROIplts = sum([l.nROI for l in self.lst])

        self.fig = plt.figure(figsize=figsize)

        # Devide images in a grid
        if self.naxes == 0:
            nrow = 1
            ncol = 1
        if self.naxes == 2:
            nrow = 1
            ncol = 2
        else:
            ncol = int(np.ceil(np.sqrt(self.naxes)))
            nrow = ncol

        # Add plot to the grid is needed
        if nROIplts != 0:
            if self.naxes == 0:
                gs = gridspec.GridSpec(nrow, ncol)
                self.ax2 = plt.subplot(gs[0, 0])
            else:
                gs = gridspec.GridSpec(nrow, 2*ncol)
                self.ax2 = plt.subplot(gs[0:nrow, ncol:])
                ncol *= 2
            self.ax2.set_xlabel('Energy (keV)')
        else:
            gs = gridspec.GridSpec(nrow, ncol)
            self.ax2 = None
        self.ax2legend = None

        self.ax = [None]*self.naxes
        cnt = np.zeros(self.naxes)
        lst = [l for l in self.lst if l.figindex is not None]
        for l in lst:
            cnt[l.figindex] += 1
            if self.ax[l.figindex] is None:
                row = l.figindex//ncol
                col = l.figindex % ncol
                self.ax[l.figindex] = plt.subplot(gs[row, col])
                self.ax[l.figindex].set_xlabel('X ($\mu$m)')
                self.ax[l.figindex].set_ylabel('Y ($\mu$m)')
        for l in lst:
            if cnt[l.figindex] == 1 and not l.notitle:
                self.ax[l.figindex].set_title(l.name)

        self.set_axes_limits(limits)

        plt.tight_layout()

    def set_axes_limits(self, limits):
        self.origin = [None]*self.naxes
        for figindex in range(self.naxes):
            if self.ax[figindex] is None:
                continue
            ylimits, xlimits = self.get_axes_limits(figindex)
            self.origin[figindex] = [ylimits[0], xlimits[0]]
            xlimits = [0, xlimits[1]-xlimits[0]]
            ylimits = [0, ylimits[1]-ylimits[0]]

            for c in self.transform:
                if c == 't':  # Transpose
                    ylimits, xlimits = xlimits, ylimits
                elif c == 'v':  # Flip vertical
                    ylimits = ylimits[::-1]
                elif c == 'h':  # Flip horizontal
                    xlimits = xlimits[::-1]

            if figindex < len(limits):
                if len(limits[figindex]) == 2:
                    xlimits2, ylimits2 = limits[figindex]
                    if len(xlimits2) == 2:
                        xlimits = xlimits2
                    if len(ylimits2) == 2:
                        ylimits = ylimits2
            self.ax[figindex].set_xlim(xlimits[0], xlimits[1])
            self.ax[figindex].set_ylim(ylimits[0], ylimits[1])
            self.ax[figindex].relim()

    def get_axes_limits(self, figindex):
        lst = [l for l in self.lst if l.figindex == figindex]
        if len(lst) == 0:
            return []

        dim1 = []
        dim2 = []
        for l in lst:
            tmp1, tmp2, _ = l.get2ddims()
            dim1 += tmp1
            dim2 += tmp2

        return [min(dim1), max(dim1)], [min(dim2), max(dim2)]

    def drawframe(self, frame):
        if frame < 0:
            frame += self.nframes

        for l in self.lst:
            if l.figindex is None:
                ax = None
                origin = [0, 0]
            else:
                ax = self.ax[l.figindex]
                origin = self.origin[l.figindex]
            l.plot(ax, self.ax2, origin, frame, transform="")
        self.update_axes()

    def update_axes(self):
        if self.ax2 is not None:
            self.ax2.relim()
            self.ax2.autoscale(True)
            self.ax2legend = self.ax2.legend(loc=self.legendloc)
        # Does not work well
        # for figindex in range(self.naxes):
        #    if self.ax[figindex] is None:
        #        continue
        #    self.ax[figindex].relim()
        #    self.ax[figindex].autoscale_view(True)


class lstAnimation(lstPlot, animation.TimedAnimation):
    def __init__(self, lst, limits=[], figsize=None, nframes=None, transform="", **kwargs):
        lstPlot.__init__(self, lst, limits=limits, figsize=figsize,
                         nframes=nframes, transform=transform)
        animation.TimedAnimation.__init__(self, self.fig, **kwargs)

    def _draw_frame(self, frame):
        print("Frame {}".format(frame))
        self.drawframe(frame)

    def new_frame_seq(self):
        return iter(range(self.nframes))

    def _init_draw(self):
        for l in self.lst:
            if l.figindex is None:
                ax = None
                origin = [0, 0]
            else:
                ax = self.ax[l.figindex]
                origin = self.origin[l.figindex]
            l.plot(ax, self.ax2, origin, 0)
        if self.ax2legend is not None:
            self.ax2legend.remove()
            del self.ax2legend
            self.ax2legend = None


def animate(lst, save="", interactive=True, dpi=300, format="wmv", **kwargs):
    ani = lstAnimation(lst, **kwargs)
    if save != "":
        # ffmpeg -codecs
        if format == "mp4":
            ani.save("{}.mp4".format(save), bitrate=1024,
                     dpi=dpi, codec="mpeg4")
        elif format == "wmv":
            ani.save("{}.wmv".format(save),
                     bitrate=1024, dpi=dpi, codec="wmv1")
        else:
            raise ValueError("Format {} not supported.".format(format))
    if interactive:
        plt.show()


def plot(lst, frame, save="", interactive=True, dpi=300, **kwargs):
    p = lstPlot(lst, **kwargs)
    p.drawframe(frame)

    if save != "":
        p.fig.savefig("{}.png".format(save), bbox_inches='tight', dpi=dpi)
    if interactive:
        plt.show()
