# -*- coding: utf-8 -*-
#
#   Copyright (C) 2017 European Synchrotron Radiation Facility, Grenoble, France
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

from . import scene
from . import scene_data
from ..io import nxfs

import collections
import pandas as pd


class Image(scene.Image):

    def updatedata(self, **params):

        data, channels, labels = self.datahandle.displaydata(index=self.index)

        params["channels"] = channels
        params["labels"] = labels
        params["axis0name"] = self.datahandle.axis0name
        params["axis1name"] = self.datahandle.axis1name
        params["lim0"] = self.datahandle.axis0values[[0, -1]]
        params["lim1"] = self.datahandle.axis1values[[0, -1]]

        super(Image, self).updatedata(data, **params)


class Text(scene.Text):

    def updatedata(self, **params):
        params["labels"] = self.datahandle.labels
        params["axis0name"] = self.datahandle.axis0name
        params["axis1name"] = self.datahandle.axis1name

        super(Text, self).updatedata(self.datahandle.coordinates0,
                                     self.datahandle.coordinates1, **params)


class ZapRoiMap(Image):

    def __init__(self, filenames, items, plotparams=None, **dataparams):
        """
        Args:
            filename(str|list(str)): list of edf file names
        """
        if plotparams is None:
            plotparams = {}
        self.datahandle = scene_data.EDFStack(filenames, items, **dataparams)

        index = plotparams.pop("channels", None)
        data, channels, labels = self.datahandle.displaydata(index=index)

        plotparams["channels"] = channels
        plotparams["labels"] = plotparams.get("labels", labels)
        plotparams["axis0name"] = plotparams.get(
            "axis0name", self.datahandle.axis0name)
        plotparams["axis1name"] = plotparams.get(
            "axis1name", self.datahandle.axis1name)

        super(ZapRoiMap, self).__init__(data, lim0=self.datahandle.axis0values[[0, -1]],
                                        lim1=self.datahandle.axis1values[[
                                            0, -1]],
                                        **plotparams)


class Nexus(Image):

    def __init__(self, nxgroup, items, plotparams=None, **dataparams):
        """
        Args:
            nxgroup(str): NXdata or NXprocess
            items: list(str)
        """
        if plotparams is None:
            plotparams = {}
        self.datahandle = scene_data.NexusStack(
            nxfs.factory(nxgroup), items, **dataparams)

        self.index = plotparams.pop("channels", None)
        data, channels, labels = self.datahandle.displaydata(index=self.index)

        plotparams["channels"] = channels
        plotparams["labels"] = plotparams.get("labels", labels)
        plotparams["axis0name"] = plotparams.get(
            "axis0name", self.datahandle.axis0name)
        plotparams["axis1name"] = plotparams.get(
            "axis1name", self.datahandle.axis1name)

        super(Nexus, self).__init__(data, lim0=self.datahandle.axis0values[[0, -1]],
                                    lim1=self.datahandle.axis1values[[0, -1]],
                                    **plotparams)


class XanesSpec(Text):

    def __init__(self, filenames, specnumbers, plotparams=None, **dataparams):
        """
        Args:
            filename(str|list(str)): list of edf file names
            specnumbers(list|list(list)): empty list of numbers => all xanes spectra
        """
        if plotparams is None:
            plotparams = {}
        self.output = dataparams.pop("output", None)

        self.datahandle = scene_data.XanesSpec(
            filenames, specnumbers, **dataparams)
        plotparams["labels"] = plotparams.get("labels", self.datahandle.labels)
        plotparams["axis0name"] = plotparams.get(
            "axis0name", self.datahandle.axis0name)
        plotparams["axis1name"] = plotparams.get(
            "axis1name", self.datahandle.axis1name)

        super(XanesSpec, self).__init__(self.datahandle.coordinates0, self.datahandle.coordinates1,
                                        **plotparams)

    def interpolate(self):
        result = collections.OrderedDict()

        k = "{}({:~})".format(self.axis0name,
                              self.datahandle.coordinates0.units)
        result[k] = self.datahandle.coordinates0.magnitude
        k = "{}({:~})".format(self.axis1name,
                              self.datahandle.coordinates1.units)
        result[k] = self.datahandle.coordinates1.magnitude

        for item in self.scene:
            try:
                result.update(item.datahandle.interpolate(
                    self.datahandle.coordinates0, self.datahandle.coordinates1))
            except AttributeError:
                pass

        return result

    def interpolatesave(self):
        if self.output is not None:
            writer = self.output.get("writer", None)
            if writer is not None:
                sheet = self.output.get("sheet", "Sheet1")
                df = pd.DataFrame(self.interpolate(), index=self.labels)
                df.to_excel(writer, sheet)
                worksheet = writer.sheets[sheet]
                worksheet.set_column(0, len(df.columns)+1, 25)
                worksheet.freeze_panes(1, 1)
