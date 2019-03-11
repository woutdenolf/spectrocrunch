# -*- coding: utf-8 -*-
#
#   Copyright (C) 2018 European Synchrotron Radiation Facility, Grenoble, France
#
#   Principal author:   Wout De Nolf (wout.de_nolf@esrf.eu)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the 'Software'), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import logging
import os

from . import nxprocess
from ..visualization import scene
from ..visualization import scene_view
from ..visualization import defaultsettings
from ..utils import instance
from ..utils import units
from ..io import xiaedf

logger = logging.getLogger(__name__)


def scene_modifyaxes(sc, parameters):
    if parameters.get('originzero', False):
        sc.originzero()
    if parameters.get('transpose', False):
        sc.transpose(True)
    if parameters.get('flipx', False):
        sc.flipx(increasing=parameters.get('incx', None))
    elif parameters.get('incx', False):
        sc.increasingx(True)
    if parameters.get('flipy', False):
        sc.flipy(increasing=parameters.get('incy', None))
    elif parameters.get('incy', False):
        sc.increasingy(True)
    if parameters.get('aspectratio', None) is not None:
        sc.aspectratio = parameters['aspectratio']


def scene_globalscaling(sc, parameters):
    rlo = parameters.get('rlo', None)
    rhi = parameters.get('rhi', None)
    if rlo is not None or rhi is not None:
        vmin, vmax = sc.vminmax
        dv = np.array(vmax)-np.array(vmin)

        if rlo is None:
            vmin = None
        else:
            rlo, func = instance.asarrayf(rlo)
            vmin = vmin + dv*func(rlo)
        if rhi is None:
            vmax = None
        else:
            rhi, func = instance.asarrayf(rhi)
            vmax = vmin + dv*func(rhi)
        sc.scale(vmin=vmin, vmax=vmax)

    alo = parameters.get('alo', None)
    ahi = parameters.get('ahi', None)
    if alo is not None or ahi is not None:
        sc.scale(vmin=alo, vmax=ahi)


def createwriter(parameters, filename):
    filename, ext = os.path.splitext(filename)
    filename = filename+'.xlsx'
    if filename not in parameters['writers']:
        writer = pd.ExcelWriter(filename)
        parameters['writers'][filename] = writer
    return parameters['writers'][filename]


def savewriters(parameters):
    for filename, writer in parameters['writers'].items():
        writer.save()
        logger.info('Saved {}'.format(filename))


def savefigures(filename, parameters):
    figs = [plt.figure(i) for i in plt.get_fignums()]
    if len(figs) > 1:
        filename, ext = os.path.splitext(filename)
        filename = '{}_{{:02d}}{}'.format(filename, ext)
    for i, fig in enumerate(figs):
        name = filename.format(i+1)
        if not os.path.exists(os.path.dirname(name)):
            os.makedirs(os.path.dirname(name))
        saveparams = parameters.get('saveparams',{})
        fig.savefig(name, **saveparams)
        logger.info('Saved {}'.format(name))


def closefigures():
    for i in plt.get_fignums():
        plt.close(i)


def specname(info):
    filename = "spec{:02d}".format(info["spec"])
    radix = "_".join((info["sample"],filename))
    filename = os.path.join(info["sample"],radix,"{}.dat".format(radix))
    return os.path.join(info["root"],filename)


def extractname(item):
    if instance.isarray(item):
        return item[0]
    else:
        return item


def ctrfilenames(info):
    radix = '_'.join((info['sample'], info['dataset']))
    root = info['root']
    sample = info['sample']
    fformat = xiaedf.xiaformat_ctr(radix, info['num'])
    return [os.path.join(root, sample, radix, 'zap',
            fformat.format(extractname(item)))
            for item in info['items']]


def specfilename(info):
    filename = 'spec{:02d}'.format(info['spec'])
    radix = '_'.join((info['sample'], filename))
    filename = os.path.join(info['sample'], radix, '{}.dat'.format(radix))
    return os.path.join(info['root'], filename)


def h5path(info):
    root = os.path.join(info["root"], info["sample"]+'.h5')
    entry = "_".join((info["sample"], info["dataset"]))+'.'+info["scan"]
    return "{}::/{}/{}".format(root, entry, info["process"])


def createscene(parameters, output=None):
    parameters = parameters.copy()
    objects = parameters.pop('objects')

    # Create scene and connect it to axes
    figsize = parameters.get('figsize', None)
    if figsize is None:
        figsize = defaultsettings.default('figure.figsize')
    defaultsettings.adapttofigsize(figsize, **parameters)
    f, ax = plt.subplots(figsize=figsize)
    if parameters.get('noaxes', False):
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        # ax.set_axis_off()
    sc = scene.Scene(unit0=parameters.get('unitfast', 'mm'),
                     unit1=parameters.get('unitslow', 'mm'),
                     title=parameters.get('title', None))
    sc.setaxes(ax)

    # Add objects to scene
    for info in objects:
        plotparams = info.get('plotparams', {}).copy()
        plotparams['scene'] = sc
        dataparams = info.get('dataparams', {}).copy()
        dataparams['instrument'] = parameters.get('instrument', None)
        pmin = plotparams.pop('lo', [])
        pmax = plotparams.pop('hi', [])
        if 'process' in info:
            uri = h5path(info)
            item = scene_view.Nexus(uri, info.get('items', None),
                                    plotparams=plotparams, **dataparams)
        elif 'uri' in info:
            item = scene_view.Nexus(info['uri'], info.get('items', None),
                                    plotparams=plotparams, **dataparams)
        elif 'spec' in info:
            filename = specname(info)
            item = scene_view.XanesSpec(filename, info.get('items', None),
                                        plotparams=plotparams, **dataparams)
        elif 'sample' in info:
            filenames = ctrfilenames(info)
            item = scene_view.ZapRoiMap(filenames, info.get('items', None),
                                        plotparams=plotparams, **dataparams)
        else:
            raise RuntimeError('Unknown object description')
        if pmin and pmax:
            item.selfscale(pmin=pmin, pmax=pmax)
        item.useaxesnames()

    # Modify axes
    scene_modifyaxes(sc, parameters)

    # Global intensity scaling
    scene_globalscaling(sc, parameters)

    # Save interpolated data
    for item in sc:
        try:
            item.interpolatesave()
        except AttributeError:
            pass

    # Plot/save figure
    sc.updateview()
    if parameters.get('tight_layout', False):
        plt.tight_layout()
    if output is not None:
        savefigures(output, parameters)

    if parameters.get('plot', True):
        # TODO: this hangs when another window is opened in the mean time:
        # sc.ax.get_figure().canvas.mpl_connect('resize_event', lambda event: sc.updateview())
        plt.show()

    closefigures()


class Task(nxprocess.Task):
    """Create scene image
    """

    def _parameters_defaults(self):
        super(Task, self)._parameters_defaults()
        self.required_parameters |= {
            'objects',
            'instrument'
        }
        self.optional_parameters |= {
            'figsize',
            'noaxes',
            'unitfast',
            'unitslow',
            'title',
            'instrument',
            'tight_layout',
            'plot',
            'originzero',
            'transpose',
            'flipx',
            'flipy',
            'incx',
            'incy',
            'aspectratio',
            'rlo',
            'rhi',
            'alo',
            'ahi',
            'saveparams'
        }

    def _execute(self):
        createscene(self.parameters, self.temp_localpath.path)

    @property
    def temp_localpath(self):
        path = super(Task, self).temp_localpath
        return path.parent[path.name+'.png']

    @property
    def output_localpath(self):
        path = super(Task, self).output_localpath
        return path.parent[path.name+'.png']
