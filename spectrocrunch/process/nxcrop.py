# -*- coding: utf-8 -*-

import numpy as np
import logging

from . import nxregulargrid
from . import basetask
from . import axis

logger = logging.getLogger(__name__)


class Task(nxregulargrid.Task):
    def _parameters_defaults(self):
        super(Task, self)._parameters_defaults()
        self.required_parameters |= {"reference"}
        self.optional_parameters |= {"roi", "nanval", "nanfull"}
        parameters = self.parameters
        if all(p not in parameters for p in ["roi", "nanval"]):
            raise basetask.MissingParameter('Specify either "nanval" or "roi"')
        if "nanval" in parameters:
            parameters["nanfull"] = parameters.get("nanfull", True)

    def _prepare_process(self):
        super(Task, self)._prepare_process()
        logger.info("Determine crop size ...")
        if "nanval" in self.parameters:
            self.roi = self.calccroproi(self.reference_signal)
        elif "roi" in self.parameters:
            self.roi = self.convertuserroi(self.reference_signal)
        self.cropped_shape = tuple([b - a for a, b in self.roi])
        self.indexin = [slice(a, b) for a, b in self.roi]
        self.indexout = [slice(None)] * len(self.roi)

    @property
    def signalout_shape(self):
        return self.cropped_shape

    def _process_axes(self):
        axes = []
        for ax, (a, b) in zip(self.signal_axes, self.roi):
            if ax.size != b - a:
                ax = self._new_axis(ax[a:b], ax)
            axes.append(ax)
        return axes

    def calccroproi(self, refgrid):
        """Determine crop ROI to remove slices other than stackdim which contain
        only nanval.

        Args:
            refgrid(RegularGrid):

        Returns:
            list(2-tuple): dimensions of refgrid
        """

        nanval = self.parameters["nanval"]
        nanfull = self.parameters["nanfull"]

        # Mask (True = valid pixel)
        if self.sliced:
            shape, indexgen = refgrid.sliceinfo(self.signal_stackdim)
            mask = np.ones(shape, dtype=np.bool)
            for i in range(refgrid.shape[self.signal_stackdim]):
                img = refgrid[indexgen(i)]
                if nanval is np.nan:
                    mask &= np.logical_not(np.isnan(img))
                else:
                    mask &= img != nanval
        else:
            if nanval is np.nan:
                mask = np.isnan(refgrid.values).sum(axis=self.signal_stackdim) == 0
            else:
                mask = (refgrid.values == nanval).sum(axis=self.signal_stackdim) == 0
            shape = mask.shape

        imask = -1
        roi = []
        for igrid in range(refgrid.ndim):
            if igrid == self.signal_stackdim:
                iroi = (0, refgrid.shape[self.signal_stackdim])
            else:
                imask += 1
                sumdims = tuple([i for i in range(refgrid.ndim - 1) if i != imask])
                indvalid = mask.sum(axis=sumdims)
                if nanfull:
                    m = np.max(indvalid)
                    if m:
                        indvalid = indvalid == m
                indvalid = np.argwhere(indvalid)[:, 0]
                if indvalid.size == 0:
                    return None
                iroi = indvalid[0], indvalid[-1] + 1
            roi.append(iroi)

        return roi

    def convertuserroi(self, refgrid):
        """Determine crop ROI to remove slices other than stackdim which contain
        only nanval.

        Args:
            refgrid(RegularGrid):

        Returns:
            list(2-tuple): dimensions of refgrid
        """
        roi = list(self.parameters["roi"])
        stackdim = self.parameters["stackdim"]
        roi.insert(stackdim, (0, refgrid.shape[stackdim]))
        return roi
