# -*- coding: utf-8 -*-

from . import nxregulargrid
from . import basetask
from ..utils import instance
from ..io import fs


class Task(nxregulargrid.Task):

    def _parameters_defaults(self):
        super(Task, self)._parameters_defaults()
        self.required_parameters |= {
            'alignmethod',
            'reference',
            'refimageindex',
            'crop',
            'roi',
            'plot',
            'pad',
            'onraw',
            "alignmultiplier",
            "alignoffset"
        }
        parameters = self.parameters
        parameters['refimageindex'] = parameters.get('refimageindex', -1)
        parameters['crop'] = parameters.get('crop', False)
        parameters['roi'] = parameters.get('roi', None)
        parameters['plot'] = parameters.get('plot', False)
        parameters['pad'] = not parameters['crop']
        parameters['onraw'] = True
        parameters['alignmultiplier'] = parameters.get('alignmultiplier', None)
        parameters['alignoffset'] = parameters.get('alignoffset', None)

        alignmethod = parameters.get('alignmethod', None)
        if alignmethod == "sift":
            from ..align.alignSift import alignSift as alignclass
        elif alignmethod == "elastix":
            from ..align.alignElastix import alignElastix as alignclass
        elif alignmethod == "fft":
            from ..align.alignFFT import alignFFT as alignclass
        elif alignmethod == "min":
            from ..align.alignSimple import alignMin as alignclass
        elif alignmethod == "centroid":
            from ..align.alignSimple import alignCentroid as alignclass
        elif alignmethod == "max":
            from ..align.alignSimple import alignMax as alignclass
        else:
            raise basetask.ParameterError(
                'Unknown alignmethod {}'.format(repr(alignmethod)))
        self.alignclass = alignclass

    def _prepare_process(self):
        # Output signal paths
        self.signalsout = []
        for signalin in self.grid.signals:
            nxdata = self.temp_nxresults[signalin.parent.name]
            if not nxdata.exists:
                nxdata = self.temp_nxresults.nxdata(name=signalin.parent.name)
            self.signalsout.append(nxdata[signalin.name])

    def _process_axes(self, o):
        axes1 = self.signal_axes
        axes2 = o.transform_axes([ax.values for ax in axes1])
        out = []
        for ax1, ax2 in zip(axes1, axes2):
            ax2 = self._new_axis(ax2, ax1)
            if ax1 == ax2:
                out.append(ax1)
            else:
                out.append(ax2)
        return out

    def rawcalc(self):
        parameters = self.parameters
        multiplier = parameters['alignmultiplier']
        offset = parameters['alignoffset']
        if multiplier is not None or offset is not None:
            if multiplier is None:
                multiplier = 1
            if offset is None:
                offset = 0
            def rawcalc(img):
                return multiplier*img + offset
            return rawcalc
        return None

    def _execute_grid(self):
        # Align image stacks
        with self.grid.open_signals() as datasets:
            parameters = self.parameters
            kwargs = {k: parameters[k] for k in ['stackdim', 'plot']}
            with self.signalsout[0].h5open() as fout:
                signalsout = [sig.path for sig in self.signalsout]
                o = self.alignclass(datasets, None, fout,
                                    signalsout, "", **kwargs)
                kwargs = {k: parameters[k] for k in [
                    'onraw', 'pad', 'crop', 'roi', 'refimageindex']}
                kwargs['rawcalc'] = self.rawcalc()
                self._prepare_reference(kwargs)
                o.align(self.reference_signal_index, **kwargs)
                axes = self._process_axes(o)
                axes = self._create_axes(axes)
                self.temp_nxresults['change-of-frame'].write(
                    data=o.absolute_cofs())

        # Set NXdata signal and axes attributes
        for signalout in self.signalsout:
            nxdata = signalout.parent
            signalout = nxdata.add_signal(signalout.name)
            if not nxdata.axes:
                nxdata.set_axes(*axes)

    def _prepare_reference(self, parameters):
        refimageindex = parameters['refimageindex']
        refgrid = self.reference_signal
        if instance.isstring(refimageindex):
            if refimageindex == 'first':
                refimageindex = 0
            elif refimageindex == 'last':
                refimageindex = -1
            elif refimageindex == 'middle':
                refimageindex = refgrid.shape[refgrid.stackdim]//2
            else:
                raise basetask.ParameterError(
                    'Unknown refimageindex {}'.format(repr(refimageindex)))
        elif instance.isinteger(refimageindex):
            pass
        elif refimageindex is None:  # pair-wise
            refimageindex = None
        else:  # fraction
            refimageindex = int(
                round(refimageindex*refgrid.shape[refgrid.stackdim]))
        parameters['refimageindex'] = refimageindex
