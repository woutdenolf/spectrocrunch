# -*- coding: utf-8 -*-

import os
import re
import itertools
import logging
import functools
import numpy as np
from glob import glob
import matplotlib.pyplot as plt
from future.utils import with_metaclass
from abc import ABCMeta, abstractmethod, abstractproperty

from ..io import localfs
from ..io import mca
from ..utils import subprocess
from ..utils import instance
from ..materials.element import Element
from ..materials.compoundfromname import compoundfromname

logger = logging.getLogger(__name__)


def installed():
    return subprocess.installed('xrmc')


def execute(*args, **kwargs):
    return subprocess.execute('xrmc', *args, **kwargs)


headerdt = np.dtype([('ninter', np.int32),  # N. of interactions
                     ('ncol', np.int32),  # N. of columns
                     ('nrow', np.int32),  # N. of rows
                     ('scol', np.float64),  # Pixel size Sx
                     ('srow', np.float64),  # Pixel size Sy
                     ('time', np.float64),  # Exposure time in sec.
                     ('type', np.int32),  # Pixel content type
                     ('nbin', np.int32),  # N. of energy bins
                     ('emin', np.float64),  # Minimum bin energy
                     ('emax', np.float64)])  # Maximum bin energy


def saveemptyresult(filename, header):
    with open(filename, "w") as f:
        header.tofile(f)
        shape = header['ninter'], header['nbin'], header['ncol'], header['nrow']
        np.zeros(np.product(shape), dtype=np.float64).tofile(f)


def loadxrmcresult(filename):
    """
    returns: Ninteractions, Nrows, Ncolumns, Nchannels
    """
    with open(filename, "rb") as f:
        header = np.fromfile(f, dtype=headerdt, count=1)[0]
        shape = header['ninter'], header['nbin'], header['ncol'], header['nrow']
        data = np.fromfile(f, dtype=np.float64).reshape(shape)
        x0 = np.arange(header['nrow'])-(header['nrow']-1)/2.
        x1 = np.arange(header['ncol'])-(header['ncol']-1)/2.
        x0 *= header['srow']
        x1 *= header['scol']
        energy = np.linspace(header['emin'], header['emax'], header['nbin'])
        zero = header['emin']
        gain = (header['emax']-header['emin'])/(header['nbin']-1.)
        info = {'time': header['time'],
                'x0': x0,
                'x1': x1,
                'xenergy': energy,
                'zero': zero,
                'gain': gain}
    return np.transpose(data, (0, 3, 2, 1)), info


def loadxrmcresults(path, basename, ext='.dat'):
    """
    returns: Ninteractions, Nstack, Nrows, Ncolumns, Nchannels
    """
    filenames = glob(os.path.join(path, basename+'*'+ext))
    data = []
    info = {}
    for filename in sorted(filenames):
        datai, infoi = loadxrmcresult(filename)
        data.append(datai)
        if not info:
            info = infoi
    data = np.stack(data, axis=1)
    return data, info


def xrmcresult_to_mca(xrmcfile, mcafile, mode='w'):
    """
    Save sum spectrum as mca
    """
    data, info = loadxrmcresult(xrmcfile)
    mcadata = data.sum(axis=tuple(range(data.ndim - 1)))
    mca.save(mcadata, mcafile, mode=mode,
             zero=info['zero'], gain=info['gain'])


def showxrmcresult(data, x0=None, x1=None, xenergy=None, time=None, ylog=False):
    print('Data shape: {}'.format(data.shape))

    if len(x0) > 1:
        x0min, x0max = x0[0], x0[-1]
        dx0 = (x0[1]-x0[0])*0.5
        x0min -= dx0
        x0max += dx0
        x1min, x1max = x1[0], x1[-1]
        dx1 = (x1[1]-x1[0])*0.5
        x1min -= dx1
        x1max += dx1
        extent = x1min, x1max, x0min, x0max
    else:
        extent = None

    axes = tuple(range(2, data.ndim))
    min = data.min(axes).T
    max = data.max(axes).T
    sum = data.sum(axes).T
    for i, (mi, ma, su) in enumerate(zip(min, max, sum)):
        print('Step {}'.format(i))
        print(' Min counts/pixel (for each order): {}'.format(mi))
        print(' Max counts/pixel (for each order): {}'.format(ma))
        print(' Total counts (for each order): {}'.format(su))

    for order, stack in enumerate(data):
        for islice, chunk in enumerate(stack):
            chunk = np.squeeze(chunk)
            if chunk.ndim == 0:
                # Diode signal
                continue
            elif chunk.ndim == 1:
                # XRF spectrum
                spectrum = chunk
                img = None
            elif chunk.ndim == 2:
                # Transmission image
                spectrum = None
                img = chunk
            else:
                # XRF image
                img = chunk.sum(axis=-1)
                spectrum = chunk.sum(axis=(0, 1))
            label = 'Order {}, Step {}'.format(order, islice)
            if spectrum is not None:
                plt.figure(1)
                if xenergy is None:
                    plt.plot(spectrum, label=label)
                    plt.xlabel('MCA (channels)')
                else:
                    plt.plot(xenergy, spectrum, label=label)
                    plt.xlabel('MCA (keV)')
                plt.title('Exposure time = {} sec'.format(time))
                plt.ylabel('Counts')
                if ylog:
                    plt.yscale('log')
                plt.legend()
            if img is not None:
                plt.figure()
                plt.imshow(chunk.sum(axis=-1), origin='lower', extent=extent)
                plt.title(label + ' (looking upstream)')
                plt.xlabel('Synchrotron plane (cm)')
                plt.ylabel('Vertical direction (cm)')


class xrmc_file(object):

    def __init__(self):
        super(xrmc_file, self).__init__()

    @property
    def filename(self):
        if self.fileprefix and self.filesuffix:
            return self.fileprefix + '_' + self.filesuffix + '.dat'
        elif self.filesuffix:
            return self.filesuffix + '.dat'
        else:
            raise RuntimeError('Needs at least a file suffix')

    @property
    def filepath(self):
        return os.path.join(self.path, self.filename)

    @abstractproperty
    def body(self):
        pass

    @abstractproperty
    def header(self):
        pass

    @property
    def content(self):
        body = self.body
        if not body:
            return []
        lines = self.header
        if not lines:
            lines = []
        lines = ['; '+s for s in lines]
        lines.append('')
        lines += body
        lines += ['', 'End']
        return lines

    @property
    def empty(self):
        return not bool(self.body)

    def save(self, mode='w'):
        filepath = self.filepath
        logger.info('Saving ' + filepath)
        localfs.Path(self.path).mkdir()
        localfs.Path(filepath).remove()
        lines = self.content
        if not lines:
            return
        with open(filepath, mode=mode) as f:
            for line in lines:
                f.write(line+'\n')

    @staticmethod
    def title(text):
        return ';%%%%%%%%%%%%%%%%%%%%% {} %%%%%%%%%%%%%%%%%%%%%'.format(text.upper())


class xrmc_main(xrmc_file):

    def __init__(self, path, fileprefix=None):
        self.path = path
        self.fileprefix = fileprefix
        self.devices = []
        self.loops = []  # first is the outer loop
        super(xrmc_main, self).__init__()

    def __repr__(self):
        s = ''
        lst = '\n '.join(map(repr, self.devices))
        if lst:
            s += 'Devices:\n ' + lst
        lst = '\n '.join([repr(loop[0]) for loop in self.loops])
        if lst:
            if s:
                s += '\n'
            s += 'Loops:\n ' + lst
        return s

    @property
    def filesuffix(self):
        return 'main'

    @property
    def header(self):
        return ['X-ray Monte Carlo (xrmc) input file',
                'Usage: xrmc {}'.format(self.filename)]

    @property
    def body(self):
        lines = []
        if self.devices:
            lines += ['', self.title('DEVICES USED FOR SIMULATION')]
            for device in self.devices:
                if not device.empty:
                    lines += ['', 'Load ' + device.filename]
        detectors = self.detectors()
        if detectors:
            lines += ['', self.title('START SIMULATIONS')]
            if self.loops:
                fmt = self.loop_suffix_format
                bfirst = True
                for idx, steps in self.loop_steps():
                    lines.append('')
                    if bfirst:
                        bfirst = False
                    else:
                        lines += ['Load ' + step.filename for step in steps]
                    for detector in detectors:
                        if detector.TYPE == 'detectorconvolute':
                            imgType = 'ConvolutedImage'
                        else:
                            imgType = 'Image'
                        lines += ['Run {}'.format(detector.name),
                                  'Save {} {} {}'.format(detector.name,
                                                         imgType,
                                                         detector.reloutput(suffix=fmt.format(*idx)))]
            else:
                for detector in detectors:
                    if detector.TYPE == 'detectorconvolute':
                        imgType = 'ConvolutedImage'
                    else:
                        imgType = 'Image'
                    lines += ['Run {}'.format(detector.name),
                              'Save {} {} {}'.format(detector.name,
                                                     imgType,
                                                     detector.reloutput())]
        return lines

    def output_suffixes(self):
        if self.loops:
            fmt = self.loop_suffix_format
            return [fmt.format(*idx) for idx, _ in self.loop_steps()]
        else:
            return ['']

    def add_device(self, cls, *args, **kwargs):
        device = cls(self, *args, **kwargs)
        self.devices.append(device)
        return device

    def save(self, **kwargs):
        super(xrmc_main, self).save(**kwargs)
        for device in self.devices:
            device.save(**kwargs)
        for transform, rtransform, _ in self.loops:
            transform.save(**kwargs)
            rtransform.save(**kwargs)

    def source(self, **kwargs):
        return self.finddevice(cls=Source, **kwargs)

    def spectrum(self, **kwargs):
        return self.finddevice(cls=Spectrum, **kwargs)

    def sample(self, **kwargs):
        return self.finddevice(cls=Sample, **kwargs)

    def quadrics(self, **kwargs):
        return self.finddevice(cls=Quadrics, **kwargs)

    def compositions(self, **kwargs):
        return self.finddevice(cls=Compositions, **kwargs)

    def objects(self, **kwargs):
        return self.finddevice(cls=Objects, **kwargs)

    def detectors(self, first=False, **kwargs):
        return self.finddevice(cls=Detector, first=first, **kwargs)

    def finddevice(self, cls=None, name=None, first=True):
        lst = []
        for device in self.devices:
            badd = True
            if cls is not None:
                badd &= isinstance(device, cls)
            if name is not None:
                badd &= device.name == name
            if badd:
                lst.append(device)
        if first:
            if lst:
                return lst[0]
            else:
                return None
        else:
            return lst

    def removedevice(self, cls=None, name=None):
        if cls is None and name is None:
            return
        lst = []
        for device in self.devices:
            bremove = True
            if cls is not None:
                bremove &= isinstance(device, cls)
            if name is not None:
                bremove &= device.name == name
            if not bremove:
                lst.append(device)
        self.devices = lst

    def simulate(self):
        if not installed():
            raise RuntimeError('xrmc is not installed')
        for detector in self.detectors():
            detector.removeoutput()
        localfs.Path(self.path)['output'].mkdir()
        _, err, returncode = execute(self.filepath,
                                     cwd=self.path, stderr=True)
        success = returncode == 0
        if not success:
            if err:
                print('XRMC errors:')
                print(err)
            if 'pulsetrain maximum reached' in err:
                self.create_empty_output()
                success = True
        return success

    def create_empty_output(self):
        ninteractions = self.sample().ninteractions
        for detector in self.detectors():
            detector.saveemptyoutput(ninteractions)

    def addloop(self, name, device, step=None, nsteps=1):
        """
        First call specifies the outer loop and last call the inner loop

        Args:
            name(str):
            device(object):
            step(list(2-tuple)): [('Rotate', (...)), ('Translate', (...)), ('TranslateAll', (...)), ...]
            nsteps(num): number of steps
        """
        transform = Transform(self, name, device, step)
        rtransform = transform.retour(nsteps)
        self.loops.append((transform, rtransform, nsteps))

    def removeloops(self):
        self.loops = []

    @property
    def loop_suffix_format(self):
        fmt = ''
        for _, _, nsteps in self.loops:
            fmt += '_{{:0{}d}}'.format(max(len(str(nsteps)), 1))
        return fmt

    def loop_steps(self):
        loops = self.loops
        loopidx = [list(range(nsteps+1)) for _, _, nsteps in loops]
        idx_prev = (-1,)*len(loops)
        for idx in itertools.product(*loopidx):
            idxsteps = []
            for (step, rstep, _), a, b in zip(loops, idx_prev, idx):
                if b > a:
                    idxsteps.append(step)
                elif b < a:
                    idxsteps.append(rstep)
            yield idx, idxsteps[::-1]
            idx_prev = idx


class xrmc_child(xrmc_file):

    def __init__(self, parent):
        self.parent = parent
        super(xrmc_child, self).__init__()

    @property
    def fileprefix(self):
        return self.parent.fileprefix

    @fileprefix.setter
    def fileprefix(self, value):
        self.parent.fileprefix = value

    @property
    def path(self):
        return self.parent.path

    @path.setter
    def path(self, value):
        self.parent.path = value

    @property
    def body(self):
        lines = []
        code = self.code
        if code:
            fmt = '{{:{}s}}  ; {{}}'.format(max(len(tpl[0]) for tpl in code if tpl))
            for tpl in code:
                if tpl:
                    var, comment, enabled = tpl
                    line = fmt.format(var, comment)
                    if not enabled:
                        line = ';' + line
                else:
                    line = ''
                lines.append(line)
        return lines


class xrmc_device(xrmc_child):

    TYPE = None

    def __init__(self, parent, name):
        super(xrmc_device, self).__init__(parent)
        self.name = name

    def __repr__(self):
        return self.name

    @property
    def code(self):
        return [('Newdevice ' + self.TYPE, 'Device type', True),
                (self.name, 'Device name', True)]

    @property
    def filesuffix(self):
        return self.name

    @filesuffix.setter
    def filesuffix(self, value):
        self.name = value.lower()


class xrmc_reference_frame(object):
    """
    The User sample reference frame:

    * Z-axis: beam direction
    * X-axis: horizontal and pointing to the right when looking upstream
    * Y-axis: vertical and pointing upwards

    The XRMC sample reference frame is related by x', y', z' = -x, z, y.

    Additionally, the User defines coordinate transformation while
    XMRC expects change-of-basis.
    """

    def __init__(self):
        super(xrmc_reference_frame, self).__init__()

    @staticmethod
    def translateUserToXrmc(x, y, z):
        return -x, z, y

    def rotationUserToXrmc(self, x0, y0, z0, ux, uy, uz, angle):
        if not angle:
            return ()
        return self.translateUserToXrmc(x0, y0, z0) + \
               self.translateUserToXrmc(ux, uy, uz) + \
               (-angle,)

    @staticmethod
    def xrmcArrayInput(name, comment, a, n):
        b = bool(a)
        if not b:
            a = [0]*n
        return ' '.join([name]+list(map(str, a))), comment, b

    def xrmcTranslationInput(self, name, *args):
        arr = self.translateUserToXrmc(*args)
        comment = 'position in the sample frame (x, y, z; cm)'
        return self.xrmcArrayInput(name, comment, arr, len(args))

    def xrmcRotationInput(self, name, *args):
        arr = self.rotationUserToXrmc(*args)
        comment = 'rotation around axis passing through (x, y, z; cm), with direction (u,v,w; cm) and rotation angle (deg; change-of-frame)'
        return self.xrmcArrayInput(name, comment, arr, len(args))


def addPositionalAttrs(names):
    """
    Inherited attributes from parent
    """
    def inner(cls):
        for name in names:
            def getAttr(self, attr_name=name):
                if self.lock_to is None:
                    return getattr(self, "_" + attr_name)
                else:
                    return getattr(self.lock_to, "_" + attr_name)
            def setAttr(self, value, attr_name=name):
                if self.lock_to is None:
                    setattr(self, "_" + attr_name, value)
                else:
                    setattr(self.lock_to, "_" + attr_name, value)
            prop = property(getAttr, setAttr)
            setattr(cls, name, prop)
        return cls
    return inner


@addPositionalAttrs(('radius', 'polar', 'azimuth', 'hoffset', 'voffset'))
class xrmc_positional_device(xrmc_device, xrmc_reference_frame):
    """
    Device reference point is positioned at spherical coordinates [r, polar, azimuth].
    Device surface goes through the reference point.
    Device surface-center conincides with reference point (unless offsets are not zero)
    """

    def __init__(self, parent, name, radius=None, polar=None, azimuth=None,
                 hoffset=0, voffset=0, lock_to=None):
        super(xrmc_positional_device, self).__init__(parent, name)
        self.lock_to = lock_to
        self._radius = radius
        self._polar = polar
        self._azimuth = azimuth
        self._hoffset = hoffset
        self._voffset = voffset

    def xrmcSphericalRotationInput(self, name):
        """
        Transformation of point [0, 0, r] (cartesian) to [r, polar, azimuth] (spherical)
        """
        # Spherical: [r, polar, azimuth]
        # Cartesian: RotZ(azimuth).RotY(polar).[0, 0, r]
        return [self.xrmcRotationInput(name, 0, 0, 0, 0, 1, 0, self.polar),
                self.xrmcRotationInput(name, 0, 0, 0, 0, 0, 1, self.azimuth)]

    def xrmcPositionInput(self):
        """
        Position device at [r, polar, azimuth] (spherical)
        """
        T = self.xrmcTranslationInput('X', self.hoffset, self.voffset, self.radius)
        # Rotate reference point
        R = self.xrmcSphericalRotationInput('Rotate')
        return T, R

    def add_rotationloop(self, suffix, x, y, z, u, v, w, step, nsteps):
        """
        Rotate around axis going through [x, y, z] with direction [u, v, w]
        """
        params = x, y, z, u, v, w, step
        if isinstance(self, Quadrics):
            step = [('RotateAll', params)]
        else:
            step = [('Rotate', params)]
        self.parent.addloop(self.name+suffix, self, step=step, nsteps=nsteps)

    def add_Xrotationloop(self, step, nsteps):
        self.add_rotationloop('_Xrot', 0, 0, 0, 1, 0, 0, step, nsteps)

    def add_Yrotationloop(self, step, nsteps):
        self.add_rotationloop('_Yrot', 0, 0, 0, 0, 1, 0, step, nsteps)

    def add_Zrotationloop(self, step, nsteps):
        self.add_rotationloop('_Zrot', 0, 0, 0, 0, 0, 1, step, nsteps)

    def add_polarrotationloop(self, step, nsteps):
        """
        Keeps azimuth constant modulo 180 degrees
        """
        angle = np.radians(self.azimuth+90)
        u = np.cos(angle)
        v = np.sin(angle)
        self.add_rotationloop('_polarrot', 0, 0, 0, u, v, 0, step, nsteps)

    def add_azimuthalrotationloop(self, step, nsteps):
        """
        Keeps polar angle constant
        """
        self.add_rotationloop('_azimuthrot', 0, 0, 0, 0, 0, 1, step, nsteps)


class Transform(xrmc_child, xrmc_reference_frame):

    def __init__(self, parent, name, device, transformations):
        super(Transform, self).__init__(parent)
        self.filesuffix = name
        self.device = device
        self.transformations = transformations
        for cmd, parameters in transformations:
            if cmd.startswith('Rotate'):
                if len(parameters) != 7:
                    raise ValueError('A rotation needs 7 parameters: x,y,z,u,v,w,angle')
            elif cmd.startswith('Translate'):
                if len(parameters) != 3:
                    raise ValueError('A translation needs 3 parameters: dx,dy,dz')
            else:
                raise ValueError('Unknown transformation {}'.format(repr(cmd)))

    def __repr__(self):
        return ','.join(cmd for cmd, _ in self.transformations) + '({})'.format(self.device)

    @property
    def header(self):
        return ['Transformation of device {}.'.format(repr(self.device.name)),
                'Executed every time this file is loaded.']

    @property
    def code(self):
        lines = [('Device {}'.format(self.device.name), 'device name', True)]
        for cmd, parameters in self.transformations:
            if cmd.startswith('Rotate'):
                lines.append(self.xrmcRotationInput(cmd, *parameters))
            elif cmd.startswith('Translate'):
                lines.append(self.xrmcTranslationInput(cmd, *parameters))
            else:
                raise NotImplementedError(cmd)
        return lines

    def retour(self, nsteps):
        rtransformations = []
        for cmd, parameters in self.transformations:
            if cmd.startswith('Rotate'):
                parameters = tuple(parameters[:-1]) + (-nsteps*parameters[-1],)
            elif cmd.startswith('Translate'):
                parameters = (-nsteps*parameters[0],
                              -nsteps*parameters[1],
                              -nsteps*parameters[2])
            else:
                raise NotImplementedError(cmd)
            rtransformations.append((cmd, parameters))
        return self.__class__(self.parent, self.filesuffix+'_retour', self.device, rtransformations)


class WithObjects(object):

    def __init__(self, parent, name, **kwargs):
        super(WithObjects, self).__init__(parent, name, **kwargs)
        self.quadrics = parent.add_device(Quadrics,
                                          name + '_quadrics',
                                          lock_to=self)
        self.objects = parent.add_device(Objects,
                                         name + '_objects')

    def add_layer(self, material=None, thickness=None,
                  dhor=None, dvert=None, ohor=0, overt=0):
        if not material:
            material = 'Vacuum'
        ilayer = self.quadrics.add_layer(thickness, dhor, dvert, ohor, overt)
        box = ['layer{}a'.format(ilayer), 'layer{}b'.format(ilayer),
               'blayer{}_xmin'.format(ilayer), 'blayer{}_xmax'.format(ilayer),
               'blayer{}_zmin'.format(ilayer), 'blayer{}_zmax'.format(ilayer)]
        name = self.parent.compositions().addmaterial(material)
        self.objects.add('layer{}'.format(ilayer), material, box)
        return name

    def clear(self):
        self.quadrics.clear()
        self.objects.clear()


class Source(WithObjects, xrmc_positional_device):

    TYPE = 'source'

    def __init__(self, parent, name, distance=None, beamsize=None):
        super(Source, self).__init__(parent, name,
                                     radius=distance,
                                     polar=180,
                                     azimuth=0,
                                     hoffset=0,
                                     voffset=0)
        self.beamsize = beamsize

    @property
    def distance(self):
        return self.radius

    @distance.setter
    def distance(self, value):
        self.radius = value

    @property
    def position(self):
        return -self.distance

    @property
    def header(self):
        return ['Source position/orientation file',
                '',
                'The sample reference frame is defined so that the beam travels along',
                'the Y-axis, the Z-axis points upwards and the X-axis is parallel with',
                'the synchrotron plane, pointing to the left when looking upstream.',
                '',
                'In the local reference from of the source, the beam',
                'travels along the Z-axis.'
                '',
                'The source is positioned at y=-distance and can have a divergence'
                '(horizontal, vertical) and size (horizontal, vertical, longitudinal)',
                'in the local reference frame.']

    @property
    def divergence(self):
        return np.arctan2(self.beamsize*0.5, self.distance)

    @property
    def code(self):
        lines = super(Source, self).code
        div = self.divergence
        position = self.xrmcTranslationInput('X', 0, 0, -self.distance)
        lines += [('SpectrumName {}'.format(self.parent.spectrum().name), 'Spectrum input device name', True),
                  position,
                  ('ui 1 0 0', 'direction of X in the sample frame (x, y, z; cm)', True),
                  ('uk 0 1 0', 'direction of Z in the sample frame (x, y, z; cm)', True),
                  ('Divergence {} {}'.format(div, div), 'beam divergency (divH, divV; rad)', True),
                  ('Size 0.0 0.0 0.0', 'source size (sH, sV, sL; cm)', True)]
        return lines


class Spectrum(xrmc_device):

    TYPE = 'spectrum'

    def __init__(self, parent, name, lines=None, polarization=True, multiplicity=1):
        super(Spectrum, self).__init__(parent, name)
        self.polarization = polarization
        self.lines = lines
        self.multiplicity = multiplicity

    @property
    def header(self):
        return ['Spectrum file',
                '',
                'Discrete X-ray source spectrum + continuum',
                '',
                'Polarization is either unpolarized or linear polarized.',
                'The orientation of the electric field vector is defined'
                'by (Ehorizontal, Evertical) in the source reference frame (see source).']

    @property
    def code(self):
        lines = super(Spectrum, self).code
        m = int(np.round(self.multiplicity))
        lines += [('PolarizedFlag {:d}'.format(self.polarization), 'unpolarized/polarized beam (0/1)', True),
                  ('LoopFlag 1', '0: extract random energies on the whole spectrum', True),
                  ('', '1: loop on all lines and sampling points', True),
                  ('ContinuousPhotonNum {:d}'.format(m), 'Multiplicity of events for each spectrum interval', True),
                  ('LinePhotonNum {:d}'.format(m), 'Multiplicity of events for each spectrum line', True),
                  ('RandomEneFlag 1', 'enable random energy on each interval (0/1)', True),
                  ('Lines', 'discrete energy lines of the spectrum', True),
                  (str(len(self.lines)), 'Number of lines in the spectrum', True)]
        for energy, sigma, mult in self.lines:
            lines.append(('{} {} {} {}'.format(energy, sigma, mult, 0), 'Energy (keV), sigma (keV), Eh (ph/s), Ev (ph/s)', True))
        return lines

    @property
    def emax(self):
        return max(energy for energy, sigma, intensity in self.lines)

    @property
    def totalflux(self):
        return sum(intensity for energy, sigma, intensity in self.lines)


class Detector(WithObjects, xrmc_positional_device):

    TYPE = 'detectorarray'

    def __init__(self, parent, name, distance=None,
                 poissonnoise=True, ebinsize=None,
                 forcedetect=False, multiplicity=None,
                 emin=None, emax=None,
                 time=1, pixelshape='elliptical', **kwargs):
        super(Detector, self).__init__(parent, name,
                                       radius=distance,
                                       **kwargs)
        self.emin = emin
        self.emax = emax
        self.ebinsize = ebinsize
        self.poissonnoise = poissonnoise
        self.forcedetect = forcedetect
        self.multiplicity = multiplicity
        self.time = time
        self.pixelshape = pixelshape.lower()

    @property
    def distance(self):
        return self.radius

    @distance.setter
    def distance(self, value):
        self.radius = value

    @property
    def position(self):
        return self.distance

    @property
    def header(self):
        return ['Detector array parameter file',
                '',
                'The detector reference frame (before transformation) is defined',
                'so that the Z-axis points upstream, the X-axis (detector image rows)',
                'points upwards and the Y-axis (detector image columns) points to the',
                'right when looking upstream.',
                '',
                'A detector can have 1 or more elliptical or square pixels.']

    @property
    def code(self):
        lines = super(Detector, self).code
        energydispersive = self.energydispersive
        m = int(np.round(self.multiplicity))
        elliptical = self.pixelshape == 'elliptical'
        if energydispersive:
            pixeltype = 2
        else:
            pixeltype = 0
        
        position, rotations = self.xrmcPositionInput()
        lines += [('SourceName {}'.format(self.parent.sample().name), 'Source input device name', True),
                  ('NPixels {} {}'.format(*self.dims[::-1]), 'Pixel number (nx x ny, ncol x nrow)', True),
                  ('PixelSize {} {}'.format(*self.pixelsize[::-1]), 'Pixel Size (sx x sy, scol x srow; cm)', True),
                  ('Shape {:d}'.format(elliptical), 'Pixel shape (0 rectangular, 1 elliptical)', True),
                  ('dOmegaLim {}'.format(2*np.pi), 'Limit solid angle (default = 2*PI)', False),
                  position,
                  ('ui 0 0 1', 'direction of X (rows) in the sample frame (x, y, z; cm)', True),
                  ('uk 0 -1 0', 'direction of Z in the sample frame (x, y, z; cm)', True),
                  ('ExpTime {}'.format(self.time), 'Exposure time (sec)', True),
                  ('PhotonNum {:d}'.format(m), 'Multiplicity of simulated events per pixel', True),
                  ('RandomPixelFlag 1', 'Enable random point on pixels (0/1)', True),
                  ('PoissonFlag {:d}'.format(self.poissonnoise), 'Enable Poisson statistic on pix. counts (0/1)', True),
                  ('RoundFlag 0', 'Round pixel counts to integer (0/1)', True),
                  ('AsciiFlag 0', 'binary(0) or ascii(1) file format', False),
                  ('ForceDetectFlag {:d}'.format(self.forcedetect), 'MC variance reduction (0/1)', True),
                  ('HeaderFlag 1', 'Use header in output file (0/1)', True),
                  ('PixelType {:d}'.format(pixeltype), 'Pixel content type:', True),
                  ('', '0: fluence,      1: energy fluence,', True),
                  ('', '2: fluence(E),   3: energy fluence(E)', True),
                  ('Emin {}'.format(self.emin), 'Emin', energydispersive),
                  ('Emax {}'.format(self.emax), 'Emax', energydispersive),
                  ('NBins {}'.format(self.nbins), 'NBins', energydispersive),
                  ('SaturateEmin 0', 'Saturate energies lower than Emin (0/1)', energydispersive),
                  ('SaturateEmax 0', 'Saturate energies greater than Emax (0/1)', energydispersive)]
        lines += rotations
        return lines

    def reloutput(self, suffix=None):
        if not suffix:
            suffix = ''
        return 'output/{}{}.dat'.format(self.name, suffix)

    def absoutput(self, suffix=None):
        if not suffix:
            suffix = ''
        return os.path.join(self.outpath, '{}{}.dat'.format(self.name, suffix))

    @property
    def outpath(self):
        return os.path.join(self.path, 'output')

    @property
    def energydispersive(self):
        return bool(self.nbins)

    @property
    def emax(self):
        if self._emax is None:
            return self.parent.spectrum().emax+1
        else:
            return self._emax

    @emax.setter
    def emax(self, value):
        self._emax = value

    @property
    def emin(self):
        if self._emin is None:
            return 0
        else:
            return self._emin

    @emin.setter
    def emin(self, value):
        self._emin = value

    @property
    def ebinsize(self):
        if self.energydispersive:
            return (self.emax-self.emin)/(self.nbins-1.)
        else:
            return 0

    @ebinsize.setter
    def ebinsize(self, value):
        if bool(value):
            self.nbins = max(int(np.round((self.emax-self.emin)/float(value)+1)), 1)
        else:
            self.nbins = 1

    @property
    def mcagain(self):
        return self.ebinsize

    @property
    def mcazero(self):
        return 0

    def result(self):
        if self.parent.loops:
            suffix = '_*'
        else:
            suffix = None
        filename = self.absoutput(suffix=suffix)
        dirname = os.path.dirname(filename)
        basename, ext = os.path.splitext(os.path.basename(filename))
        return loadxrmcresults(dirname, basename, ext)

    def removeoutput(self):
        filename = self.absoutput()
        if os.path.isfile(filename):
            os.remove(filename)

    def saveemptyoutput(self, ninteractions):
        energydispersive = self.energydispersive
        if energydispersive:
            pixeltype = 2
        else:
            pixeltype = 0
        nrow, ncol = self.dims
        srow, scol = self.pixelsize
        header = np.array((ninteractions, ncol, nrow, scol, srow, self.time,
                          pixeltype, self.nbins, self.emin, self.emax),
                          dtype=headerdt)
        suffixes = self.parent.output_suffixes()
        for suffix in suffixes:
            saveemptyresult(self.absoutput(suffix=suffix), header)


class AreaDetector(Detector):

    TYPE = 'detectorarray'

    def __init__(self, parent, name, pixelsize=None, dims=None,
                 hpoffset=None, vpoffset=None, **kwargs):
        self.dims = dims
        self.pixelsize = pixelsize
        super(AreaDetector, self).__init__(parent, name, **kwargs)
        if hpoffset is not None:
            self.hpoffset = hpoffset
        if vpoffset is not None:
            self.vpoffset = vpoffset

    @property
    def hpoffset(self):
        return self.hoffset/float(self.pixelsize[0])

    @hpoffset.setter
    def hpoffset(self, value):
        self.hoffset = value*self.pixelsize[0]

    @property
    def vpoffset(self):
        return self.voffset/float(self.pixelsize[1])

    @vpoffset.setter
    def vpoffset(self, value):
        self.voffset = value*self.pixelsize[1]


class SingleElementDetector(Detector):

    TYPE = 'detectorarray'

    def __init__(self, parent, name, activearea=None, **kwargs):
        self.dims = 1, 1
        self.activearea = activearea
        super(SingleElementDetector, self).__init__(parent, name, **kwargs)

    @property
    def pixelsize(self):
        a = self.activearea**0.5
        return a, a


class SDD(SingleElementDetector):

    TYPE = 'detectorconvolute'

    def __init__(self, parent, name, material=None, thickness=None,
                 windowmaterial=None, windowthickness=0, dtfrac=None,
                 countrate=None, pulseproctime=0, noise=None, fano=None, **kwargs):
        self.material = material
        self.thickness = thickness
        if not windowmaterial:
            windowmaterial = "Vacuum"
        self.windowmaterial = windowmaterial
        self.windowthickness = windowthickness
        self.noise = noise
        self.fano = fano
        self.set_pulseproctime(dtfrac, countrate, pulseproctime=pulseproctime)
        super(SDD, self).__init__(parent, name, **kwargs)

    @property
    def code(self):
        lines = super(SDD, self).code
        lines += [('CompositionName {}'.format(self.parent.compositions().name), 'Composition input device name', True),
                  ('CrystalPhase {}'.format(self.material), 'Detector crystal material', True),
                  ('CrystalThickness {}'.format(self.thickness), 'Detector crystal thickness', True),
                  ('WindowPhase {}'.format(self.windowmaterial), 'Window material', True),
                  ('WindowThickness {}'.format(self.windowthickness), 'Window thickness', True),
                  ('Noise {}'.format(self.noise), 'Detector energy noise (keV)', True),
                  ('FanoFactor {}'.format(self.fano), 'Detector fano noise (dimensionless)', True),
                  ('PulseWidth {}'.format(self.pulseproctime), 'Time to process one pulse (sec)', bool(self.pulseproctime))]
        return lines

    def set_pulseproctime(self, dtfrac, countrate, pulseproctime=0):
        # Assume extendable deadtime
        if dtfrac and countrate:
            self.pulseproctime = -np.log(1-dtfrac)/countrate
        else:
            self.pulseproctime = pulseproctime


class Sample(WithObjects, xrmc_positional_device):

    TYPE = 'sample'

    def __init__(self, parent, name, multiplicity=(1, 1), **kwargs):
        super(Sample, self).__init__(parent, name, **kwargs)
        self.multiplicity = multiplicity

    @property
    def ninteractions(self):
        return len(self.multiplicity)

    @property
    def header(self):
        return ['Sample parameters file',
                '',
                'maximum scattering order',
                '    0: transmission',
                '    1: first order scattering or fluorescence emission',
                '    2: second order scattering or fluorescence emission']

    @property
    def code(self):
        lines = super(Sample, self).code
        lines += [('SourceName {}'.format(self.parent.source().name), 'Source input device name', True),
                  ('Geom3DName {}'.format(self.parent.objects().name), 'Geom3d input device name', True),
                  ('CompName {}'.format(self.parent.compositions().name), 'Composition input device name', True),
                  ('WeightedStepLength 1', 'Weighted step length (0/1)', True),
                  ('FluorFlag {:d}'.format(len(self.multiplicity)>1), 'Activate Fluorescence (0/1)', True),
                  ('ScattOrderNum {:d}'.format(len(self.multiplicity)-1), 'Maximum scattering order', True)]
        for i, mi in enumerate(self.multiplicity):
            m = str(int(np.round(mi)))
            lines.append((m, 'Multiplicity of simulated events for order {}'.format(i), True))
        return lines


class Objects(xrmc_device):

    TYPE = 'geom3d'

    def __init__(self, parent, name):
        super(Objects, self).__init__(parent, name)
        self.clear()

    @property
    def header(self):
        return ['3D object geometric description file',
                '',
                'Makes objects with quadrics and compositions.']

    @property
    def code(self):
        if not self._dict:
            return []
        lines = super(Objects, self).code
        lines += [('QArrName {}'.format(self.parent.quadrics().name), 'Quadric array input device name', True), 
                  ('CompName {}'.format(self.parent.compositions().name), 'Composition input device name', True)]
        for name, (compoundin, compoundout, quadricnames) in self._dict.items():
            compoundin = materialname(compoundin)
            compoundout = materialname(compoundout)
            lines += [None,
                      ('Object {}'.format(name), '', True),
                      ('{} {}'.format(compoundin, compoundout), 'Phase in, phase out', True),
                      (str(len(quadricnames)), 'Number of quadrics', True),
                      (' '.join(quadricnames), 'Quadric names', True)]
        return lines

    def add(self, name, compoundin, quadricnames):
        try:
            compoundout = self.parent.compositions().atmosphere
        except AttributeError:
            compoundout = None
        self._dict[name] = [compoundin, compoundout, quadricnames]

    def change_atmosphere(self, compoundout):
        for lst in self._dict.values():
            lst[1] = compoundout

    def clear(self):
        self._dict = {}


class Quadrics(xrmc_positional_device):

    TYPE = 'quadricarray'

    def __init__(self, parent, name, **kwargs):
        """
        Args:
            parent
            name
            position: position of the first layer
        """
        super(Quadrics, self).__init__(parent, name, **kwargs)
        self.clear()
 
    @property
    def position(self):
        try:
            return [0, self.lock_to.position, 0]
        except AttributeError:
            return [0, 0, 0]

    @property
    def header(self):
        return ['Quadric array file',
                '',
                'A quadric (a.k.a. quadratic surface) is a parametric surface (plane, cylinder, sphere, ...):',
                '  x.A.x^T = 0   (where x homogeneous coordinates)',
                'For internal space x.A.x^T < 0 and for external space x.A.x^T > 0']

    @property
    def code(self):
        if not self._lst:
            return []
        lines = super(Quadrics, self).code
        for quadric in self._lst:
            name, typ, params = quadric
            params = list(params)
            lines.append(('{} {}'.format(typ, name), '', True))
            if typ == 'Plane':
                for i in range(3):
                    params[i] += self.position[i]
                comment = 'Plane through (x0,y0,z0) with normal vector (nx,ny,nz)'
            elif typ.startswith('Cylinder'):
                axis = typ[-1]
                iaxis = ['X', 'Y', 'Z'].index(axis)
                idx = [j for j in range(3) if j != iaxis]
                for i, j in enumerate(idx):
                    params[i] += self.position[j]
                comment = 'Cylinder with main axis along {} and with coordinates (xa, xb) in the perpendicular plane. (Ra, Rb) are the elliptical semi-axes.'.format(axis)
            else:
                comment = ''
            lines.append((' '.join(list(map(str, params))), comment, True))
        lines += self.xrmcSphericalRotationInput('RotateAll')
        return lines

    def add_plane(self, name, x0, y0, z0, nx, ny, nz, fixed=False):
        if name in self.names:
            raise ValueError('Plane {} already exists'.format(repr(name)))
        params = [x0, y0, z0, nx, ny, nz]
        if fixed:
            params.append('BlockTransformAll')
        self._lst.append((name, 'Plane', params))

    def add_cylinder(self, name, axis, xa, xb, Ra, Rb, fixed=False):
        if name in self.names:
            raise ValueError('Cylinder {} already exists'.format(repr(name)))
        axis = axis.upper()
        if axis not in ['X', 'Y', 'Z']:
            raise ValueError("Axis must be 'X', 'Y' or 'Z'")
        params = [xa, xb, Ra, Rb]
        if fixed:
            params.append('BlockTransformAll')
        self._lst.append((name, 'Cylinder'+axis, params))

    def add_layer(self, thickness, dhor, dvert, ohor, overt, dthickness=1e-10,
                  **kwargs):
        """
        Norm parallel with the beam direction (Y-axis)
        Beam runs though the geometric center

        Args:
            thickness: can be negative
            dhor: length along the horizontal axis
            dvert: length along the vertical axis
            ohor: shift along horizontal axis
            overt: shift along vertical axis
            dthickness: layers cannot be touching
        """
        surface = 0
        ilayer = -1
        if thickness*dthickness < 0:
            dthickness = -dthickness
        for name, clsname, params in self._lst:
            if clsname != 'Plane':
                continue
            m = re.match(r'layer(\d+)b', name)
            if m:
                i = int(m.groups()[0])
                if i > ilayer:
                    ilayer = i
                    surface = params[1]  # end of previous layer
        ilayer += 1
        ta = surface + dthickness
        tb = ta + thickness
        self.add_plane('layer{}a'.format(ilayer), 0, ta, 0, 0, -1, 0, **kwargs)
        self.add_plane('layer{}b'.format(ilayer), 0, tb, 0, 0, 1, 0, **kwargs)
        self.add_plane('blayer{}_xmin'.format(ilayer), ohor-dhor/2, 0, 0, -1, 0, 0, **kwargs)
        self.add_plane('blayer{}_xmax'.format(ilayer), ohor+dhor/2, 0, 0, 1, 0, 0, **kwargs)
        self.add_plane('blayer{}_zmin'.format(ilayer), 0, 0, overt-dvert/2, 0, 0, -1, **kwargs)
        self.add_plane('blayer{}_zmax'.format(ilayer), 0, 0, overt+dvert/2, 0, 0, 1, **kwargs)
        return ilayer

    @property
    def names(self):
        return [tpl[0] for tpl in self._lst]

    def clear(self):
        self._lst = []


def materialname(material):
    try:
        name = material.name
    except AttributeError:
        name = material
    try:
        if 'vacuum' in name.lower():
            name = 'Vacuum'
    except AttributeError:
        name = 'Vacuum'
    name = name.replace(' ', '_')
    return name


class Compositions(xrmc_device):

    TYPE = 'composition'

    def __init__(self, parent, name, atmosphere=None):
        super(Compositions, self).__init__(parent, name)
        self.clear()
        self.atmosphere = atmosphere

    @property
    def header(self):
        return ['Composition file',
                '',
                'List of compounds, defined by their density and elemental composition (weight fractions)']

    @property
    def code(self):
        if not self._dict:
            return []
        lines = super(Compositions, self).code
        for name, (elements, massfractions, density) in self._dict.items():
            lines += [None,
                      ('Phase {}'.format(name), 'Name', True),
                      ('NElem {}'.format(len(elements)), ' # elements', True)]
            for formula, wfrac in zip(elements, massfractions):
                lines.append(('{} {}'.format(formula, wfrac*100), 'mass fraction', True))
            lines += [('Rho {}'.format(density), 'Mass density (g/cm3)', True)]
        return lines

    def add(self, name, elements, massfractions, density):
        self._dict[name] = list(map(str, elements)), massfractions, density

    def clear(self):
        self._dict = {}

    def addmaterial(self, material, name=None):
        """
        Args:
            material(str or Element/Compound/Mixture)
            name(str)
        Returns:
            str: xrmc composition name
        """
        if material is None:
            material = 'vacuum'
        if name is None:
            name = materialname(material)
        if instance.isstring(material):
            if 'vacuum' in material.lower():
                return 'Vacuum'  # No need to add this
            material = compoundfromname(material)
        wfrac = material.elemental_massfractions()
        self.add(name, list(wfrac.keys()), list(wfrac.values()), material.density)
        return name

    @property
    def atmosphere(self):
        return self._atmosphere

    @atmosphere.setter
    def atmosphere(self, material):
        self._atmosphere = self.addmaterial(material)
        for objects in self.parent.objects(first=False):
            objects.change_atmosphere(material)


class XrmcWorldBuilder(object):

    def __init__(self, path, atmosphere=None):
        self.main = main = xrmc_main(path)
        self.sample = main.add_device(Sample, 'sample')
        main.add_device(Compositions, 'compoundlib', atmosphere=atmosphere)

    def __repr__(self):
        return repr(self.main)

    def define_source(self, flux=None, energy=None, distance=None, beamsize=None, multiplicity=1):
        self.main.removedevice(cls=Source)
        self.source = self.main.add_device(Source, 'synchrotron',
                                           beamsize=beamsize, distance=distance)
        self.spectrum = self.main.add_device(Spectrum, 'dcmspectrum',
                                             lines=[[energy, 0, flux]],
                                             multiplicity=multiplicity)

    def add_diode(self, distance=None, activearea=None, ebinsize=None,
                  polar=0, azimuth=0, poissonnoise=False, forcedetect=False,
                  multiplicity=1, time=1):
        self.main.removedevice(cls=Detector)
        self.detector = self.main.add_device(SingleElementDetector, 'detector',
                                    distance=distance, activearea=activearea,
                                    polar=polar, azimuth=azimuth,
                                    ebinsize=ebinsize, poissonnoise=poissonnoise,
                                    forcedetect=forcedetect, multiplicity=multiplicity, time=time)

    def add_xrfdetector(self, distance=None, activearea=None, ebinsize=None,
                        polar=0, azimuth=0, hoffset=0, voffset=0,
                        poissonnoise=False, forcedetect=True,
                        multiplicity=1, time=1, response=None,
                        emin=None, emax=None):
        self.main.removedevice(cls=Detector)
        if response:
            cls = SDD
        else:
            cls = SingleElementDetector
            response = {}
        self.detector = self.main.add_device(cls, 'detector',
                                    distance=distance, activearea=activearea,
                                    polar=polar, azimuth=azimuth,
                                    hoffset=hoffset, voffset=voffset,
                                    ebinsize=ebinsize, poissonnoise=poissonnoise,
                                    forcedetect=forcedetect, multiplicity=multiplicity,
                                    emin=emin, emax=emax,
                                    time=time, **response)

    def add_areadetector(self, distance=None, activearea=None, ebinsize=None,
                         polar=0, azimuth=0, hpoffset=0, vpoffset=0,
                         poissonnoise=False, forcedetect=True,
                         multiplicity=1, time=1, pixelsize=None, dims=None):
        self.main.removedevice(cls=Detector)
        self.detector = self.main.add_device(AreaDetector, 'detector',
                                    distance=distance,
                                    pixelsize=pixelsize, dims=dims,
                                    polar=polar, azimuth=azimuth,
                                    hpoffset=hpoffset, vpoffset=vpoffset,
                                    ebinsize=ebinsize, poissonnoise=poissonnoise,
                                    forcedetect=forcedetect, multiplicity=multiplicity, time=time)

    def addmaterial(self, material, name=None):
        return self.main.compositions().addmaterial(material, name=name)

    def finalize(self, interactions=None):
        sample = self.sample
        sample.multiplicity = interactions
        self.main.save()

    def simulate(self, interactions=None):
        return self.main.simulate()
