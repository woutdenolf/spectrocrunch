import os
from abc import ABCMeta, abstractmethod, abstractproperty
from future.utils import with_metaclass
import numpy as np
from ..utils import instance
import matplotlib.pyplot as plt


class xrmc_file(object):

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
        lines = self.header
        lines = ['; '+s for s in lines]
        lines.append('')
        lines += self.body
        lines += ['', 'End']
        return lines

    def save(self, mode='w'):
        lines = self.content
        with open(self.filepath, mode=mode) as f:
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
                lines += ['', 'Load ' + device.filename]
        detectors = self.detectors(first=False)
        if detectors:
            lines += ['', self.title('START SIMULATIONS')]
            for detector in detectors:
                lines += ['',
                          'Run {}'.format(detector.name),
                          'Save {} Image {}'.format(detector.name, detector.reloutput)]
        return lines

    def add_device(self, typ, *args, **kwargs):
        cls = DEVICETOCLASS[typ]
        device = cls(self, *args, **kwargs)
        self.devices.append(device)
        return device

    def save(self, **kwargs):
        super(xrmc_main,self).save(**kwargs)
        for device in self.devices:
            device.save(**kwargs)

    def source(self, **kwargs):
        return self._finddevice(cls=Source, **kwargs)
    
    def spectrum(self, **kwargs):
        return self._finddevice(cls=Spectrum, **kwargs)
        
    def sample(self, **kwargs):
        return self._finddevice(cls=Sample, **kwargs)

    def quadrics(self, **kwargs):
        return self._finddevice(cls=Quadrics, **kwargs)
    
    def compositions(self, **kwargs):
        return self._finddevice(cls=Compositions, **kwargs)
    
    def objects(self, **kwargs):
        return self._finddevice(cls=Objects, **kwargs)

    def detectors(self, **kwargs):
        return self._finddevice(cls=Detector, **kwargs)

    def _finddevice(self, cls=None, name=None, first=True):
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

    def simulate(self):
        print('TODO')


class xrmc_child(xrmc_file):

    def __init__(self, parent):
        self.parent = parent

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
    

class xrmc_device(xrmc_child):

    TYPE = None

    def __init__(self, parent, name):
        super(xrmc_device, self).__init__(parent)
        self.name = name

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


class Source(xrmc_device):

    TYPE = 'source'

    def __init__(self, parent, name, distance=None, beamsize=None):
        super(Source, self).__init__(parent, name)
        self.distance = distance
        self.beamsize = beamsize

    @property
    def header(self):
        return ['Source position/orientation file'
                ''
                'The source reference frame (a.k.a. local coordinate system)'
                'is defined by unit vectors {is, js, ks} and origin rs.'
                'The spherical coordinates in the source frame are {R, phi, theta}'
                '(theta is the polar angle, i.e. the angle with ks).'
                ''
                'The central axis of the source runs along ks.'
                'The source divergence is defined as an elliptical cone (thetax, thetay).'
                'The source can be a point or a 3D Gaussian (sigmax sigmay, sigmaz).']

    @property
    def code(self):
        lines = super(Source, self).code
        div = np.arctan2(self.beamsize, self.distance)
        lines += [('SpectrumName {}'.format(self.parent.spectrum().name), 'Spectrum input device name', True),
                  ('X 0 0 -{}'.format(self.distance), 'position of the source in the world frame (x, y, z; cm)', True),
                  ('uk 0 0 1', 'direction of ks in the world frame (x, y, z; cm)', True),
                  ('ui 1 0 0', 'direction of is in the world frame (x, y, z; cm)', True),
                  ('Divergence {} {}'.format(div, div), 'beam divergency (thetax, thetay)', True),
                  ('Size 0.0 0.0 0.0', 'source size (sigmax sigmay, sigmaz; cm)', True)]
        return lines


class Spectrum(xrmc_device):

    TYPE = 'spectrum'

    def __init__(self, parent, name, lines=None, polarization=True):
        super(Spectrum, self).__init__(parent, name)
        self.polarization = polarization
        self.lines = lines

    @property
    def header(self):
        return ['Spectrum file',
                '',
                'Discrete X-ray source spectrum + continuum']

    @property
    def code(self):
        lines = super(Spectrum, self).code
        lines += [('PolarizedFlag {:d}'.format(self.polarization), 'unpolarized/polarized beam (0/1)', True),
                  ('LoopFlag 1', '0: extract random energies on the whole spectrum', True),
                  ('', '1: loop on all lines and sampling points', True),
                  ('ContinuousPhotonNum 1', 'Multiplicity of events for each spectrum interval', True),
                  ('LinePhotonNum 1' ,'Multiplicity of events for each spectrum line', True),
                  ('RandomEneFlag 1', 'enable random energy on each interval (0/1)', True),
                  ('Lines', 'discrete energy lines of the spectrum', True),
                  (str(len(self.lines)), 'Number of lines in the spectrum', True)]
        for energy, sigma, mult in self.lines:
            lines.append(('{} {} {} {}'.format(energy, sigma, mult, 0), 'Energy (keV), sigma (keV), X multiplicity, Y multiplicity', True))
        return lines

    @property
    def emax(self):
        return max(energy for energy, sigma, intensity in self.lines)

    @property
    def multiplicity(self):
        return sum(intensity for energy, sigma, intensity in self.lines)


class Detector(xrmc_device):

    TYPE = 'detectorarray'

    def __init__(self, parent, name, distance=None, activearea=None, 
                 energydispersive=True, poissonnoise=True, ebinsize=None,
                 orientation_inplane=None, orientation_outplane=None,
                 boutputflux=False, flux=None):
        super(Detector, self).__init__(parent, name)
        self.orientation_inplane = orientation_inplane
        self.orientation_outplane = orientation_outplane
        self.distance = distance
        self.activearea = activearea
        self.energydispersive = energydispersive
        self.ebinsize = ebinsize
        self.poissonnoise = poissonnoise
        self.boutputflux = boutputflux
        self.flux = flux

    @property
    def header(self):
        return ['Detector array parameter file',
                '',
                'The detector reference frame (a.k.a. local coordinate system)',
                'is defined by unit vectors {id, jd, kd} and origin rd.',
                '',
                'The origin origin lies on the geometric center of the detector.',
                'The normal to the detector runs along kd.',
                'A detector can have 1 or more pixels (elliptical pixel useful for',
                'single element detector) with rows and columns parallel to',
                'id and jd respectively.']

    @property
    def code(self):
        lines = super(Detector, self).code
        spectrum = self.parent.spectrum()
        energydispersive = self.energydispersive
        detsize = self.activearea**0.5
        lines += [('SourceName {}'.format(self.parent.sample().name), 'Source input device name', True),
                  ('NPixels 1 1', 'Pixel number (NX x NY)', True),
                  ('PixelSize {} {}'.format(detsize, detsize), 'Pixel Size (cm)', True),
                  ('Shape 1', 'Pixel shape (0 rectangular, 1 elliptical)', True),
                  ('dOmegaLim {}'.format(2*np.pi), 'Limit solid angle (default = 2*PI)', False),
                  ('X 0 0 {}'.format(self.distance), 'origin of the detector in the world frame (x, y, z; cm)', True),
                  ('uk 0 0 1', 'direction of kd in the world frame (x, y, z)', True),
                  ('ui 1 0 0', 'direction of id in the world frame (x, y, z)', True),
                  ('ExpTime 1', 'Exposure time (sec)', True),
                  ('PhotonNum {:d}'.format(int(np.round(self.flux))), 'Multiplicity of simulated events per pixel', True),
                  ('RandomPixelFlag 1', 'Enable random point on pixels (0/1)', True),
                  ('PoissonFlag {:d}'.format(self.poissonnoise), 'Enable Poisson statistic on pix. counts (0/1)', True),
                  ('RoundFlag 1', 'Round pixel counts to integer (0/1)', True),
                  ('ForceDetectFlag {:d}'.format(self.boutputflux), 'MC variance reduction (0/1)', True),
                  ('HeaderFlag 0','Use header in output file (0/1)', True),
                  ('PixelType 0', 'Pixel content type:', True),
                  ('', '0: fluence,      1: energy fluence,', True),
                  ('', '2: fluence(E),   3: energy fluence(E)', True),
                  ('Emin 0', 'Emin', energydispersive),
                  ('Emax {}'.format(spectrum.emax), 'Emax', energydispersive),
                  ('NBins {}'.format(self.nbins), 'NBins', energydispersive),
                  ('SaturateEmin 0', 'Saturate energies lower than Emin (0/1)', energydispersive),
                  ('SaturateEmax 0', 'Saturate energies greater than Emax (0/1)', energydispersive),
                  ('Rotate 0 0 0 1 0 0 {}'.format(self.orientation_inplane), 'rotation around axis passing through (x0,y0,z0; cm), with direction (u,v,w; cm) and rotation angle (deg)', energydispersive),
                  ('Rotate 0 0 0 0 1 0 {}'.format(self.orientation_outplane), 'rotation around axis passing through (x0,y0,z0; cm), with direction (u,v,w; cm) and rotation angle (deg)', energydispersive)]
        return lines

    @property
    def reloutput(self):
        return 'output/{}.dat'.format(self.name)
    
    @property
    def nbins(self):
        if self.energydispersive:
            Emax = self.parent.spectrum().emax
            return int(np.round(Emax/self.ebinsize))
        else:
            return 1

    def result(self):
        """
        N. of scattering orders × N. of energy bins × N. of columns × N. of rows
        """
        filename = os.path.join(self.path, 'output', '{}.dat'.format(self.name))
        shape = self.parent.sample().order + 1, self.nbins, 1, 1
        return np.fromfile(filename, dtype=np.double).reshape(shape)


class Sample(xrmc_device):

    TYPE = 'sample'

    def __init__(self, parent, name, order=1):
        super(Sample, self).__init__(parent, name)
        self.order = order

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
                  ('FluorFlag 1', 'Activate Fluorescence (0/1)', True), 
                  ('ScattOrderNum {}'.format(self.order), 'Maximum scattering order', True)]
        for i in range(self.order+1):
            lines.append(('1', 'Multiplicity of simulated events for order {}'.format(i), True))
        return lines


class Objects(xrmc_device):

    TYPE = 'geom3d'

    def __init__(self, parent, name):
        super(Objects, self).__init__(parent, name)
        self._dict = {}

    @property
    def header(self):
        return ['3D object geometric description file',
                '',
                'Makes objects with quadrics and compositions.']

    @property
    def code(self):
        lines = super(Objects, self).code
        lines += [('QArrName {}'.format(self.parent.quadrics().name), 'Quadric array input device name', True), 
                  ('CompName {}'.format(self.parent.compositions().name), 'Composition input device name', True)]
        for name, (compoundin, compoundout, quadricnames) in self._dict.items():
            lines += [None,
                      ('Object {}'.format(name), '', True),
                      ('{} {}'.format(compoundin, compoundout), 'Phase in, phase out', True),
                      (str(len(quadricnames)), 'Number of quadrics', True),
                      (' '.join(quadricnames), 'Quadric names', True)]
        return lines

    def add(self, name, compoundin, compoundout, quadricnames):
        self._dict[name] = compoundin, compoundout, quadricnames


class Quadrics(xrmc_device):

    TYPE = 'quadricarray'

    def __init__(self, parent, name):
        super(Quadrics, self).__init__(parent, name)
        self._lst = []

    @property
    def header(self):
        return ['Quadric array file',
                '',
                'A quadric (a.k.a. quadratic surface) is a parametric surface (plane, cylinder, sphere, ...):',
                '  x.A.x^T = 0   (where x homogeneous coordinates)',
                'For internal space x.A.x^T < 0 and for external space x.A.x^T > 0']

    @property
    def code(self):
        lines = super(Quadrics, self).code
        for quadric in self._lst:
            name, typ, params = quadric
            lines.append(('{} {}'.format(typ, name), '', True))
            if typ == 'Plane':
                comment = 'Plane through (x0,y0,z0) with normal vector (nx,ny,nz)'
            elif typ.startswith('Cylinder'):
                comment = 'Cylinder with main axis along {} and with coordinates (d0, d1) in the perpendicular plane. (R0, R1) are the elliptical semi-axes.'.format(typ[-1])
            else:
                comment = ''
            lines.append((' '.join(list(map(str, params))), comment, True))
        return lines

    def add_plane(self, name, x0, y0, z0, nx, ny, nz):
        self._lst.append((name, 'Plane', [x0, y0, z0, nx, ny, nz]))

    def add_cylinder(self, name, axis, d0, d1, R0, R1):
        axis = axis.upper()
        if axis not in ['X', 'Y', 'Z']:
            raise ValueError("Axis must be 'X', 'Y' or 'Z'")
        self._lst.append((name, 'Cylinder'+axis, [d0, d1, R0, R1]))


class Compositions(xrmc_device):

    TYPE = 'composition'

    def __init__(self, parent, name):
        super(Compositions, self).__init__(parent, name)
        self._dict = {}

    @property
    def header(self):
        return ['Composition file',
                '',
                'List of compounds, defined by their density and elemental composition (weight fractions)']

    @property
    def code(self):
        lines = super(Compositions, self).code
        for name, (elements, massfractions, density) in self._dict.items():
            lines += [None,
                      ('Phase {}'.format(name), '', True),
                      ('NElem {}'.format(len(elements)), '', True)]
            for formula, wfrac in zip(elements, massfractions):
                lines.append(('{} {}'.format(formula, wfrac), 'mass fraction', True))
            lines += [('Rho {}'.format(density), 'Mass density (g/cm3)', True)]
        return lines

    def add(self, name, elements, massfractions, density):
        self._dict[name] = list(map(str, elements)), massfractions, density


DEVICETOCLASS = {
    Detector.TYPE: Detector,
    Source.TYPE: Source,
    Spectrum.TYPE: Spectrum,
    Sample.TYPE: Sample,
    Objects.TYPE: Objects,
    Quadrics.TYPE: Quadrics,
    Compositions.TYPE: Compositions,
}


if __name__=='__main__':
    # X: synchrotron plane
    # Z: beam direction
    flux = 1e10

    main = xrmc_main(r'Z:\inhouse\wout\tmp\xrmctest')
    main.add_device('source', 'synchrotron', beamsize=1e-4, distance=10)
    spectrum = main.add_device('spectrum', 'dcmspectrum', lines=[(7, 0, 5e4)])

    if True:
        detector = main.add_device('detectorarray', 'leia', distance=0.1,
                        activearea=0.8, energydispersive=False, ebinsize=5e-3,
                        orientation_inplane=120, orientation_outplane=0,
                        flux=flux/spectrum.multiplicity, boutputflux=True)
    else:
        detector = main.add_device('detectorarray', 'idet', distance=0.1,
                        activearea=0.8, energydispersive=False,
                        orientation_inplane=0, orientation_outplane=0,
                        flux=flux/spectrum.multiplicity, boutputflux=False)

    compositions = main.add_device('composition', 'compositions')
    compositions.add('calcite', ['Ca', 'C', 'O'], [0.5, 0.2, 0.3], 2.71)

    quadrics = main.add_device('quadricarray', 'quadrics')
    quadrics.add_plane('layer1_top', 0, 0, 0, 0, 0, -1)
    quadrics.add_plane('layer1_bottom', 0, 0, 10e-4, 0, 0, 1)
    size = 10
    quadrics.add_plane('sample_xmin', -size, 0, 0, -1, 0, 0)
    quadrics.add_plane('sample_xmax', size, 0, 0, 1, 0, 0)
    quadrics.add_plane('sample_ymin', 0, -size, 0, 0, -1, 0)
    quadrics.add_plane('sample_ymax', 0, size, 0, 0, 1, 0)
    
    objects = main.add_device('geom3d', 'objects')
    objects.add('layer1', 'calcite', 'Vacuum', ['layer1_top', 'layer1_bottom', 'sample_xmin', 'sample_xmax', 'sample_ymin', 'sample_ymax'])

    main.add_device('sample', 'flatsample')
    main.save()
    main.simulate()

    data = detector.result()
    print(data.sum())
    #for spectrum in data:
    #    plt.plot(spectrum[:,0,0])
    #plt.show()
