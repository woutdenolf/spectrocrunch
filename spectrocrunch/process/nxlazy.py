import numpy as np

from ..io import edf
from ..io import xiaedf


class LazyFunction(object):
    def __init__(self, samemerge=False):
        self.samemerge = samemerge

    def __str__(self):
        return self._func.__class__.__name__

    def __eq__(self, other):
        return str(self) == str(other)

    def __ne__(self, other):
        return not self.__eq__(other)

    def merge(self, other):
        if self == other:
            return self.samemerge
        else:
            return False


class lazy_transmission(LazyFunction):
    def __call__(self, fluxt, flux0):
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.divide(fluxt, flux0)

    def __str__(self):
        return "transmission"


transmission_func = lazy_transmission()


class lazy_absorbance(LazyFunction):
    def __call__(self, transmission):
        with np.errstate(divide="ignore", invalid="ignore"):
            return -np.log(np.clip(transmission, 0, 1))

    def __str__(self):
        return "absorbance"


absorbance_func = lazy_absorbance()


class lazy_xrfnorm(LazyFunction):
    def __call__(self, xrf, flux, fluxref, xiaimage, detnr):
        if fluxref:
            norm = fluxref.to("Hz").magnitude / xiaedf.normalizer(flux)
        else:
            norm = 1
        if xiaimage:
            xiaimage.onlyicrocr(True)
            xiaimage.exclude_detectors = []
            xiaimage.include_detectors = [detnr]
            stats = xiaimage.stats
            dtcor = xiaedf.deadtimecorrector(stats[..., 0, 0], stats[..., 1, 0])
            dtcor = dtcor.reshape(xrf.shape)
        else:
            dtcor = 1

        return xrf * norm * dtcor

    def __str__(self):
        return "xrfnorm"


xrfnorm_func = lazy_xrfnorm()


class lazy_nanmean(LazyFunction):
    def __init__(self):
        super(lazy_nanmean, self).__init__(samemerge=True)

    def __call__(self, x):
        return np.nanmean(list(x), axis=0)

    def __str__(self):
        return "nanmean"


nanmean_func = lazy_nanmean()


class lazy_nansum(LazyFunction):
    def __init__(self):
        super(lazy_nansum, self).__init__(samemerge=True)

    def __call__(self, x):
        return np.nansum(list(x), axis=0)

    def __str__(self):
        return "nansum"


nansum_func = lazy_nansum()


class lazy_nanmax(LazyFunction):
    def __init__(self):
        super(lazy_nanmax, self).__init__(samemerge=True)

    def __call__(self, x):
        return np.nanmax(list(x), axis=0)

    def __str__(self):
        return "nanmax"


nanmax_func = lazy_nanmax()


class lazy_sum(LazyFunction):
    def __init__(self):
        super(lazy_sum, self).__init__(samemerge=True)

    def __call__(self, x):
        return sum(x)

    def __str__(self):
        return "sum"


sum_func = lazy_sum()


class lazy_readedf(LazyFunction):
    def __init__(self):
        super(lazy_readedf, self).__init__(samemerge=True)

    def __call__(self, x):
        return x

    def __str__(self):
        return "readedf"


readedf_func = lazy_readedf()


class LazyArgument(object):
    def __init__(self, arg):
        self._arg = arg

    def data(self, *args):
        return self._arg

    def __repr__(self):
        return self._arg

    def __str__(self):
        return self.__repr__()


class LazyArgumentEdf(LazyArgument):
    def __init__(self, filename):
        self._filename = filename

    def data(self, *args):
        return edf.edfimage(self._filename).data

    def __repr__(self):
        return self._filename


class LazyArgumentH5Dataset(LazyArgument):
    def __init__(self, path):
        self._path = path

    def data(self, islice, stackdim):
        with self._path.open(mode="r") as dset:

            if stackdim == 0:
                data = dset[islice, ...]
            elif stackdim == 1:
                data = dset[:, islice, :]
            else:
                data = dset[..., islice]

            return data

    def __repr__(self):
        return self._path.__repr__()

    def __str__(self):
        return self.__repr__()


class LazyStackSlice(LazyArgument):
    def __init__(self, func=None, unpackargs=True):
        if func is None:
            self._func = readedf_func
        else:
            self._func = func
        self._args = []
        self._unpackargs = unpackargs

    def data(self, *info):
        if self._unpackargs:
            return self._func(*list(self._arggen(*info)))
        else:
            return self._func(self._arggen(*info))

    def _arggen(self, *info):
        for x in self._args:
            if isinstance(x, LazyArgument):
                yield x.data(*info)
            else:
                yield x

    def appendarg(self, arg):
        if isinstance(arg, self.__class__):
            if self._func.merge(arg._func):
                self._args.extend(arg._args)
                return
        self._args.append(arg)

    def appendarg_edf(self, filename):
        self.appendarg(LazyArgumentEdf(filename))

    def appendarg_h5dataset(self, path):
        self.appendarg(LazyArgumentH5Dataset(path))

    def __repr__(self):
        return "{}({})".format(self._func, ",".join([str(arg) for arg in self._args]))
