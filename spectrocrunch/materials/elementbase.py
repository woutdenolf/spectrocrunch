from ..utils.hashable import CompHashable
from ..utils.copyable import Copyable
from ..math.utils import weightedsum
from . import xrayspectrum
from ..utils import instance
from ..patch.pint import ureg
import numpy as np


def refractive_index_factor(energy, density):
    """Factor in g/mol"""
    return ureg.Quantity(energy, "keV").to("cm", "spectroscopy") ** 2 * (
        ureg.classical_electron_radius
        * ureg.particles_per_mol
        * ureg.Quantity(density, "g/cm^3")
        / (2 * np.pi)
    )


def refractive_index_delta_calc(energy, e_wfrac, density, **kwargs):
    delta = sum(e_wfrac[e] / e.MM * e.scatfact_real(energy, **kwargs) for e in e_wfrac)
    delta = ureg.Quantity(delta, "mol/g") * refractive_index_factor(energy, density)
    return delta.to("dimensionless").magnitude


def refractive_index_beta_calc(energy, e_wfrac, density, **kwargs):
    # TODO: modify sign in scatfact_imag?
    beta = -sum(e_wfrac[e] / e.MM * e.scatfact_imag(energy, **kwargs) for e in e_wfrac)
    beta = ureg.Quantity(beta, "mol/g") * refractive_index_factor(energy, density)
    return beta.to("dimensionless").magnitude


class ElementBase(Copyable, CompHashable):
    def refractive_index_delta(self, E, fine=False, decomposed=False, **kwargs):
        """n = 1-delta-i*beta"""
        if hasattr(self, "structure") and fine:
            environ = self
        else:
            environ = None
        return refractive_index_delta_calc(
            E, self.elemental_massfractions(), self.density, environ=environ, **kwargs
        )

    def refractive_index_beta(self, E, fine=False, decomposed=False, **kwargs):
        """n = 1-delta-i*beta"""
        if hasattr(self, "structure") and fine:
            environ = self
        else:
            environ = None
        return refractive_index_beta_calc(
            E, self.elemental_massfractions(), self.density, environ=environ, **kwargs
        )

    def refractive_index_real(self, E, **kwargs):
        """Real part of the refractive index"""
        return 1 - self.refractive_index_delta(E)

    def refractive_index_imag(self, E, **kwargs):
        """Imaginary part of the refractive index"""
        return -self.refractive_index_beta(E)

    def xrayspectrum(self, E, source=None, weights=None, emin=0, emax=None, **kwargs):
        E = instance.asarray(E)
        if emax is None:
            emax = np.max(E)
        self.markabsorber(energybounds=(emin, emax))

        spectrum = xrayspectrum.Spectrum()

        if source is None:
            spectrum.update(
                self.fluorescence_cross_section_lines(E, decomposed=False, **kwargs)
            )
            spectrum[xrayspectrum.RayleighLine(E)] = self.rayleigh_cross_section(
                E, decomposed=False, **kwargs
            )
            spectrum[xrayspectrum.ComptonLine(E)] = self.compton_cross_section(
                E, decomposed=False, **kwargs
            )
            spectrum.type = spectrum.TYPES.crosssection
        else:
            spectrum.density = self.density
            spectrum.update(
                self.diff_fluorescence_cross_section(E, decomposed=False, **kwargs)
            )
            spectrum[xrayspectrum.RayleighLine(E)] = self.diff_rayleigh_cross_section(
                E, source=source, decomposed=False, **kwargs
            )
            spectrum[xrayspectrum.ComptonLine(E)] = self.diff_compton_cross_section(
                E, source=source, decomposed=False, **kwargs
            )
            spectrum.type = spectrum.TYPES.diffcrosssection

        spectrum.density = self.density
        spectrum.xlim = [emin, emax]
        spectrum.title = str(self)
        spectrum.geomkwargs = kwargs
        spectrum.apply_weights(weights)

        return spectrum

    def fisxgroups(self, emin=0, emax=np.inf):
        self.markabsorber(energybounds=[emin, emax])
        return {el: el.shells for el in self.elements}

    pymcamaterial_prefix = "Material_"
    pymcacomment_prefix = "From spectrocrunch: "

    @property
    def pymcaname(self):
        return self.pymcamaterial_prefix + self.name

    @property
    def pymcacomment(self):
        return self.pymcacomment_prefix + self.name

    @classmethod
    def namefrompymca(cls, string):
        for prefix in [cls.pymcamaterial_prefix, cls.pymcacomment_prefix]:
            if string.startswith(prefix):
                return string[len(prefix) :]
        return string

    def absorbance(self, energy, thickness, weights=None, decomposed=False, **kwargs):
        muL = self.mass_att_coeff(energy, decomposed=decomposed, **kwargs) * (
            thickness * self.density
        )
        # TODO: decomposed -> apply recursively
        if weights is None:
            return muL
        else:
            return weightedsum(muL, weights=weights)

    def transmission(self, energy, thickness, decomposed=False, **kwargs):
        A = self.absorbance(energy, thickness, decomposed=decomposed, **kwargs)
        if decomposed:
            return A  # TODO: apply recursively
        else:
            return np.exp(-A)
