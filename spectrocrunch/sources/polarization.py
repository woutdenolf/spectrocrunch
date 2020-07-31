# -*- coding: utf-8 -*-

"""Polarization state of transverse waves
"""

import cmath
import logging
import collections

from ..utils import units
from ..utils import instance
from ..patch.pint import ureg

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

logger = logging.getLogger(__name__)


def jonesnormsq_to_intensity(x):
    """Convert from squared Jones vector norm (V^2/m^2) to intensity (W/m^2)
    """
    return (
        (
            units.Quantity(x / 2.0, "V^2/m^2")
            * (ureg.vacuum_permittivity * ureg.speed_of_light)
        )
        .to("W/m^2")
        .magnitude
    )


def intensity_to_jonesnormsq(x):
    """Convert from intensity (W/m^2) to squared Jones vector norm (V^2/m^2)
    """
    return (
        (
            units.Quantity(x * 2.0, "W/m^2")
            / (ureg.vacuum_permittivity * ureg.speed_of_light)
        )
        .to("V^2/m^2")
        .magnitude
    )


def JonesMatrixRotation(angle):
    """Coordinate transformation matrix

    Args:
        angle(num): vector rotation (degree)

    Returns:
        array: Coordinate transformation (2x2)
    """
    ph = np.radians(angle)
    cosph = np.cos(ph)
    sinph = np.sin(ph)
    return np.array([[cosph, -sinph], [sinph, cosph]])


def MuellerMatrixRotation(angle):
    """Mueller matrix for rotation (equivalent of JonesMatrixRotation)

    Args:
        angle(num): vector rotation (degree)

    Returns:
        array: Mueller matrix (4x4)
    """
    ph = np.radians(2 * angle)
    cosph = np.cos(ph)
    sinph = np.sin(ph)
    return np.array(
        [[1, 0, 0, 0], [0, cosph, -sinph, 0], [0, sinph, cosph, 0], [0, 0, 0, 1]]
    )


def JonesMatrixThomson(polar):
    """Jones matrix for Thomson scattering (assume azimuth=90)

    Args:
        polar(num): scattering angle (degree)

    Returns:
        array: Coordinate transformation (2x2)
    """
    return np.array([[1, 0], [0, np.cos(np.radians(polar))]])


def ThomsonRotationAngle(azimuth):
    """Reference frame rotation so that JonesMatrixThomson holds

    Args:
        azimuth(num): scattering azimuth (degree)

    Returns:
        angle: rotation of the axes in degrees
    """
    # First component axis perpendicular to the scattering plane
    return azimuth - 90


def ComptonRotationAngle(azimuth):
    """Reference frame rotation so that MuellerMatrixCompton holds

    Args:
        azimuth(num): scattering azimuth (degree)

    Returns:
        angle: rotation of the axes in degrees
    """
    # First component axis perpendicular to the scattering plane
    return azimuth - 90


def MuellerMatrixThomson(azimuth, polar):
    """Mueller matrix for Thomson scattering

    Args:
        azimuth(num): scattering azimuth (degree)
        polar(num): scattering angle (degree)

    Returns:
        array: Mueller matrix (4x4)
    """

    costh = np.cos(np.radians(polar))
    cossq_polar = costh ** 2
    a = (1 + cossq_polar) / 2.0
    b = 1 - a  # = sinsq_polar/2

    ph = -np.radians(2 * ThomsonRotationAngle(azimuth))
    cosph = np.cos(ph)
    sinph = np.sin(ph)

    return np.array(
        [
            [a, b * cosph, -b * sinph, 0],
            [b, a * cosph, -a * sinph, 0],
            [0, costh * sinph, costh * cosph, 0],
            [0, 0, 0, costh],
        ]
    )


def kev_to_mc2(x):
    return ureg.Quantity(x, "keV").to("m_e*c^2").magnitude


def MuellerMatrixCompton(azimuth, polar, energy):
    """Mueller matrix for Compton scattering

    Args:
        azimuth(num): scattering azimuth (degree)
        polar(num): scattering angle (degree)
        energy(num): incident photon energy (keV)

    Returns:
        array: Mueller matrix (4x4)
    """
    # Same incident and scattered reference frame as MuellerMatrixThomson
    costh = np.cos(np.radians(polar))
    cossq_polar = costh ** 2
    a = (1 + cossq_polar) / 2.0
    b = 1 - a  # = sinsq_polar/2

    ph = -np.radians(2 * ComptonRotationAngle(azimuth))
    cosph = np.cos(ph)
    sinph = np.sin(ph)

    k = kev_to_mc2(energy)
    ksc = k / (1 + k * (1 - costh))
    c = (k - ksc) * (1 - costh) / 2.0

    return np.array(
        [
            [a + c, b * cosph, -b * sinph, 0],
            [b, a * cosph, -a * sinph, 0],
            [0, costh * sinph, costh * cosph, 0],
            [0, 0, 0, costh * (1 + c)],
        ]
    ) * (ksc ** 2 / k ** 2)


class Jones(object):
    """Parameterization of second order statistics of a transverse wave when fully polarized
    """

    def __init__(self, V0, V1):
        """Two perpendicular components of a transverse wave:

        .. math::
            \vec{E}(t,\vec{x}) = (V_0.\hat{e}_0 + V_1.\hat{e}_1).e^{i(\vec{k}\cdot\vec{x} - \omega t)}

            \hat{e}_0\times\hat{e}_1 = \hat{e}_2
            \hat{e}_1\times\hat{e}_2 = \hat{e}_0
            \vec{k}\cdot\hat{e}_2 = 0

        The components $V_0$ and $V_1$ are complex numbers who's moduli have units $V/m$.

        """
        self.V = [V0, V1]
        if self.V.size != 2:
            raise ValueError("Expected a Jones vector (2 components)")

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return np.array_equal(self.V, other.V)
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __str__(self):
        params = (
            ("intensity", "W/m^2"),
            ("dolp", ""),
            ("polangle", "deg"),
            ("handedness", ""),
            ("phase0", "deg"),
        )

        p = self.to_params()
        out = []
        for k, u in params:
            v = p[k]
            if u:
                out.append("{} = {} ({})".format(k, v, u))
            else:
                out.append("{} = {}".format(k, v))
        return "\n".join(out)

    @property
    def V(self):
        return self._V

    @V.setter
    def V(self, value):
        self._V = np.array(value, dtype=complex)

    def __getitem__(self, index):
        return self.V[index]

    @classmethod
    def _hermitian_inner_product(cls, U, V):
        """The standard Hermitian inner product in $\mathbb{C}^n$
        """
        return U.dot(V.conjugate())

    @classmethod
    def _hermitian_outer_product(cls, U, V):
        """The standard Hermitian outer product in $\mathbb{C}^n$
        """
        return np.outer(U, V.conjugate())

    @classmethod
    def _hermitian_sqnorm(cls, V):
        """The standard Hermitian squared norm in $\mathbb{C}^n$
        """
        return cls._hermitian_inner_product(V, V).real  # imag == 0

    @property
    def sqnorm(self):
        """Squared norm of the electric field vector in $V^2/m^2$
        """
        return self._hermitian_sqnorm(self.V)

    @property
    def norm(self):
        """Norm of the electric field vector in $V/m$
        """
        return np.sqrt(self.sqnorm)

    def efield(self, phase):
        """Electric field vector (norm in $V/m$)

        ..math::

            \vec{E}(t,\vec{x}) = (V_0.\hat{e}_0 + V_1.\hat{e}_1).e^{-2\pi\nu i)}

        Args:
            phase(num or array): $360.\nu$ in degrees

        Returns:
            tuple: E0,E1
        """
        return np.outer(self.V, np.exp(-1j * np.radians(phase))).real

    @property
    def intensity(self):
        """Expectation of the Poynting vector norm in $W/m^2$
        """
        return jonesnormsq_to_intensity(self.sqnorm)

    @intensity.setter
    def intensity(self, value):
        self.V *= np.sqrt(value / self.intensity)

    @property
    def coherency_matrix(self):
        # https://cel.archives-ouvertes.fr/cel-00521501v4/file/Polarization_handout.pdf
        # C = <J.J^T*>  (J as column vector)
        # C00 = V0.V0* = |V0|^2
        # C11 = V1.V1* = |V1|^2
        # C01 = V0.V1* = |V0|.|V1|.exp(i(ph0-ph1)) = |V0|.|V1|.exp(-i.delta)
        # C10 = V1.V0* = |V0|.|V1|.exp(i(ph1-ph0)) = |V0|.|V1|.exp(i.delta)
        # delta = ph1-ph0
        return self._hermitian_outer_product(self.V, self.V)

    @staticmethod
    def _coherency_matrix_to_jones(C, phase0=0):
        if not np.allclose(C, C.T.conjugate()):
            raise RuntimeError("The coherency matrix must be Hermitian")

        C00, C01, C10, C11 = C.flatten()
        ph0 = np.radians(phase0)
        V0 = np.sqrt(C00) * np.exp(ph0 * 1j)
        dph = cmath.phase(C10)
        ph1 = dph + ph0
        V1 = np.sqrt(C11) * np.exp(ph1 * 1j)
        return [V0, V1]

    @coherency_matrix.setter
    def coherency_matrix(self, C):
        self.V = self._coherency_matrix_to_jones(C)

    @classmethod
    def from_coherency_matrix(cls, C, phase0=0):
        return cls(*cls._coherency_matrix_to_jones(C, phase0=phase0))

    def to_numpy(self):
        return self.V

    def to_stokes(self, silent=False):
        if not silent:
            logger.warning("Absolute phase is lost")
        return Stokes.from_coherency_matrix(self.coherency_matrix)

        # Equivalent:
        # V0,V1 = abs(self.V)
        # delta = self.phase_difference
        # S0 = V0**2+V1**2
        # S1 = V0**2-V1**2
        # S2 = 2*V0*V1*np.cos(delta)
        # S3 = 2*V0*V1*np.sin(delta)
        # return Stokes(S0,S1,S2,S3)

    @property
    def dop(self):
        """Degree of polarization = Ipol / Itot = 1 (use Stokes for partial or unpolarized)
        """
        # sqrt(S1^2 + S2^2 + S3^2)/S0 = 1
        return 1

    @property
    def dolp(self):
        """Degree of linear polarization
        """
        return self.to_stokes(silent=True).dolp

    @property
    def hdolp(self):
        """Horizontal degree of linear polarization
        """
        return self.to_stokes(silent=True).hdolp

    @property
    def docp(self):
        """Degree of circular polarization
        """
        return self.to_stokes(silent=True).docp

    @property
    def polangle(self):
        """Polarization angle
        """
        return self.to_stokes(silent=True).polangle

    @property
    def handedness(self):
        """Handedness
        """
        return self.to_stokes(silent=True).handedness

    @property
    def phase0(self):
        """Phase of the first component
        """
        return np.degrees(cmath.phase(self[0]))

    @property
    def phase1(self):
        """Phase of the second component
        """
        return np.degrees(cmath.phase(self.V[1]))

    @property
    def phase_difference(self):
        """Phase difference phase1-phase0
        """
        return self.phase1 - self.phase0

    def to_params(self):
        S = self.to_stokes(silent=True)
        return {
            "intensity": self.intensity,
            "dolp": S.dolp,
            "polangle": S.polangle,
            "handedness": S.handedness,
            "phase0": self.phase0,
        }

    @classmethod
    def from_params(cls, **kwargs):
        kwargs["dop"] = 1
        phase0 = kwargs.pop("phase0", 0)
        return Stokes.from_params(**kwargs).to_jones(phase0=phase0)

    def plot_efield(self, animate=False, **kwargs):
        """Upstream view of the electric field vector in the plane perpendicular to the propagation direction
        """
        # Axes with labels
        ax = plt.gca()
        fig = ax.get_figure()
        ax.set_xlabel("$E_0$ (V/m)")
        ax.set_ylabel("$E_1$ (V/m)")
        p = self.to_params()
        ax.set_title(
            "dolp={:.02f}, polangle={:.02f}$^\circ$, docp={:.02f} ({}), $\phi_0$={:.02f}$^\circ$".format(
                p["dolp"], p["polangle"], self.docp, p["handedness"], p["phase0"]
            )
        )

        # Set limits
        mV0, mV1 = abs(self.V)
        if abs(mV0) < abs(mV1) * 1e-2:
            mV0 = mV1 * 0.1
        elif abs(mV1) < abs(mV0) * 1e-2:
            mV1 = mV0 * 0.1
        ax.set_xlim(-mV0, mV0)
        ax.set_ylim(-mV1, mV1)
        ax.set_aspect("equal")

        # Plot ellipse
        def phase(ph):
            return np.exp(1j * ph)

        x, y = self.efield(np.linspace(0, 360, 100))
        lines = ax.plot(x, y, **kwargs)
        color = lines[0].get_color()

        # Add handedness arrows
        head_width = max(mV0, mV1) * 0.1
        for a in [0, 45, 90, 135, 180, 225, 270, 315]:
            x, y = self.efield(np.array([a, a + 5]))
            ax.arrow(
                x[0],
                y[0],
                x[1] - x[0],
                y[1] - y[0],
                lw=0,
                length_includes_head=True,
                head_width=head_width,
                color=color,
            )

        # Add rotation electric field vector
        if animate:
            annotation = ax.annotate(
                "",
                xy=(0, 0),
                xytext=(0, 0),
                xycoords="data",
                textcoords="data",
                arrowprops={"arrowstyle": "<-"},
                horizontalalignment="left",
                verticalalignment="bottom",
                color=color,
            )

            def xy(ph):
                return (self.V * phase(ph)).real

            def init():
                return (annotation,)

            def drawframe(ph):
                annotation.set_position(xy(ph))
                return (annotation,)

            ani = animation.FuncAnimation(
                fig,
                drawframe,
                frames=np.linspace(0, 2 * np.pi, 100),
                blit=True,
                interval=100,
                repeat=False,
                init_func=init,
            )
        else:
            ani = None

        return ani

    def rotate(self, angle):
        """Rotate in the plane perpendicular to the propagation direction

        Args:
            angle(num): rotation of the axes in degrees

        Returns:
            Jones
        """
        R = JonesMatrixRotation(-angle)
        return self.__class__(*R.dot(self.V))

    def thomson_scattering(self, azimuth, polar):
        """Thomson scattering (incident wave along the Z-axis)

        Args:
            azimuth(num): scattering direction in spherical coordinates (deg)
            polar(num): scattering direction in spherical coordinates (deg)

        Returns:
            Stokes
        """
        # We could do C' = (M.R).C.(M.R)^H but then we loose the phase
        J = self.rotate(ThomsonRotationAngle(azimuth))
        Mth = JonesMatrixThomson(polar)
        C = Mth.dot(J.coherency_matrix).dot(Mth.T.conjugate())
        return self.from_coherency_matrix(C, phase0=J.phase0)


class Stokes(object):
    """Parameterization of second order statistics of a transverse wave when partially polarized
    """

    # Convention: propagation in the direction of e_2 (Z-axis)
    pauli422 = np.array(
        [[[1, 0], [0, 1]], [[1, 0], [0, -1]], [[0, 1], [1, 0]], [[0, -1j], [1j, 0]]]
    )
    # inverse = pauli44.T.conjugate()/2
    pauli44 = np.array([[1, 0, 0, 1], [1, 0, 0, -1], [0, 1, 1, 0], [0, 1j, -1j, 0]])

    def __init__(self, S0, S1, S2, S3):
        """S0 in V^2/m^2
        """
        self.S = [S0, S1, S2, S3]
        if self.S.size != 4:
            raise ValueError("Expected a 2D Stokes vector (4 components)")

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return np.array_equal(self.S, other.S)
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __str__(self):
        params = (
            ("intensity", "W/m^2"),
            ("dop", ""),
            ("dolp", ""),
            ("polangle", "deg"),
            ("handedness", ""),
        )

        p = self.to_params()
        out = []
        for k, u in params:
            v = p[k]
            if u:
                out.append("{} = {} ({})".format(k, v, u))
            else:
                out.append("{} = {}".format(k, v))
        return "\n".join(out)

    @property
    def S(self):
        return self._S

    @S.setter
    def S(self, value):
        self._S = np.array(value)

    def __getitem__(self, index):
        return self.S[index]

    @property
    def intensity(self):
        """Expectation of the Poynting vector norm in $W/m^2$
        """
        return jonesnormsq_to_intensity(self.S[0])

    @intensity.setter
    def intensity(self, value):
        self.S *= value / self.intensity

    @property
    def intensity_polarized(self):
        # proportional to sqrt(S1^2 + S2^2 + S3^2)
        return jonesnormsq_to_intensity(np.sqrt((self.S[1:] ** 2).sum()))

    @property
    def intensity_unpolarized(self):
        return jonesnormsq_to_intensity(self.S[0] - np.sqrt((self.S[1:] ** 2).sum()))

    @property
    def dop(self):
        """Degree of polarization = Ipol / Itot
        """
        # sqrt(S1^2 + S2^2 + S3^2)/S0
        return np.sqrt((self.S[1:] ** 2).sum()) / self.S[0]

    @property
    def dolp(self):
        """Degree of linear polarization
        """
        # sqrt(S1^2 + S2^2)/S0
        return np.sqrt((self.S[1:3] ** 2).sum()) / self.S[0]

    @property
    def hdolp(self):
        """Horizontal degree of linear polarization
        """
        # S1/S0 = (|V0|^2 - |V1|^2)/(|V0|^2 + |V1|^2)
        return self.S[1] / self.S[0]

    @property
    def docp(self):
        """Degree of circular polarization
        """
        # S3/S0
        return self.S[3] / self.S[0]

    @property
    def polangle(self):
        """Polarization angle
        """
        return np.degrees(np.arctan2(self.S[2], self.S[1]) / 2.0)

    @property
    def ellipangle(self):
        """Ellipticity angle
        """
        return np.arctan2(self.S[3], np.sqrt((self.S[1:3] ** 2).sum())) / 2.0

    @property
    def handedness(self):
        if self.S[3] > 0:
            return "left"
        elif self.S[3] < 0:
            return "right"
        else:
            return "none"

    @property
    def phase_difference(self):
        """Phase difference between Jones components
        """
        # delta = ph1-ph0
        return np.degrees(np.arctan2(self.S[3], self.S[2]))

    @property
    def coherency_matrix(self):
        # a.k.a. density matrix
        # C00 = (S0+S1)/2
        # C01 = (S2-i.S3)/2
        # C10 = (S2+i.S3)/2
        # C11 = (S0-S1)/2
        C00, C01, C10, C11 = self.pauli44.T.conjugate().dot(self.S) / 2.0
        return np.array([[C00, C01], [C10, C11]])

        # Equivalent:
        # return self.pauli422.T.dot(self.S).T/2.

    @classmethod
    def _coherency_matrix_to_stokes(cls, C):
        if not np.allclose(C, C.T.conjugate()):
            raise RuntimeError("The coherency matrix must be Hermitian")

        # delta = ph1-ph0
        # S0 = C00+C11 = |V0|^2 + |V1|^2
        # S1 = C00-C11 = |V0|^2 - |V1|^2
        # S2 = C01+C10 = |V0|.|V1|.[exp(-i.delta)+exp(i.delta)] = 2.|V0|.|V1|.cos(delta)
        # S3 = i.(C01-C10) = i.|V0|.|V1|.[exp(-i.delta)-exp(i.delta)] = 2.|V0|.|V1|.sin(delta)
        # Sometimes you see another convention of handedless: S3 = i.(C10-C01)
        return cls.pauli44.dot(C.flatten()).real

        # Equivalent:
        # return cls.pauli422.dot(C).trace(axis1=1,axis2=2).real

    @coherency_matrix.setter
    def coherency_matrix(self, C):
        self.S = self._coherency_matrix_to_stokes(C)

    @classmethod
    def from_coherency_matrix(cls, C):
        return cls(*cls._coherency_matrix_to_stokes(C))

    def to_numpy(self):
        return self.S

    def decompose(self):
        """Decompose in polarized and unpolarized fraction
        """
        # unpolarized = [(1-dop)*S0,0,0,0]
        # polarized = [dop*S0,S1,S2,S3]
        tmp = np.sqrt((self.S[1:] ** 2).sum())
        unpol = self.__class__(self.S[0] - tmp, 0, 0, 0)
        pol = self.__class__(tmp, self.S[1], self.S[2], self.S[3])
        return {"pol": pol, "unpol": unpol}

    def to_jones(self, allowloss=False, atol=1e-7, phase0=0):
        d = self.dop - 1
        if d == 0:
            pol = self
        else:
            if abs(d) > atol:
                if allowloss:
                    logger.warning("Unpolarized fraction is lost")
                else:
                    raise RuntimeError("Unpolarized fraction is lost")
            pol = self.decompose()["pol"]
        return Jones.from_coherency_matrix(pol.coherency_matrix, phase0=phase0)

    @classmethod
    def from_params(
        cls, intensity=None, dop=None, dolp=None, polangle=None, handedness=None
    ):
        if dolp > dop:
            raise RuntimeError(
                "Degree of linear polarization must be smaller or equal to the degree of polarization."
            )
        S0 = intensity_to_jonesnormsq(intensity)
        S12 = (dolp * S0) ** 2  # S12 = S1^2 + S2^2
        S123 = (dop * S0) ** 2  # S123 = S1^2 + S2^2 + S3^2
        polangle2 = np.radians(2 * polangle)  # in [-180,180]
        S2S1 = np.tan(polangle2)  # S2S1 = S2/S1

        S3 = np.sqrt(S123 - S12)
        if handedness == "right":
            S3 = -S3

        # S12 = S1^2 + S2^2 = S1^2 *(1+S2S1^2)
        S1sq = S12 / (1.0 + S2S1 ** 2)
        S1 = np.sqrt(S1sq) * np.sign(np.cos(polangle2))
        S2 = np.sqrt(S12 - S1sq) * np.sign(np.sin(polangle2))

        return cls(S0, S1, S2, S3)

    def to_params(self):
        return {
            "intensity": self.intensity,
            "dop": self.dop,
            "dolp": self.dolp,
            "polangle": self.polangle,
            "handedness": self.handedness,
        }

    @classmethod
    def jones_to_mueller(cls, M):
        """Convert Jones matrix to Mueller matrix

        Args:
            M(array): Jones matrix

        Returns:
            array: Mueller matrix
        """

        return (
            cls.pauli44.dot(np.kron(M, M.T.conjugate())).dot(cls.pauli44.T.conjugate())
            / 2.0
        )

        # Equivalent
        # return cls.pauli422.dot(M).dot(cls.pauli422).dot(M.T.conjugate()).trace(axis1=1,axis2=3)/2.

    def rotate(self, azimuth):
        """Rotate in the plane perpendicular to the propagation direction

        Args:
            azimuth(num): rotation of the axes in degrees

        Returns:
            Stokes
        """
        M = MuellerMatrixRotation(-azimuth)
        return self.__class__(*M.dot(self.S))

    def thomson_scattering(self, azimuth, polar):
        """Thomson scattering (incident wave along the Z-axis)

        Args:
            azimuth(num): scattering direction in spherical coordinates (deg)
            polar(num): scattering direction in spherical coordinates (deg)

        Returns:
            Stokes
        """

        M = MuellerMatrixThomson(azimuth, polar)
        return self.__class__(*M.dot(self.S))

        # Equivalent
        # M1 = MuellerMatrixRotation(-ThomsonRotationAngle(azimuth)) # rotate so that azimuth = 90 deg
        # M2 = self.jones_to_mueller(JonesMatrixThomson(polar)).real
        # return self.__class__(*M2.dot(M1).dot(self.S))

    def compton_scattering(self, azimuth, polar, energy):
        """Compton scattering (incident wave along the Z-axis)

        Args:
            azimuth(num): scattering direction in spherical coordinates (deg)
            polar(num): scattering direction in spherical coordinates (deg)
            energy(num): incident photon energy in keV

        Returns:
            Stokes
        """

        M = MuellerMatrixCompton(azimuth, polar, energy)
        return self.__class__(*M.dot(self.S))

        # Equivalent
        # M1 = MuellerMatrixRotation(-ComptonRotation(azimuth)) # rotate so that azimuth = 90 deg
        # M2 = self.jones_to_mueller(JonesMatrixCompton(polar)).real
        # return self.__class__(*M2.dot(M1).dot(self.S))

    @property
    def thomson_K(self):
        """Directional dependent part of the Thomson scattering intensity (angles in radians)
        """
        if self.S[1] == 0 and self.S[2] == 0:

            def K(azimuth=None, polar=None):
                return (1 + np.cos(polar) ** 2) / 2.0

        else:
            S10 = self.S[1] / self.S[0]
            S20 = self.S[2] / self.S[0]

            def K(azimuth=None, polar=None):
                a = (1 + np.cos(polar) ** 2) / 2.0
                ph = 2 * azimuth
                cosph = np.cos(ph)
                sinph = np.sin(ph)
                return a - (1 - a) * (S10 * cosph + S20 * sinph)

        # Equivalent
        # b = np.tan(2*np.radians(self.polangle))
        # P = self.hdolp
        # ph = -np.radians(2*ThomsonRotationAngle(azimuth))
        # cosph = np.cos(ph)
        # sinph = np.sin(ph)
        # K = a + (1-a)*P*(cosph - b*sinph)
        # K = a + np.sign(self.S[1])*(1-a)*self.dolp/np.sqrt(b**2 + 1)*(cosph - b*sinph)
        # K = a + (1-a)*(P*cosph - np.sqrt(1-P**2)*np.cos(np.radians(self.phase_difference))*sinph)

        return K

    def compton_K(self, energy):
        """Directional dependent part of the Compton scattering intensity (angles in radians)
        """
        k = kev_to_mc2(energy)

        if self.S[1] == 0 and self.S[2] == 0:

            def K(azimuth=None, polar=None):
                costh = np.cos(polar)
                a = (1 + costh ** 2) / 2.0
                ph = 2 * azimuth
                cosph = np.cos(ph)
                sinph = np.sin(ph)
                ksc = k / (1 + k * (1 - costh))
                c = (k - ksc) * (1 - costh) / 2.0
                return ksc ** 2 / k ** 2 * (a + c)

        else:
            S10 = self.S[1] / self.S[0]
            S20 = self.S[2] / self.S[0]

            def K(azimuth=None, polar=None):
                costh = np.cos(polar)
                a = (1 + costh ** 2) / 2.0
                ph = 2 * azimuth
                cosph = np.cos(ph)
                sinph = np.sin(ph)
                ksc = k / (1 + k * (1 - costh))
                c = (k - ksc) * (1 - costh) / 2.0

                return (
                    ksc ** 2 / k ** 2 * (a + c - (1 - a) * (S10 * cosph + S20 * sinph))
                )

        return K

    def thomson_intensity(self, azimuth, polar):
        """
        Args:
            azimuth(num): scattering direction in spherical coordinates (deg)
            polar(num): scattering direction in spherical coordinates (deg)

        Returns:
            num
        """
        return self.intensity * self.thomson_K(np.radians(azimuth), np.radians(polar))

    def compton_intensity(self, azimuth, polar, energy):
        """
        Args:
            azimuth(num): scattering direction in spherical coordinates (deg)
            polar(num): scattering direction in spherical coordinates (deg)
            energy(num): incident photon energy in keV

        Returns:
            num
        """
        return self.intensity * self.compton_K(energy)(
            np.radians(azimuth), np.radians(polar)
        )

    def plot_efield(self, **kwargs):
        self.to_jones().plot_efield(**kwargs)
