# -*- coding: utf-8 -*-

from sage.all_cmdline import *  # import sage library

_sage_const_2 = Integer(2)
_sage_const_1 = Integer(1)
_sage_const_0 = Integer(0)
_sage_const_1j = ComplexNumber(0, "1")
_sage_const_4 = Integer(4)  #!/usr/bin/env sage

from sage.all import *
from sage.symbolic.integration.integral import definite_integral


def MuellerMatrixRotation(angle):
    return matrix(
        SR,
        _sage_const_4,
        _sage_const_4,
        [
            _sage_const_1,
            _sage_const_0,
            _sage_const_0,
            _sage_const_0,
            _sage_const_0,
            cos(_sage_const_2 * angle),
            -sin(_sage_const_2 * angle),
            _sage_const_0,
            _sage_const_0,
            sin(_sage_const_2 * angle),
            cos(_sage_const_2 * angle),
            _sage_const_0,
            _sage_const_0,
            _sage_const_0,
            _sage_const_0,
            _sage_const_1,
        ],
    )


def JonesMatrixThomson(angle):
    return matrix(
        SR,
        _sage_const_2,
        _sage_const_2,
        [_sage_const_1, _sage_const_0, _sage_const_0, cos(angle)],
    )


pauli44 = matrix(
    SR,
    _sage_const_4,
    _sage_const_4,
    [
        _sage_const_1,
        _sage_const_0,
        _sage_const_0,
        _sage_const_1,
        _sage_const_1,
        _sage_const_0,
        _sage_const_0,
        -_sage_const_1,
        _sage_const_0,
        _sage_const_1,
        _sage_const_1,
        _sage_const_0,
        _sage_const_0,
        _sage_const_1j,
        -_sage_const_1j,
        _sage_const_0,
    ],
)
pauli44i = pauli44.H / _sage_const_2


def jones_to_mueller(M):
    return pauli44 * M.tensor_product(M.H) * pauli44i


def mueller_to_jones(M):
    return (
        pauli44i * M * pauli44
    )  # decompose in M.tensor_product(M.H) not generally possible


def thomson():
    var("theta,phi,beta,S0,S1,S2,S3,a", domain="real")

    S = vector([S0, S1, S2, S3])

    M1 = MuellerMatrixRotation(beta)
    M2 = jones_to_mueller(JonesMatrixThomson(theta))
    # print mueller_to_jones(M2)

    M2 = M2.substitute(theta=acos(sqrt(_sage_const_2 * a - _sage_const_1))).simplify()
    M = (M2 * M1).simplify()
    print "\nMueller matrix:"
    print M
    S_scat = (M * S).simplify()

    I_scat = S_scat[_sage_const_0].substitute(beta=pi / _sage_const_2 - phi).simplify()
    print "\nScattered intensity:"
    print I_scat

    It_scat = I_scat.substitute(
        a=(_sage_const_1 + cos(theta) ** _sage_const_2) / _sage_const_2
    )

    print "\nScattered intensity (unpol):"
    print It_scat.substitute(S1=_sage_const_0).substitute(
        S2=_sage_const_0
    ).simplify_full()

    print "\nScattered intensity (horizontal linear pol.):"
    print It_scat.substitute(S1=S0).substitute(S2=_sage_const_0).simplify_full()

    It_scat = definite_integral(It_scat * sin(theta), theta, _sage_const_0, pi)
    It_scat = definite_integral(It_scat, phi, _sage_const_0, _sage_const_2 * pi)
    print "\nThomson cross-section (units of r_e^2):"
    print It_scat / S0


def compton():
    var("theta,phi,beta,S0,S1,S2,S3,a,c,costheta,s,E,Esc", domain="real")
    assume(E > _sage_const_0)
    assume(Esc > _sage_const_0)

    S = vector([S0, S1, S2, S3])

    M1 = MuellerMatrixRotation(beta)
    M2 = (
        matrix(
            SR,
            _sage_const_4,
            _sage_const_4,
            [
                a + c,
                _sage_const_1 - a,
                _sage_const_0,
                _sage_const_0,
                _sage_const_1 - a,
                a,
                _sage_const_0,
                _sage_const_0,
                _sage_const_0,
                _sage_const_0,
                cos(theta),
                _sage_const_0,
                _sage_const_0,
                _sage_const_0,
                _sage_const_0,
                cos(theta) * (_sage_const_1 + c),
            ],
        )
        * s
    )
    # print mueller_to_jones(M2)
    M2 = M2.substitute(theta=acos(sqrt(_sage_const_2 * a - _sage_const_1))).simplify()

    M = (M2 * M1).simplify()
    print "\nMueller matrix:"
    print M
    S_scat = (M * S).simplify()

    I_scat = S_scat[_sage_const_0].substitute(beta=pi / _sage_const_2 - phi).simplify()
    print "\nScattered intensity:"
    print I_scat

    # Energy in units of m_e*c^2
    It_scat = I_scat.substitute(
        a=(_sage_const_1 + cos(theta) ** _sage_const_2) / _sage_const_2
    )
    It_scat = It_scat.substitute(s=Esc ** _sage_const_2 / E ** _sage_const_2)
    It_scat = It_scat.substitute(
        c=(E - Esc) / _sage_const_2 * (_sage_const_1 - cos(theta))
    )
    It_scat = It_scat.simplify_full()

    print "\nScattered intensity (unpol):"
    print It_scat.substitute(S1=_sage_const_0).substitute(
        S2=_sage_const_0
    ).simplify_full()

    It_scatu1 = It_scat.substitute(S1=_sage_const_0).substitute(S2=_sage_const_0)
    It_scatu2 = (
        Esc ** _sage_const_2
        / E ** _sage_const_2
        * S0
        * (E / Esc + Esc / E - sin(theta) ** _sage_const_2)
        / _sage_const_2
    )
    It_scatu1 = (
        It_scatu1.substitute(Esc=E / (_sage_const_1 + E * (_sage_const_1 - cos(theta))))
        .simplify_full()
        .simplify_trig()
    )
    It_scatu2 = (
        It_scatu2.substitute(Esc=E / (_sage_const_1 + E * (_sage_const_1 - cos(theta))))
        .simplify_full()
        .simplify_trig()
    )
    print bool(It_scatu1 == It_scatu2)

    It_scatu1 = It_scat.substitute(S1=S0).substitute(S2=_sage_const_0)
    It_scatu2 = (
        Esc ** _sage_const_2
        / E ** _sage_const_2
        * S0
        * (
            E / Esc
            + Esc / E
            - _sage_const_2 * sin(theta) ** _sage_const_2 * cos(phi) ** _sage_const_2
        )
        / _sage_const_2
    )
    It_scatu1 = (
        It_scatu1.substitute(Esc=E / (_sage_const_1 + E * (_sage_const_1 - cos(theta))))
        .simplify_full()
        .simplify_trig()
    )
    It_scatu2 = (
        It_scatu2.substitute(Esc=E / (_sage_const_1 + E * (_sage_const_1 - cos(theta))))
        .simplify_full()
        .simplify_trig()
    )
    print bool(It_scatu1 == It_scatu2)

    print "\nScattered intensity (horizontal linear pol.):"
    print It_scat.substitute(S1=S0).substitute(S2=_sage_const_0).simplify_full()

    It_scat = It_scat.substitute(
        Esc=E / (_sage_const_1 + E * (_sage_const_1 - cos(theta)))
    )
    It_scat = It_scat.simplify_full()
    It_scat = definite_integral(It_scat * sin(theta), theta, _sage_const_0, pi)
    It_scat = definite_integral(It_scat, phi, _sage_const_0, _sage_const_2 * pi)
    print "\nKlein-Nishina cross-section (units of r_e^2):"
    print (It_scat / S0).simplify_full()


thomson()
compton()
