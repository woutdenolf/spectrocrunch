"""Quadric surface in homogeneous coordinates x^T.A.x = 0

    A[0,0].x^2 + A[1,1].y^2 + A[2,2].z^2 + 
    (A[0,1]+A[1,0]).xy + (A[0,2]+A[2,0]).xz + (A[1,2]+A[2,1]).yz +
    (A[0,3]+A[3,0]).x + (A[1,3]+A[3,1]).y + (A[2,3]+A[3,2]).z + A[3,3] = 0
"""

import numpy as np


def plane(x0, u):
    """
    u.(x - x0) = 0

    Args:
        x0(array): point on plane
        u(array): plane normal (norm can be anything)
    Returns:
        array: n+1 x n+1
    """
    n1 = len(x0) + 1
    A = np.zeros((n1, n1))
    x0 = np.asarray(x0)
    u = np.asarray(u)
    u2 = u / 2
    A[-1, :-1] = u2
    A[:-1, -1] = u2
    A[-1, -1] = -x0.dot(u)
    return A


def ellipsoid(x0, a):
    """
    sum_i (yi-x0)^2/(ai*ai) = 1

    Args:
        x0(array): center
        a(array): semi-axes parallel to the basis
    Returns:
        array: n+1 x n+1
    """
    n = len(x0)
    ind = list(range(n))
    return _quadsum1(x0, a, n, ind)


def cylinder(x0, a, axis):
    """
    sum_i (yi-x0)^2/(ai*ai) = 1  with i != axis

    Args:
        x0(array): point on central axis
        a(array): semi-axes
        axis(int): central axis parallel to this axis
    Returns:
        array: n+1 x n+1
    """
    n = len(x0) + 1
    ind = list(range(n))
    ind.pop(axis)
    return _quadsum1(x0, a, n, ind)


def _quadsum1(x0, a, n, ind):
    """
    sum_i (yi-x0)^2/(ai*ai) = 1  with i in ind

    Args:
        x0(array):
        a(array):
        n(int): number of dimensions
        ind(list(int)): number of dimensions involved
    Returns:
        array: n+1 x n+1
    """
    x0 = np.asarray(x0)
    a = np.asarray(a)
    A = np.zeros((n + 1, n + 1))
    a2 = a * a
    p = x0 / a2
    A[-1, ind] = -p
    A[ind, -1] = -p
    A[ind, ind] = 1 / a2
    A[-1, -1] = p.dot(x0) - 1
    return A


def transform(A, C):
    """Quadric matrix under change-of-frame:

        x = C.x'
        x^T.A.x = x'^T.A'.x' = 0
        A' = C^T.A.C

    Args:
        A(array): quardic
        C(array): change-of-frame matrix
    """
    return C.T.dot(A.dot(C))
