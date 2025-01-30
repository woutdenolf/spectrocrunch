# See http://maggotroot.blogspot.ch/2013/11/constrained-linear-least-squares-in.html for more info
"""
A simple library to solve constrained linear least squares problems
with sparse and dense matrices. Uses cvxopt library for
optimization
"""

__author__ = "Valeriy Vishnevskiy"
__email__ = "valera.vishnevskiy@yandex.ru"
__version__ = "1.0"
__date__ = "22.11.2013"
__license__ = "WTFPL"


import itertools
import numbers
from scipy import sparse
import numpy as np
from cvxopt import solvers, matrix, spmatrix, mul


def scipy_sparse_to_spmatrix(A):
    coo = A.tocoo()
    return spmatrix(coo.data, coo.row.tolist(), coo.col.tolist())


def spmatrix_sparse_to_scipy(A):
    data = np.array(A.V).squeeze()
    rows = np.array(A.I).squeeze()
    cols = np.array(A.J).squeeze()
    return sparse.coo_matrix((data, (rows, cols)))


def sparse_None_vstack(A1, A2):
    if A1 is None:
        return A2
    else:
        return sparse.vstack([A1, A2])


def numpy_None_vstack(A1, A2):
    if A1 is None:
        return A2
    else:
        return np.vstack([A1, A2])


def numpy_None_concatenate(A1, A2):
    if A1 is None:
        return A2
    else:
        return np.concatenate([A1, A2])


def get_shape(A):
    if isinstance(A, spmatrix):
        return A.size
    else:
        return A.shape


def numpy_to_cvxopt_matrix(A):
    if A is None:
        return A
    elif is_sparse(A):
        if isinstance(A, sparse.spmatrix):
            return scipy_sparse_to_spmatrix(A)
        else:
            return A
    else:
        if isinstance(A, np.ndarray):
            if A.ndim == 1:
                return matrix(A, (A.shape[0], 1), "d")
            else:
                return matrix(A, A.shape, "d")
        else:
            return A


def cvxopt_to_numpy_matrix(A):
    if A is None:
        return A
    elif isinstance(A, spmatrix):
        return spmatrix_sparse_to_scipy(A)
    elif isinstance(A, matrix):
        return np.array(A).squeeze()
    else:
        return np.array(A).squeeze()


def as_double(A):
    if isinstance(A, np.ndarray):
        return A.astype(np.double)
    elif isinstance(A, numbers.Number):
        return np.array([A], dtype=np.double)
    else:
        return A


def is_sparse(A):
    return sparse.issparse(A) or isinstance(A, spmatrix)


def lsqlin(
    C,
    d,
    reg=0,
    A=None,
    b=None,
    Aeq=None,
    beq=None,
    lb=None,
    ub=None,
    x0=None,
    opts=None,
):
    """
    Solve linear constrained l2-regularized least squares. Can
    handle both dense and sparse matrices. Matlab's lsqlin
    equivalent. It is actually wrapper around CVXOPT QP solver.

        min_x ||C*x  - d||^2_2 + reg * ||x||^2_2
        s.t.  A * x <= b
              Aeq * x = beq
              lb <= x <= ub

    Input arguments:
        C   is n x m dense or sparse matrix
        d   is n x 1 dense matrix
        reg is regularization parameter
        A   is p x n dense or sparse matrix
        b   is p x 1 dense matrix
        Aeq is q x n dense or sparse matrix
        beq is q x 1 dense matrix
        lb  is n x 1 matrix or scalar
        ub  is n x 1 matrix or scalar

    Output arguments:
        Return dictionary, the output of CVXOPT QP.
    """
    C = as_double(C)
    d = as_double(d)
    A = as_double(A)
    b = as_double(b)
    Aeq = as_double(Aeq)
    beq = as_double(beq)
    lb = as_double(lb)
    ub = as_double(ub)
    x0 = as_double(x0)

    sparse_case = False
    if is_sparse(A):
        sparse_case = True
        # We need A to be scipy sparse, as I couldn't find how
        # CVXOPT spmatrix can be vstacked
        if isinstance(A, spmatrix):
            A = spmatrix_sparse_to_scipy(A)

    C = numpy_to_cvxopt_matrix(C)
    d = numpy_to_cvxopt_matrix(d)
    Q = C.T * C
    q = -d.T * C
    nvars = C.size[1]

    if reg > 0:
        if sparse_case:
            I = scipy_sparse_to_spmatrix(sparse.eye(nvars, nvars, format="coo"))
        else:
            I = matrix(np.eye(nvars), (nvars, nvars), "d")
        Q = Q + reg * I

    lb = cvxopt_to_numpy_matrix(lb)
    ub = cvxopt_to_numpy_matrix(ub)
    b = cvxopt_to_numpy_matrix(b)

    if lb is not None:  # Modify 'A' and 'b' to add lb inequalities
        if lb.size == 1:
            lb = np.repeat(lb, nvars)

        if sparse_case:
            lb_A = -sparse.eye(nvars, nvars, format="coo")
            A = sparse_None_vstack(A, lb_A)
        else:
            lb_A = -np.eye(nvars)
            A = numpy_None_vstack(A, lb_A)
        b = numpy_None_concatenate(b, -lb)

    if ub is not None:  # Modify 'A' and 'b' to add ub inequalities
        if ub.size == 1:
            ub = np.repeat(ub, nvars)
        if sparse_case:
            ub_A = sparse.eye(nvars, nvars, format="coo")
            A = sparse_None_vstack(A, ub_A)
        else:
            ub_A = np.eye(nvars)
            A = numpy_None_vstack(A, ub_A)
        b = numpy_None_concatenate(b, ub)

    # Convert data to CVXOPT format
    A = numpy_to_cvxopt_matrix(A)
    Aeq = numpy_to_cvxopt_matrix(Aeq)
    b = numpy_to_cvxopt_matrix(b)
    beq = numpy_to_cvxopt_matrix(beq)

    # Set up options
    if opts is not None:
        for k, v in opts.items():
            solvers.options[k] = v

    # Run CVXOPT.SQP solver
    sol = solvers.qp(Q, q.T, G=A, h=b, A=Aeq, b=beq, init_vals=x0)
    return sol


def lsqnonneg(C, d):
    """
    Solves nonnegative linear least-squares problem:

    min_x ||C*x - d||_2^2,  where x >= 0
    """
    return lsqlin(C, d, lb=0, opts={"show_progress": False})


def lsqnorm(C, d):
    """
    Solves normalized linear least-squares problem:

    min_x ||C*x - d||_2^2,  where 0 <= x <= 1 and sum(x) == 1
    """
    Aeq = np.ones((1, get_shape(C)[1]))
    return lsqlin(C, d, Aeq=Aeq, beq=1, lb=0, ub=1, opts={"show_progress": False})
