# -*- coding: utf-8 -*-

from .types import transformationType
from ..math.fit1d import lstsq

import numpy as np
from scipy import stats
import scipy.ndimage.interpolation as ndtransform
import skimage.transform as sktransform


class Mapping(object):
    def transformimage(self, img):
        raise NotImplementedError()

    def fromkeypoints(self, xsrc, ysrc, xdest, ydest):
        raise NotImplementedError()

    def fromtransform(self, transformobj):
        raise NotImplementedError()


class LinearMapping(Mapping):
    """
    Transformation:
        L.Xold = Xnew (coordinate transformation)
        Xold = C.Xnew (change-of-frame)
    Combine transformations:
        L2.L1.Xold = Xnew (coordinate transformation, order = multiply from the left)
        Xold = C1.C2.Xnew (change-of-frame, order = multiply from the right)
    Change of reference frame:
        Xold = C1.Xnew
        Xold = C2.Xold'
        Xnew = C2.Xnew'
        Xold' = C2^(-1).C1.C2.Xnew'
        Xnew = L1.Xold
        Xnew' = C2^(-1).L1.C2.Xold'
    """

    def __init__(self, transfotype, dtype=float, cval=np.nan, **interpolationargs):
        self.transfotype = transfotype
        self.dtype = dtype
        self.cval = cval

        if "order" not in interpolationargs:
            interpolationargs["order"] = 1
        if "mode" not in interpolationargs:
            interpolationargs["mode"] = "constant"
        self.interpolationargs = interpolationargs

        # (change of frame matrices, not change of coordinates!)
        self.cof = self.getidentity()

    def __copy__(self):
        m = LinearMapping(
            self.transfotype, dtype=self.dtype, cval=self.cval, **self.interpolationargs
        )
        m.cof = self.cof
        return m

    def copy(self):
        return self.__copy__()

    def transformimage(self, img):
        # self.cof: change-of-frame matrix for coordinates (x,y)
        if self.isidentity():
            return img
        if self.transfotype == transformationType.translation:
            # shift: takes transformation vector for coordinates (y,x)
            return ndtransform.shift(
                img, -self.cof[1::-1], cval=self.cval, **self.interpolationargs
            )
        else:
            if self.islinearidentity():
                # shift: takes transformation vector for coordinates (y,x)
                return ndtransform.shift(
                    img, -self.cof[1::-1, 2], cval=self.cval, **self.interpolationargs
                )
            elif self.isprojidentity():
                # affine_transform: takes change-of-frame matrix
                #                   for coordinates (y,x)
                return ndtransform.affine_transform(
                    img,
                    self.cof[0:2, 0:2].T,
                    offset=self.cof[1::-1, 2],
                    cval=self.cval,
                    **self.interpolationargs
                )
            else:
                raise NotImplementedError()

    def __repr__(self):
        if self.transfotype == transformationType.translation:
            return "tx = {}, ty = {}".format(self.cof[0], self.cof[1])
        else:
            return "tx = {}, ty = {}\nR={}\npx = {}, py = {}".format(
                self.cof[0, 2],
                self.cof[1, 2],
                self.cof[0:2, 0:2],
                self.cof[2, 0],
                self.cof[2, 1],
            )

    def __str__(self):
        return self.__repr__()

    def getnumpy(self):
        return self.cof[:]

    def getnumpyhomography(self):
        if self.transfotype == transformationType.translation:
            return np.array(
                [[1, 0, self.cof[0]], [0, 1, self.cof[1]], [0, 0, 1]], dtype=self.dtype
            )
        else:
            return self.cof[:]

    def settranslation(self, trn, ty=None):
        """
        Args:
            trn(number or 2-array-like): x and optionally y translation
            ty(Optional(number)): y translation
        """
        if ty is not None:
            trn = [trn, ty]
        if self.transfotype == transformationType.translation:
            self.cof[:] = trn
        elif self.transfotype == transformationType.rigid:
            self.cof[0:2, 2] = trn
        elif self.transfotype == transformationType.similarity:
            self.cof[0:2, 2] = trn
        elif self.transfotype == transformationType.affine:
            self.cof[0:2, 2] = trn
        elif self.transfotype == transformationType.projective:
            self.cof[0:2, 2] = trn
        else:
            raise ValueError("Transformation does not have an translation part")

    def gettranslation(self):
        if self.transfotype == transformationType.translation:
            return self.cof[:]
        elif self.transfotype == transformationType.rigid:
            return self.cof[0:2, 2]
        elif self.transfotype == transformationType.similarity:
            return self.cof[0:2, 2]
        elif self.transfotype == transformationType.affine:
            return self.cof[0:2, 2]
        elif self.transfotype == transformationType.projective:
            return self.cof[0:2, 2]
        else:
            raise ValueError("Transformation does not have an translation part")

    def setrigid(self, theta, tx, ty):
        if self.transfotype == transformationType.translation:
            # trn = self.gettranslation()
            self.transfotype == transformationType.affine
            self.cof = self.getidentity()
        elif self.transfotype == transformationType.rigid:
            pass
        elif self.transfotype == transformationType.similarity:
            self.transfotype = transformationType.rigid
        elif self.transfotype == transformationType.affine:
            self.transfotype = transformationType.rigid
        elif self.transfotype == transformationType.projective:
            pass
        else:
            raise ValueError("Transformation does not have an linear part")
        costheta = np.cos(theta)
        sintheta = np.sin(theta)
        self.cof[0:2, :] = np.array(
            [[costheta, -sintheta, tx], [sintheta, costheta, ty]], self.dtype
        )

    def getrigid(self):
        if self.transfotype == transformationType.translation:
            return 0, self.cof[0], self.cof[1]
        elif self.transfotype == transformationType.rigid:
            return np.arccos(self.cof[0, 0]), self.cof[0, 2], self.cof[1, 2]
        else:
            raise ValueError(
                "Transformation should be a translation or a rigid transformation"
            )

    def setlinear(self, M):
        if self.transfotype == transformationType.translation:
            trn = self.gettranslation()
            self.transfotype = transformationType.affine
            self.cof = self.getidentity()
            self.settranslation(trn)
        elif self.transfotype == transformationType.rigid:
            pass
        elif self.transfotype == transformationType.similarity:
            pass
        elif self.transfotype == transformationType.affine:
            pass
        elif self.transfotype == transformationType.projective:
            pass
        else:
            raise ValueError("Transformation does not have an linear part")
        self.cof[0:2, 0:2] = M

    def getlinear(self):
        if self.transfotype == transformationType.translation:
            return np.identity(2, dtype=self.dtype)
        elif self.transfotype == transformationType.rigid:
            return self.cof[0:2, 0:2]
        elif self.transfotype == transformationType.similarity:
            return self.cof[0:2, 0:2]
        elif self.transfotype == transformationType.affine:
            return self.cof[0:2, 0:2]
        elif self.transfotype == transformationType.projective:
            return self.cof[0:2, 0:2]
        else:
            raise ValueError("Transformation does not have an linear part")

    def setaffine(self, M):
        if self.transfotype == transformationType.translation:
            trn = self.gettranslation()
            self.transfotype = transformationType.affine
            self.cof = self.getidentity()
            self.settranslation(trn)
        elif self.transfotype == transformationType.rigid:
            pass
        elif self.transfotype == transformationType.similarity:
            pass
        elif self.transfotype == transformationType.affine:
            pass
        elif self.transfotype == transformationType.projective:
            pass
        else:
            raise ValueError("Transformation does not have an affine part")
        self.cof[0:2, :] = M

    def getaffine(self):
        if self.transfotype == transformationType.translation:
            return np.array(
                [[1, 0, self.cof[0]], [0, 1, self.cof[1]]], dtype=self.dtype
            )
        elif self.transfotype == transformationType.rigid:
            return self.cof[:, 0:2]
        elif self.transfotype == transformationType.similarity:
            return self.cof[:, 0:2]
        elif self.transfotype == transformationType.affine:
            return self.cof[:, 0:2]
        elif self.transfotype == transformationType.projective:
            return self.cof[:, 0:2]
        else:
            raise ValueError("Transformation does not have an affine part")

    def setprojective(self, M):
        if self.transfotype == transformationType.translation:
            trn = self.gettranslation()
            self.transfotype = transformationType.projective
            self.cof = self.getidentity()
            self.settranslation(trn)
        elif self.transfotype == transformationType.rigid:
            self.transfotype == transformationType.projective
        elif self.transfotype == transformationType.similarity:
            self.transfotype == transformationType.projective
        elif self.transfotype == transformationType.affine:
            self.transfotype == transformationType.projective
        elif self.transfotype == transformationType.projective:
            pass
        else:
            raise ValueError("Transformation does not have an affine part")
        self.cof[:] = M

    def getprojective(self):
        if self.transfotype == transformationType.translation:
            return np.array(
                [[1, 0, self.cof[0]], [0, 1, self.cof[1]], [0, 0, 1]], dtype=self.dtype
            )
        elif self.transfotype == transformationType.rigid:
            return self.cof[:]
        elif self.transfotype == transformationType.similarity:
            return self.cof[:]
        elif self.transfotype == transformationType.affine:
            return self.cof[:]
        elif self.transfotype == transformationType.projective:
            return self.cof[:]
        else:
            raise ValueError("Transformation does not have an affine part")

    def setidentity(self):
        if self.transfotype == transformationType.translation:
            self.cof[:] = np.zeros(2, dtype=self.dtype)
        else:
            self.cof[:] = np.identity(3, dtype=self.dtype)

    def getidentity(self):
        if self.transfotype == transformationType.translation:
            return np.zeros(2, dtype=self.dtype)
        else:
            return np.identity(3, dtype=self.dtype)

    def isidentity(self):
        if self.transfotype == transformationType.translation:
            return np.array_equal(self.cof, np.zeros(2, dtype=self.dtype))
        else:
            return np.array_equal(self.cof, np.identity(3, dtype=self.dtype))

    def islinearidentity(self):
        if self.transfotype == transformationType.translation:
            return True
        else:
            return np.array_equal(self.cof[0:2, 0:2], np.identity(2, dtype=self.dtype))

    def isprojidentity(self):
        if self.transfotype == transformationType.translation:
            return True
        else:
            return np.array_equal(self.cof[2, 0:2], np.zeros(2, dtype=self.dtype))

    def _combine(self, C, right=True):
        """
        :param LinearMapping C:
        :param bool right: multiply self from the right (self is the first transformation)
        """
        if (
            self.transfotype == transformationType.translation
            and C.transfotype == transformationType.translation
        ):
            cof = self.cof + C.cof
            return cof, transformationType.translation

        if self.transfotype == transformationType.translation:
            C1 = np.identity(3, dtype=self.dtype)
            C1[0:2, 2] = self.cof
        else:
            C1 = self.cof

        if C.transfotype == transformationType.translation:
            C2 = np.identity(3, dtype=self.dtype)
            C2[0:2, 2] = C.cof
        else:
            C2 = C.cof

        if right:
            cof = np.dot(C1, C2)
        else:
            cof = np.dot(C2, C1)

        if (
            self.transfotype == transformationType.projective
            or C.transfotype == transformationType.projective
        ):
            transfotype = transformationType.projective
        elif (
            self.transfotype == transformationType.affine
            or C.transfotype == transformationType.affine
        ):
            transfotype = transformationType.affine
        elif (
            self.transfotype == transformationType.similarity
            or C.transfotype == transformationType.similarity
        ):
            transfotype = transformationType.similarity
        else:
            transfotype = transformationType.rigid

        return cof, transfotype

    def _dot(self, C, right=True):
        cof, transfotype = self._combine(C, right=right)
        ret = transform(transfotype, dtype=cof.dtype, cval=self.cval)
        ret.cof[:] = cof
        return ret

    def _dotinplace(self, C, right=True):
        cof, transfotype = self._combine(C, right=right)
        if len(cof) == len(self.cof):
            self.cof[:] = cof
        else:
            self.cof = cof
        self.transfotype = transfotype
        self.dtype = cof.dtype

    def after(self, C):
        # C is the first transformation
        return self._dot(C, right=False)

    def before(self, C):
        # self is the first transformation
        return self._dot(C, right=True)

    def after_inplace(self, C):
        return self._dotinplace(C, right=False)

    def before_inplace(self, C):
        return self._dotinplace(C, right=True)

    def new_frame(self, C):
        return C.inverse().before(self).before(C)

    def new_frame_inplace(self, C):
        self.after_inplace(C.inverse())
        self.before_inplace(C)

    def inverse(self):
        ret = transform(self.transfotype, dtype=self.cof.dtype, cval=self.cval)
        if self.transfotype == transformationType.translation:
            ret.cof[:] = -self.cof
        else:
            ret.cof[:] = np.linalg.inv(self.cof)
        return ret

    def inverse_inplace(self):
        if self.transfotype == transformationType.translation:
            self.cof[:] = -self.cof
        else:
            self.cof[:] = np.linalg.inv(self.cof)

    def transformcoordinates(self, xy):
        # Anew = C^-1.Aold
        if self.transfotype == transformationType.translation:
            Ci = -self.cof
        else:
            Ci = np.linalg.inv(self.cof)
        return self._transformcoordinates(Ci, xy)

    def transformcoordinatesinverse(self, xy):
        # Aold = C.Anew
        return self._transformcoordinates(self.cof, xy)

    def _transformcoordinates(self, C, xy):
        # xy' = C.xy
        if self.transfotype == transformationType.translation:
            if len(xy.shape) == 1:
                if xy.shape[0] == 2:
                    return xy + C
                elif xy.shape[0] == 3:
                    return np.append(xy[0:2] + C, xy[2])
                else:
                    raise ValueError("Shape of coordinates cannot be handled")
            else:
                if xy.shape[0] == 2:
                    return xy + C[:, np.newaxis]
                elif xy.shape[0] == 3:
                    return np.vstack((xy[0:2, :] + C[:, np.newaxis], xy[2, :]))
                else:
                    raise ValueError("Shape of coordinates cannot be handled")
        else:
            return np.dot(C, xy)

    def fromtransform(self, transformobj):
        self.transfotype = transformobj.transfotype
        self.dtype = transformobj.dtype
        self.cval = transformobj.cval
        self.cof[:] = transformobj.cof[:]

    def _solvelinearsystem(self, A, b):
        sol = lstsq(A, b)
        if sol is not None:
            sol = sol.flatten()
        return sol

        ##### Using pseudo-inverse #####
        # A-1* = (A^T.A)^(-1).A^T
        # try:
        #    S = np.dot(A.T, A)
        #    sol = np.dot(np.linalg.inv(S), np.dot(A.T, b))
        # except np.linalg.LinAlgError as err:
        #    logger = logging.getLogger(__name__)
        #    logger.error("Singular matrix in calculating a transformation from SIFT keypoints")
        #    sol = None

        # sol = np.dot(numpy.linalg.pinv(A),b) #slower?

        ##### Using SVD #####
        # sol = np.linalg.lstsq(A,b)[0] # computing the numpy solution

        ##### Using QR #####
        # Q,R = np.linalg.qr(A) # qr decomposition of A
        # Qb = np.dot(Q.T,b) # computing Q^T*b (project b onto the range of A)
        # sol = np.linalg.solve(R,Qb) # solving R*x = Q^T*b

        # result
        # MSE = np.linalg.norm(b - np.dot(A,sol))**2/N #Mean Squared Error

        # return sol

    def centroid(self, x):
        # return np.mean(x)
        # return np.median(x)
        return stats.trim_mean(x, 0.1)  # trimmed mean (trim 10% at both sides)

    def translation_fromkeypoints(self, xsrc, ysrc, xdest, ydest):
        # REMARK: we need the change-of-frame matrix, so we need the map dest -> src
        self.setidentity()
        self.settranslation([self.centroid(xsrc - xdest), self.centroid(ysrc - ydest)])

    def rigid_fromkeypoints(self, xsrc, ysrc, xdest, ydest):
        # REMARK: we need the change-of-frame matrix, so we need the map dest -> src
        self.setidentity()

        # https://igl.ethz.ch/projects/ARAP/svd_rot.pdf
        #
        # R.(X-Xcen) = Y-Ycen
        # R.X + T = Y   and   T = Ycen - R.Xcen

        censrc = np.asarray([self.centroid(xdest), self.centroid(ydest)])[:, np.newaxis]
        cendest = np.asarray([self.centroid(xsrc), self.centroid(ysrc)])[:, np.newaxis]

        X = np.vstack([xdest, ydest])
        Y = np.vstack([xsrc, ysrc])

        S = np.dot(
            np.dot(X - censrc, np.identity(len(xdest))), np.transpose(Y - cendest)
        )
        U, s, V = np.linalg.svd(S, full_matrices=False)
        C = np.dot(V, np.transpose(U))
        self.setlinear(C)

        # self.settranslation(cendest - C.dot(censrc))
        # This seems to be more accurate:
        Ynoshift = Y - C.dot(X)
        self.settranslation(
            [self.centroid(Ynoshift[0, :]), self.centroid(Ynoshift[1, :])]
        )

    def similarity_fromkeypoints(self, xsrc, ysrc, xdest, ydest):
        # REMARK: we need the change-of-frame matrix, so we need the map dest -> src
        self.setidentity()

        # Similarity transformation:
        #    x' = a.x - b.y + t0
        #    y' = b.x + a.y + t1
        #    sol = [a,b,t0,t1]
        #
        # xdest = [x1,x2,...]
        # ydest = [y1,y2,...]
        # xsrc = [x1',x2',...]
        # ysrc = [y1',y2',...]
        #
        # X = x1 -y1  1  0
        #     y1  x1  0  1
        #     x2 -y2  1  0
        #     y2  x2  0  1
        #     ...
        #
        # Y = x1'
        #     y1'
        #     x2'
        #     y2'
        #     ...

        N = len(xdest)

        X = np.zeros((2 * N, 4))
        X[::2, 0] = xdest
        X[1::2, 0] = ydest
        X[::2, 1] = -ydest
        X[1::2, 1] = xdest
        X[::2, 2] = 1
        X[1::2, 3] = 1

        Y = np.zeros((2 * N, 1))
        Y[::2, 0] = xsrc
        Y[1::2, 0] = ysrc

        sol = self._solvelinearsystem(X, Y)
        if sol is not None:
            self.setlinear([[sol[0], -sol[1]], [sol[1], sol[0]]])
            self.settranslation(sol[2:].flatten())

    def affine_fromkeypoints(self, xsrc, ysrc, xdest, ydest):
        # REMARK: we need the change-of-frame matrix, so we need the map dest -> src
        self.setidentity()

        # Affine transformation:
        #    x' = a.x + b.y + t0
        #    y' = c.x + d.y + t1
        #    sol = [a,b,t0,c,d,t1]
        #
        # xdest = [x1,x2,...]
        # ydest = [y1,y2,...]
        # xsrc = [x1',x2',...]
        # ysrc = [y1',y2',...]
        #
        # X = x1 y1  1  0  0  0
        #      0  0  0 x1 y1  1
        #     x2 y2  1  0  0  0
        #      0  0  0 x2 y2  1
        #     ...
        #
        # Y = x1'
        #     y1'
        #     x2'
        #     y2'
        #     ...

        # raise NotImplementedError("Sift doesn't support affine transformations (keypoints not invariant).")

        N = len(xdest)

        X = np.zeros((2 * N, 6))
        X[::2, 0] = xdest
        X[::2, 1] = ydest
        X[::2, 2] = 1
        X[1::2, 3] = xdest
        X[1::2, 4] = ydest
        X[1::2, 5] = 1

        Y = np.zeros((2 * N, 1))
        Y[::2, 0] = xsrc
        Y[1::2, 0] = ysrc

        sol = self._solvelinearsystem(X, Y)
        if sol is not None:
            self.setaffine(sol.reshape((2, 3)))

    def projective_fromkeypoints(self, xsrc, ysrc, xdest, ydest):
        # REMARK: we need the change-of-frame matrix, so we need the map dest -> src
        self.setidentity()

        # Projective transformation:
        #    x' = (a.x + b.y + t0)/(px.x+py.y+1)
        #    y' = (c.x + d.y + t1)/(px.x+py.y+1)
        #     x' = a.x + b.y + t0 - px.x.x' - py.y.x'
        #     y' = c.x + d.y + t1 - px.x.y' - py.y.y'
        #    sol = [a,b,t0,c,d,t1,px,py]
        #
        # xdest = [x1,x2,...]
        # ydest = [y1,y2,...]
        # xsrc = [x1',x2',...]
        # ysrc = [y1',y2',...]
        #
        # X = x1 y1  1  0  0  0 -x1.x1' -y1.x1' -x1'
        #      0  0  0 x1 y1  1 -x1.y1' -y1.y1' -y1'
        #     x2 y2  1  0  0  0 -x2.x2' -y2.x2' -x1'
        #      0  0  0 x2 y2  1 -x2.y2' -y2.y2' -y1'
        #     ...
        #
        # Y = 0
        #     0
        #     0
        #     0
        #     ...

        raise NotImplementedError(
            "Sift doesn't support homographies (keypoints not invariant)."
        )
        # if self.usekernel:
        #    raise NotImplementedError("Sift doesn't support this type of transformation.")

        N = len(xdest)

        X = np.zeros((2 * N, 8))
        X[::2, 0] = xdest
        X[::2, 1] = ydest
        X[::2, 2] = 1
        X[1::2, 3] = xdest
        X[1::2, 4] = ydest
        X[1::2, 5] = 1
        X[::2, 6] = -xdest * xsrc
        X[1::2, 6] = -xdest * ysrc
        X[::2, 7] = -ydest * xsrc
        X[1::2, 7] = -ydest * ysrc

        Y = np.zeros((2 * N, 1))
        Y[::2, 0] = xsrc
        Y[1::2, 0] = ysrc

        sol = self._solvelinearsystem(X, Y)
        if sol is not None:
            self.setprojective(np.append(sol, 1).reshape((3, 3)))

    def fromkeypoints(self, xsrc, ysrc, xdest, ydest):
        """Least-squares transformation parameters to map src to dest

        Remark: the rigid transformation is the most problematic (cfr. test_transform)
        """

        if self.transfotype == transformationType.translation:
            self.translation_fromkeypoints(xsrc, ysrc, xdest, ydest)
        elif self.transfotype == transformationType.rigid:
            self.rigid_fromkeypoints(xsrc, ysrc, xdest, ydest)
        elif self.transfotype == transformationType.similarity:
            self.similarity_fromkeypoints(xsrc, ysrc, xdest, ydest)
        elif self.transfotype == transformationType.affine:
            self.affine_fromkeypoints(xsrc, ysrc, xdest, ydest)
        elif self.transfotype == transformationType.projective:
            self.projective_fromkeypoints(xsrc, ysrc, xdest, ydest)
        else:
            raise NotImplementedError(
                "Sift doesn't support this type of transformation."
            )

        # xsrc2,ysrc2,norm = self.transformcoordinatesinverse(np.vstack([xdest,ydest,ydest*0+1]))
        # print("Maximum keypoint residual: [{}, {}]".format(np.max(np.abs(xsrc-xsrc2/norm)),np.max(np.abs(ysrc-ysrc2/norm))))


def transform(transfotype, **kwargs):
    if transfotype in [
        transformationType.translation,
        transformationType.rigid,
        transformationType.similarity,
        transformationType.affine,
        transformationType.projective,
    ]:
        return LinearMapping(transfotype, **kwargs)
