import numpy as np
from matplotlib import pyplot as plt
from shapely.geometry import Polygon
from shapely.geometry import LineString
import logging

from .alignSource import alignSource
from .alignDest import alignDest
from .types import alignType
from .types import transformationType
from .transform import transform
from ..utils.roi import cliproi
from ..utils import instance


logger = logging.getLogger(__name__)


class align(object):
    """
    Allows for the alignment of several stacks based on one stack.
    "Alignment" is the process of determining a transformation between
    two images that should represent the same thing but transformed/deformed.
    """

    def __init__(
        self,
        source,
        sourcelist,
        dest,
        destlist,
        extension,
        stackdim=None,
        overwrite=False,
        cval=np.nan,
        plot=False,
        transfotype=transformationType.translation,
    ):
        """"""
        # Data IO
        self.source = alignSource(source, sourcelist, stackdim=stackdim)
        self.dest = alignDest(
            dest,
            destlist,
            extension,
            stackdim=self.source.stackdim,
            overwrite=overwrite,
        )

        # Missing data
        self.cval = cval

        # Transformation settings (set before actual transformation)
        try:
            one = self.source.dtype.type(1)
        except AttributeError:
            one = self.source.dtype(1)
        self.dtype = (np.float32(1) * one).dtype.type
        self.alignonraw = True
        self.usekernel = False  # Doesn't work well for Elastix!
        self.pre_align = {"roi": None, "func": None}
        self.pre_transform = {"pad": False}
        self.post_transform = {"crop": False}
        self.pre_transform_requested = {"pad": False}
        self.post_transform_requested = {"crop": False}
        self.extend = ((0, 0), (0, 0))  # negative: crop, positive: pad

        # Transformation (change of frame matrices, not change of coordinates!)
        self.transfotype = transfotype
        self.transfos = [self.defaulttransform() for _ in range(self.source.nimages)]
        self.pre_transfos = [
            self.defaulttransform() for _ in range(self.source.nimages)
        ]
        self.prealign_to_raw = [
            self.defaulttransform() for _ in range(self.source.nimages)
        ]
        self.prealign_to_pretransform = [
            self.defaulttransform() for _ in range(self.source.nimages)
        ]

        # Plot
        self.plotinfo = {"ON": plot, "fig": None, "axes": None}

    def defaulttransform(self, ttype=None, dtype=None):
        if ttype is None:
            ttype = self.transfotype
        if dtype is None:
            dtype = self.dtype
        return transform(ttype, dtype=dtype, cval=self.cval)

    def enableplot(self):
        self.plotinfo["ON"] = True

    def disableplot(self):
        self.plotinfo["ON"] = False

    def plot(self, img, index, title):
        """Visualize alignment in progress"""
        if not self.plotinfo["ON"]:
            return

        if self.plotinfo["fig"] is None:
            self.plotinfo["fig"], self.plotinfo["axes"] = plt.subplots(1, 3)
        ax = self.plotinfo["axes"][index]
        ax.cla()

        if img.size in img.shape:
            ax.plot(img.flatten())
        else:
            # img2 = img.copy()
            # img2[np.isnan(img2)] = 0
            ax.imshow(img, origin="lower", interpolation="nearest", cmap="jet")
        ax.set_title(title)

        plt.pause(0.01)

    def padfromextend(self):
        return (
            (max(self.extend[0][0], 0), max(self.extend[0][1], 0)),
            (max(self.extend[1][0], 0), max(self.extend[1][1], 0)),
        )

    def cropfromextend(self, dim1, dim2):
        return (
            (max(-self.extend[0][0], 0), dim1 - max(-self.extend[0][1], 0)),
            (max(-self.extend[1][0], 0), dim2 - max(-self.extend[1][1], 0)),
        )

    def pad(self, img):
        """Apply padding"""
        pad = self.padfromextend()
        if np.count_nonzero(pad) != 0:
            return np.pad(img, pad, "constant", constant_values=(self.cval, self.cval))
        else:
            return img

    def crop(self, img):
        """Apply cropping"""
        dim1, dim2 = img.shape
        crop = self.cropfromextend(dim1, dim2)
        if (
            crop[0][0] != 0
            or crop[1][0] != 0
            or crop[0][1] != dim1
            or crop[1][1] != dim2
        ):
            return img[crop[0][0] : crop[0][1], crop[1][0] : crop[1][1]]
        else:
            return img

    def roi(self, img, roi):
        """Extract ROI"""
        [[ya, yb], [xa, xb]] = cliproi(img.shape, roi)
        if xb <= xa or yb <= ya:
            raise ValueError(
                "ROI reduces image size to zero: [{}:{},{}:{}]".format(ya, yb, xa, xb)
            )
        return img[ya:yb, xa:xb]

    def writeimg(self, img, datasetindex, imageindex):
        """Save 1 image in 1 stack."""
        self.dest.writeimg(img, datasetindex, imageindex)

    def copyimg(self, datasetindex, imageindex):
        """Copy 1 image in 1 stack."""
        img = self.readimgraw(datasetindex, imageindex)
        self.writeimg(img, datasetindex, imageindex)

    def readimgraw(self, datasetindex, imageindex):
        """Read 1 image in 1 stack."""
        return self.source.readimgas(datasetindex, imageindex, self.dtype)

    def readimgrawprep(self, datasetindex, imageindex):
        """Get raw image, preprocessed for alignment"""
        img = self.readimgraw(datasetindex, imageindex)
        img = self.dopre_align(img, imageindex)
        if 0 in img.shape or len(img.shape) != 2:
            raise ValueError(
                "Image preprocessed for alignment has shape {}".format(img.shape)
            )
        return img

    def nopre_align(self):
        """
        Returns:
            bool: True when alignment is done on the raw image
        """
        return self.pre_align["roi"] is None

    def dopre_align(self, img, i):
        if callable(self.pre_align["func"]):
            img = self.pre_align["func"](img)
        transfo = self.pre_transfos[i]
        if not transfo.isidentity():
            img = self.execute_transform(img, i, transfo)
        if self.pre_align["roi"] is not None:
            img = self.roi(img, self.pre_align["roi"])
        return img

    def nopre_transform(self, i):
        """
        Returns:
            bool: True when transformation is done on the raw image
        """
        return (
            np.all(np.asarray(self.extend) <= 0) and self.pre_transfos[i].isidentity()
        )

    def dopre_transform(self, img, i):
        """Manual transformation before the real transformation (not used in alignment)"""
        transfo = self.pre_transfos[i]
        if not transfo.isidentity():
            img = self.execute_transform(img, i, transfo)
        if self.pre_transform["pad"]:
            img = self.pad(img)
        return img

    def nopost_transform(self):
        """
        Returns:
            bool: True when transformation doesn't have any post processing
        """
        return not self.post_transform["crop"]

    def dopost_transform(self, img):
        """Manual transformation after the real transformation (not used in alignment)"""
        if self.post_transform["crop"]:
            img = self.crop(img)
        return img

    def execute_alignkernel(self, img):
        raise NotImplementedError()

    def execute_transformkernel(self, img):
        raise NotImplementedError()

    def absolute_cofs(self, homography=False, include_pre=False):
        """Change-of-frame matrix (i.e. inverse of coordinate transformation matrix)
        to convert each raw image to the its aligned version.

        The columns of the COF matrix are the coordinates of the new basis vectors with
        respect to the old basis vectors. The last column is the origin of the new frame
        with respect to the old reference frame.
        """
        if self.source.nimages == 0:
            return None
        if include_pre:
            transfos = [
                ptr.before(tr) for ptr, tr in zip(self.pre_transfos, self.transfos)
            ]
        else:
            transfos = self.transfos
        if homography:
            return np.array([t.getnumpyhomography() for t in transfos])
        else:
            return np.array([t.getnumpy() for t in transfos])

    def cof_in_raw_frame(self, i):
        C = self.prealign_to_raw[i]
        C2 = self.transfos[i]  # defined in pre-align frame
        if C.isidentity():
            return C2
        else:
            return C2.new_frame(C)

    def cof_in_pretransform_frame(self, i):
        C = self.prealign_to_pretransform[i]
        C2 = self.transfos[i]  # defined in pre-align frame
        if C.isidentity():
            return C2
        else:
            return C2.new_frame(C)

    def calccof_prealign_to_raw(self):
        """
        Calculate transformation (image alignment):
            fixed image:   raw -> dopre_align (cof=C1) -> img1
            moving image:  raw -> dopre_align (cof=C1) -> img2 -> execute_alignkernel(img1,img2) -> img3, C2

            C1: raw to pre-align (in order: pretransform, roi)
            C2: cof in pre-align frame
        """
        CroiInv = self.defaulttransform(ttype=transformationType.translation)
        if self.pre_align["roi"] is not None:
            # cof raw to pre-align
            ty, tx = self.pre_align["roi"][:][0]
            # cof pre-align to raw
            CroiInv.settranslation(-tx, -ty)
        for prealign_to_raw, prealign_to_pretransform, pretransfo in zip(
            self.prealign_to_raw, self.prealign_to_pretransform, self.pre_transfos
        ):
            prealign_to_raw.fromtransform(CroiInv)
            if not pretransfo.isidentity():
                # raw to pre-transform
                prealign_to_raw.before_inplace(pretransfo.inverse())
            prealign_to_pretransform.fromtransform(prealign_to_raw)

    def calccof_prealign_to_pretransform(self):
        """
        Apply transformation (image transformation):
            aligned image: raw -> dopre_transform (cof=C3) -> img1 -> apply_transform (cof=C2) -> img4 -> dopost_transform (cof=C4) -> img5

            C1: raw to pre-align
            C2: cof in pre-align frame
            C3: raw to pre-transform (in order: pretransform, pad)
            C4: post transformation (crop)
        """
        Cpad = self.defaulttransform(ttype=transformationType.translation)
        if self.pre_transform["pad"]:
            # raw to pre-transform
            o2min = -max(self.extend[1][0], 0)
            o1min = -max(self.extend[0][0], 0)
            Cpad.settranslation(o2min, o1min)
        for prealign_to_raw, prealign_to_pretransform, pretransfo in zip(
            self.prealign_to_raw, self.prealign_to_pretransform, self.pre_transfos
        ):
            prealign_to_pretransform.fromtransform(prealign_to_raw)
            if not pretransfo.isidentity():
                # raw to pre-transform
                prealign_to_pretransform.before(pretransfo)
            if self.pre_transform["pad"]:
                prealign_to_pretransform.before(Cpad)

    def execute_transform(self, img, i, transfo=None):
        """
        Transform according to the transformation extracted from
        the transformation kernel (see store_transformation).

        :param np.ndarray img: in the pre-transform frame
        :param int or LinearMapping transfo:
        :returns np.ndarray:
        """
        if transfo is None:
            transfo = self.transfos[i]  # defined in pre-align frame
        C = self.prealign_to_pretransform[i]
        if not C.isidentity():
            transfo.new_frame_inplace(C)  # now defined in pre-transform frame
        if not transfo.isidentity():
            if self.usekernel:
                self.set_transformkernel(transfo)
                return self.execute_transformkernel(img)
            else:
                return transfo.transformimage(img)
        return img

    def transformidentity(self, transfo):
        """Is the transformation the identity"""
        if instance.isnumber(transfo):
            transfo = self.transfos[transfo]
        return transfo.isidentity()

    def pureidentity(self, i):
        """Is the transformation the identity, including the changes applied
        before (padding) and after (cropping)
        """
        return (
            self.nopre_transform(i)
            and self.nopost_transform()
            and self.transformidentity(i)
        )

    def transform(self, img, i):
        """Apply image transformation"""
        # Return when transformation is the identity
        if self.pureidentity(i):
            return img

        # Apply initial transformation (not used in alignment)
        imgtransformed = self.dopre_transform(img, i)

        # Apply transformation
        imgtransformed = self.execute_transform(imgtransformed, i)

        # plt.figure(5)
        # plt.imshow(imgtransformed,origin='lower',interpolation='nearest',cmap='jet')
        # plt.pause(1)

        # Apply final transformation (not used in alignment)
        imgtransformed = self.dopost_transform(imgtransformed)

        return imgtransformed

    def get_alignkernel(self):
        """Get transformation from align kernel."""
        raise NotImplementedError()

    def set_transformkernel(self, transfo):
        """Set transformation in transform kernel"""
        raise NotImplementedError()

    def store_transformation(self, i, pairwise):
        """
        Store transformation from align kernel for usage
        in transform kernel
        """
        transfo = self.get_alignkernel()
        # pairwise: transfo relative to previous
        # not pairwise: transfo relative to i=iref
        if pairwise and i != 0 and self.alignonraw:
            # make transfo relative to i=0
            self.transfos[i].fromtransform(self.transfos[i - 1].before(transfo))
        else:
            self.transfos[i].fromtransform(transfo)

    def settransformidentity(self, i):
        """Make this transformation the identity"""
        self.transfos[i].setidentity()

    def genpolygon(self, lst):
        p = Polygon(lst)
        if p.area == 0:
            p = LineString(lst)
        return p

    def polygonempty(self, p):
        if isinstance(p, Polygon):
            return p.area == 0
        else:
            return p.length == 0

    def untransformedimagepolygon(self):
        add0 = 0.0
        add1 = 0.0
        return self.genpolygon(
            [
                (add0, add0),
                (self.source.imgsize[1] - 1 + add1, add0),
                (self.source.imgsize[1] - 1 + add1, self.source.imgsize[0] - 1 + add1),
                (add0, self.source.imgsize[0] - 1 + add1),
            ]
        )

    def transformedimagepolygons(self):
        add0 = 0.0
        add1 = 0.0

        # Corners of the image in the transformed frame: A'
        xy = np.empty((3, 4))
        xy[0, :] = [
            add0,
            self.source.imgsize[1] - 1 + add1,
            self.source.imgsize[1] - 1 + add1,
            add0,
        ]  # x
        xy[1, :] = [
            add0,
            add0,
            self.source.imgsize[0] - 1 + add1,
            self.source.imgsize[0] - 1 + add1,
        ]  # y
        xy[2, :] = [1, 1, 1, 1]

        # Corners of the image in the raw frame: A = C.A'
        ret = [None] * self.source.nimages
        for i in range(self.source.nimages):
            xy2 = self.cof_in_raw_frame(i).transformcoordinates(
                xy
            )  # C1^(-1).C2^(-1).C1.XY
            xy2[0, :] /= xy2[2, :]
            xy2[1, :] /= xy2[2, :]
            ret[i] = self.genpolygon(xy2[0:2, :].T)

        return ret

    def polygoncollectionbounds(self, ps, pad=True):
        p = ps[0]
        if pad:
            for i in range(1, len(ps)):
                p = p.union(ps[i])
        else:
            for i in range(1, len(ps)):
                p = p.intersection(ps[i])

            if self.polygonempty(p):
                logger.warning("Cropping skipped because there would be nothing left.")
                return ()

        xmin, ymin, xmax, ymax = p.bounds
        xmin = int(np.floor(xmin))
        ymin = int(np.floor(ymin))
        xmax = int(np.ceil(xmax))
        ymax = int(np.ceil(ymax))

        return xmin, ymin, xmax, ymax

    def minimaltransformation(self, p0, ps, centroids=False):
        """If all transformations are known, they can be reduced to minimize
        the difference with the original image
        """
        # return

        if centroids:
            # Put average centroids in the middle
            x0, y0 = p0.centroid.coords.xy
            xy = np.array([p.centroid.coords.xy for p in ps])
            x0 = x0[0]
            y0 = y0[0]
            x = np.mean(xy[:, 0])
            y = np.mean(xy[:, 1])
        else:
            # Center total boundary
            xmin0, ymin0, xmax0, ymax0 = p0.bounds
            xmin, ymin, xmax, ymax = self.polygoncollectionbounds(ps)
            x0 = np.mean([xmin0, xmax0])
            y0 = np.mean([ymin0, ymax0])
            x = np.mean([xmin, xmax])
            y = np.mean([ymin, ymax])

        # Center
        trn = self.defaulttransform(ttype=transformationType.translation)
        trn.settranslation(x - x0, y - y0)
        for t in self.transfos:
            t.before_inplace(trn)

    def setextend(self, xmin, ymin, xmax, ymax):
        # Padding/cropping <> positive/negative
        o1min = -ymin
        o1max = ymax - self.source.imgsize[0] + 1
        o2min = -xmin
        o2max = xmax - self.source.imgsize[1] + 1
        self.extend = ((o1min, o1max), (o2min, o2max))
        self.pre_transform["pad"] = np.any(np.asarray(self.extend) > 0)
        self.post_transform["crop"] = np.any(np.asarray(self.extend) < 0)

    def extendfromtransformation(self, p0, ps):
        """If all transformations are known, padding/cropping can be calculated"""
        self.extend = ((0, 0), (0, 0))
        # Smallest rectangle that contains the union (pad)
        # or intersection (crop) of all polygons
        tmp = self.polygoncollectionbounds(ps, pad=self.pre_transform_requested["pad"])
        if len(tmp) != 4:
            self.pre_transform["pad"] = False
            self.post_transform["crop"] = False
            return
        self.setextend(*tmp)

    def bextendfrommask(self, pairwise):
        return (
            self.post_transform_requested["crop"]
            and self.cval is np.nan
            and (
                self.pre_align["roi"] is None
                or self.transfotype == transformationType.translation
            )
            and not pairwise
        )

    def setextendmask(self, img, reset=False):
        mask = np.logical_not(np.isnan(img))
        # if self.cval is np.nan:
        #    mask = np.logical_not(np.isnan(img))
        # else:
        #    mask = img != self.cval
        if reset:
            self._extendmask = mask
        else:
            self._extendmask &= mask

    def extendfrommask(self):
        """If all transformations are applied, padding/cropping can be calculated"""
        indvalidrow = np.argwhere(self._extendmask.sum(axis=1))
        indvalidcol = np.argwhere(self._extendmask.sum(axis=0))
        # When pre_align["roi"]: only valid for translations
        ymin = indvalidrow[0][0]
        ymax = indvalidrow[-1][0] - self._extendmask.shape[0] + self.source.imgsize[0]
        xmin = indvalidcol[0][0]
        xmax = indvalidcol[-1][0] - self._extendmask.shape[1] + self.source.imgsize[1]
        self.setextend(xmin, ymin, xmax, ymax)

    def parsetransformation_beforeapplication(self, pairwise):
        """Adapt transformations before applying them"""
        if self.bextendfrommask(pairwise):
            self.extendfrommask()
        else:
            # Corners of the image in the transformed frame
            p0 = self.untransformedimagepolygon()

            # Corners of the transformed image in the raw frame
            ps = self.transformedimagepolygons()

            # Adapt transformation
            if (
                self.pre_transform_requested["pad"]
                or self.post_transform_requested["crop"]
            ):
                # adapt self.extend to either crop or pad
                self.extendfromtransformation(p0, ps)
            else:
                # try to fit as much data in the original image size as possible
                self.minimaltransformation(p0, ps)

        self.calccof_prealign_to_pretransform()

    def transform_axes(self, axes):
        """Image X and Y axes after transformation

        Args:
            list(array)
        Returns:
            list(array)
        """
        if not self.pre_transform["pad"] and not self.post_transform["crop"]:
            return axes

        if self.source.stackdim == 2:
            ind = [0, 1]
        elif self.source.stackdim == 1:
            ind = [0, 2]
        else:
            ind = [1, 2]

        out = list(axes)
        for i in range(len(ind)):
            j = ind[i]
            nleft = self.extend[i][0]
            nright = self.extend[i][1]
            naxis = len(axes[j])
            newaxis = np.empty(naxis + nleft + nright, dtype=axes[j].dtype)

            off0 = 0
            axis0 = 0
            axisn = naxis
            if nleft < 0:
                axis0 -= nleft
            else:
                off0 += nleft
            if nright < 0:
                axisn += nright
            offn = off0 + axisn - axis0

            if nleft > 0:
                delta = axes[j][1] - axes[j][0]
                newaxis[0:nleft] = (axes[j][0] - delta * nleft) + delta * np.arange(
                    nleft
                )

            newaxis[off0:offn] = axes[j][axis0:axisn]

            if nright > 0:
                delta = axes[j][-1] - axes[j][-2]
                newaxis[offn:] = (axes[j][-1] + delta) + delta * np.arange(nright)

            out[j] = newaxis
        return out

    def getaxesaftertransformation(self, axes):
        """Image X and Y axes after transformation

        Args:
            list(array)
        """
        if not self.pre_transform["pad"] and not self.post_transform["crop"]:
            return

        if self.source.stackdim == 2:
            ind = [0, 1]
        elif self.source.stackdim == 1:
            ind = [0, 2]
        else:
            ind = [1, 2]

        for i in range(len(ind)):
            j = ind[i]
            nleft = self.extend[i][0]
            nright = self.extend[i][1]
            naxis = len(axes[j])
            newaxis = np.empty(naxis + nleft + nright, dtype=axes[j].dtype)

            off0 = 0
            axis0 = 0
            axisn = naxis
            if nleft < 0:
                axis0 -= nleft
            else:
                off0 += nleft
            if nright < 0:
                axisn += nright
            offn = off0 + axisn - axis0

            if nleft > 0:
                delta = axes[j][1] - axes[j][0]
                newaxis[0:nleft] = (axes[j][0] - delta * nleft) + delta * np.arange(
                    nleft
                )

            newaxis[off0:offn] = axes[j][axis0:axisn]

            if nright > 0:
                delta = axes[j][-1] - axes[j][-2]
                newaxis[offn:] = (axes[j][-1] + delta) + delta * np.arange(nright)

            axes[j] = newaxis

    def dest_imgsize(self, nopost=False):
        imgsize = self.source.imgsize
        if self.pre_transform["pad"] or (self.post_transform["crop"] and not nopost):
            imgsize = (
                imgsize[0] + self.extend[0][0] + self.extend[0][1],
                imgsize[1] + self.extend[1][0] + self.extend[1][1],
            )
        return imgsize

    def preparedestination(self, img=None):
        """Allocate space for saving results"""
        if img is not None:
            self.setup_post_transform(img)

        nimages = self.source.nimages
        imgsize = self.dest_imgsize()
        self.dest.prepare(nimages, imgsize, self.dtype)

    def set_reference(self, img, previous=False):
        raise NotImplementedError()

    def doalign(self, refdatasetindex, refimageindex=None, aligntype=alignType.full):
        """Align datasets and save the result."""
        pairwise = refimageindex is None

        # Prepare destination (will be done after alignment)
        if aligntype != alignType.calctransfo:
            self.preparedestination()

        # First reference image
        if aligntype != alignType.usetransfo:
            if pairwise:
                # Pair-wise alignment: first image is the first reference
                imgref = self.readimgrawprep(refdatasetindex, 0)
                iref = 0
            else:
                # Fixed-reference alignment
                rawprep = self.readimgrawprep(refdatasetindex, refimageindex)
                iref = refimageindex
                self.plot(rawprep, 0, "Reference %d (fixed)" % iref)
                self.set_reference(rawprep)

        # from pympler import tracker
        # tr = tracker.SummaryTracker()
        # s1 = None

        # Loop over the images
        bfirst = True
        for i in range(self.source.nimages):
            if aligntype != alignType.usetransfo:
                # Image i
                rawprep = self.readimgrawprep(refdatasetindex, i)
                # np.save("img{}.npy".format(i),rawprep)

                # Update fixed image
                if pairwise:
                    self.set_reference(imgref)

                # Get align transformation
                logger.debug("Align image %d on %d" % (i, iref))

                if i == iref:
                    self.settransformidentity(i)
                    imgaligned = rawprep
                else:
                    # Align image i to reference

                    # if s1 is None:
                    #    s1 = tr.create_summary()

                    imgaligned = self.execute_alignkernel(rawprep)
                    # s2 = tr.create_summary()

                    # tr.print_diff(summary1=s1,summary2=s2)
                    if pairwise:
                        self.plot(imgref, 0, "Reference %d (pair-wise)" % (iref))
                    self.plot(rawprep, 2, "To align %d" % i)
                    self.plot(imgaligned, 1, "Aligned %d" % i)
                    self.store_transformation(i, pairwise)

                logger.debug("Change-of-frame {}".format(self.transfos[i]))

                # Reference for the next image
                if pairwise:
                    if self.alignonraw:
                        imgref = rawprep
                    else:
                        imgref = imgaligned
                    iref = i

                # Only calculation
                if aligntype == alignType.calctransfo:
                    # TODO: This is still not good enough because
                    # align and transformation kernels could be different kernels
                    if self.bextendfrommask(pairwise):
                        self.setextendmask(imgaligned, reset=bfirst)
                    bfirst = False
                    continue  # no results needed

            # Save the transformed image i of all datasets
            if self.pureidentity(i):
                for j in range(self.source.nsets):
                    self.copyimg(j, i)
            else:
                for j in range(self.source.nsets):
                    usealignedimage = (
                        aligntype != alignType.usetransfo
                        and j == refdatasetindex
                        and self.nopre_align()
                        and self.nopre_transform(i)
                        and self.nopost_transform()
                    )
                    if usealignedimage:
                        img = imgaligned
                    else:
                        img = self.readimgraw(j, i)
                        img = self.transform(img, i)
                    self.writeimg(img, j, i)

    def prepare_pre_align(self, roi=None, transfo=None):
        if not isinstance(transfo, (list, tuple)):
            if transfo is None:
                transfo = self.defaulttransform()
            transfo = [transfo] * self.source.nimages
        self.pre_transfos = transfo
        if roi is None:
            self.pre_align["roi"] = None
        else:

            def bdefault(a):
                if a is None:
                    return 0
                else:
                    return a

            def edefault(a):
                if a is None:
                    return -1
                else:
                    return a

            self.pre_align["roi"] = (
                (bdefault(roi[0][0]), edefault(roi[0][1])),
                (bdefault(roi[1][0]), edefault(roi[1][1])),
            )
        self.calccof_prealign_to_raw()

    def align(
        self,
        refdatasetindex,
        refimageindex=None,
        onraw=False,
        pad=True,
        crop=False,
        redo=False,
        roi=None,
        rawcalc=None,
        prealigntransfo=None,
    ):
        """Alignment function that needs to be called

        Args:
            refdatasetindex(int): stack to be used for alignment
            refimageindex(Optional(int)): fixed reference to align on
                                          pairwise alignment when None
            onraw(Optional(bool)): when doing pairwise alignment, use the
                                   previous raw or aligned images to align
                                   the next one on
            pad(Optional(bool)): make sure nothing is transformed outside
                                 the field of view (has priority over crop)
            crop(Optional(bool)): make sure no missing data is added
            redo(Optional(bool)): apply transformations without recalculating
                                  them (all other keywords are ignored)
            roi(Optional(array-like)): use ony part of the image to align
            rawcalc(Optional(callable)): apply function to raw data before
                                         alignment (not used in transformation)
            prealigntransfo(Optional(list(LinearMapping))): apply to raw data before alignment
                                         (propogates to transformation)
        """
        if redo:
            self.doalign(refdatasetindex, aligntype=alignType.usetransfo)
        else:
            pairwise = refimageindex is None
            if pairwise:
                self.alignonraw = onraw
            else:
                self.alignonraw = True
            self.prepare_pre_align(roi=roi, transfo=prealigntransfo)
            self.pre_transform_requested["pad"] = pad
            self.post_transform_requested["crop"] = crop
            self.pre_transform["pad"] = pad
            self.post_transform["crop"] = crop
            # center = False
            if roi or pad or crop or callable(rawcalc) or prealigntransfo:
                # Do not use transformed image of alignment procedure
                self.pre_align["func"] = rawcalc
                self.doalign(
                    refdatasetindex,
                    refimageindex=refimageindex,
                    aligntype=alignType.calctransfo,
                )
                self.pre_align["func"] = None
                self.parsetransformation_beforeapplication(pairwise)
                self.doalign(refdatasetindex, aligntype=alignType.usetransfo)
            else:
                # Use transformed image of alignment procedure
                self.doalign(refdatasetindex, refimageindex=refimageindex)
