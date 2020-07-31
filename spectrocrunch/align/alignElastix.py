# -*- coding: utf-8 -*-

import logging
import warnings
import numpy as np

try:
    import SimpleITK as sitk
except ImportError:
    sitk = None
    warnings.warn("SimpleElastix is not installed", ImportWarning)

from .align import align
from ..utils.stdout import stdout_redirect
from .types import transformationType


class alignElastix(align):
    def __init__(self, *args, **kwargs):
        super(alignElastix, self).__init__(*args, **kwargs)

        # No 1D
        if 1 in self.source.imgsize:
            raise ValueError("Elastix can only be applied on images, not 1D vectors.")

        # Prepare alignment kernel
        try:
            self.elastix = sitk.ElastixImageFilter()
        except:
            self.elastix = sitk.SimpleElastix()

        # self.elastix.LogToFolder("")
        # self.elastix.LogToFolderOff()
        self.elastix.LogToConsoleOff()
        self.SetParameterMap()
        self.fixed = None
        self.moving = None

        # Prepare transformation kernel
        try:
            self.transformix = sitk.TransformixImageFilter()
        except:
            self.transformix = sitk.SimpleTransformix()

    def replacecval(self, img):
        if self.cval is not np.nan or self.defaultvalue is not np.nan:
            if self.cval != self.defaultvalue:
                if self.defaultvalue is np.nan:
                    img[np.isnan(img)] = self.cval
                else:
                    img[img == self.defaultvalue] = self.cval

    def SetParameterMap(self):
        if self.transfotype == transformationType.translation:
            parameterMap = sitk.GetDefaultParameterMap("translation")
        elif self.transfotype == transformationType.rigid:
            parameterMap = sitk.GetDefaultParameterMap("rigid")
        # elif self.transfotype==transformationType.similarity:
        #    parameterMap = sitk.GetDefaultParameterMap('similarity') # in Elastix but not SimpleElastix
        elif self.transfotype == transformationType.affine:
            parameterMap = sitk.GetDefaultParameterMap("affine")
        else:
            raise NotImplementedError(
                "Elastix doesn't support this type transformation."
            )

        self.defaultvalue = -999  # Elastix cannot handle NaN's!!!
        parameterMap["DefaultPixelValue"] = (str(self.defaultvalue),)
        self.elastix.SetParameterMap(parameterMap)

    def execute_transformkernel(self, img):
        """Transform image according with the transformation kernel
        """
        self.changefortransform(img.shape)
        img[np.isnan(img)] = 0
        self.transformix.SetInputImage(sitk.GetImageFromArray(img))
        with stdout_redirect():
            self.transformix.Execute()
        aligned = sitk.GetArrayFromImage(self.transformix.GetResultImage())
        self.replacecval(aligned)
        return aligned

    def execute_alignkernel(self, img):
        """Align image on reference
        """
        img[np.isnan(img)] = 0
        self.moving = sitk.GetImageFromArray(img)
        self.elastix.SetMovingImage(self.moving)

        try:
            with stdout_redirect():
                self.elastix.Execute()
            aligned = sitk.GetArrayFromImage(self.elastix.GetResultImage())
            self.replacecval(aligned)
        except:
            aligned = img.copy()

        return aligned

    def set_reference(self, img, previous=False):
        """Reference for alignment
        """
        if previous:
            self.fixed = self.moving
        else:
            img[np.isnan(img)] = 0
            self.fixed = sitk.GetImageFromArray(img)
        self.elastix.SetFixedImage(self.fixed)

    def elastix_GetTransformParameterMap(self):
        if hasattr(self.elastix, "GetNumberOfParameterMaps"):
            n = self.elastix.GetNumberOfParameterMaps()
        else:
            n = 1
        if n == 0:
            logger = logging.getLogger(__name__)
            logger.info("Elastix couldn't align images")
            return []
        else:
            try:
                return self.elastix.GetTransformParameterMap()
            except:
                logger = logging.getLogger(__name__)
                # import traceback
                # logger.debug(traceback.format_exc())
                logger.info("Elastix couldn't align images")
                return []

    def get_alignkernel(self):
        """Get transformation from alignment kernel.
        """
        transform = self.defaulttransform()

        transformParameterMap = self.elastix_GetTransformParameterMap()
        if len(transformParameterMap) >= 1:
            # for k in transformParameterMap[0]:
            #    print(k,transformParameterMap[0][k])

            if self.transfotype == transformationType.translation:
                params = np.array(
                    transformParameterMap[0]["TransformParameters"], self.dtype
                )
                transform.settranslation(params[0:2])
            elif self.transfotype == transformationType.rigid:
                theta, tx, ty = np.array(
                    transformParameterMap[0]["TransformParameters"], self.dtype
                )
                cx, cy = np.array(
                    transformParameterMap[0]["CenterOfRotationPoint"], self.dtype
                )

                L = self.defaulttransform(ttype=transformationType.rigid)
                L.setrigid(theta, tx, ty)

                C = self.defaulttransform(ttype=transformationType.translation)
                C.settranslation(-cx, -cy)

                # cof = C^-1.L^-1.C

                transform.fromtransform(L.dot(C).inverse().dot(C))

                # TODO: not verified!
            else:
                raise NotImplementedError(
                    "Elastix doesn't support this type of transformation."
                )

        return transform

    def set_transformkernel(self, transform):
        """Set the transformation kernel according to the alignment kernel and adapted transformation
        """
        transformParameterMap = self.elastix_GetTransformParameterMap()
        if transform.transfotype != self.transfotype:
            raise ValueError("Transformations must have the same type")
        if self.transfotype == transformationType.translation:
            tmp = transform.gettranslation()
            transformParameterMap[0]["TransformParameters"] = (str(tmp[0]), str(tmp[1]))
        elif self.transfotype == transformationType.rigid:
            transformParameterMap[0]["CenterOfRotationPoint"] = ("0", "0")
            theta, tx, ty = transform.getrigid()
            transformParameterMap[0]["TransformParameters"] = (
                str(theta),
                str(tx),
                str(ty),
            )
        elif self.transfotype == transformationType.affine:
            raise NotImplementedError
        else:
            raise NotImplementedError(
                "Elastix doesn't support this type of transformation."
            )
        self.transformix.SetTransformParameterMap(transformParameterMap)

    def changefortransform(self, shape):
        transformParameterMap = self.elastix_GetTransformParameterMap()

        oldshape = (
            int(transformParameterMap[0]["Size"][1]),
            int(transformParameterMap[0]["Size"][0]),
        )
        if shape == oldshape:
            return

        # Not sure about this
        transformParameterMap[0]["Size"] = (str(shape[1]), str(shape[0]))
        transformParameterMap[0]["Origin"] = (str(self.origin[0]), str(self.origin[1]))

        self.transformix.SetTransformParameterMap(transformParameterMap)
