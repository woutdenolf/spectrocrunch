import SimpleITK as sitk
import numpy as np
import h5py
import matplotlib.pyplot as plt
import fabio

def align(fixedImage,movingImage):
    o = sitk.SimpleElastix()
    #o.LogToConsoleOn()

    o.SetFixedImage(sitk.GetImageFromArray(fixedImage))
    o.SetMovingImage(sitk.GetImageFromArray(movingImage))
    o.SetParameterMap(sitk.GetDefaultParameterMap("translation"))
    o.Execute()
    alignedImage = sitk.GetArrayFromImage(o.GetResultImage())

    return o.GetTransformParameterMap()[0]["TransformParameters"],alignedImage

def test():
    fixedImage = np.zeros((5,5))
    movingImage = np.zeros((5,5))
    fixedImage[2,2] = 1
    movingImage[3,3] = 1

    trn,alignedImage = align(fixedImage,movingImage)

    np.testing.assert_array_almost_equal(fixedImage,alignedImage,decimal=1)

def test2():
    h5name = "/data/visitor/hc3130/id21/spectrocrunch/crunched/bife_AG_fmap4/bife_AG_fluoXAS_4.norm.h5"

    plt.figure()

    with h5py.File(h5name) as hdf5FileObject:
        
        fixedImage = hdf5FileObject["counters"]["arr_fdet"]["data"][...,0]

        for k in range(80,85):
            print k
            movingImage = hdf5FileObject["counters"]["arr_fdet"]["data"][...,k]
            trn,alignedImage = align(fixedImage,movingImage)

            plt.imshow(movingImage)
            plt.pause(0.01)

def test3():

    dtype = float

    fh = fabio.open("/mntdirect/_data_visitor/hc3130/id21/sampleV/sampleV_V1015_small/zap/sampleV_V1015_small_arr_absorp2_0001_0000.edf")  
    #fh = fabio.open("/mntdirect/_data_visitor/hc3130/id21/sampleV/sampleV_V1015_small_macro2/zap/sampleV_V1015_small_macro2_arr_absorp2_0024_0000.edf")
    fixedImage =  fh.data.astype(dtype)       

    fh = fabio.open("/mntdirect/_data_visitor/hc3130/id21/sampleV/sampleV_V1015_small_macro2/zap/sampleV_V1015_small_macro2_arr_absorp2_0054_0000.edf")
    #fh = fabio.open("/mntdirect/_data_visitor/hc3130/id21/sampleV/sampleV_V1015_small_macro2/zap/sampleV_V1015_small_macro2_arr_absorp2_0063_0000.edf")
    movingImage =  fh.data.astype(dtype)       
    
    plt.figure(1)
    plt.imshow(fixedImage)
    plt.figure(2)
    plt.imshow(movingImage)

    trn,alignedImage = align(fixedImage,movingImage)
    plt.figure(3)
    plt.imshow(alignedImage)

    print trn

    plt.show()

if __name__ == '__main__':
    test3()
    


