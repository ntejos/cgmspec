from mpdaf.obj import Cube, WCS, WaveCoord, Spectrum
from astropy.convolution import convolve, Gaussian1DKernel
import numpy as np

datacube = Cube('ccdatanorm.fits')



modelcube = Cube('ex4.fits')



'''def filtrogauss(R,tau,dispersion=1.25,lam=4850):
    delta_lambda = lam/R
    fw = delta_lambda/dispersion

    gauss_kernel = Gaussian1DKernel(fw)
    gausflux = convolve(tau, gauss_kernel)
    return(gausflux)

for i in range(modelcube.shape[1]):
    for j in range(modelcube.shape[2]):
        spei = modelcube[:,i,j].data
        gausfluxi = filtrogauss(2100, spei)
        modelcube[:,i,j] = gausfluxi'''


compcube = datacube.clone(data_init=np.empty, var_init=np.zeros)

for i in range(datacube.shape[1]):
    for j in range(datacube.shape[2]):
        for k in range(len(compcube[:,i,j].data)):
            compcube[k,i,j] = (datacube[k,i,j] - modelcube[k,i,j])/np.sqrt(datacube.var[k,i,j])
