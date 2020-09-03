import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from mpdaf.obj import Cube, WCS, WaveCoord, Spectrum
from cgmspec.disco import Disco
from astropy.coordinates import SkyCoord
from astropy.cosmology import FlatLambdaCDM
from mpdaf.obj import iter_spe
from astropy.convolution import convolve, Gaussian1DKernel




lam1 = np.arange(4825.12, 4886.37+1.25, 0.1) #wavelenght to calculate the spectra

#cosmological parameters to calculate the distance
z = 0.73379
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
dist = cosmo.luminosity_distance(z)

#galaxy imputs
galPA = 55
galcen = SkyCoord(ra=237.52059628552735*u.degree, dec=-78.188149277705151*u.degree, frame='icrs')



datacube  = Cube('ccdatatrue.fits')
xlen = len(datacube[0,:,0].data)
ylen = len(datacube[0,0,:].data[0])

#create a model cube to put all the model spectras
head = datacube.primary_header
wcs1 = WCS(head)
wave1 = WaveCoord(cdelt=0.5, crval=4825.12, cunit=u.angstrom)
rawdat = np.ones((50, xlen, ylen))
modelcube = Cube(data=rawdat, wcs=wcs1, wave=wave1)

#Model parameters

incli = 80
col_dens = 10**15
h = 5
b = 5
v_max = 200
h_v = 10
csize = 0.1
r_0= 1000

model = Disco(h, incli, Rcore=0.1)

for i in range(xlen):
    for j in range(ylen):
        coord_sky = datacube.wcs.pix2sky([i-1, j-1], unit=u.deg)
        dec = coord_sky[0][0]
        ra = coord_sky[0][1]
        scale = 7.28 # z=0.73379  # plat scale (kpc/") for Planck
        c1 = SkyCoord(ra*u.degree,dec*u.degree,frame='icrs')
        D = scale * galcen.separation(c1).arcsec
        pa = galcen.position_angle(c1).to(u.deg)
        alpha = galPA - pa.value
        wave1 = WaveCoord(cdelt=0.1, crval=4751.37, cunit= u.angstrom, shape=247.5)
        spec = model.averagelosspec(D, alpha, lam1, 100,12,z,csize, col_dens, b, r_0, v_max, h_v, 0)
        spe = Spectrum(wave=wave1, data=spec[1])
        rspe = spe.resample(1.25)
        modelcube[:,i-1,j-1] = rspe.data

modelcube.write('modelcube_i%s' %incli  +'_N%s' %col_dens + '_h%s' %h+'_b%s' %b + '_vmax%s'%v_max +'_hv%s' %h_v+'_csize%s'%csize+'_r0%s'%r_0+'.fits')
