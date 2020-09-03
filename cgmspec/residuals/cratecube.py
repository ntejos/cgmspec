import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from mpdaf.obj import Cube, WCS, WaveCoord, Spectrum
from cgmspec.disco import Disco
from astropy.coordinates import SkyCoord
from astropy.cosmology import FlatLambdaCDM
from mpdaf.obj import iter_spe
from astropy.convolution import convolve, Gaussian1DKernel






lam1 = np.arange(4825.12, 4886.37+1.25, 0.1)

z = 0.73379
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
dist = cosmo.luminosity_distance(z)

#galaxy imputs

galPA = 55

R =  50
h = 5
galcen = SkyCoord(ra=237.52059628552735*u.degree, dec=-78.188149277705151*u.degree, frame='icrs')



datacube  = Cube('ccdatanorm.fits')
xlen = len(datacube[0,:,0].data)
ylen = len(datacube[0,0,:].data[0])

head = datacube.primary_header
wcs1 = WCS(head)
wave1 = WaveCoord(cdelt=0.5, crval=4825.12, cunit=u.angstrom)
rawdat = np.ones((50, xlen, ylen))
modelcube = Cube(data=rawdat, wcs=wcs1, wave=wave1)
centroid = []


#calcular para lista de parametros
'''inclis= [80, 80, 80, 80, 80, 80, 80, 80, 80]
csize = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
Nes = [10**15, 10**15, 10**15, 10**15, 10**15, 10**15, 10**15, 10**15, 10**15]
bes = [5, 5,  5, 5, 5,  5, 5,5, 5]
[2, 4, 6, 8, 10, 12, 14, 16, 18]
pmax = [100, 100, 100, 100, 100, 100, 100, 100, 100]
pmin = [20, 20, 20, 20, 20, 20, 20, 20, 20]
hs = [2, 4, 6, 8, 10, 12, 14, 16, 18]

for k in range(len(inclis)):
            model = Disco(R, hs[k], inclis[k], Rcore=0.1)
            coord_sky = datacube.wcs.pix2sky([6, 19], unit=u.deg)
            dec = coord_sky[0][0]
            ra = coord_sky[0][1]
            scale = 7.28 # z=0.73379  # plat scale (kpc/") for Planck
            c1 = SkyCoord(ra*u.degree,dec*u.degree,frame='icrs')
            D = scale * galcen.separation(c1).arcsec
            pa = galcen.position_angle(c1).to(u.deg)
            alpha = galPA - pa.value
            print(D,alpha)
            wave1 = WaveCoord(cdelt=0.1, crval=4751.37, cunit= u.angstrom, shape=247.5)

            spec = model.averagelosspec(D, alpha, lam1, 10,12,z,csize[k], Nes[k], bes[k], pmax[k], pmin[k])
            centi = spec[3]
            centroid.append(centi)
            spe = Spectrum(wave=wave1, data=spec[1])
            rspe = spe.resample(1.25)
            modelcube[:,15,8] = rspe.data
            numero = 5 + k
            modelcube.write('paramc_h_p3_%s' %k  +'.fits')'''


inclis= 60
csize = 0.1
Nes = 10**15
bes = 5
pmax = 100
pmin = 20
hs = 5

model = Disco(R, hs, inclis, Rcore=0.1)

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
        spec = model.averagelosspec(D, alpha, lam1, 10,12,z,csize, Nes, bes, pmax, pmin, 0)
        spe = Spectrum(wave=wave1, data=spec[1])
        rspe = spe.resample(1.25)
        modelcube[:,i-1,j-1] = rspe.data
