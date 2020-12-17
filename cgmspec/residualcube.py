import matplotlib.pyplot as plt
from cgmspec import utils as csu
from scipy.optimize import curve_fit
from mpdaf.obj import Cube, WCS, WaveCoord, Spectrum, Image
import astropy.units as u
from cgmspec.disco import Disco
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.cosmology import FlatLambdaCDM
from lmfit import Model
from scipy.optimize import least_squares
from lmfit import Minimizer, Parameters, fit_report
from matplotlib.colors import LogNorm

params = Parameters()
params.add_many(
        ('incli', 85, True),
        ('col_dens', 14, True),
        ('height', 11.5, True),
        ('dop_param', 7.0, True),
        ('csize', 0.1, False),
        ('r_0', 1, True),
        ('vel_max', 195, True),
        ('h_v', 4.5, True),
        )

z = 0.73379
scale = 7.28
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
dist = cosmo.luminosity_distance(z)

#galaxy imputs
galPA = 55
galcen = SkyCoord(ra=237.52059628552735*u.degree, dec=-78.188149277705151*u.degree, frame='icrs')

specwidth =  [14,36]

s_n_mask = np.load('out_r4_3.npy')
all_pixels = []
for i in range(len(s_n_mask)):
    posi = [s_n_mask[i][0], s_n_mask[i][1]]
    all_pixels.append(posi)

all_pixels = np.asarray(all_pixels)

datacube  = Cube('nor_cut_cube.fits')


xlen = len(datacube[0,:,0].data)
ylen = len(datacube[0,0,:].data[0])

wcs1 = WCS(crval=0, cdelt=0.2)
MyData = np.ones((xlen,ylen))
ima = Image(data=MyData, wcs=wcs1)




for i in range(xlen):
    for j in range(ylen):
        par = params.valuesdict()
        incli = par['incli']
        col_denst = par['col_dens']
        h = par['height']
        bes = par['dop_param']
        v_max = par['vel_max']
        h_vt = par['h_v']

        csize = par['csize']
        r_0t= par['r_0']

        h_v = 10**h_vt
        col_dens = 10**col_denst
        r_0 = 10**r_0t

        pixi = np.array([[i,j]])

        if (pixi[:, None] == all_pixels).all(-1).any(-1):
            print(i,j)
            coord_sky = datacube.wcs.pix2sky([i, j], unit=u.deg)
            dec = coord_sky[0][0]
            ra = coord_sky[0][1]
            scale = 7.28 # z=0.73379  # plat scale (kpc/") for Planck
            c1 = SkyCoord(ra*u.degree,dec*u.degree,frame='icrs')
            alpha = csu.get_alpha(c1, galcen,galPA)
            D = csu.get_impactparam(c1, galcen, scale)

            flux = datacube[:,i,j]
            abso = flux.data[specwidth[0]:specwidth[1]+1]
            sigma = np.sqrt(flux.var[specwidth[0]:specwidth[1]+1])
            lam2 = np.arange(4825.12, 4886.37+1.5, 0.1)
            wave1 = WaveCoord(cdelt=0.1, crval=4825.12, cunit= u.angstrom)

            model = Disco(h, incli, Rcore=0.1)
            spec = model.averagelos(D, alpha, lam2, 100,12,z,csize, col_dens, bes, r_0, v_max,h_v, 0)
            spe = Spectrum(wave=wave1, data=spec)
            rspe = spe.resample(1.25)
            fluxmodel = rspe.data
            absomodel = fluxmodel[specwidth[0]:specwidth[1]+1]

            #ymodelt = np.concatenate((ymodel[0][14:36],ymodel[1][16:36],ymodel[2][16:36]))
            dif = np.sum(((abso-absomodel)/sigma)**2)
            ima.data[i,j] = dif
            ima.mask[i,j] = False

        else:
            ima.data[i,j] = 0
            ima.mask[i,j] = True


ima.plot(colorbar='v')
plt.show()
