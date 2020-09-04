
from mpdaf.obj import Cube, WCS, WaveCoord, Spectrum
from mpdaf.obj import iter_spe
import numpy as np
import astropy.units as u

datacube  = Cube('PSZ1GA311_cube_absplane.fits')

norcube = datacube.clone(data_init=np.empty, var_init=np.zeros)

wavep = 4855

for i in range(datacube.shape[1]):
    for j in range(datacube.shape[2]):
        spe = datacube[:,i,j]
        var = spe.var
        sigma = np.sqrt(var)
#arrays
        flux = spe.data
        sigma = sigma.data
        wave = spe.wave.coord()
        #vel = (wave/(wave0*(1+zabs))-1)*299792.458
#mask pixels for continuum definition
        p1 = np.int(spe.wave.pixel(wavep -60*1.25))
        p2= np.int(spe.wave.pixel(wavep -30*1.25))
        p3 = np.int(spe.wave.pixel(wavep +30*1.25))
        p4 = np.int(spe.wave.pixel(wavep +100*1.25))
        flux1  = flux[np.r_[p1:p2,p3:p4]]
        wave1  = wave[np.r_[p1:p2,p3:p4]]
#continuum
        z = np.polyfit(wave1,flux1,2)
        p = np.poly1d(z)
        cont = p(wave)
        nflux = flux/cont
        stdev = sigma/cont   #for displaying
#mpdaf normalized continuum
        wave_aux = WaveCoord(cdelt=1.25, crval=wave[0], cunit= u.angstrom)
        nspe = Spectrum(wave=wave_aux, data=nflux, var=(stdev**2))
        norcube[:,i,j]  = nspe
        #co[:] = sp-ai
