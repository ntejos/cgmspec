from cgmspec.disco import Disco
from lmfit import Minimizer, Parameters, fit_report
import numpy as np
import pylab as plt
from mpdaf.obj import Cube, WCS, WaveCoord, Spectrum
import astropy.units as u

params = Parameters()
params.add_many(
('incli', 60, True),
('col_dens', 14, True),
('height', 11.5, True),
('dop_param', 7.0, True),
('csize', 0.1, False),
('r_0', 1, True),
('vel_max', 195, True),
('h_v', 4.5, True),
)

par = params.valuesdict()
incli = par['incli']
col_denst = par['col_dens']
h = par['height']
bes = par['dop_param']
v_max = par['vel_max']
h_vt = par['h_v']
csize = par['csize']
r_0t= par['r_0']
r_0 = 10**r_0t
h_v = 10**h_vt

model = Disco(h, incli, Rcore=0.1)
lam2 = np.arange(4825.12, 4886.37+1.5, 0.1)
wave1 = WaveCoord(cdelt=0.1, crval=4825.12, cunit= u.angstrom)
spec = model.averagelos(5, 25, lam2, 100,12,0.73379,csize, bes, r_0, v_max,h_v, 0)
spe = Spectrum(wave=wave1, data=spec)
rspec = spe.resample(1.25)
wave = np.arange(rspec.get_start(), rspec.get_end()+rspec.get_step(), rspec.get_step())
print(rspec)
plt.plot(lam2, spec,  'b')
plt.plot(wave, rspec.data, 'g')
plt.show()
