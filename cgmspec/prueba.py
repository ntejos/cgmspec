from pylab import *
import numpy as np
import math
from math import log
from astropy import constants as const
from scipy.special import wofz
from astropy import units as u
from astropy.coordinates import SkyCoord

xgal_l=237.52059628552735
ygal_l=-78.188149277705151

posgal = SkyCoord(xgal_l*u.deg, ygal_l*u.deg, distance=4512.83311699*u.mpc, frame='icrs')

Ds1 = [-3.6, 1.4, 3.6,6.8,9.9, 13.2, 16.5, 19.8, 13.0, 26.3,29.7]
alphas1 = [10,10,10,10,10,10,10,10,10,10,10]
xabs_l = np.load('xabs_l.npy')
yabs_l = np.load('yabs_l.npy')

alphas = []
for i in  range(len(xabs_l)):
     locals()["pos"+str(i)] = SkyCoord(xabs_l[i]*u.deg, yabs_l[i]*u.deg, distance=4512.83311699*u.mpc, frame='icrs')
     locals()["pa"+str(i)] = posgal.position_angle(locals()["pos"+str(i)]).to(u.deg)
     alphas.append(90 - locals()["pa"+str(i)].value)

Ds = []
for i in range(len(xabs_l)):
    sep  = posgal.separation_3d(locals()["pos"+str(i)]).value  *1000
    Ds.append(sep)

x =  []
y = []

x1 = []
y1 =[]

for i in range(len(Ds)):
     xi = Ds[i] * np.cos(np.radians(alphas[i]))
     yi = -Ds[i] * np.sin(np.radians(alphas[i]))

     x.append(xi)
     y.append(yi)

     xi1 = Ds1[i] * np.cos(np.radians(alphas1[i]))
     yi1 = -Ds1[i] * np.sin(np.radians(alphas1[i]))

     x1.append(xi1)
     y1.append(yi1)

plt.figure()
plt.scatter(x, y)
plt.scatter(x1, y1)
