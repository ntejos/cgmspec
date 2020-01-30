import numpy as np
import math
from math import log
from astropy import constants as const
from scipy.special import wofz

"""Utilities for cgmspec"""

def prob_hit(r, rmin, rmax, prob_rmin=100., prob_rmax=20.):
    """
    Probability of hitting a cloud at distance r in the plane xy of a disc of radius rmax

    :param r: np.array, array with distances to the center of disc in plane xy in kpc
    :param rmin: float, minimum radius of the disc in plane xy in kpc, below this value probability is set by prob_rmin
    :param rmax: float, radius of the disc in plane xy in kpc, above this value prob is set to zero
    :param prob_rmin: float, probability at rmin or below
    :param prob_rmax: float, probability at rmax
    :return: float, probability of hitting a cloud
    """

    ind = np.log(prob_rmin/prob_rmax)/np.log(rmin/rmax)
    A = prob_rmax/(rmax**(ind))
    prob = A*(r**(ind))
    prob = np.where(r>rmax, 0., prob)
    prob = np.where(r<=rmin, prob_rmin, prob)

    '''b = np.log10(prob_rmax / prob_rmin) / np.log10(rmax)
    prob = prob_rmin * (r ** b)
    prob = np.where(r>rmax, 0., prob)
    prob = np.where(r<=rmin, prob_rmin, prob)'''
    return prob


def los_disc_intersect(D, i, alpha, R, h):
    """For a given set of (i, D, alpha) + (R, h), this function return (x0, yt1, yt2, zt1, zt2)
    If the sightline does not intersect the disk, it returns False
    """

    radio = R
    hdis = h
    incli = i

    al_rad = np.radians(alpha)
    x0 = D * np.cos(al_rad)

    # check if sightline intersects the disk
    if x0 > radio:
        return False

    if incli == 90:  # edge-on
        z0 = D * np.sin(al_rad)
        if np.fabs(z0) > hdis / 2.:
            return False  # outside projected disk
        else:
            yt = [-np.sqrt((radio ** 2) - x0 ** 2), np.sqrt((radio ** 2) - x0 ** 2)]
            zt = [D * np.sin(al_rad), D * np.sin(al_rad)]

    elif incli == 0.0:  # face-on
        yt = [D * np.sin(al_rad), D * np.sin(al_rad)]
        zt = [-hdis / 2, hdis / 2]
        incli_rad = np.radians(incli)
        # print(yt, zt)

    else:
        incli_rad = np.radians(incli)
        y0 = D * np.sin(al_rad) / np.cos(incli_rad)
        # print(y0)

        z0 = 0

        radio1 = np.sqrt((radio ** 2) - x0 ** 2)
        # print(x0,y0)
        ymin = y0 - (hdis / 2) * np.tan(incli_rad)
        ymax = y0 + (hdis / 2) * np.tan(incli_rad)
        zmin = -(np.sqrt((-radio1 - y0) ** 2) / np.tan(incli_rad))
        zmax = (np.sqrt((radio1 - y0) ** 2) / np.tan(incli_rad))
        ys = [ymin, ymax, -radio1, radio1]
        zs = [-hdis / 2, hdis / 2, zmin, zmax]
        yt = []
        zt = []
        #  print(ymin, ymax, zmin, zmax)

        for i in range(len(ys)):
            if abs(ys[i]) < radio1:
                #           print(ys[i], zs[i])
                yt.append(ys[i])
                zt.append(zs[i])
            else:
                pass
        for i in range(len(zs)):
            if abs(zs[i]) < hdis / 2:
                 #           print(ys[i], zs[i])
                yt.append(ys[i])
                zt.append(zs[i])
            else:
                pass

        yt = np.sort(yt)
        zt = np.sort(zt)

    #  print(yt,zt)
    return x0, yt[0], yt[1], zt[0], zt[1]



def Tau(lam,vel):
    N = 10**(14)
    redshift = 0.7
    b = 5

    lam0, f, gamma, mass = [2796, 0.6155, 2.68e8, 24.305]
    c  = const.c.to('cm/s').value
    sigma0 = 0.0263

    lamc = ((vel/const.c.to('km/s').value)+1)*((lam0))
    nu = c/(lam*1e-8)
    nu0 = c/(lamc*1e-8)

    dnu = nu - (nu0/(1+redshift))
    dnud = (b*100000)*nu/c

    x = dnu/dnud
    y = gamma/(4*np.pi*dnud)
    z = x + 1j*y
    v = np.asarray(np.real(wofz(z)/(np.sqrt(np.pi)*dnud)))

    taut = N * sigma0 * f  *v
    vt = np.exp(-taut)
    obs_spec_wave = lam / (1+redshift)
    dv = (const.c.to('km/s').value * ((lam / (lamc *(1 + redshift)))-1))

    velt = []
    for i in range(len(dv)):
        velt.append(dv[i]+vel)

    return(taut,dv)


def normflux(taun):
    taun = np.asarray(taun)
    return(np.exp(-taun))


def sumtau(taus):
    tautosum =  []
    for i in range(len(taus)):
        tautosum.append(taus[i][0])
    tautosum = np.asarray(tautosum)
    #print(tautosum)

    tausum = []
    for i in range(len(tautosum[0])):
        tot = 0
        for j in range(len(tautosum)):
            tot = tot + tautosum[j][i]
        tausum.append(tot)
    #print(tausum)

    return(tausum)
