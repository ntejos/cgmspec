import numpy as np
import math
from math import log
from astropy import constants as const
from scipy.special import wofz

"""Utilities for cgmspec"""

def prob_hit_new(r, rmin, rmax, prob_rmin=100., prob_rmax=20.):
    """
    Probability of hitting a cloud at distance r in the plane xy of a disc of radius rmax

    :param r: np.array, array with distances to the center of disc in plane xy in kpc
    :param rmin: float, minimum radius of the disc in plane xy in kpc, below this value probability is set by prob_rmin
    :param rmax: float, radius of the disc in plane xy in kpc, above this value prob is set to zero
    :param prob_rmin: float, probability at rmin or below
    :param prob_rmax: float, probability at rmax
    :return: float, probability of hitting a cloud
    """

    b = np.log10(prob_rmax / prob_rmin) / np.log10(rmax)
    prob = prob_rmin * (r ** b)
    prob = np.where(r>rmax, 0., prob)
    prob = np.where(r<rmin, prob_rmin, prob)
    return prob


def prob_hit(r, rmax):
    A = 100.
    #b = np.log10(20. / A) / np.log10(rmax)
    b = log(10/A, rmax)
    return (A * (r ** b))


def Tau(lam,vel):
    N = 10**(14)
    redshift = 0.7
    b = 10

    lam0, f, gamma, mass = [2852, 0.6155, 2.68e8, 24.305]
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