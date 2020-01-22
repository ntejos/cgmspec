import numpy as np
import math

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

    b = np.log10(prob_rmax / prob_rmin) / np.log10(rmax)
    prob = prob_rmin * (r ** b)
    prob = np.where(r>rmax, 0., prob)
    prob = np.where(r<rmin, prob_rmin, prob)
    return prob

