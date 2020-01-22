import numpy as np
from cgmspec import utils as csu

"""This Class Disco represents the CGM of a galaxy from a disc model"""

class Disco:
    """Represents the CGM of a idealized disc projected in the sky"""

    def __init__(self, R, h, incl):
        """

        :param R: float, Radius of disc in kpc
        :param h: float, height of disc in kpc
        :param incl: float, inclination angle of disc in degrees
        """

        self.R = R
        self.h = h
        self.incl = incl

    def prob_hit(self, r):
        """
        :param r: float, distance to the center of disc in plane xy in kpc
        :return: float, covering fraction at distance r
        """

        return csu.prob_hit(r, self.R)



    def spec_clouds(self):
        pass

