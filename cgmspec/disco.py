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

    def los_vel(self, D, al, y, vR=200):
        """
        :param D: float, impact parameter in kpc
        :param al: float, angle between the mayor axis and the line-of-sight, clockwise, in degrees
        :param y: distance of the line-of-sight to the y axis of the disk
        :param vR: velocity of rotation of the disk
        :return: line-of-sight velocity for a cloud in the given position
        """

        inc = np.radians(self.incl)
        al = np.radians(al)

        x0 = (D*1000)*np.cos(al)
        y0 = (D*1000)*np.sin(al)/np.cos(inc)
        a = np.sin(inc)/np.sqrt(1+(y/x0)**2)
        vr = vR*a
        return(vr)



    def spec_clouds(self):
        pass
