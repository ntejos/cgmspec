import numpy as np
from cgmspec import utils as csu
from matplotlib import pyplot as plt
from matplotlib import patches
from astropy import constants as const

"""This Class Disco represents the CGM of a galaxy from a disc model"""

class Disco:
    """Represents the CGM of a idealized disc projected in the sky"""

    def __init__(self, R, h, incl, Rcore=0.1):
        """

        :param R: float, Radius of disc in kpc
        :param h: float, height of disc in kpc
        :param incl: float, inclination angle of disc in degrees
        :param Rcore: float, radius of disk core in kpc
        """

        self.R = R
        self.Rcore = Rcore
        self.h = h
        self.incl = incl
        self.incl_rad = np.radians(self.incl)

    def prob_hit(self, r):
        """
        :param r: float, distance to the center of disc in plane xy in kpc
        :return: float, covering fraction at distance r
        """

        return csu.prob_hit(r, self.Rcore, self.R)

    def los_vel(self, y, D, alpha, vR=180, hv=5000):
        """
        See Ho et al. 2017

        :param y: np.array, distance of the line-of-sight to the semi-major axis (x-axis) along the y-axis
        :param D: float, impact parameter in kpc
        :param alpha: float, angle between the major axis and the line-of-sight, clockwise, in degrees
        :param vR: maximum velocity of rotation of the disk in km/s
        :param hv: velocity scale height in kpc
        :return: line-of-sight velocity for a cloud in the given position

        """

        al_rad = np.radians(alpha)

        x0 = D * np.cos(al_rad)  # this is p in Ho et al. 2019, fig 10.
        y0 = D * np.sin(al_rad) / np.cos(self.incl_rad)  # this is y0 in the same fig.
        a = np.sin(self.incl_rad) / np.sqrt(1 + (y/x0)**2)
        b = np.exp(-np.fabs(y - y0) / hv * np.tan(self.incl_rad))
        vr = vR*a*b
        return(vr)

    def get_clouds(self, D, alpha, grid_size=1):
        """
        Obtain main parameters to describe clouds intersected by the sightline along the disc

        :param D:
        :param alpha:
        :param radio:
        :param hdis:
        :param grid_size: size of grid in kpc to model a hit or no hit of a cloud
        :return: (nclouds, vel_array, r_array) : (float, array, array)
        """

        incli = self.incl
        radio = self.R
        hdis = self.h
        dnub = grid_size


        answer = csu.los_disc_intersect(D, incli, alpha, radio, hdis)
        if answer is False:
            return 0, [], []
        else:
            x0, yt1, yt2, zt1, zt2 = answer
        yt = [yt1, yt2]
        zt = [zt1, zt2]
        radio1 = np.sqrt((radio ** 2) - x0 ** 2)

        if len(yt) == 0:
            # print('no entra al disco')
            return (0, [], [])
        else:
            # print(radio)
            yminc = yt[0]
            ymaxc = yt[1]
            zminc = zt[0]
            zmaxc = zt[1]
            # print(yminc, ymaxc)
            dm = dnub / 10.
            dz = zminc
            dy = yminc
            nxy = int(radio1 * 2 / dnub)

            nh = int(hdis / dnub)
            # print('nh', nh)
            xy = np.linspace(-radio1, radio1, nxy + 1)
            h = np.linspace(-hdis / 2, hdis / 2, nh + 1)
            # print(h)
            dlos = np.sqrt(((ymaxc - yminc) ** 2) + (zmaxc - zminc) ** 2)
            # print(dlos)
            # print(dlos/dm)

            # print(dz + len(los)*dm*sin(incli))
            n = 0
            ngrill = 0
            velos = []
            radios = []
            ygrillmint = 0
            zgrillmint = 0

            # print('aaa',dy, xy)
            # print(dz, h)

            for i in range(int(dlos / dm)):
                dz = dz + (dm * np.cos(self.incl_rad))
                dy = dy + (dm * np.sin(self.incl_rad))

                for j in range(len(xy) - 1):
                    if xy[j] <= dy < xy[j + 1]:
                        ygrillmin = xy[j]
                        ygrillmax = xy[j + 1]
                    else:
                        pass
                        # ygrillmin = dy
                        # ygrillmax = dy
                for k in range(len(h)):
                    if h[k] < dz < h[k + 1]:
                        zgrillmin = h[k]
                        zgrillmax = h[k + 1]
                    else:
                        pass
                        # zgrillmin = dz
                        # zgrillmax = dz
                # print('y', ygrillmin, ygrillmax)
                # print('z', zgrillmin, zgrillmax)

                if ygrillmin == ygrillmint and zgrillmin == zgrillmint:
                    pass
                elif ygrillmin != ygrillmint:
                    # print(ygrillmint, ygrillmin)
                    ygrillmint = ygrillmin
                    ngrill = ngrill + 1
                    rc = np.sqrt((((ygrillmax + ygrillmin) / 2) ** 2) + x0 ** 2)
                    yp = (ygrillmax + ygrillmin) / 2
                    # print(yp)
                    veli = self.los_vel(D, alpha, yp)
                    # print(veli)

                    prob = self.prob_hit(rc)

                    selec = np.random.uniform(0, 100)
                    if selec < prob:
                        n = n + 1
                        velos.append(veli)
                        radios.append(rc)
                    else:
                        pass
                elif zgrillmin != zgrillmint:
                    # print(zgrillmint, zgrillmin)
                    zgrillmint = zgrillmin
                    ngrill = ngrill + 1

                    # print(veli)

                    prob = self.prob_hit(rc)

                    selec = np.random.uniform(0, 100)
                    if selec < prob:
                        rc = np.sqrt((((ygrillmax + ygrillmin) / 2) ** 2) + x0 ** 2)
                        yp = (ygrillmax + ygrillmin) / 2
                        # print(yp)
                        veli = self.los_vel(D, alpha, yp)
                        n = n + 1
                        velos.append(veli)
                        radios.append(rc)

                        # print(1,dy,dz)
            print(ngrill)
            return (n, velos, radios)

    def losspec(self, D, alp, lam):

        Ns = self.get_clouds(D, alp)
        print(Ns)
        taus = []
        if Ns[0] == 0:
            return ([], [])
        else:
            for i in range(len(Ns[1])):
                vel = Ns[1][i]
                tau = csu.Tau(lam, vel)
                taus.append(tau)
            # print(taus[0])
            vele = (const.c.to('km/s').value * ((lam / (2852 * (1 + 0.7))) - 1))
            # print(len(taus))
            # print(len(taus[0]))
            sumataus = csu.sumtau(taus)
            flux = csu.normflux(sumataus)
            return (vele, flux, Ns[0])

    def plotspecandelipse(self, impar, alf, lam):

        flux = self.losspec(impar, alf, lam)
        # import pdb; pdb.set_trace()

        fig = plt.figure(figsize=(15, 5))
        grid = plt.GridSpec(1, 3, wspace=0.4, hspace=0.3)
        elipse = fig.add_subplot(grid[0, 0])
        spectro = fig.add_subplot(grid[0, 1:])
        b = self.R* 2 * np.cos(self.incl_rad)
        e1 = patches.Ellipse((0, 0), self.R * 2, b, alpha=0.5)
        elipse.axis('equal')
        elipse.add_patch(e1)
        # elipse.set_xlim((-radio)-1,(radio)+1)
        # elipse.set_ylim((-radio)-1,(radio)+1)
        x = impar * np.cos(np.radians(alf))
        y = impar * np.sin(np.radians(alf))
        elipse.plot(x, y, 'r*')
        eyi = -b / 2
        eyf = b / 2
        exi = -self.R
        exf = self.R
        elipse.plot((0, 0), (eyi, eyf), 'k--')
        elipse.plot((exi, exf), (0, 0), 'k--')
        elipse.set_title('incl:%s,' % self.incl + ' D:%s,' % impar + ' alf:%s' % alf)
        elipse.set_ylabel('kpc')
        elipse.set_xlabel('kpc')
        spectro.plot(flux[0], flux[1])
        spectro.set_title('Number of cluds:%s' % flux[2])
        spectro.set_ylim(-0.05, 1.05)
        spectro.set_xlim(-300, 300)
        spectro.axvline(ls='--', lw=1)
        spectro.axhline(ls='--', lw=1)
        spectro.set_xlabel('LOS vel [km/s]')
        spectro.set_ylabel('Norm. Flux')
        plt.show()

    def spec_clouds(self):
        pass
