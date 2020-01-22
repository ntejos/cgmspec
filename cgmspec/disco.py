import numpy as np
from cgmspec import utils as csu
from matplotlib import pyplot as plt
from matplotlib import patches
from astropy import constants as const

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
        self.incl_rad = np.radians(self.incl)

    def prob_hit(self, r):
        """
        :param r: float, distance to the center of disc in plane xy in kpc
        :return: float, covering fraction at distance r
        """

        return csu.prob_hit(r, self.R)

    def los_vel(self, D, al, y, vR=-180):
        """
        :param D: float, impact parameter in kpc
        :param al: float, angle between the major axis and the line-of-sight, clockwise, in degrees
        :param y: distance of the line-of-sight to the semi-major axis along the y-axis
        :param vR: maximum velocity of rotation of the disk in km/s
        :return: line-of-sight velocity for a cloud in the given position
        """

        al_rad = np.radians(al)

        x0 = (D*1000) * np.cos(al_rad)
        y0 = (D*1000) * np.sin(al_rad) / np.cos(self.incl_rad)
        a = np.sin(self.incl_rad) / np.sqrt(1 + (y/x0)**2)
        vr = vR*a
        return(vr)

    def numnubestrue(self, D, al, incli, radio, hdis, dnub=1000):

        incli = self.incl
        radio = self.R
        hdis = self.h

        # print('start numnb')
        radio = radio * 1000
        al_rad = np.radians(al)
        hdis = hdis * 1000
        x0 = (D * 1000) * np.cos(al_rad)

        if x0 > radio:
            yt = []
            zt = []
        else:
            # print(incli)
            if incli == 90:
                z0 = D * 1000 * np.sin(al_rad)
                #   print(z0, hdis/2)
                if abs(z0) > hdis / 2:
                    yt = []
                    zt = []
                else:
                    yt = [-np.sqrt((radio ** 2) - x0 ** 2), np.sqrt((radio ** 2) - x0 ** 2)]
                    zt = [D * 1000 * np.sin(al_rad), D * 1000 * np.sin(al_rad)]
                    #       print(yt,zt)
                incli_rad = np.radians(incli)
                radio1 = radio
            elif incli == 0.0:
                yt = [D * 1000 * np.sin(al_rad), D * 1000 * np.sin(al_rad)]
                zt = [-hdis / 2, hdis / 2]
                incli_rad = np.radians(incli)
                radio1 = radio
                # print(yt, zt)
            else:
                incli_rad = np.radians(incli)
                y0 = (D * 1000) * np.sin(al_rad) / np.cos(incli_rad)
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
            dm = dnub / 10
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
            los = np.linspace(0, dlos, int(dlos / dm))

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
                dz = dz + (dm * np.cos(incli_rad))
                dy = dy + (dm * np.sin(incli_rad))

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
                    veli = self.los_vel(D, al, yp)
                    # print(veli)

                    prob = csu.prob_hit(rc, 1000*self.R)
                    # import pdb; pdb.set_trace()

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
                        veli = self.los_vel(D, al, yp)
                        n = n + 1
                        velos.append(veli)
                        radios.append(rc)

                        # print(1,dy,dz)
            print(ngrill)
            return (n, velos, radios)

    def losspec(self, D, alp, incli, rad, h, lam):

        Ns = self.numnubestrue(D, alp, incli, rad, h)
        import pdb; pdb.set_trace()
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

    def plotspecandelipse(self, impar, alf, inc, radio, h, lam):

        flux = self.losspec(impar, alf, inc, radio, h, lam)
        # import pdb; pdb.set_trace()

        fig = plt.figure(figsize=(15, 5))
        grid = plt.GridSpec(1, 3, wspace=0.4, hspace=0.3)
        elipse = fig.add_subplot(grid[0, 0])
        spectro = fig.add_subplot(grid[0, 1:])
        b = radio * 2 * np.cos(self.incl_rad)
        e1 = patches.Ellipse((0, 0), radio * 2, b, alpha=0.5)
        elipse.axis('equal')
        elipse.add_patch(e1)
        # elipse.set_xlim((-radio)-1,(radio)+1)
        # elipse.set_ylim((-radio)-1,(radio)+1)
        x = impar * np.cos(np.radians(alf))
        y = impar * np.sin(np.radians(alf))
        elipse.plot(x, y, 'r*')
        eyi = -b / 2
        eyf = b / 2
        exi = -radio
        exf = radio
        elipse.plot((0, 0), (eyi, eyf), 'k--')
        elipse.plot((exi, exf), (0, 0), 'k--')
        elipse.set_title('incl:%s,' % inc + ' D:%s,' % impar + ' alf:%s' % alf)
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
