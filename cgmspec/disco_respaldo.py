import numpy as np
from cgmspec import utils as csu
from matplotlib import pyplot as plt
from matplotlib import patches
from astropy import constants as const
from matplotlib.pyplot import cm

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

    def prob_hit(self, r, r_0):
        """
        :param r: float, distance to the center of disc in plane xy in kpc
        :return: float, covering fraction at distance r
        """

        return csu.prob_hit(r, self.Rcore, r_0)

    def los_vel(self, y, D, alpha, v_inf,vR=196, hv=10):
        """
        See Ho et al. 2017

        :param y: np.array, distance of the line-of-sight to the semi-major axis (x-axis) along the y-axis
        :param D: float, impact parameter in kpc
        :param alpha: float, angle between the major axis and the line-of-sight, clockwise, in degrees
        :param vR: maximum velocity of rotation of the disk in km/s
        :param hv: velocity scale height in kpc
        :return: line-of-sight velocity for a cloud in the given position

        """

        v_los_inf = (v_inf * np.sin(self.incl_rad)) * (y/(np.sqrt((y**2) + D**2)))
        al_rad = np.radians(alpha)

        R = D * np.sqrt(1+(np.sin(al_rad)**2)*np.tan(self.incl_rad)**2)
        vrot = (2/np.pi)*np.arctan2(R,1)

        x0 = D * np.cos(al_rad)  # this is p in Ho et al. 2019, fig 10.
        y0 = D * np.sin(al_rad) / np.cos(self.incl_rad)  # this is y0 in the same fig.
        if x0>=0:
              a = np.sin(self.incl_rad) / np.sqrt(1 + (y/x0)**2)
        else:
            a = np.sin(self.incl_rad) / np.sqrt(1 + (y/x0)**2)

        b = np.exp(-np.fabs(y - y0) / hv * np.tan(self.incl_rad))
        #print(b)
        vr = (vR*vrot*a*b) + v_los_inf


        return(vr)

    def get_clouds(self, D, alpha, grid_size, r_0, losinter, v_max, h_v, v_inf):
        """
        Obtain main parameters to describe clouds intersected by the sightline along the disc

        :param D: float, impact parameter in kpc
        :param alpha: float, angle between the major axis and the line-of-sight, clockwise, in degrees
        :param grid_size: size of grid in kpc to model a hit or no hit of a cloud
        :return: (nclouds, vel_array, r_array) : (float, array, array)
        """

        incli = self.incl
        radio = self.R
        hdis = self.h
        dnub = grid_size


        #answer = csu.los_disc_intersect(D, incli, alpha, radio, hdis)
        answer = losinter
        if answer is False:
            return 0, [], []
        else:
            x0, yt1, yt2, zt1, zt2 = answer
        yt = [yt1, yt2]
        zt = [zt1, zt2]
        radio1 = np.sqrt((radio ** 2) - x0 ** 2)
        #print('aaaaa', yt, zt)

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
            #ygrillmin = 'a'
            #zgrillmin = 'a'
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
                for k in range(len(h)-1):
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
                    veli = self.los_vel(yp, D, alpha, v_inf, v_max, h_v)
                    # print(veli)

                    prob = self.prob_hit(rc, r_0 )

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
                    rc = np.sqrt((((ygrillmax + ygrillmin) / 2) ** 2) + x0 ** 2)
                    # print(veli)

                    prob = self.prob_hit(rc, r_0)

                    selec = np.random.uniform(0, 100)
                    if selec < prob:
                        yp = (ygrillmax + ygrillmin) / 2
                        # print(yp)
                        veli = self.los_vel(yp, D, alpha, v_inf, v_max, h_v)
                        n = n + 1
                        velos.append(veli)
                        radios.append(rc)

                #elif zgrillmin == 'a' or ygrillmin == 'a':
                #    n  = 0
                #    velos = []
                #    radios = []
                         # print(1,dy,dz)
            #print(ngrill)
            return (n, velos, radios)

    def losspec(self, D, alpha, lam, grid_size, X=1, z=0.73379):
        """
        Obtain the spectra for a LOS crossing the disk in velocity scale

        :param D: float, impact parameter in kpc
        :param alpha: float, angle between the major axis and the line-of-sight, clockwise, in degrees
        :lam: array, wavelenghts where the spectra is calculated
        :return: (vele, flux, nclouds) : (array, array, float)
        """
        if X == 1:
            lam0 = 2796
        if X ==2:
            lam0 = 2803


        Ns = self.get_clouds(D, alpha, grid_size)
        #print(Ns)
        taus = []
        if Ns[0] == 0:
            vele = (const.c.to('km/s').value * ((lam / (lam0 * (1 + z))) - 1))
            return (vele, np.ones(len(lam)), 0)
        else:
            for i in range(len(Ns[1])):
                vel = Ns[1][i]
                tau = csu.Tau(lam, vel, X)
                taus.append(tau)
            # print(taus[0])
            vele = (const.c.to('km/s').value * ((lam / (lam0 * (1 + z))) - 1))
            # print(len(taus))
            # print(len(taus[0]))
            sumataus = csu.sumtau(taus)
            flux = csu.normflux(sumataus)
            return (vele, flux, Ns[0])


    def averagelosspec(self, D, alpha, lam, iter,X, z, grid_size, N, b, r_0, v_max, h_v, v_inf):
        """
        Obtain the spectra for a LOS crossing the disk in velocity scale

        :param D: float, impact parameter in kpc
        :param alpha: float, angle between the major axis and the line-of-sight, clockwise, in degrees
        :lam: array, wavelenghts where the spectra is calculated
        :iter: numer of iterations to do the average
        :return: (vele, flux, nclouds) : (array, array, float)
        """

        print('incl',self.incl)

        if X == 1:
            lam0 = 2796.35
        if X ==2:
            lam0 = 2803.53
        if X == 12:
            lam0 = 2796.35

        vele1 = (const.c.to('km/s').value * ((lam / (lam0 * (1 + z))) - 1))

        answer = csu.los_disc_intersect(D, self.incl, alpha, self.R, self.h)

        if answer is False:
            return (vele1, np.ones(len(vele1)), 0)
        else:
            pass

        flux_to_average1 = []
        flux_to_average2 = []
        nclouds_to_average = []
        vele = []
        average_flux = []






        for i in range(iter):
             Ns = self.get_clouds(D, alpha, grid_size, r_0, answer,v_max, h_v, v_inf)

             taus1 = []

             if Ns[0] == 0:
                 nclouds_to_average.append(0)
                 flux = np.ones(len(lam))
                 flux_to_average1.append(flux)


             else:
                 for i in range(len(Ns[1])):
                     vel = Ns[1][i]
                     tau1 = csu.Tau(lam, vel, X, N, b)
                     taus1.append(tau1)

            # print(taus[0])
            # print(len(taus))
            # print(len(taus[0]))
                 sumataus1 = csu.sumtau(taus1)
                 flux1 = csu.normflux(sumataus1)
                 flux_to_average1.append(flux1)
                 nclouds_to_average.append(Ns[0])


        vele.append(vele1)
        vele = np.asarray(vele)
        flux_to_average1 = np.asarray(flux_to_average1)
        #print('aa', flux_to_average1)
        nclouds_to_average = np.asarray(nclouds_to_average)
        average_flux1 = np.median(flux_to_average1, axis=0)
    #    print('bb', len(average_flux1), len(average_flux2))
        average_nclouds = np.median(nclouds_to_average)
        return (vele1, average_flux1, average_nclouds)


    def plotspecandelipse(self, D, alpha, lam, iter,X, z, grid_size, N, b, prob_rmax, prob_rmin):
        """
        Returns a figure with two plots. On the left is the disk proyected in the sky plane and the LOS position, in the right the spectra generated.

        :param D: float, impact parameter in kpc
        :param alpha: float, angle between the major axis and the line-of-sight, clockwise, in degrees
        :lam: array, wavelenghts where the spectra is calculated
        :return: Plot
        """

        flux = self.averagelosspec(D, alpha, lam, iter,X, z, grid_size, N, b, prob_rmax, prob_rmin)[1]

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
        x = D * np.cos(np.radians(alpha))
        y = -D * np.sin(np.radians(alpha))
        elipse.plot(x, y, 'r*')
        eyi = -b / 2
        eyf = b / 2
        exi = -self.R
        exf = self.R
        elipse.plot((0, 0), (eyi, eyf), 'k--')
        elipse.plot((exi, exf), (0, 0), 'k--')
        elipse.set_title('incl:%s,' % self.incl + ' D:%s,' % D + ' alf:%s' % alpha)
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

    def plotmanylos(self, D, alpha, lam):
        """
        Returns a figure with many LOS in the disk
        :param D: array, impact parameters in kpc
        :param alpha: array, angles between the major axis and the line-of-sight, clockwise, in degrees
        :lam: array, wavelenghts where the spectra is calculated
        :return: Plot
        """

        fluxes = []
        for i in range(len(D)):
            flux = self.losspec(D[i], alpha[i],lam)
            fluxes.append(flux)

        fig = plt.figure(figsize=(15, 5))
        grid = plt.GridSpec(1, 3, wspace=0.4, hspace=0.3)
        elipse = fig.add_subplot(grid[0, 0])
        spectro = fig.add_subplot(grid[0, 1:])
        b = self.R* 2 * np.cos(self.incl_rad)
        e1 = patches.Ellipse((0, 0), self.R * 2, b, alpha=0.5)
        elipse.axis('equal')
        elipse.add_patch(e1)
        color=cm.rainbow(np.linspace(0,1,len(fluxes)))

        for i in range(len(fluxes)):
            x = D[i] * np.cos(np.radians(alpha[i]))
            y = -D[i] * np.sin(np.radians(alpha[i]))
            elipse.plot(x, y, color=color[i], marker='*')
            spectro.plot(fluxes[i][0], fluxes[i][1], color=color[i])
        eyi = -b / 2
        eyf = b / 2
        exi = -self.R
        exf = self.R
        elipse.plot((0, 0), (eyi, eyf), 'k--')
        elipse.plot((exi, exf), (0, 0), 'k--')
        elipse.set_title('incl:%s,' % self.incl + ' D:%s,' % D + ' alf:%s' % alpha)
        elipse.set_ylabel('kpc')
        elipse.set_xlabel('kpc')
        spectro.set_ylim(-0.05, 1.05)
        spectro.set_xlim(-300, 300)
        spectro.axvline(ls='--', lw=1)
        spectro.axhline(ls='--', lw=1)
        spectro.set_xlabel('LOS vel [km/s]')
        spectro.set_ylabel('Norm. Flux')
        plt.show()


    def plot_aver_spax(self, D, alpha, lam, sizex, sizey):

        alpha_rad = np.radians(alpha)
        xc = D*np.sin(alpha_rad)
        yc = -D*np.cos(alpha_rad)
        p0 = (xc, yc)
        p1 = (xc, yc+sizey)
        p2 = (xc+sizex, yc)
        p3 = (xc, yc-sizey)
        p4 =  (xc-sizex, yc)
        ps = np.asarray((p0,p1,p2,p3,p4))
        #print(ps)

        Ds =  []
        alphas =  []
        #print('ln',len(ps))
        for i in range(len(ps)):
            #print('jajaja', ps[i][0], ps[i][1])
            #print(ps[i][1]/ps[i][0])
            alpha = -np.arctan2(ps[i][1],ps[i][0])
            #print('al',alpha)
            D = ps[i][0]/np.cos(alpha)
            Ds.append(D)
            alphas.append(alpha)
        self.plotmanylos(Ds, alphas,lam)

    def plotelipse(self,D, alpha):

        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        x = D * np.cos(np.radians(alpha))
        y = -D * np.sin(np.radians(alpha))
        plt.plot(x, y, 'r*')
        b = self.R* 2 * np.cos(self.incl_rad)
        e1 = patches.Ellipse((0, 0), self.R * 2, b, alpha=0.5)
        eyi = -b / 2
        eyf = b / 2
        exi = -self.R
        exf = self.R
        plt.plot((0, 0), (eyi, eyf), 'k--')
        plt.plot((exi, exf), (0, 0), 'k--')
        plt.title('incl:%s,' % self.incl + ' D:%s,' % D + ' alf:%s' % alpha)
        plt.ylabel('kpc')
        plt.xlabel('kpc')
        plt.axis('equal')
        ax.add_patch(e1)
        return(fig)
