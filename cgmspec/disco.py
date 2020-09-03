import numpy as np
from cgmspec import utils as csu
from matplotlib import pyplot as plt
from matplotlib import patches
from astropy import constants as const
from matplotlib.pyplot import cm

"""This Class Disco represents the CGM of a galaxy from a disc model"""

class Disco:
    """Represents the CGM of a idealized disc projected in the sky"""

    def __init__(self,h, incl, Rcore=0.1):
        """
        :param h: float, height of disc in kpc
        :param incl: float, inclination angle of disc in degrees
        :param Rcore: float, radius of disk core in kpc, where the probability is maximum
        """

        self.Rcore = Rcore
        self.h = h
        self.incl = incl
        self.incl_rad = np.radians(self.incl)



    def prob_hit(self, r, r_0, prob_rmin=100):
        """
        Probability of hitting a cloud at distance r in the plane xy of a disc of radius rmax. For the moment is a power law.

        :param r: np.array, array with distances to the center of disc in plane xy in kpc
        :param r_0: float, the characteristicradius of the power law
        :param prob_rmin: float, probability at Rcore or below, default 100% probability of crossing a cloud
        :return: float, probability of hitting a cloud
        """
        rmin = self.Rcore
        ind = np.log(prob_rmin)/(np.log(rmin)-np.log(r_0))
        A = r/r_0
        prob = A**(ind)
        prob = np.where(r<=rmin, prob_rmin, prob)
        return prob

    '''def prob_hit(self, r, r_0):
        """
        :param r: float, distance to the center of disc in plane xy in kpc
        :param r_0: float, characteristic
        :return: float, covering fraction at distance r
        """

        return csu.prob_hit(r, self.Rcore, r_0)'''

    def los_vel(self, y, D, alpha, vR, hv, v_inf=0):
        """
        line of sight velocity of a cloud in a disc. See Ho et al. 2017

        :param y: np.array, distance of the line-of-sight to the semi-major axis (x-axis) along the y-axis
        :param D: float, impact parameter in kpc
        :param alpha: float, angle between the major axis and the line-of-sight, clockwise, in degrees
        :param vR: maximum velocity of rotation of the disk in km/s
        :param hv: velocity scale height in kpc
        :param v_inf: the infall velocity, for the moment i'm not using it.
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

    def get_cells(self,D,alpha,size):
        """
        which cells does the line of sight crosses given the parameters.

        :param D: float, impact parameter in kpc
        :param alpha: float, angle between the major axis and the line-of-sight, clockwise, in degrees
        :size: float, size of grid where a cloud could be in kpc, to model a hit or no hit of a cloud
        :return: (number of cells, positions(y,z)): (float, array of tuples)
        """


        h = self.h
        incli = self.incl

        m = -np.tan(np.radians(incli))
        y0 = D*np.sin(np.radians(alpha))/np.cos(np.radians(incli))
        n = -m*y0
        y1 = ((h/2)-n)/m
        y2 = (-(h/2)-n)/m

        z1 = h/2
        z2 = -h/2
        #print('y1,y2', y1, y2)
        t = 0
        y = y1
        z = z1
        loslenght = np.sqrt(((y2-y1)**2)+(z2-z1)**2)
        if y1<0:
            nygrill = int(y1/size)
        else:
            nygrill = int((y1/size)-1)

        nzgrill = int((h/2)/size)
        #print('grillas', nygrill, nzgrill)

        cellnum = 0
        cellpos = []
        ylos = np.linspace(y1,y2, 100)
        zlos = n + m*ylos

        while t<loslenght:
              #print('grillas', nygrill, nzgrill)
              ypos = (size*nygrill)+(size/2)
              zpos = (size*nzgrill)-(size/2)
              #print('posicion', ypos,zpos)
              cellpos.append([ypos,zpos])
              cellnum = cellnum+1
              nexty = ((nygrill+1)*size, n+(m*(nygrill+1)*size))
              dnexty = np.sqrt(((nexty[0]-y)**2)+ (nexty[1]-z)**2)
              nextz = ((((nzgrill-1)*size)-n)/m, (nzgrill-1)*size)
              dnextz = np.sqrt(((nextz[0]-y)**2)+ (nextz[1]-z)**2)

              if dnexty < dnextz:
                #  print(0)
                  t = t + dnexty
                  y = nexty[0]
                  z = nexty[1]
                  nygrill = nygrill+1

              else:
                 # print(1)
                  t = t + dnextz
                  y = nextz[0]
                  z = nextz[1]
                  nzgrill = nzgrill-1

        return(cellnum, cellpos)

    def get_clouds(self, grids, D, alpha, grid_size, r_0, v_max, h_v, v_inf):
        """
        Obtain main parameters to describe clouds intersected by the sightline along the disc

        :param D: float, impact parameter in kpc
        :param alpha: float, angle between the major axis and the line-of-sight, clockwise, in degrees
        :param grid_size: size of grid in kpc to model a hit or no hit of a cloud
        :return: (nclouds, vel_array, r_array) : (float, array, array)
        """

        incli = self.incl
        hdis = self.h

        al_rad = np.radians(alpha)
        x0 = D * np.cos(al_rad)
        n = 0
        ngrill = 0
        velos = []
        radios = []

        #answer = csu.los_disc_intersect(D, incli, alpha, radio, hdis)

        #print('aaaaa', yt, zt)
        for i in range(len(grids[1])):
            rc = np.sqrt((grids[1][i][0] ** 2) + x0 ** 2)
            prob = self.prob_hit(rc, r_0 )
            selec = np.random.uniform(0, 100)

            if selec < prob:
                # print(yp)
                veli = self.los_vel(grids[1][i][0], D, alpha, v_inf, v_max, h_v)
                n = n + 1
                velos.append(veli)
                radios.append(rc)
            else:
                pass

        return(n, velos, radios)


    def averagelosspec(self, D, alpha, lam, iter,X, z, grid_size, N, b, r_0, v_max, h_v, v_inf):
        """
        Obtain the spectra for a LOS crossing the disk in velocity scale

        :param D: float, impact parameter in kpc
        :param alpha: float, angle between the major axis and the line-of-sight, clockwise, in degrees
        :lam: array, wavelenghts where the spectra is calculated
        :iter: numer of iterations to do the average
        :param X: can be 1,2 or 12, depending if you want MgII 2796.35 , MII 2803.53 or both transitions
        :param z: float, the redshift of  the absorbing system
        :param N: column density of the absorption produced by one cloud
        :param b: doppler parameter of the absorption
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

        grids = self.get_cells(D,alpha,grid_size)
        #print(grids)
        flux_to_average1 = []
        flux_to_average2 = []
        nclouds_to_average = []
        vele = []
        average_flux = []

        for i in range(iter):

            Ns = self.get_clouds(grids, D, alpha, grid_size, r_0, v_max, h_v, v_inf)

            #print(Ns)

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
