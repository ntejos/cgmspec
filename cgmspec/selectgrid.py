
import numpy as np
import pylab as plt
from timeit import default_timer as timer
import math
from math import log
from astropy import constants as const
from scipy.special import wofz





def los_vel(incli, y, D, alpha, vR, hv, v_inf=0):
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
    incli=np.radians(incli)
    start = timer()
    v_los_inf = (v_inf * np.sin(incli)) * (y/(np.sqrt((y**2) + D**2)))
    al_rad = np.radians(alpha)

    R = D * np.sqrt(1+(np.sin(al_rad)**2)*np.tan(incli)**2)
    vrot = (2/np.pi)*np.arctan2(R,1)

    x0 = D * np.cos(al_rad)  # this is p in Ho et al. 2019, fig 10.
    y0 = D * np.sin(al_rad) / np.cos(incli)  # this is y0 in the same fig.
    if x0>=0:
          a = np.sin(incli) / np.sqrt(1 + (y/x0)**2)
    else:
        a = np.sin(incli) / np.sqrt(1 + (y/x0)**2)

    b = np.exp(-np.fabs(y - y0) / hv * np.tan(incli))
    #print(b)
    vr = (vR*vrot*a*b) + v_los_inf


    return(vr)

def prob_hit(r, r_0, prob_rmin=100):
    """
    Probability of hitting a cloud at distance r in the plane xy of a disc of radius rmax. For the moment is a power law.

    :param r: np.array, array with distances to the center of disc in plane xy in kpc
    :param r_0: float, the characteristicradius of the power law
    :param prob_rmin: float, probability at Rcore or below, default 100% probability of crossing a cloud
    :return: float, probability of hitting a cloud
    """
    start = timer()
    rmin = 0.1
    ind = np.log(prob_rmin)/(np.log(rmin)-np.log(r_0))
    A = r/r_0
    prob = A**(ind)
    prob = np.where(r<=rmin, prob_rmin, prob)
    print('calcular prob', timer() - start)
    return prob


def get_cells(h,incli,D,alpha,size):
    start = timer()
    m = -np.tan(np.radians(90-incli))
    y0 = D*np.sin(np.radians(alpha))/np.cos(np.radians(incli))
    n = -m*y0
    y1 = ((h/2)-n)/m
    y2 = (-(h/2)-n)/m

    z1 = h/2
    z2 = -h/2
    print('y1,y2', y1, y2)
    t = 0
    y = y1
    z = z1
    loslenght = np.sqrt(((y2-y1)**2)+(z2-z1)**2)
    if y1<0:
        nygrill = int(y1/size)
    else:
        nygrill = int((y1/size)-1)

    nzgrill = int((h/2)/size)

    print('grillas', nygrill, nzgrill)

    cellnum = 0
    cellpos = []



    while t<loslenght:
          print('grillas', nygrill, nzgrill)
          ypos = (size*nygrill)+(size/2)
          zpos = (size*nzgrill)-(size/2)
          print('posicion', ypos,zpos)
          cellpos.append([ypos,zpos])
          cellnum = cellnum+1
          nexty = ((nygrill+1)*size, n+(m*(nygrill+1)*size))
          dnexty = np.sqrt(((nexty[0]-y)**2)+ (nexty[1]-z)**2)
          nextz = ((((nzgrill-1)*size)-n)/m, (nzgrill-1)*size)
          dnextz = np.sqrt(((nextz[0]-y)**2)+ (nextz[1]-z)**2)

          if dnexty < dnextz:
              print(0)
              t = t + dnexty
              y = nexty[0]
              z = nexty[1]
              nygrill = nygrill+1

          else:
              print(1)
              t = t + dnextz
              y = nextz[0]
              z = nextz[1]
              nzgrill = nzgrill-1



    return(cellnum, cellpos)


def dist_line_poin(p1x,p1y,p2x,p2y,p3x,p3y):
    p1=np.array([p1x,p1y])
    p2=np.array([p2x,p2y])
    p3=np.array([p3x,p3y])
    return(np.cross(p2-p1,p3-p1)/np.linalg.norm(p2-p1))


def get_cells2(h,incli,D,alpha,size,r_0, vR,hv):
    start = timer()
    m = -np.tan(np.radians(90-incli))

    x0 = D * np.cos(np.radians(alpha))
    y0 = D*np.sin(np.radians(alpha))/np.cos(np.radians(incli))
    n = -m*y0

    y1 = ((h/2)-n)/m
    y2 = (-(h/2)-n)/m

    mindis = np.sqrt(2*(size**2))/2
    z1 = h/2
    z2 = -h/2
    b = -1
    zgrid = np.arange((-h/2) + (size/2), (h/2) + (size/2), size)
    ymin = int(y1/size) * size + (size/2)
    ymax = int(y2/size)*size +(size/2)
    ygrid = np.arange(ymin,ymax,size)
    points = abs((m * ygrid + b * zgrid[:,None] + n)) / (np.sqrt(m * m + b * b))
    selected = points <= mindis
    yv, zv = np.meshgrid(ygrid, zgrid)
    ypos = yv[selected]
    zpos = zv[selected]

    radios = np.sqrt((x0**2)+ypos**2)
    probs = prob_hit(radios,r_0)
    velos = los_vel(incli, ypos, D, alpha, vR, hv)
    return(ypos,zpos, probs, velos)


def get_clouds(ypos,zpos,probs,velos):
    randomnum = np.random.uniform(0, 100, len(probs))
    selected = probs >= randomnum
    return(velos[selected])


def Tau(lam,vel,X,N, b,z=0.73379):
    if X ==1:
        lam0 = [2796.35]
        f = [0.6155]
    if X ==2:
        lam0 = [2803.53]
        f = [0.3054]
    if X ==12:
        lam0 = [2796.35, 2803.53]
        f = [0.6155, 0.3054]


    taus = []
    for i in range(len(lam0)):
        gamma, mass = [2.68e8, 24.305]
        c  = const.c.to('cm/s').value
        sigma0 = 0.0263
        lamc = ((vel[:,None]/const.c.to('km/s').value)+1)*((lam0[i]))
        #print('lamc', lamc)
        nu = c/(lam*1e-8)
        nu0 = c/(lamc*1e-8)

        dnu = nu - (nu0/(1+z))
        dnud = (b*100000)*nu/c

        x = dnu/dnud
        y = gamma/(4*np.pi*dnud)
        zi = x + 1j*y
        v = np.asarray(np.real(wofz(zi)/(np.sqrt(np.pi)*dnud)))

        taut = N * sigma0 * f[i] * v
        taus.append(taut)
        dv = (const.c.to('km/s').value * ((lam / (lamc *(1 + z)))-1))


    taus = np.asarray(taus)
    taust = taus.sum(axis=0)
    return(taust)


def losspec(lam,velos,X,N, b,z=0.73379):
    taus = Tau(lam,velos,X,N, b,z=0.73379)
    tottau = np.sum(taus,axis=0)
    return(np.exp(-tottau))


def averagelos(incli,h, D, alpha, lam, iter,X, z, grid_size, N, b, r_0, v_max, h_v, v_inf):
    start = timer()
    cells = get_cells2(h,incli,D,alpha,grid_size, r_0,v_max,h_v)

    results = [0]*iter
    results = [get_clouds(cells[0],cells[1],cells[2],cells[3]) for x in results]
    results = np.asarray(results)

    fluxes = [0]*iter
    fluxtoaver = [losspec(lam,results[x],X,N,b) for x in fluxes]
    fluxtoaver = np.asarray(fluxtoaver)
    totflux = np.median(fluxtoaver, axis=0)
    print('tiempo', timer() - start)
    return(totflux)







def Tau1(lam,vel,X,N, b,z=0.73379):
    if X ==1:
        lam0 = [2796.35]
        f = [0.6155]
    if X ==2:
        lam0 = [2803.53]
        f = [0.3054]
    if X ==12:
        lam0 = [2796.35, 2803.53]
        f = [0.6155, 0.3054]


    taus = []
    for i in range(len(lam0)):
        gamma, mass = [2.68e8, 24.305]
        c  = const.c.to('cm/s').value
        sigma0 = 0.0263
        lamc = ((vel/const.c.to('km/s').value)+1)*((lam0[i]))
        #print('lamc', lamc)
        nu = c/(lam*1e-8)
        nu0 = c/(lamc*1e-8)

        dnu = nu - (nu0/(1+z))
        dnud = (b*100000)*nu/c

        x = dnu/dnud
        y = gamma/(4*np.pi*dnud)
        zi = x + 1j*y
        v = np.asarray(np.real(wofz(zi)/(np.sqrt(np.pi)*dnud)))

        taut = N * sigma0 * f[i] * v
        taus.append(taut)
        dv = (const.c.to('km/s').value * ((lam / (lamc *(1 + z)))-1))

    #print(taus)
    taus = np.asarray(taus)
    taust = taus.sum(axis=0)
    return(taust,dv)


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
