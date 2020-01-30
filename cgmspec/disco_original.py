import numpy as np
from math import pi, sqrt, sin, cos, radians, tan, asin, exp, log
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from random import randrange
import tqdm
import random
from numpy.linalg import norm
from scipy import special
from astropy import constants as const
from scipy.special import wofz


#Parametros para el modelo

dnub = 1000 #diametro de las nubes en pc
#hdis = 5000 #altura del disco pc



#nxy = radis/dnub
#nh = hdis/dnub

#xy = np.linspace(-radis, radis, nxy)
#h = np.linspace(-hdis/2, hdis/2, nh)

#volnub = (4/3)*pi*(dnub/2)**3

#De = np.load('r_mage_kpc.npy')

#Velocidades

hv = 5000 #kpc  -----> parametro libre
vR = 200 #km/s ----> parametro libre

#velobs = np.load('vel1.npy')

def lamda(lami, lamf, sampling):
    dlam = (lamf - lami)/sampling
    return(np.linspace(lami, lamf, dlam))

def probabilidad(r,rmax):
        A = 100
        b = log(20/A,rmax)
        return(A*(r**b))

def vel(D, incl, al, y):

    print('start vel')
    x0 = (D*1000)*cos(al)
    y0 = (D*1000)*sin(al)/cos(incl)
    print(sin(incl), y/x0)
    a = sin(incl)/sqrt(1+(y/x0)**2)
    b = exp(-abs(y-y0)/hv*tan(incl))

    vr = vR*a*b
    print('end vel')
    return(vr)

def numnubestrue(D, al, incli,radio, hdis):
    #print('start numnb')
    radio = radio*1000
    al  = radians(al)
    hdis = hdis*1000
    x0 = (D*1000)*cos(al)
    if x0>radio:
        yt = []
        zt = []
    else:
        #print(incli)
        if incli == 90:
           z0 = D*1000*sin(al)
        #   print(z0, hdis/2)
           if abs(z0)>hdis/2:
               yt = []
               zt = []
           else:
               yt = [-sqrt((radio**2)-x0**2),sqrt((radio**2)-x0**2)]
               zt = [D*1000*sin(al), D*1000*sin(al)]
        #       print(yt,zt)
           incli = radians(incli)
           radio1 = radio
        elif incli == 0.0:
           yt = [D*1000*sin(al), D*1000*sin(al)]
           zt = [-hdis/2, hdis/2]
           incli = radians(incli)
           radio1 = radio
           #print(yt, zt)
        else:
           incli = radians(incli)
           y0 = (D*1000)*sin(al)/cos(incli)
    #print(y0)

           z0 = 0

           radio1 = sqrt((radio**2)-x0**2)
          # print(x0,y0)
           ymin = y0 - (hdis/2)*tan(incli)
           ymax = y0 + (hdis/2)*tan(incli)
           zmin = -(sqrt((-radio1-y0)**2)/tan(incli))
           zmax = (sqrt((radio1-y0)**2)/tan(incli))
           ys = [ymin, ymax, -radio1, radio1]
           zs = [-hdis/2, hdis/2, zmin, zmax]
           yt = []
           zt = []
         #  print(ymin, ymax, zmin, zmax)

           for i in range(len(ys)):
              if abs(ys[i])<radio1:
        #           print(ys[i], zs[i])
                   yt.append(ys[i])
                   zt.append(zs[i])
              else:
                   pass
           for i in range(len(zs)):
               if abs(zs[i])<hdis/2:
        #           print(ys[i], zs[i])
                   yt.append(ys[i])
                   zt.append(zs[i])
               else:
                   pass
           yt = np.sort(yt)
           zt = np.sort(zt)

         #  print(yt,zt)

    if len(yt)==0:
        #print('no entra al disco')
        return (0,[], [])
    else:
        #print(radio)
        yminc = yt[0]
        ymaxc = yt[1]
        zminc = zt[0]
        zmaxc = zt[1]
        #print(yminc, ymaxc)
        dm = dnub/10
        dz = zminc
        dy = yminc
        nxy = int(radio1*2/dnub)

        nh = int(hdis/dnub)
        #print('nh', nh)
        xy = np.linspace(-radio1, radio1, nxy+1)
        h = np.linspace(-hdis/2, hdis/2, nh+1)
        #print(h)
        dlos = sqrt(((ymaxc-yminc)**2)+(zmaxc-zminc)**2)
        #print(dlos)
        #print(dlos/dm)
        los = np.linspace(0, dlos, int(dlos/dm))

        #print(dz + len(los)*dm*sin(incli))
        n = 0
        ngrill = 0
        velos=[]
        radios=[]
        ygrillmint = 0
        zgrillmint = 0

        #print('aaa',dy, xy)
        #print(dz, h)
        for i in range(int(dlos/dm)):
            dz = dz +(dm*cos(incli))
            dy = dy +(dm*sin(incli))

            for j in range(len(xy)-1):
                if xy[j]<=dy<xy[j+1]:
                    ygrillmin = xy[j]
                    ygrillmax = xy[j+1]
                else:
                    pass
                    #ygrillmin = dy
                    #ygrillmax = dy
            for k in range(len(h)):
                if h[k]<dz<h[k+1]:
                    zgrillmin = h[k]
                    zgrillmax = h[k+1]
                else:
                    pass
                    #zgrillmin = dz
                    #zgrillmax = dz
            #print('y', ygrillmin, ygrillmax)
            #print('z', zgrillmin, zgrillmax)

            if ygrillmin == ygrillmint and zgrillmin == zgrillmint:
                pass
            elif ygrillmin != ygrillmint:
                #print(ygrillmint, ygrillmin)
                ygrillmint = ygrillmin
                ngrill = ngrill+1
                rc = sqrt((((ygrillmax+ygrillmin)/2)**2)+x0**2)
                yp = (ygrillmax+ygrillmin)/2
                #print(yp)
                veli = vel(D,incli,al,yp)
                #print(veli)

                prob = probabilidad(rc,radio)

                selec = random.uniform(0, 100)
                if selec<prob:
                   n = n+1
                   velos.append(veli)
                   radios.append(rc)
                else:
                    pass
            elif zgrillmin != zgrillmint:
                #print(zgrillmint, zgrillmin)
                zgrillmint = zgrillmin
                ngrill = ngrill+1

                #print(veli)

                prob = probabilidad(rc, radio)

                selec = random.uniform(0, 100)
                if selec<prob:
                   rc = sqrt((((ygrillmax+ygrillmin)/2)**2)+x0**2)
                   yp = (ygrillmax+ygrillmin)/2
                   #print(yp)
                   veli = vel(D,incli,al,yp)
                   n = n+1
                   velos.append(veli)
                   radios.append(rc)

            #print(1,dy,dz)
        print(ngrill)
        return(n, velos, radios)


'''def numnubestrue(D, al, incli,radio):
    radio = radio*1000
    incli = radians(incli)
    oincli = radians(90-incli)
    al  = radians(al)
    y0 = (D*1000)*sin(al)/cos(incli)
    #print(y0)
    x0 = (D*1000)*cos(al)
    z0 = 0
    ymin = y0 - (hdis/2)*tan(incli)
    ymax = y0 + (hdis/2)*tan(incli)
    #print(ymin, ymax)
    zmin = -hdis/2
    zmax = hdis/2
    dm = dnub/10
    dz = zmin
    dy = ymin
    nxy = radio/dnub
    xy = np.linspace(-radio, radio, nxy)
    dlos = sqrt(((ymax-ymin)**2)+(zmax-zmin)**2)
    #print('dlos',dlos)
    los = np.linspace(0,  dlos, dlos/dm)
    n = 0
    ngrill = 0
    velos=[]
    radios=[]

    for i in range(len(los)):
        if abs(dz)>hdis/2 or abs(dy)>radio or abs(x0)>radio:
            dz = dz +(dm*cos(incli))
            dy = dy +(dm*sin(incli))
            print(0,dy,dz)



        else:
            ygrill = min(xy, key=lambda x:abs(x-dy))
            zgrill = min(h, key=lambda x:abs(x-dz))
            if ygrill<dy:
                ygrillmin = ygrill
                ygrillmax = ygrill + dnub
                print(ygrillmin,ygrillmax)

            else:
                ygrillmin = ygrill - dnub
                ygrillmax = ygrill
                print(ygrillmin,ygrillmax)

            if zgrill<dz:
                zgrillmin = zgrill
                zgrillmax = zgrill + dnub
                print(zgrillmin,zgrillmax)
            else:
                zgrillmin = zgrill - dnub
                zgrillmax = zgrill
                print(zgrillmin,zgrillmax)

            if (abs(dy - ygrillmin))<(dm*cos(oincli)) or (abs(dz - zgrillmin))< (dm*sin(oincli)):
                ngrill = ngrill +1
                rc = sqrt((((ygrillmax+ygrillmin)/2)**2)+x0**2)
                yp = (ygrillmax+ygrillmin)/2
                #print(yp)
                veli = vel(D,incli,al,yp)
                #print(veli)

                prob = probabilidad(rc)

                selec = random.uniform(0, 100)
                if selec<prob:
                   n = n+1
                   velos.append(veli)
                   radios.append(rc)

                else:
                    pass



            else:
                pass

            dz = dz +(dm*sin(oincli))
            dy = dy +(dm*cos(oincli))
            print(1,dy,dz)

    #print(velos)

    b = radis * cos(incli)
    x = D*1000 * cos(al)
    y = - D*1000 * sin(al)
    fig = plt.figure()
    ax = plt.axes()
    e1 = patches.Ellipse((0,0), radis*2, b*2, alpha=0.5)
    ax.scatter(x,y)
    ax.add_patch(e1)
    plt.plot()
    plt.ylim(-100000, 100000)
    plt.xlim(-100000, 100000)

    plt.show()
    print(ngrill)
    return(n, velos, radios)

def numnubestrue(impactparam):
    D1 = impactparam*1000
    xc = D1 * cos(alpha)
    y0 = D1 * sin(alpha)/cos(inc)
    z = np.linspace(-hdis/2, hdis/2, nh+1)
    n = 0
    velos = []
    radios = []
    for i in range(len(z)):
        zc = z[i]
        if zc<=0:
            yc = (-1) * (((zc+1000)/(2*tan(inc)))-y0)
        else:
            yc = (((zc+1000)/(2*tan(inc)))+y0)

        rc = sqrt(xc**2 + yc**2)

        if rc < radis:
           veli = vel(rc, impactparam)
           prob = probabilidad(rc)

           selec = random.uniform(0, 100)
           if selec<prob:
               n = n+1
               velos.append(veli)
               radios.append(rc)
           else:
               pass
        else:
            pass
    return(n,velos,radios)'''

#calcular tau



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
    y = gamma/(4*pi*dnud)
    z = x + 1j*y
    v = np.asarray(np.real(wofz(z)/(sqrt(pi)*dnud)))

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

def losspec(D, alp, incli, rad, h, lam):
    Ns = numnubestrue(D, alp, incli, rad, h)
    print(Ns)
    taus = []
    if Ns[0] == 0:
        return ([],[])
    else:

        for i in range(len(Ns[1])):
            vel = Ns[1][i]
            tau = Tau(lam,vel)
            taus.append(tau)
        #print(taus[0])
        vele = (const.c.to('km/s').value * ((lam / (2852 *(1 + 0.7)))-1))
        #print(len(taus))
        #print(len(taus[0]))
        sumataus = sumtau(taus)
        flux = normflux(sumataus)
        return(vele, flux, Ns[0])


from astropy.convolution import convolve, Gaussian1DKernel

def filtrogauss(fw, tau):
    gauss_kernel = Gaussian1DKernel(fw)
    gausflux = convolve(tau[1], gauss_kernel)
    return(tau[0],gausflux)


#lam1 = lamda(4800, 5000, 0.1)
#caca = losspec(10, 14, 45, 50, 10, lam1)
#plt.plot(caca[0],caca[1])
#plt.show()
