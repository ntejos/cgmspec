import disco
import pylab as plt
from math import pi, sqrt, sin, cos, radians
from matplotlib import patches

def plotlos(flux):
    plt.plot(flux[0], flux[1])
    plt.ylim(-0.05, 1.05)
    plt.show()

def plotallclouds(fluxes):
    for i in range(len(fluxes[1])):
        plt.plot(fluxes[0], fluxes[1][i], 'b')
        plt.ylim(-0.05, 1.05)
    plt.show()

def plotelipse(inc, radio):
    b = radio * cos(radians(inc))
    fig = plt.figure()
    ax = plt.axes()
    e1 = patches.Ellipse((0,0), radio, b)
    ax.add_patch(e1)
    plt.ylim((-radio/2)-1,(radio/2)+1)
    plt.xlim((-radio/2)-1,(radio/2)+1)

    plt.show()

def plotspecandelipse(impar, alf, inc, radio,h, lam):
    flux = disco.losspec(impar,alf, inc,radio,h,lam)
    fig = plt.figure(figsize=(15,5))
    grid = plt.GridSpec(1, 3, wspace=0.4, hspace=0.3)
    elipse =fig.add_subplot(grid[0, 0])
    spectro = fig.add_subplot(grid[0, 1:])
    b = radio*2 * cos(radians(inc))
    e1 = patches.Ellipse((0,0), radio*2, b, alpha=0.5)
    elipse.axis('equal')
    elipse.add_patch(e1)
    #elipse.set_xlim((-radio)-1,(radio)+1)
    #elipse.set_ylim((-radio)-1,(radio)+1)
    x = impar * cos(radians(81))
    y = impar * sin(radians(81))
    elipse.plot(x,y,'r*')
    eyi = -b/2
    eyf = b/2
    exi = -radio
    exf = radio
    elipse.plot((0,0),(eyi,eyf),'k--')
    elipse.plot((exi,exf),(0,0),'k--')
    elipse.set_title('incl:%s,'%inc + ' D:%s,'%impar + ' alf:%s' %alf)
    elipse.set_ylabel('kpc')
    elipse.set_xlabel('kpc')
    spectro.plot(flux[0],flux[1])
    spectro.set_title('Number of cluds:%s'%flux[2])
    spectro.set_ylim(-0.05, 1.05)
    spectro.set_xlim(-300,300)
    spectro.axvline(ls='--',lw=1)
    spectro.axhline(ls='--',lw=1)
    spectro.set_xlabel('LOS vel [km/s]')
    spectro.set_ylabel('Norm. Flux')
    plt.show()


lam1 = disco.lamda(4000, 5500, 0.1)
plotspecandelipse(16.8,99,75,100,10,lam1)
