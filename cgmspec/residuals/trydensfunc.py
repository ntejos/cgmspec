import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from mpdaf.obj import Cube, WCS, WaveCoord, Spectrum
from cgmspec.disco import Disco
from astropy.coordinates import SkyCoord
from astropy.cosmology import FlatLambdaCDM
from mpdaf.obj import iter_spe
from astropy.convolution import convolve, Gaussian1DKernel
from scipy.stats import  moment
from astropy import constants as const

lam1 = np.arange(4825.12, 4886.37+1.25, 1.25)
pixs = [(8,15), (12,12), (15,10), (18,6), (20,4), (23,2)]

z = 0.73379
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
dist = cosmo.luminosity_distance(z)
galPA = 55

galcen = SkyCoord(ra=237.52059628552735*u.degree, dec=-78.188149277705151*u.degree, frame='icrs')

datacube  = Cube('ccdatatrue.fits')

def modelcoord(galcen, pixl):
            coord_sky = datacube.wcs.pix2sky([pixl[1],pixl[0]], unit=u.deg)
            dec = coord_sky[0][0]
            ra = coord_sky[0][1]
            scale = 7.28 # z=0.73379  # plat scale (kpc/") for Planck
            c1 = SkyCoord(ra*u.degree,dec*u.degree,frame='icrs')
            D = scale * galcen.separation(c1).arcsec
            pa = galcen.position_angle(c1).to(u.deg)
            alpha = galPA - pa.value
            return(D, alpha)

'''def eqw(spe):
    flux = spe.data[3:46]
    tosum = []
    for i in range(len(flux)):
       if flux[i] < 0.9:
           sw = 1.25*(1-flux[i])
           tosum.append(sw)
       else:
           pass
    tosum = np.array(tosum)

    return(np.sum(tosum))

def centroid(spe):
    if spe[0] == 1.:
        return(0)

    flux = spe.data[3:46]
    vels = []
    weight =  []
    for i in range(len(flux)):
        if flux[i] < 0.9:
            weighti = 1 - flux[i]
            vel = (const.c.to('km/s').value * ((lam1[i] / (2796 * (1 + 0.73379))) - 1))
            vels.append(vel)
            weight.append(weighti)
    aver  = np.average(vels,weights=weight)
    return(aver)'''

def eqw(spe, p1, p2):
    flux = spe[p1:p2]
    tosum = []
    for i in range(len(flux)):
        sw = 1.25*(1-flux[i])
        tosum.append(sw)
    tosum = np.array(tosum)

    return(np.sum(tosum))

def centroid(spe, p1, p2):
    flux = spe[p1:p2]
    vels = []
    weight =  []
    for i in range(len(flux)):
        weighti = 1 - flux[i]
        vel = (const.c.to('km/s').value * ((lam1[i] / (2796 * (1 + 0.73379))) - 1))
        vels.append(vel)
        weight.append(weighti)
    aver  = np.average(vels,weights=weight)
    return(aver)

Ds = []
alphas = []

for i in range(len(pixs)):
    coords = modelcoord(galcen, pixs[i])
    alphai = coords[1]
    Di = coords[0]
    Ds.append(Di)
    alphas.append(alphai)


#plot data  and var
'''for i in range(len(pixs)):
    ydata = datacube[:,pixs[i][1]-1,pixs[i][0]-1].data
    varst = datacube[:,pixs[i][1]-1,pixs[i][0]-1].var
    sigma =  np.sqrt(varst)
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    plt.step(lam1, ydata, label='D = %s' %Ds[i])
    ax.fill_between(lam1, 0.,varst,color='blue',alpha=.1)
    plt.legend()
    plt.show()'''


ews = []
mins = []
cent = []
variance = []
skewness = []

for i in range(len(pixs)):
    spe = datacube[:,pixs[i][1]-1,pixs[i][0]-1].data
    flux = spe[18:25]
    ewi = eqw(spe,18,25)
    minsi = min(flux)
    centi = centroid(spe,18,36)
    vari = moment(flux,2)
    skei = moment(flux,3)
    ews.append(ewi)
    mins.append(minsi)
    cent.append(centi)
    variance.append(vari)
    skewness.append(skei)

prob = [(100,10), (100,50) , (80, 10), (80,50)]
R = [10, 30, 50]

plt.figure()
plt.plot(Ds, ews)
plt.show()

plt.figure()
plt.plot(Ds, mins)
plt.show()

plt.figure()
plt.plot(Ds, cent)
plt.show()
plt.figure()
plt.plot(Ds, variance)
plt.show()
plt.figure()
plt.plot(Ds, skewness)
plt.show()

'''for i in range(len(prob)):
    for  j in range(len(R)):
        for k in range(len(Ds)):
            model = Disco(R[j], 10, 70, Rcore=0.1)
            spec = model.averagelosspec(Ds[k], alphas[k], lam1, 10,12,z,0.1, 10**14, 5, prob[i][0], prob[i][1], 0)
            np.save('pix%s' %k +'_R%s' %j + '_prob%s' %i, spec[1])'''

'''ews = []
mins = []
cent = []
variance = []
skewness = []

for i in range(len(R)):
    ewsj = []
    minsj = []
    centj = []
    variancej = []
    skewnessj = []
    for j in range(len(prob)):
        ewsk = []
        minsk = []
        centk = []
        variancek = []
        skewnessk = []
        for k in range(len(Ds)):
            print(i,j,k)
            spe = np.load('pix%s' %k +'_R%s' %i + '_prob%s' %j + '.npy')
            fluxii = spe[18:36]
            ewski = eqw(spe)
            minski = min(fluxii)
            centki = centroid(spe)
            varianceki = moment(fluxii,2)
            skewnesski = moment(fluxii,3)
            ewsk.append(ewski)
            minsk.append(minski)
            centk.append(centki)
            variancek.append(varianceki)
            skewnessk.append(skewnesski)
        ewsj.append(ewsk)
        minsj.append(minsk)
        centj.append(centk)
        variancej.append(variancek)
        skewnessj.append(skewnessk)
    ews.append(ewsj)
    mins.append(minsj)
    cent.append(centj)
    variance.append(variancej)
    skewness.append(skewnessj)'''




'''fig = plt.figure()
plt.suptitle('Prob(80-10)', fontsize=14)
ax1 = fig.add_subplot(2, 3, 1)
ax1.set_title('equivalent widths')
ax1.plot(Ds, ews[0][3][:], 'r', label='R = 10 kpc')
ax1.plot(Ds, ews[1][3][:], 'b', label='R = 30 kpc')
ax1.plot(Ds, ews[2][3][:], 'g', label='R = 50 kpc')
#ax1.set_xscale('log')
ax1.set_xlabel('impact parameter')
ax1.set_ylabel('eq.width')
ax1.legend()
ax2 = fig.add_subplot(2, 3, 2)
ax2.set_title('minimum absorption')
ax2.plot(Ds, mins[0][3][:], 'r')
ax2.plot(Ds, mins[1][3][:], 'b')
ax2.plot(Ds, mins[2][3][:], 'g',)
#ax2.set_xscale('log')
ax2.set_xlabel('impact parameter')
ax2.set_ylabel('min. flux')
ax2 = fig.add_subplot(2, 3, 4)
ax2.set_title('centroid')
ax2.plot(Ds, cent[0][3][:], 'r')
ax2.plot(Ds, cent[1][3][:], 'b')
ax2.plot(Ds, cent[2][3][:], 'g')
#ax2.set_xscale('log')
ax2.set_xlabel('impact parameter')
ax2.set_ylabel('centroid in velocity')
ax2 = fig.add_subplot(2, 3, 5)
ax2.set_title('variance')
ax2.plot(Ds, variance[0][3][:], 'r')
ax2.plot(Ds, variance[1][3][:], 'b')
ax2.plot(Ds, variance[2][3][:], 'g')
#ax2.set_xscale('log')
ax2.set_xlabel('impact parameter')
ax2.set_ylabel('variance')
ax2 = fig.add_subplot(2, 3, 6)
ax2.set_title('skewness')
ax2.plot(Ds, skewness[0][3][:], 'r')
ax2.plot(Ds, skewness[1][3][:], 'b')
ax2.plot(Ds, skewness[2][3][:], 'g')
#ax2.set_xscale('log')
ax2.set_xlabel('impact parameter')
ax2.set_ylabel('skewness')
fig.tight_layout()
plt.show()'''
