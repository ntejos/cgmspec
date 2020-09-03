from cgmspec.disco import Disco
from pylab import *
from astropy import units as u
from astropy.coordinates import SkyCoord
from matplotlib.ticker import NullFormatter
from astropy.convolution import convolve, Gaussian1DKernel
from mpdaf.obj import Cube
from mpdaf.obj import Spectrum, WaveCoord

props = dict(boxstyle='round', facecolor='papayawhip')

xabs_l = np.load('xabs_l.npy')
yabs_l = np.load('yabs_l.npy')

z = 0.73379
wave0 = 2796
cube=Cube('SW_nor_cube_12apr19.fits')

fluxobs = []
velobs = []

for i in range(len(xabs_l)):
    spe = cube[:,(10-i),0]
    wave = spe.wave.coord()

    fluxobsi = spe.data
    nflux = spe.data
    wave = spe.wave.coord()
    vel = (wave/(wave0*(1+z))-1)*299792.458
    velobs.append(vel)
    fluxobs.append(fluxobsi)

lam1 = np.arange(4000,5500,0.1)
museR = 4000
XshooterR = 8000

xgal_l=237.52059628552735
ygal_l=-78.188149277705151

posgal = SkyCoord(xgal_l*u.deg, ygal_l*u.deg, distance=4512.83311699*u.mpc, frame='icrs')



'''alphas = []
for i in  range(len(xabs_l)):
     locals()["pos"+str(i)] = SkyCoord(xabs_l[i]*u.deg, yabs_l[i]*u.deg, distance=4512.83311699*u.mpc, frame='icrs')
     locals()["pa"+str(i)] = posgal.position_angle(locals()["pos"+str(i)]).to(u.deg)
     alphas.append(70 - locals()["pa"+str(i)].value)

Ds = []
for i in range(len(xabs_l)):
    sep  = posgal.separation_3d(locals()["pos"+str(i)]).value  *1000
    Ds.append(sep)'''




galPA = 55
incli = 49
R = 100
h  = 10
slitPA = 42

D = [-3.6, 1.4, 3.6, 6.8, 9.9, 13.2, 16.5, 19.8, 23.0, 26.3, 29.7]
alphas = [10, 10,  10, 10, 10 ,10 ,10 ,10, 10, 10 ,10]
#alphas = [28, 28,  28, 28, 28 ,28 ,28 ,28, 28, 28 ,28]
model = Disco(R, h, incli, Rcore=0.1)

fig = plt.figure()  # 2 panels figsize=(7.5, 10.4)
ax = fig.add_subplot(1, 1, 1)
ax.set_xlabel(r'Velocity [km s$^{-1}$] from $z=$'+str(z), labelpad=16)
ax.set_ylabel(r'Normalized flux', labelpad=16)
ax.yaxis.set_tick_params(labelleft=False)
#ax.yaxis.set_tick_params(labelbottom=False)
#ax.xaxis.set_tick_params(labelleft=False)
ax.xaxis.set_tick_params(labelbottom=False)

def filtrogauss(fw, tau):
    gauss_kernel = Gaussian1DKernel(fw)
    gausflux = convolve(tau[1], gauss_kernel)
    return(tau[0],gausflux)

#ax.text(0.35,1.03,'MgII 2796,2803',fontsize=15)
'''for i in range(len(D)):

    Di = D[i]
    alphai = alphas[i]
    speci = model.averagelosspec(Di, alphai, lam1, 10,1,2,z)
    speci1 = [speci[0][0], speci[1][0]]
    speci2 = [speci[0][0], speci[1][1]]
    tspeci1m = filtrogauss(6.2,speci1)
    tspeci2m = filtrogauss(6.2,speci2)
    tspeci1x = filtrogauss(2.8,speci1)
    tspeci2x = filtrogauss(2.8,speci2)
    ax1 = fig.add_subplot(11, 2, i*2 +1)
    ax2 = fig.add_subplot(11, 2, i*2 +2, sharey=ax1 )

    ax1.plot(tspeci1m[0],tspeci1m[1], 'g')
    ax1.plot(tspeci2m[0],tspeci2m[1], 'g')
    ax1.plot(velobs[i], fluxobs[i], alpha=0.5, color='k',linewidth=1,marker=None, ls='-',drawstyle='steps-mid')
    #ax1.spines['top'].set_color('none')
    #ax2.spines['left'].set_color('none')
    #ax1.spines['right'].set_color('none')
    ax1.set_xticks([-800,-500,0,500, 1000])
    ax2.set_xticks([-1000,-500,0,500,1000,1500])
    if i < (len(D)-1):
        ax1.set_xticks([])
        ax2.set_xticks([])

    ax2.yaxis.set_tick_params(labelleft=False)
    ax1.set_yticks([0,1])
    ax2.plot(tspeci1x[0],tspeci1x[1], 'g')
    ax2.plot(tspeci2x[0],tspeci2x[1], 'g')
    plt.subplots_adjust(wspace=0,hspace=0)
    ax1.set_ylim(0,1.4)
    ax2.set_ylim(0,1.4)
    ax1.set_xlim(-1000,1500)
    ax2.set_xlim(-1000,1500)
plt.show()'''

for i in range(len(D)):

    Di = D[i]
    alphai = alphas[i]
    speci = model.averagelosspec(Di, alphai, lam1, 10,1,2,z)
    speci1 = [speci[0][0], speci[1][0]]
    speci2 = [speci[0][0], speci[1][1]]
    tspeci1m = filtrogauss(6.2,speci1)
    tspeci2m = filtrogauss(6.2,speci2)
    ax1 = fig.add_subplot(11, 1, i+1)

    ax1.plot(tspeci1m[0],tspeci1m[1], 'g')
    ax1.plot(tspeci2m[0],tspeci2m[1], 'g')
    ax1.plot(velobs[i], fluxobs[i], alpha=0.5, color='k',linewidth=1,marker=None, ls='-',drawstyle='steps-mid')
    #ax1.spines['top'].set_color('none')
    #ax2.spines['left'].set_color('none')
    #ax1.spines['right'].set_color('none')
    ax1.set_xticks([-500,0,500, 1000])
    if i < (len(D)-1):
        ax1.set_xticks([])

    ax1.set_yticks([0,1])

    plt.subplots_adjust(wspace=0,hspace=0)
    ax1.set_ylim(0,1.4)
    ax1.set_xlim(-800,1500)
    ax1.text(0.03, 0.90, 'D = %s' % D[i] + ' kpc', transform=ax1.transAxes, fontsize=8,
            verticalalignment='top', bbox=props)
plt.show()
