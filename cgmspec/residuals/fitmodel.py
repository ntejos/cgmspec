import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from mpdaf.obj import Cube, WCS, WaveCoord, Spectrum
import astropy.units as u
from cgmspec.disco import Disco
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.cosmology import FlatLambdaCDM
from lmfit import Model
from scipy.optimize import least_squares
from lmfit import Minimizer, Parameters, fit_report
from matplotlib.colors import LogNorm

params = Parameters()
params.add_many(
        ('incli', 85.0, True),
        ('col_dens', 15.0, True),
        ('height', 5.0, True),
        ('dop_param', 5.0, False),
        ('csize', 0.1, False),
        ('r_0', 1000, False),
        ('vel_max', 200, False),
        ('h_v', 10, False),
        )


truepix = [[11,17],[10,17],[9,17],[8,17],[8,16],[9,16],[10,16],[11,16],[10,15],[9,15],[8,15],[9,14],[10,14],[11,14],[13,13],[14,10],[14,9],[16,9],[12,12],[12,9],[15,9],[17,9]]
mabypix = [[8,18], [7,16], [12,16], [11,15], [7,15],[6,15], [6,14],[7,14],[8,14],[12,14],[13,14],[15,14],[14,13],[12,13],[11,13],[7,13],[14,12],[18,9],[17,7], [19,12],[16,11],[14,11],[18,9],[19,6]]
#pixel = [[8,17]]
pixel = [[8,17],[8,15],[12,12],[16,10],[18,9],[20,6]]

spectralpixel =  [14,36]

lam1 = np.arange(4825.12, 4886.37+1.25, 1.25)
lam2 = np.arange(4825.12, 4886.37+1.5, 0.1)



z = 0.73379
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
dist = cosmo.luminosity_distance(z)

#galaxy imputs

galPA = 55


galcen = SkyCoord(ra=237.52059628552735*u.degree, dec=-78.188149277705151*u.degree, frame='icrs')



datacube  = Cube('ccdatatrue.fits')
head = datacube.primary_header
wcs1 = WCS(head)
wave1 = WaveCoord(cdelt=0.5, crval=4825.12, cunit=u.angstrom)

ydata = []
varst = []
sigma = []
sigmat = []
for i in range(len(pixel)):
    ydatai = datacube[:,pixel[i][1]-1,pixel[i][0]-1].data
    varsti = datacube[:,pixel[i][1]-1,pixel[i][0]-1].var
    sigmai  = np.sqrt(varsti)
    ydata.append(ydatai)
    varst.append(varsti)
    sigma.append(sigmai)

for i in range(len(sigma)):
    sigma[i] = sigma[i][spectralpixel[0]:spectralpixel[1]]

sigmatt=np.asarray(sigma)
sigmat = sigmatt.flatten()

'''ydata = [datacube[:,pixel[0][1]-1,pixel[0][0]-1].data, datacube[:,pixel[1][1]-1,pixel[1][0]-1].data, datacube[:,pixel[2][1]-1,pixel[2][0]-1].data]
    varst = [datacube[:,pixel[0][1]-1,pixel[0][0]-1].var, datacube[:,pixel[1][1]-1,pixel[1][1]-1].var, datacube[:,pixel[2][1]-1,pixel[2][0]-1].var]
    sigma =  [np.sqrt(varst[0]), np.sqrt(varst[1]), np.sqrt(varst[2])]
    sigmat = np.concatenate((sigma[0][14:36],sigma[1][16:36], sigma[2][16:36]))'''

#calcular d y alfa para un pixel
'''
coord_sky = datacube.wcs.pix2sky([pixel[1]-1,pixel[0]-1], unit=u.deg)
dec = coord_sky[0][0]
ra = coord_sky[0][1]
scale = 7.28 # z=0.73379  # plat scale (kpc/") for Planck
c1 = SkyCoord(ra*u.degree,dec*u.degree,frame='icrs')
D = scale * galcen.separation(c1).arcsec
pa = galcen.position_angle(c1).to(u.deg)
alpha = galPA - pa.value
print(D,alpha)
wave1 = WaveCoord(cdelt=0.1, crval=4751.37, cunit= u.angstrom, shape=247.5)'''

#calcular d  y  alfa para lista
alphas=[]
Ds = []
for i in range(len(pixel)):
    coord_sky = datacube.wcs.pix2sky([pixel[i][1]-1,pixel[i][0]-1], unit=u.deg)
    dec = coord_sky[0][0]
    ra = coord_sky[0][1]
    scale = 7.28 # z=0.73379  # plat scale (kpc/") for Planck
    c1 = SkyCoord(ra*u.degree,dec*u.degree,frame='icrs')
    D = scale * galcen.separation(c1).arcsec
    pa = galcen.position_angle(c1).to(u.deg)
    alpha = galPA - pa.value
    wave1 = WaveCoord(cdelt=0.1, crval=4751.37, cunit= u.angstrom, shape=247.5)
    alphas.append(alpha)
    Ds.append(D)

'''
for i in range(len(ydata)):
    plt.figure()
    redlam = lam1[14:36]
    redy = ydata[i][14:36]
    bluelam1 = lam1[0:15]
    bluey1 =  ydata[i][0:15]
    bluelam2 = lam1[35:]
    bluey2 = ydata[i][35:]
    plt.step(redlam, redy, color='r')
    plt.step(bluelam1, bluey1, color='b')
    plt.step(bluelam2, bluey2, color='b')
    plt.ylim(0,1.5)
    plt.show()'''


#plot data and var
'''incli = 45
R =  50
#h = 5
csize = 0.1
Nes = 10**14.5
bes = 5
pmax = 100
pmin = 50
h = 7

model = Disco(R, h, incli, Rcore=0.1)
spec = model.averagelosspec(D, alpha, lam2, 10,12,z,csize, Nes, bes, pmax, pmin,0)
spe = Spectrum(wave=wave1, data=spec[1])
rspe = spe.resample(1.25)

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
plt.step(lam1, ydata)
plt.plot()
ax.fill_between(lam1, 0.,varst,color='blue',alpha=.1)
plt.plot(lam1,rspe.data, 'r')
plt.ylabel('Norm. Flux')
plt.xlabel('wavelenght')'''

def chisq(p):
    par = p.valuesdict()
    incli = par['incli']
    col_dens = par['col_dens']
    h = par['height']
    bes = par['dop_param']
    v_max = par['vel_max']
    h_v = par['h_v']
    csize = par['csize']
    r_0= par['r_0']
    print(len(ydata), len(ydata[0]))
    ydatatt = []
    for i in range(len(ydata)):
        ydatatt.append(ydata[i][spectralpixel[0]:spectralpixel[1]])
    ydatat=np.array(ydatatt)
    y=ydatat.flatten()
    #y = np.concatenate((ydata[0][14:36],ydata[1][16:36], ydata[2][16:36]))
    N = 10**(col_dens)
    ymodel = []
    lamt = []
    for i in range(len(pixel)):
        lam2 = np.arange(4825.12, 4886.37+1.5, 0.1)
        model = Disco(h, incli, Rcore=0.1)
        spec = model.averagelosspec(Ds[i], alphas[i], lam2, 100,12,z,csize, N, bes, r_0, v_max,h_v, 0)
        spe = Spectrum(wave=wave1, data=spec[1])
        rspe = spe.resample(1.25)
        ym = rspe.data
        ymodel.append(ym)
    for i in range(len(ymodel)):
        ymodel[i] = ymodel[i][spectralpixel[0]:spectralpixel[1]]
    ymodeltt=np.asarray(ymodel)
    ymodelt = ymodeltt.flatten()
    #ymodelt = np.concatenate((ymodel[0][14:36],ymodel[1][16:36],ymodel[2][16:36]))
    print(ymodelt)
    tosum = np.empty(len(ymodelt))
    for i in range(len(ymodelt)):
        print('y',len(y))
        print('ymodelt', len(ymodelt))
        print('sigmat', len(sigmat))
        rest = ((y[i]-ymodelt[i])/sigmat[i])**2
        print('rest', rest)
        tosum[i]=rest
    return(np.sum(tosum))

#hacer la minimizacion
params['incli'].set(min=40, max=80, brute_step=5)
params['col_dens'].set(min=13, max=15, brute_step=0.5)
params['height'].set(min=2, max=10, brute_step=2)
#params['dop_param'].set(min=2, max=10, brute_step=1)
#params['vel_max'].set(min=195, max=205, brute_step=2)
#params['h_v'].set(min=10, max=100, brute_step=20)



fitter = Minimizer(chisq, params)
result = fitter.minimize(method='brute')


#plot many absorptions

'''incli = 45
R =  60
#h = 5
csize = 0.1
N = 10**14.5
bes = 5
r_0 = 1000
h = 8
v_max=200
h_v  =10
ymodel=[]
for i in range(len(pixel)):
    lam2 = np.arange(4825.12, 4886.37+1.5, 0.1)
    model = Disco(R, h, incli, Rcore=0.1)
    spec = model.averagelosspec(Ds[i], alphas[i], lam2, 10,12,z,csize, N, bes, r_0, v_max,h_v, 0)
    spe = Spectrum(wave=wave1, data=spec[1])
    rspe = spe.resample(1.25)
    ym = rspe.data
    ymodel.append(ym)

fig = plt.figure()
for i in range(len(pixel)):
    ax1 = fig.add_subplot((len(pixel)), 1, i+1)
    ax1.step(lam1,ydata[i])
    ax1.plot(lam1,ymodel[i],'r')
    ax1.set_ylim(0,1.5)'''


'''ax2 = fig.add_subplot(3, 1, 2)
ax2.step(lam1,ydata[1])
ax2.plot(lam1,ymodel[1],'r')
ax2.set_ylim(0,1.5)
ax3 = fig.add_subplot(3, 1, 3)
ax3.step(lam1,ydata[2])
ax3.plot(lam1,ymodel[2],'r')
ax3.set_ylim(0,1.5)
plt.show()'''





'''def model(x, incli):
    print('icli', incli)
    lam2 = np.arange(4825.12, 4886.37+1.5, 0.1)
    csize = 0.1
    Nes = 10**15
    bes = 5
    pmax = 100
    pmin = 80

    model = Disco(R, h, incli, Rcore=0.1)
    spec = model.averagelosspec(D, alpha, lam2, 10,12,z,csize, Nes, bes, pmax, pmin,0)
    spe = Spectrum(wave=wave1, data=spec[1])
    rspe = spe.resample(1.25)
    return(rspe.data)'''


'''
plt.figure()
plt.plot(lam1, ydata, 'b')
plt.plot(lam1, model(lam1,48.34880853), 'r')
plt.show()'''


#plot results
def plot_results_brute(result, best_vals=True, varlabels=None,
                       output=None):
    """Visualize the result of the brute force grid search.

    The output file will display the chi-square value per parameter and contour
    plots for all combination of two parameters.

    Inspired by the `corner` package (https://github.com/dfm/corner.py).

    Parameters
    ----------
    result : :class:`~lmfit.minimizer.MinimizerResult`
        Contains the results from the :meth:`brute` method.

    best_vals : bool, optional
        Whether to show the best values from the grid search (default is True).

    varlabels : list, optional
        If None (default), use `result.var_names` as axis labels, otherwise
        use the names specified in `varlabels`.

    output : str, optional
        Name of the output PDF file (default is 'None')
    """
    npars = len(result.var_names)
    _fig, axes = plt.subplots(npars, npars)

    if not varlabels:
        varlabels = result.var_names
    if best_vals and isinstance(best_vals, bool):
        best_vals = result.params

    for i, par1 in enumerate(result.var_names):
        for j, par2 in enumerate(result.var_names):

            # parameter vs chi2 in case of only one parameter
            if npars == 1:
                axes.plot(result.brute_grid, result.brute_Jout, 'o', ms=3)
                axes.set_ylabel(r'$\chi^{2}$')
                axes.set_xlabel(varlabels[i])
                if best_vals:
                    axes.axvline(best_vals[par1].value, ls='dashed', color='r')

            # parameter vs chi2 profile on top
            elif i == j and j < npars-1:
                if i == 0:
                    axes[0, 0].axis('off')
                ax = axes[i, j+1]
                red_axis = tuple([a for a in range(npars) if a != i])
                ax.plot(np.unique(result.brute_grid[i]),
                        np.minimum.reduce(result.brute_Jout, axis=red_axis),
                        'o', ms=3)
                ax.set_ylabel(r'$\chi^{2}$')
                ax.yaxis.set_label_position("right")
                ax.yaxis.set_ticks_position('right')
                ax.set_xticks([])
                if best_vals:
                    ax.axvline(best_vals[par1].value, ls='dashed', color='r')

            # parameter vs chi2 profile on the left
            elif j == 0 and i > 0:
                ax = axes[i, j]
                red_axis = tuple([a for a in range(npars) if a != i])
                ax.plot(np.minimum.reduce(result.brute_Jout, axis=red_axis),
                        np.unique(result.brute_grid[i]), 'o', ms=3)
                ax.invert_xaxis()
                ax.set_ylabel(varlabels[i])
                if i != npars-1:
                    ax.set_xticks([])
                elif i == npars-1:
                    ax.set_xlabel(r'$\chi^{2}$')
                if best_vals:
                    ax.axhline(best_vals[par1].value, ls='dashed', color='r')

            # contour plots for all combinations of two parameters
            elif j > i:
                ax = axes[j, i+1]
                red_axis = tuple([a for a in range(npars) if a not in (i, j)])
                X, Y = np.meshgrid(np.unique(result.brute_grid[i]),
                                   np.unique(result.brute_grid[j]))
                lvls1 = np.linspace(result.brute_Jout.min(),
                                    np.median(result.brute_Jout)/2.0, 7, dtype='int')
                lvls2 = np.linspace(np.median(result.brute_Jout)/2.0,
                                    np.median(result.brute_Jout), 3, dtype='int')
                lvls = np.unique(np.concatenate((lvls1, lvls2)))
                ax.contourf(X.T, Y.T, np.minimum.reduce(result.brute_Jout, axis=red_axis),
                            lvls, norm=LogNorm())
                ax.set_yticks([])
                if best_vals:
                    ax.axvline(best_vals[par1].value, ls='dashed', color='r')
                    ax.axhline(best_vals[par2].value, ls='dashed', color='r')
                    ax.plot(best_vals[par1].value, best_vals[par2].value, 'rs', ms=3)
                if j != npars-1:
                    ax.set_xticks([])
                elif j == npars-1:
                    ax.set_xlabel(varlabels[i])
                if j - i >= 2:
                    axes[i, j].axis('off')

    if output is not None:
        plt.savefig(output)


#plot density function
'''incli = 10
R =  50
h = 5
csize = 0.1
N = 10**13
bes = 5
r_0 = 1000
h = 5
v_max=200
h_v  =10
ymodel=[]
for i in range(len(pixel)):
    lam2 = np.arange(4825.12, 4886.37+1.5, 0.1)
    model = Disco(R, h, incli, Rcore=0.1)

rs = np.linspace(0,30,1000)
ps = []
for i in range(len(rs)):
    p = model.prob_hit(rs[i], r_0)
    ps.append(p)
plt.figure()
plt.plot(rs,ps)'''
