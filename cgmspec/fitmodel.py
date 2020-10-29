
import matplotlib.pyplot as plt
from cgmspec import utils as csu
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


#select spaxels
s_n_mask = np.load('out_r4_3.npy')
all_pixels = []
for i in range(len(s_n_mask)):
    posi = [s_n_mask[i][0], s_n_mask[i][1]]
    all_pixels.append(posi)

all_pixels = np.asarray(all_pixels)

all_pixels1 = all_pixels[0:25]
all_pixels2 = all_pixels[25:55]
all_pixels3 = all_pixels[55:-1]

pixel1 = all_pixels1[np.random.choice(all_pixels1.shape[0],6,replace=False)]
pixel2 = all_pixels2[np.random.choice(all_pixels2.shape[0],6,replace=False)]
pixel3 = all_pixels3[np.random.choice(all_pixels3.shape[0],6,replace=False)]

pixel = np.concatenate((pixel1, pixel2,pixel3))

print(pixel)
#All parameters in my model
params = Parameters()
params.add_many(
        ('incli', 80, True),
        ('col_dens', 14, True),
        ('height', 12.0, True),
        ('dop_param', 8.0, True),
        ('csize', 0.1, False),
        ('r_0', 1, True),
        ('vel_max', 200, True),
        ('h_v', 4.5, True),
        )



spectralpixel =  [14,36] #spectralpixels to "count" an absorption

lam1 = np.arange(4825.12, 4886.37+1.25, 1.25) #wavelenghts to plot the rebinned data
lam2 = np.arange(4825.12, 4886.37+1.5, 0.1) # wavelenghts to calculatethe model, it defines how fine is my model going to be



#cosmological parameters to calculate the distance
z = 0.73379
scale = 7.28
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
dist = cosmo.luminosity_distance(z)

#galaxy imputs
galPA = 55
galcen = SkyCoord(ra=237.52059628552735*u.degree, dec=-78.188149277705151*u.degree, frame='icrs')

datacube  = Cube('nor_cut_cube.fits')
head = datacube.primary_header
wcs1 = WCS(head)
wave1 = WaveCoord(cdelt=0.1, crval=4825.12, cunit=u.angstrom)







def residual(p, datacube=datacube, pos=None, specwidth=spectralpixel, galcen=galcen, galPA=galPA):
        par = p.valuesdict()
        incli = par['incli']
        col_denst = par['col_dens']
        h = par['height']
        bes = par['dop_param']
        v_max = par['vel_max']
        h_vt = par['h_v']

        csize = par['csize']
        r_0t= par['r_0']

        print(par)
        h_v = 10**h_vt
        col_dens = 10**col_denst
        r_0 = 10**r_0t
        coord_sky = datacube.wcs.pix2sky([int(pos[0]),int(pos[1])], unit=u.deg)
        dec = coord_sky[0][0]
        ra = coord_sky[0][1]
        c1 = SkyCoord(ra*u.degree,dec*u.degree,frame='icrs')

        flux = datacube[:,int(pos[0]),int(pos[1])]
        abso = flux.data[specwidth[0]:specwidth[1]+1]
        sigma = np.sqrt(flux.var[specwidth[0]:specwidth[1]+1])

        scale = 7.28 # z=0.73379  # plat scale (kpc/") for Planck
        c1 = SkyCoord(ra*u.degree,dec*u.degree,frame='icrs')
        #D = scale * galcen.separation(c1).arcsec
        #pa = galcen.position_angle(c1).to(u.deg)
        #alpha = galPA - pa.value

        alpha = csu.get_alpha(c1, galcen,galPA)
        D = csu.get_impactparam(c1, galcen, scale)
        lam2 = np.arange(4825.12, 4886.37+1.5, 0.1)
        wave1 = WaveCoord(cdelt=0.1, crval=4825.12, cunit= u.angstrom)

        model = Disco(h, incli, Rcore=0.1)
        spec = model.averagelos(D, alpha, lam2, 100,12,z,csize, col_dens, bes, r_0, v_max,h_v, 0)
        spe = Spectrum(wave=wave1, data=spec)
        rspe = spe.resample(1.25)
        fluxmodel = rspe.data
        absomodel = fluxmodel[specwidth[0]:specwidth[1]+1]

        #ymodelt = np.concatenate((ymodel[0][14:36],ymodel[1][16:36],ymodel[2][16:36]))
        rest = ((abso-absomodel)/sigma)**2
        dif = abso-absomodel
        return(fluxmodel, rest, datacube[:,int(pos[0]),int(pos[1])].data, sigma, dif)




def chisq(p):
    tosum = []
    absos = []
    data = []
    sigma = []
    for i in range(len(pixel)):
        res = residual(p, pos=pixel[i])
        tosum.append(res[1])

    tosum = np.asarray(tosum)

    tosum = tosum.flatten()
    chi = np.sum(tosum)
    return(chi)


#Do the minimization with the parameters that i want

'''params['incli'].set(min=75, max=90, brute_step=5)
params['col_dens'].set(min=12, max=18, brute_step=2)
params['height'].set(min=11, max=13, brute_step=0.5)
params['r_0'].set(min=0.5, max=3, brute_step=0.5)
params['dop_param'].set(min=5, max=9, brute_step=1)
params['vel_max'].set(min=190, max=205, brute_step=5)
params['h_v'].set(min=3.5, max=5.5, brute_step=0.5)
#params['csize'].set(min=-2, max=1, brute_step=1)


fitter = Minimizer(chisq, params)
result = fitter.minimize(method='brute')

print('best fit parameters:', result.brute_x0)
print('result.brute_fval:', result.brute_fval)'''


#plot many absorptions with model

def plot_results(p, pixel):
    par = p.valuesdict()
    incli = par['incli']
    col_denst = par['col_dens']
    h = par['height']
    bes = par['dop_param']
    v_max = par['vel_max']
    h_vt = par['h_v']

    csize = par['csize']
    r_0t= par['r_0']

    col_dens = 10**col_denst
    h_v = 10**h_vt
    r_0 = 10**r_0t
    alphas=[]
    Ds = []
    for i in range(len(pixel)):
        coord_sky = datacube.wcs.pix2sky([pixel[i][0],pixel[i][1]], unit=u.deg)
        dec = coord_sky[0][0]
        ra = coord_sky[0][1]
        c1 = SkyCoord(ra*u.degree,dec*u.degree,frame='icrs')

        alpha = csu.get_alpha(c1, galcen,galPA)
        D = csu.get_impactparam(c1, galcen, scale)
        alphas.append(alpha)
        Ds.append(D)

    Ds  = np.asarray(Ds)
    alphas = np.asarray(alphas)
    pixel = np.asarray(pixel)
    alphas =  alphas[Ds.argsort()]
    pixel = pixel[Ds.argsort()]
    Ds = np.sort(Ds)

    print(Ds)
    print(alphas)


    #modelflux = np.asarray(modelflux)

    fluxes = []
    chis = []
    modelflux = []
    sigma = []
    diff = []
    for i in range(len(pixel)):


         chii = residual(p, datacube=datacube, pos=pixel[i], specwidth=spectralpixel, galcen=galcen, galPA=galPA)
         chis.append(chii[1])
         modelflux.append(chii[0])
         fluxes.append(chii[2])
         sigma.append(chii[3])
         diff.append(chii[4])
         '''lam2 = np.arange(4825.12, 4886.37+1.5, 0.1)
         wave1 = WaveCoord(cdelt=0.1, crval=4825.12, cunit= u.angstrom)

         model = Disco(h, incli, Rcore=0.1)
         spec = model.averagelos(Ds[i], alphas[i], lam2, 100,12,z,csize, col_dens, bes, r_0, v_max,h_v, 0)
         spe = Spectrum(wave=wave1, data=spec)
         rspe = spe.resample(1.25)
         fluxmodel = rspe.data
         modelflux.append(fluxmodel)'''


    '''alphas = alphas[Ds.argsort()]
    chis = chis[Ds.argsort()]
    model = model[Ds.argsort()]'''


    fig, axs = plt.subplots(int((len(pixel)/2)), 2, sharex=True, sharey=True)
    bbox = dict(boxstyle="round", fc="0.8")
    fig.text(0.5, 0.04, 'Wavelenght angstrom', ha='center', va='center')
    fig.text(0.06, 0.5, 'Normalized Flux', ha='center', va='center', rotation='vertical')

    for i in range(int(len(pixel))):

        if i < int(len(pixel))/2:
            ax1 = axs[i, 0]
        else:
            ax1 = axs[int(i-len(pixel)/2), 1]

        print(lam1)
        print(fluxes[i])
        ax1.step(lam1,fluxes[i])
        ax1.plot(lam1,modelflux[i],'r')
        #ax1.fill_between(lam1[spectralpixel[0]:spectralpixel[1]+1], diff[i], step='pre', alpha=0.5)
        ax1.fill_between(lam1[spectralpixel[0]:spectralpixel[1]+1], sigma[i], step='pre', alpha = 0.5)
        #ax1.annotate('D=%s' %str(round(Ds[i],2)) + 'Kpc', (lam1[0], 1.5), bbox=bbox)
        ax1.annotate('D=%s' %str(round(Ds[i],2)) + 'Kpc, chi=%s' %str(round(np.sum(chis[i]),2)), (lam1[0], 1.5), bbox=bbox)
        ax1.set_ylim(0,1.5)









#plot results this function is taken from the lmfit page
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

#--------------------------------------------------------------------------------------------------------------------------


'''alphas=[]
Ds = []
for i in range(len(pixel)):
    coord_sky = datacube.wcs.pix2sky([pixel[i][1]-1,pixel[i][0]-1], unit=u.deg)
    dec = coord_sky[0][0]
    ra = coord_sky[0][1]
    c1 = SkyCoord(ra*u.degree,dec*u.degree,frame='icrs')

    alpha = csu.get_alpha(c1, galcen,galPA)
    D = csu.get_impactparam(c1, galcen, scale)
    alphas.append(alpha)
    Ds.append(D)



# cut data to match the size of the model spectra
ydata = []
varst = []
sigma = []


for i in range(len(pixel)):
    flux = datacube[:,pixel[i][1]-1,pixel[i][0]-1]
    ydatai = flux.data[spectralpixel[0]:spectralpixel[1]+1]
    sigmai = np.sqrt(flux.var[spectralpixel[0]:spectralpixel[1]+1])
    ydata.append(ydatai)
    sigma.append(sigmai)

#for i in range(len(sigma)):
    #sigma[i] = sigma[i][spectralpixel[0]:spectralpixel[1]+1]

sigma=np.asarray(sigma)
ydata = np.asarray(ydata)
sigma = sigma.flatten()
ydata = ydata.flatten()'''



#calculate D and alpha for each position in pixel
'''alphas=[]
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
    Ds.append(D)'''



'''Ds = np.asarray(Ds)
alphas = np.asarray(alphas)


chis = chisq(params)

alphas = alphas[Ds.argsort()]
chis = np.asarray(chis[0])
chis = chis[Ds.argsort()]
ydata = ydata[Ds.argsort()]
Ds = np.sort(Ds)
incli = 70

#h = 5
csize = 0.1
N = 10**14
bes = 5
r_0 = 100
h = 11
v_max=200
h_v  =10**4
ymodel=[]
for i in range(len(pixel)):
    lam2 = np.arange(4825.12, 4886.37+1.5, 0.1)
    model = Disco(h, incli, Rcore=0.1)
    spec = model.averagelos(Ds[i], alphas[i], lam2, 100,12,z,csize, N, bes, r_0, v_max,h_v, 0)
    spe = Spectrum(wave=wave1, data=spec)
    rspe = spe.resample(1.25)
    ym = rspe.data
    ymodel.append(ym)

fig, axs = plt.subplots(int((len(pixel)/2)), 2, sharex=True, sharey=True)
bbox = dict(boxstyle="round", fc="0.8")
for i in range(int(len(pixel))):

    print(axs)
    if i < int(len(pixel))/2:
        ax1 = axs[i, 0]
    else:
        ax1 = axs[int(i-len(pixel)/2), 1]

    ax1.step(lam1,ydata[i])
    ax1.plot(lam1,ymodel[i],'r')
    ax1.annotate('D=%s' %str(round(Ds[i],2)) + 'Kpc, chi=%s' %str(round(chis[i],2)), (lam1[0], 1.5), bbox=bbox)
    ax1.set_ylim(0,1.5)'''



'''def chisq(p):
    par = p.valuesdict()
    incli = par['incli']
    col_denst = par['col_dens']
    h = par['height']
    bes = par['dop_param']
    v_max = par['vel_max']
    h_vt = par['h_v']

    csize = par['csize']
    r_0= par['r_0']

    print(incli)
    h_v = 10**h_vt
    col_dens = 10**col_denst
    #for i in range(len(ydata)):
#        ydatatt.append(ydata[i][spectralpixel[0]:spectralpixel[1]])
#    ydatat=np.array(ydatatt)
#    y=ydatat.flatten()
    #y = np.concatenate((ydata[0][14:36],ydata[1][16:36], ydata[2][16:36]))


    ymodelt = []
    lamt = []
    ymodel =[]
    for i in range(len(pixel)):
        lam2 = np.arange(4825.12, 4886.37+1.5, 0.1)
        wave1 = WaveCoord(cdelt=0.1, crval=4825.12, cunit= u.angstrom)

        model = Disco(h, incli, Rcore=0.1)
        spec = model.averagelos(Ds[i], alphas[i], lam2, 100,12,z,csize, col_dens, bes, r_0, v_max,h_v, 0)
        spe = Spectrum(wave=wave1, data=spec)
        rspe = spe.resample(1.25)
        ym = rspe.data

        ymodelt.append(ym)

    for i in range(len(ymodelt)):
        ymodel.append(ymodelt[i][spectralpixel[0]:spectralpixel[1]+1])
    ymodel=np.asarray(ymodel)
    ymodel = ymodel.flatten()
    #modelflux.append(ymodel)
    #ymodelt = np.concatenate((ymodel[0][14:36],ymodel[1][16:36],ymodel[2][16:36]))

    tosum = []

    rest = ((ydata-ymodel)/sigma)**2
    #for i in range(len(ymodel)):
    #    rest = ((ydata[i]-ymodel[i])/sigma[i])**2

    #    tosum.append(rest)
    #chipix = []
    #for i in range(len(pixel)):
    #    j = i * 22
    #    tosumpix = tosum[j:j+22]
    #    chii = np.sum(tosumpix)
    #    chipix.append(chii)
    #restos.append(tosum)
    tosum = np.asarray(rest)
    return(np.sum(tosum))'''



#pixel = [[7,19], [5,17], [9,18], [9,16], [7,14], [12,16], [12,14], [11,13], [15,12], [13,9], [13,8], [16,10], [18,9], [17,5], [19,6],[20,6],[17,3],[22,4]]
#pixel = [[11,17],[10,17],[9,17],[8,17],[8,16],[9,16],[10,16],[11,16],[10,15],[9,15],[8,15],[9,14],[10,14],[11,14],[13,13],[14,10],[14,9],[16,9],[12,12],[12,9],[15,9],[17,9]]
#pixel = [[8,17]]


#pixel =  [[11,17],[9,17],[8,16],[10,16],[8,15],[10,14],[13,13],[14,9],[12,12],[15,9]]#position in datacube to do the fit
