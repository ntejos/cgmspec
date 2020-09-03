import pylab as plt
from mpdaf.obj import Cube, WCS, WaveCoord, Spectrum
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u

lam1 = np.arange(4825.12, 4886.37+1.25, 1.25)
lam2 = np.arange(4825.12, 4886.37+0.1, 0.1)

datacube  = Cube('ccdatanorm.fits')
modelcube11 = Cube('param_i_p1_0.fits')
modelcube12 = Cube('param_i_p1_1.fits')
modelcube13 = Cube('param_i_p1_2.fits')
modelcube14 = Cube('param_i_p1_3.fits')
modelcube15 = Cube('param_i_p1_4.fits')
modelcube16 = Cube('param_i_p1_5.fits')
modelcube17 = Cube('param_i_p1_6.fits')
modelcube18 = Cube('param_i_p1_7.fits')
modelcube19 = Cube('param_i_p1_8.fits')
modelcube110 = Cube('param_i_p1_9.fits')

modelcube21 = Cube('param_i_p2_0.fits')
modelcube22 = Cube('param_i_p2_1.fits')
modelcube23 = Cube('param_i_p2_2.fits')
modelcube24 = Cube('param_i_p2_3.fits')
modelcube25 = Cube('param_i_p2_4.fits')
modelcube26 = Cube('param_i_p2_5.fits')
modelcube27 = Cube('param_i_p2_6.fits')
modelcube28 = Cube('param_i_p2_7.fits')
modelcube29 = Cube('param_i_p2_8.fits')
modelcube210 = Cube('param_i_p2_9.fits')

modelcube31 = Cube('param_i_p3_0.fits')
modelcube32 = Cube('param_i_p3_1.fits')
modelcube33 = Cube('param_i_p3_2.fits')
modelcube34 = Cube('param_i_p3_3.fits')
modelcube35 = Cube('param_i_p3_4.fits')
modelcube36 = Cube('param_i_p3_5.fits')
modelcube37 = Cube('param_i_p3_6.fits')
modelcube38 = Cube('param_i_p3_7.fits')
modelcube39 = Cube('param_i_p3_8.fits')
modelcube310 = Cube('param_i_p3_9.fits')


'''modelcube21 = Cube('resaex_p2_0.fits')
modelcube22 = Cube('resaex_p2_1.fits')
modelcube23 = Cube('resaex_p2_2.fits')

modelcube31 = Cube('resaex_p3_0.fits')
modelcube32 = Cube('resaex_p3_1.fits')
modelcube33 = Cube('resaex_p3_2.fits')'''

models1 = [[modelcube11, modelcube12,  modelcube13, modelcube14, modelcube15,  modelcube16, modelcube17, modelcube18,  modelcube19, modelcube110],[modelcube21, modelcube22,  modelcube23, modelcube24, modelcube25,  modelcube26, modelcube27, modelcube28,  modelcube29, modelcube210],[modelcube31, modelcube32,  modelcube33, modelcube34, modelcube35,  modelcube36, modelcube37, modelcube38,  modelcube39, modelcube310]]

#models1 = [[modelcube11, modelcube12,  modelcube13, modelcube14, modelcube15,  modelcube16, modelcube17, modelcube18,  modelcube19, modelcube110], [modelcube21, modelcube22,modelcube23], [modelcube31, modelcube32, modelcube33]]
#models1 = [modelcube21, modelcube22,modelcube23]
#models1 = [modelcube31, modelcube32, modelcube33]
galcen = SkyCoord(ra=237.52059628552735*u.degree, dec=-78.188149277705151*u.degree, frame='icrs')

hs = ['13', '15', '17']

pixels = [(9,15), (14,11), (19,6)]

props = dict(boxstyle='round', facecolor='papayawhip')
fig = plt.figure()  # 2 panels figsize=(7.5, 10.4)
ax = fig.add_subplot(1, 1, 1)
ax.set_title('b=5, h=5, vel=196km/s,  i=85', y=1.08)
ax.set_ylabel(r'Normalized flux', labelpad=16)
ax.yaxis.set_tick_params(labelleft=False)
#ax.yaxis.set_tick_params(labelbottom=False)
#ax.xaxis.set_tick_params(labelleft=False)
ax.xaxis.set_tick_params(labelbottom=False)

for i in range(3):
    for  j in range(10):


      # R = D* np.sqrt(1+(np.sin(al_rad)**2)*np.tan(self.incl_rad)**2)


       dspe = datacube[:,pixels[i][1], pixels[i][0]]
       mspe = models1[i][j][:,15,8]
       ax1 = fig.add_subplot(10, 3, j*3+1+i)

       '''if j == 0:
           ax1 = fig.add_subplot(10, 3, i*5+1+j)
       if j == 1:
           ax1 = fig.add_subplot(10, 3, i*5+1+j)
       if j == 2:
           ax1 = fig.add_subplot(10, 3, i*5+1+j)'''
       ax1.step(lam1,dspe.data, 'g')
       ax1.plot(lam1,mspe.data, 'b')

       if j == 0:
           coord_sky = datacube.wcs.pix2sky([pixels[i][1], pixels[i][0]], unit=u.deg)
           dec = coord_sky[0][0]
           ra = coord_sky[0][1]
           scale = 7.28 # z=0.73379  # plat scale (kpc/") for Planck
           c1 = SkyCoord(ra*u.degree,dec*u.degree,frame='icrs')
           D = round(scale * galcen.separation(c1).arcsec,3)
           ax1.set_title('D = %s' %D)

       tosum = []
       for k in range(len(dspe.data)):
            jaj = ((dspe[k] - mspe[k])/  np.sqrt(dspe.var[k]))**2
            '''if mspe[k] > 0.000001:
               jaj = ((dspe[k] - mspe[k])/  np.sqrt(dspe.var[k]))**2
            else:
               jaj = 0'''
            print('s%i' %i, dspe[k], mspe[k])
            tosum.append(jaj)

       tosum = np.array(tosum)
       suma = np.sum(tosum)
       print('t%i' %i, suma)
      # chi = suma
       chi = round(suma/len(dspe.data),3)

       print('aaaa', j)


       ax1.set_yticks([0,1])

       plt.subplots_adjust(wspace=0,hspace=0)
       ax1.set_ylim(0,1.4)

       ax1.text(0.03, 0.90, 'N = 10^%s' % hs[i] + ' / ' + 'chi_sq = %s' % chi , transform=ax1.transAxes, fontsize=8,
            verticalalignment='top', bbox=props)

plt.show()
