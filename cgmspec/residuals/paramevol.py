from mpdaf.obj import Cube
import numpy as np
import pylab as plt
from astropy import constants as const
from scipy.stats import  moment

lam1 = np.arange(4825.12, 4886.37+1.25, 1.25)

modelcube11 = Cube('paramc_i_p1_0.fits')
modelcube12 = Cube('paramc_i_p1_1.fits')
modelcube13 = Cube('paramc_i_p1_2.fits')
modelcube14 = Cube('paramc_i_p1_3.fits')
modelcube15 = Cube('paramc_i_p1_4.fits')
modelcube16 = Cube('paramc_i_p1_5.fits')
modelcube17 = Cube('paramc_i_p1_6.fits')
modelcube18 = Cube('paramc_i_p1_7.fits')
modelcube19 = Cube('paramc_i_p1_8.fits')


modelcube21 = Cube('paramc_i_p2_0.fits')
modelcube22 = Cube('paramc_i_p2_1.fits')
modelcube23 = Cube('paramc_i_p2_2.fits')
modelcube24 = Cube('paramc_i_p2_3.fits')
modelcube25 = Cube('paramc_i_p2_4.fits')
modelcube26 = Cube('paramc_i_p2_5.fits')
modelcube27 = Cube('paramc_i_p2_6.fits')
modelcube28 = Cube('paramc_i_p2_7.fits')
modelcube29 = Cube('paramc_i_p2_8.fits')


modelcube31 = Cube('paramc_i_p3_0.fits')
modelcube32 = Cube('paramc_i_p3_1.fits')
modelcube33 = Cube('paramc_i_p3_2.fits')
modelcube34 = Cube('paramc_i_p3_3.fits')
modelcube35 = Cube('paramc_i_p3_4.fits')
modelcube36 = Cube('paramc_i_p3_5.fits')
modelcube37 = Cube('paramc_i_p3_6.fits')
modelcube38 = Cube('paramc_i_p3_7.fits')
modelcube39 = Cube('paramc_i_p3_8.fits')

#cents = [np.load('cent_h_p1.npy'),np.load('cent_h_p2.npy') ,np.load('cent_h_p3.npy')  ]

models1 = [[modelcube11, modelcube12,  modelcube13, modelcube14, modelcube15,  modelcube16, modelcube17, modelcube18,  modelcube19],[modelcube21, modelcube22,  modelcube23, modelcube24, modelcube25,  modelcube26, modelcube27, modelcube28,  modelcube29],[modelcube31, modelcube32,  modelcube33, modelcube34, modelcube35,  modelcube36, modelcube37, modelcube38,  modelcube39]]
inclis  = [0,10,20,30,40,50,60,70,80]

def eqw(spe):
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
    return(aver)

ews = []
mins = []
cent = []
variance = []
skewness = []

for i in range(len(models1)):
    ewsi = []
    minsi = []
    centi = []
    varii = []
    skei = []

    for j in range(len(models1[0])):
        spec = models1[i][j][:,15,8]
        flux = spec.data[3:46]
        ew = eqw(spec)
        mini = min(spec.data)
        cen = centroid(spec)
        var = moment(flux,2)
        ske = moment(flux,3)
        ewsi.append(ew)
        minsi.append(mini)
        centi.append(cen)
        varii.append(var)
        skei.append(ske)
    ews.append(ewsi)
    mins.append(minsi)
    cent.append(centi)
    variance.append(varii)
    skewness.append(skei)

plt.figure()
plt.plot(lam1, models1[0][1][:,15,8].data, 'r')
plt.plot(lam1, models1[1][1][:,15,8].data, 'b')
plt.plot(lam1, models1[2][1][:,15,8].data, 'g')
plt.show()

plt.figure()
plt.plot(lam1, models1[0][4][:,15,8].data, 'r')
plt.plot(lam1, models1[1][4][:,15,8].data, 'b')
plt.plot(lam1, models1[2][4][:,15,8].data, 'g')
plt.show()

plt.figure()
plt.plot(lam1, models1[0][8][:,15,8].data, 'r')
plt.plot(lam1, models1[1][8][:,15,8].data, 'b')
plt.plot(lam1, models1[2][8][:,15,8].data, 'g')
plt.show()

'''spe = modelcube35[:,15,8]

green = []
lgreen = []
blue = []
lblue = []
spec = spe.data[3:46]

for i in range(len(spec)):
    spec = spe.data[3:46]
    if spec[i] < 0.994:
        green.append(spec[i])
        lgreen.append(lam1[i])
    else:
        blue.append(spec[i])
        lblue.append(lam1[i])

plt.figure()
plt.scatter(lgreen,green,c='g')
plt.scatter(lblue, blue, c='b')
plt.show()'''


fig = plt.figure()
plt.suptitle('Inclination', fontsize=14)
ax1 = fig.add_subplot(2, 3, 1)
ax1.set_title('equivalent widths')
ax1.plot(inclis, ews[0], 'r', label='D = 3.847 kpc')
ax1.plot(inclis, ews[1], 'b', label='D = 12.928 kpc')
ax1.plot(inclis, ews[2], 'g', label='D = 23.223 kpc')
#ax1.set_xscale('log')
ax1.set_xlabel('inclination')
ax1.set_ylabel('eq.width')
ax1.legend()
ax2 = fig.add_subplot(2, 3, 2)
ax2.set_title('minimum absorption')
ax2.plot(inclis, mins[0], 'r')
ax2.plot(inclis, mins[1], 'b')
ax2.plot(inclis, mins[2], 'g')
#ax2.set_xscale('log')
ax2.set_xlabel('inclination')
ax2.set_ylabel('min. flux')
ax2 = fig.add_subplot(2, 3, 4)
ax2.set_title('centroid')
ax2.plot(inclis, cent[0], 'r')
ax2.plot(inclis, cent[1], 'b')
ax2.plot(inclis, cent[2], 'g')
#ax2.set_xscale('log')
ax2.set_xlabel('inclination')
ax2.set_ylabel('centroid in velocity')
ax2 = fig.add_subplot(2, 3, 5)
ax2.set_title('variance')
ax2.plot(inclis, variance[0], 'r')
ax2.plot(inclis, variance[1], 'b')
ax2.plot(inclis, variance[2], 'g')
#ax2.set_xscale('log')
ax2.set_xlabel('inclination')
ax2.set_ylabel('variance')
ax2 = fig.add_subplot(2, 3, 6)
ax2.set_title('skewness')
ax2.plot(inclis, skewness[0], 'r')
ax2.plot(inclis, skewness[1], 'b')
ax2.plot(inclis, skewness[2], 'g')
#ax2.set_xscale('log')
ax2.set_xlabel('inclination')
ax2.set_ylabel('skewness')
fig.tight_layout()
plt.show()
