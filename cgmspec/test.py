"""Module for testing cgmspec"""
import numpy as np
from cgmspec.disco import Disco
from cgmspec import utils as csu
from matplotlib import pyplot as plt


lam1 = np.arange(4000,5500,0.1)

D = 10  # in kpc
alpha = 5.  # in degrees
incl = 10.  # in degrees
R = 30.  # disk radius in kpc
h = 5.  # disk height in kpc

model_1 = Disco(R, h, incl, Rcore=0.1)

'''model_1.plotspecandelipse(D, alpha, lam1)
plt.show()

 test prob_hit'''
r = np.arange(0,200, 1)
prob_util = csu.prob_hit(r, rmin=1, rmax=R)
prob_disc = model_1.prob_hit(r)
plt.plot(r,prob_util, 'b')
plt.plot(r, prob_disc+0.1, 'r')
plt.ylabel('Prob.')
plt.xlabel('r (kpc)')
plt.show()


# test los_vel
'''y = np.linspace(-1.1*R, 1.1*R, 1000)
vel_1 = model_1.los_vel(y, D, alpha, vR=200., hv = 5000)
plt.plot(y, vel_1)
plt.show()

Ds  = np.asarray((20,20))
alphas = np.asarray((90,270))
model_1.plotmanylos(Ds,alphas,lam1)'''
