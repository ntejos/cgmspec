"""Module for testing cgmspec"""
import numpy as np
from cgmspec.disco import Disco
from cgmspec import utils as csu
from matplotlib import pyplot as plt


lam1 = np.arange(4000,5500,0.1)

D = 16.8  # in kpc
alpha = 15.  # in degrees
incl = 75.  # in degrees
R = 100.  # disk radius in kpc
h = 10.  # disk height in kpc

model_1 = Disco(R, h, incl)
model_1.plotspecandelipse(D, alpha, incl, R, h, lam1)
plt.show()