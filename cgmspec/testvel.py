import numpy as np
from matplotlib import pyplot as plt

def los_vel(inc, y, D, alpha, vR=196, hv=5000):
    """
    See Ho et al. 2017

    :param y: np.array, distance of the line-of-sight to the semi-major axis (x-axis) along the y-axis
    :param D: float, impact parameter in kpc
    :param alpha: float, angle between the major axis and the line-of-sight, clockwise, in degrees
    :param vR: maximum velocity of rotation of the disk in km/s
    :param hv: velocity scale height in kpc
    :return: line-of-sight velocity for a cloud in the given position

    """
    incl_rad = np.radians(inc)
    al_rad = np.radians(alpha)

    R = D * np.sqrt(1+(np.sin(al_rad)**2)*np.tan(incl_rad)**2)
    vrot = (2/np.pi)*np.arctan2(R,3)

    x0 = D * np.cos(al_rad)  # this is p in Ho et al. 2019, fig 10.
    y0 = D * np.sin(al_rad) / np.cos(incl_rad)  # this is y0 in the same fig.
    if x0>=0:
          a = np.sin(incl_rad) / np.sqrt(1 + (y/x0)**2)
    else:
        a = - np.sin(incl_rad) / np.sqrt(1 + (y/x0)**2)

    b = np.exp(-np.fabs(y - y0) / hv * np.tan(incl_rad))
    #print(b)
    vr = vR*vrot*a*b


    return(vr)

x = np.linspace(0,30, 1000)

y = (2/np.pi)*np.arctan2(x,3)

plt.plot(x,y)
plt.show()

Ds = np.linspace(0,30,1000)
alphas = np.zeros(1000)
vels = []

for i in range(len(Ds)):
    veli = los_vel(90, 0, Ds[i], 0)
    vels.append(veli)

plt.plot(Ds, vels)
plt.show()
