"""Module for testing cgmspec"""
import numpy as np
from cgmspec.disco import Disco
from cgmspec import utils as csu
from matplotlib import pyplot as plt


'''x = np.linspace(1,100,10)
y = np.linspace(2,20,10)

fig = plt.figure()
plot = plt.plot(x,y)

fig.savefig('test1/foo.png', bbox_inches='tight')'''


'''t1  = np.linspace(0,360, 100)
t2  = np.radians(t1)
t3 = np.cos(t2)
t4 = -np.sin(t2)

plt.plot(t1,t3,'b')
plt.plot(t1,t4,'r')
plt.show()'''

lam1 = np.arange(4000,5500,0.1)

D1 = 10
D2 = 10
D3 = 5 # in kpc
alpha1 = 20  # in degrees
incl1 = 10.  # in degrees
R = 30.  # disk radius in kpc
h = 5.  # disk height in kpc

incli2 = 45
incli3 = 80

alpha2 = 110
alpha3 = 200

model1 = Disco(R, h, incl1, Rcore=0.1)
model2 = Disco(R, h, incli2, Rcore=0.1)
model3 = Disco(R, h, incli3, Rcore=0.1)

'''fig1 = model1.plotelipse(D1,alpha1)
fig2 = model2.plotelipse(D1,alpha2)
fig3 = model3.plotelipse(D3,alpha3)'''

#plot the elipses
'''fig1.savefig('test1/model1.png', bbox_inches='tight')
fig2.savefig('test1/model2.png', bbox_inches='tight')
fig3.savefig('test1/model3.png', bbox_inches='tight')'''

#Test for number of iterations
def plotiter(model,D,alpha,niter):
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    fig = plt.figure()
    grid = plt.GridSpec(2,2)
    p1 = fig.add_subplot(grid[0, 0])
    p2 = fig.add_subplot(grid[1, 0])
    p3 = fig.add_subplot(grid[0, 1])
    p4 = fig.add_subplot(grid[1, 1])
    s1 = model.averagelosspec(D, alpha, lam1, niter)
    s2 = model.averagelosspec(D, alpha, lam1, niter)
    s3 = model.averagelosspec(D, alpha, lam1,  niter)
    s4 = model.averagelosspec(D, alpha, lam1, niter)
    p1.plot(s1[0][0],s1[1][0])
    p1.text(0.05, 0.95, 'aver.clouds:%s' % s1[2], transform=p1.transAxes, fontsize=10,
        verticalalignment='top', bbox=props)
    p1.set_ylim(-0.05, 1.05)
    p1.set_xlim(-300, 300)
    p1.axvline(ls='--', lw=1)
    p1.axhline(ls='--', lw=1)
    p1.set_xlabel('LOS vel [km/s]')
    p1.set_ylabel('Norm. Flux')
    p2.plot(s2[0][0],s2[1][0])
    p2.text(0.05, 0.95, 'aver.clouds:%s' % s2[2], transform=p2.transAxes, fontsize=10,
        verticalalignment='top', bbox=props)
    p2.set_ylim(-0.05, 1.05)
    p2.set_xlim(-300, 300)
    p2.axvline(ls='--', lw=1)
    p2.axhline(ls='--', lw=1)
    p2.set_xlabel('LOS vel [km/s]')
    p2.set_ylabel('Norm. Flux')
    p3.plot(s3[0][0],s3[1][0])
    p3.text(0.05, 0.95, 'aver.clouds:%s' % s3[2], transform=p3.transAxes, fontsize=10,
        verticalalignment='top', bbox=props)
    p3.set_ylim(-0.05, 1.05)
    p3.set_xlim(-300, 300)
    p3.axvline(ls='--', lw=1)
    p3.axhline(ls='--', lw=1)
    p3.set_xlabel('LOS vel [km/s]')
    p3.set_ylabel('Norm. Flux')
    p4.plot(s4[0][0],s4[1][0])
    p4.text(0.05, 0.95, 'aver.clouds:%s' % s4[2], transform=p4.transAxes, fontsize=10,
        verticalalignment='top', bbox=props)
    p4.set_ylim(-0.05, 1.05)
    p4.set_xlim(-300, 300)
    p4.axvline(ls='--', lw=1)
    p4.axhline(ls='--', lw=1)
    p4.set_xlabel('LOS vel [km/s]')
    p4.set_ylabel('Norm. Flux')
    return(fig)

e1 = plotiter(model1, D1, alpha1, 1)

iters = [1,10,100,1000]

for i in range(len(iters)):
    '''e1 = plotiter(model1, D1, alpha1, iters[i])
    e2 = plotiter(model2, D2, alpha2, iters[i])'''
    e3 = plotiter(model3, D3, alpha3, iters[i])
    '''e1.savefig('test1/itermod%s' %'mod1'+ '%s' %iters[i] + 'median' + '.png', bbox_inches='tight')
    e2.savefig('test1/itermod%s' %'mod2'+ '%s' %iters[i] + 'median' + '.png', bbox_inches='tight')''
    e3.savefig('test1/itermod%s' %'mod3'+ '%s' %iters[i] + 'median' + '.png', bbox_inches='tight')


#test other parameters
'''def testplot(model, D, alpha):
     props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
     fig = plt.figure()
     ax1 = plt.subplot(2,1,1)
     ax2 = plt.subplot(2,1,2, sharey= ax1, sharex=ax1)
     ax1.set_ylim(-0.05, 1.05)
     ax1.set_xlim(-300, 300)
     fig.subplots_adjust(wspace=0)
     s = model.averagelosspec(D, alpha, lam1, 10, 1, 2)
     ax1.plot(s[0][0], s[1][0])
     ax2.plot(s[0][1], s[1][1], 'g')
     ax1.text(0.05, 0.95, 'aver.clouds:%s' % s[2], transform=ax1.transAxes, fontsize=10,
         verticalalignment='top', bbox=props)

     return(fig)


e1 = testplot(model1, D1, alpha1)
e2 = testplot(model2, D2, alpha2)
e3 = testplot(model3, D3, alpha3)
e2.savefig('test3/mod%s' %'2'+ 'b%s' %'5' + '.png', bbox_inches='tight')
e3.savefig('test3/mod%s' %'3'+ 'b%s' %'5' + '.png', bbox_inches='tight')
e1.savefig('test3/mod%s' %'1'+ 'b%s' %'5' + '.png', bbox_inches='tight')
#test prob_hit'''




''''prob_util = csu.prob_hit(r, rmin=1, rmax=R)
prob_disc = model_1.prob_hit(r)
plt.plot(r,prob_util, 'b')
plt.plot(r, prob_disc+0.1, 'r')
plt.ylabel('Prob.')
plt.xlabel('r (kpc)')
plt.show()'''


# test los_vel
'''y = np.linspace(-1.1*R, 1.1*R, 1000)
vel_1 = model_1.los_vel(y, D, alpha, vR=200., hv = 5000)
plt.plot(y, vel_1)
plt.show()'''

#Many LOS
'''Ds  = np.asarray((20,15,10,5,5,10,15,20))
alphas = np.asarray((45,45,45,45,-135,-135,-135,-135))
model_1.plotmanylos(Ds,alphas,lam1)'''

#model_1.plot_aver_spax(10, 40, lam1, 3, 6)
'''flu = model_1.averagelosspec(D, alpha, lam1, 1000)
plt.xlim(-200,200)
plt.plot(flu[0],flu[1])
plt.show()'''
