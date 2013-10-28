import numpy as np
import matplotlib.pyplot as plt
from itertools import product

# multiporosity flow data
mp = np.load('r.npz')['arr_0']

print 'multiporosity results',mp.shape

# warren root flow data
wr = np.load('r2.npz')['arr_0']

print 'warren-root results',wr.shape

# times
t = np.logspace(-6,7,100)

nvals = 4
mpnvars = 4
wrnvars = 2

# parameter values
lamVec =   np.logspace(-2,10,nvals)
DmVec =    np.logspace(-10,2,nvals)
LDVec =    np.logspace(0,  4,nvals)
omegaVec = np.logspace(-8,-1,nvals)

mpvars = [{'vals':lamVec,'var':'$\\lambda$', 0:3},
          {'vals':DmVec,'var':'$D_m$',       0:2},
          {'vals':LDVec,'var':'$L_D$',       0:3},
          {'vals':omegaVec,'var':'$\\omega$',0:2}] 

wrvars = [{'vals':lamVec,'var':'$\\lambda$', 0:2},
          {'vals':omegaVec,'var':'$\\omega$',0:2}]

# first plot everything

fig = plt.figure(1)
ax = fig.add_subplot(111)
for i,j,k,l in product(range(nvals),repeat=mpnvars):
    ax.loglog(t,mp[i,j,k,l,:100],'-')

for i,j in product(range(nvals),repeat=wrnvars):
    ax.loglog(t,wr[i,j,:100],'--')

ax.set_xlabel('$t_D$')
ax.set_ylabel('$q$ at specified head')
plt.savefig('all-pars-compare.png',dpi=300)
plt.close(1)

# plot ranges, holding all other parameters constant


fig = plt.figure(1)

ii = mpvars[0][0]
jj = mpvars[1][0]
kk = mpvars[2][0]
ll = mpvars[3][0]

ax = fig.add_subplot(221)
for i in range(nvals):
    ax.loglog(t,mp[i,jj,kk,ll,:100])
    ax.set_title(mpvars[0]['var'])

ax = fig.add_subplot(222)
for j in range(nvals):
    ax.loglog(t,mp[ii,j,kk,ll,:100])
    ax.set_title(mpvars[1]['var'])

ax = fig.add_subplot(223)
for k in range(nvals):
    ax.loglog(t,mp[ii,jj,k,ll,:100])
    ax.set_title(mpvars[2]['var'])

ax = fig.add_subplot(224)
for l in range(nvals):
    ax.loglog(t,mp[ii,jj,kk,l,:100])
    ax.set_title(mpvars[3]['var'])

plt.savefig('individual-variability-%i-%i-%i-%i.png' % (ii,jj,kk,ll))
plt.close(1)

