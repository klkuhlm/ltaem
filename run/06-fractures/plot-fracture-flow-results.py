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
lamVec =   np.logspace(-10,2,nvals)
DmVec =    np.logspace(-10,2,nvals)
LDVec =    np.logspace(-2, 2,nvals)
omegaVec = np.logspace(-8,-1,nvals)

mpvars = [{'vals':lamVec,'var':'$\\lambda$',0:2},
          {'vals':DmVec,'var':'$D_m$',0:2},
          {'vals':omegaVec,'var':'$\omega$',0:2},
          {'vals':LDVec,'var':'$L_D$',0:2}] 

wrvars = [{'vals':lamVec,'var':'$\\lambda$',0:2},
          {'vals':omegaVec,'var':'$\omega$',0:2}]

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
