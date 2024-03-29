import numpy as np
import matplotlib.pyplot as plt
import subprocess 
from itertools import product
from shutil import move

PIPE = subprocess.PIPE

input_string = """T  F  F  T  T  10  ::  calc?, particle?, contour?, log(t)_deriv?, Qcalc?, output_flag
%s
1.0D0   %.5g   %.5g       ::  por, k, Ss 
0   1.0D0   1.0D-4   1.0D0   :: LEAKY:      leakFlag, K2, Ss2, b2
False 1.5D-1  2.0D0  1.0D0   :: UNCONFINED: unconfinedFlag?, Sy, Kz, b
True  %.5g   %.5g   1   %.5g   %i  :: MultiPORO  Flag?, matrixSs, lambda, diffusion idx, Dm, Nterms
0  0  %i :: nx, ny, nt
FRAC01
 0.0  0.0  0.0  2.0  :: x 
-0.25 0.25 0.75 2.0  :: y
LOGVEC -6.0 7.0 :: t
1.0D-16  1.0D-15  20  :: alpha, tolerance, M
1  fractures_circles.in    :: number of circular elements, circle data file
0  20 fractures_ellipses.in    :: number of elliptical elements, ellipse MS, ellipse data file
not_used   ::  particle data file

"""
nlines = 1
ntobs = 100
nvals = 5

lamVec =   np.logspace(-5,10,nvals)
DmVec =    np.logspace(-10,3,nvals)
LDVec =    np.logspace(0, 3,nvals)
#Dm = 1.0
#LD = 1.0

omegaVec = np.logspace(-8,-1,nvals)

r =  np.empty((nvals,nvals,nvals,nvals,nlines*ntobs))
#r2 =  np.empty((nvals,nvals,nlines*ntobs))

fk = 1.0E-3
fSs = 1.0E-5

infn = 'fractures_timeseries-points.in'
exe = ['./ltaem',infn]
outfn = 'fractures_timeseries-points.dat'
qfn = infn.replace('.in','.Q')

for i,j,k,l in product(range(nvals),repeat=4):
    
    lam = lamVec[i]
    #omega = omegaVec[j]
    Dm = DmVec[j]
    omega = omegaVec[k]
    LD = LDVec[l]
    
    print i,j,k,l,lam,Dm,omega,LD
    #print i,j,lam,omega

    # compute matrix Ss from omega and fracture Ss
    matrixSs = fSs*(1/omega - 1)

    # variables for substituting into input file
    v = (outfn,fk,fSs,matrixSs,lam,Dm,LD,ntobs)

    fh = open(infn,'w')
    fh.write(input_string % v)
    fh.close()

    p = subprocess.Popen(exe,stdout=PIPE,stderr=PIPE)
    stdout,stderr = p.communicate()

    fh = open('debug.stdout','w')
    fh.write(stdout)
    fh.close()

    fh = open('debug.stderr','w')
    fh.write(stderr)
    fh.close()

    r[i,j,k,l,:] = np.loadtxt(qfn,usecols=(1,))

    uniquename = '-%i-%i-%i-%i' % (i,j,k,l)
    #uniquename = '-%i-%i' % (i,j)

    move(outfn,outfn + uniquename)
    move(qfn,qfn + uniquename)
    
print 'saving ...'
np.savez('r4.npz',r)
