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
True  %.5g   %.5g   1   %.5g   %.5g  :: MultiPORO  Flag?, matrixSs, lambda, diffusion idx, Dm, LD
1  1  %i :: nx, ny, nt
FRAC01
 0.0  0.0  0.0  2.0  :: x 
-0.25 0.25 0.75 2.0  :: y
LOGVEC -6.0 7.0 :: t
1.0D-6  1.0D-8  10  :: alpha, tolerance, M
0  fractures_circles.in    :: number of circular elements, circle data file
%i  20 fractures_ellipses.in    :: number of elliptical elements, ellipse MS, ellipse data file
not_used   ::  particle data file

"""
nlines = 3
ntobs = 100
n = 4

lamVec =   np.logspace(-10,2,n)
DmVec =    np.logspace(-8, 2,n)
omegaVec = np.logspace(-6,-1,n)

r =  np.empty((n,n,n,nlines*ntobs))

fk = 1.0E-3
fSs = 1.0E-5
LD = 1.0

infn = 'fractures_timeseries-pydrive.in'
exe = ['./ltaem',infn]
outfn = 'fractures_timeseries-pydrive.dat'
qfn = infn.replace('.in','.Q')

for i,j,k in product(range(n),repeat=3):
    
    lam = lamVec[i]
    Dm = DmVec[j]
    omega = omegaVec[k]
    
    print i,j,k,lam,Dm,omega

    # compute matrix Ss from omega and fracture Ss
    matrixSs = fSs*(1/omega - 1)

    # variables for substituting into input file
    v = (outfn,fk,fSs,matrixSs,lam,Dm,LD,ntobs,nlines)

    fh = open(infn,'w')
    fh.write(input_string % v)
    fh.close()

    p = subprocess.Popen(exe,stdout=PIPE,stderr=PIPE)
    stdout,stderr = p.communicate()

    r[i,j,k,:] = np.loadtxt(qfn,usecols=(1,))

    uniquename = '-%i-%i-%i' % (i,j,k)

    move(outfn,outfn + uniquename)
    move(qfn,qfn + uniquename)
    
print 'saving ...'
np.savez('r.npz',r)
