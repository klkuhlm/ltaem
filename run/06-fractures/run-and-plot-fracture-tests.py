import numpy as np
import matplotlib.pyplot as plt
import subprocess 
from itertools import product

PIPE = subprocess.PIPE

input_string = """T  F  F  T  T  10  ::  calc?, particle?, contour?, log(t)_deriv?, Qcalc?, output_flag
%s
1.0D0   %.5g   %.5g       ::  por, k, Ss 
0   1.0D0   1.0D-4   1.0D0   :: LEAKY:      leakFlag, K2, Ss2, b2
False 1.5D-1  2.0D0  1.0D0   :: UNCONFINED: unconfinedFlag?, Sy, Kz, b
True  %.5g   %.5g   1   %.5g   %.5g  :: MultiPORO  Flag?, matrixSs, lambda, diffusion idx, Dm, LD
1  1  100 :: nx, ny, nt
FRAC01
 0.0  0.0  0.0  2.0  :: x 
-0.25 0.25 0.75 2.0  :: y
LOGVEC -6.0 8.0 :: t
1.0D-6  1.0D-8  12  :: alpha, tolerance, M
0  fractures_circles.in    :: number of circular elements, circle data file
3  20 fractures_ellipses.in    :: number of elliptical elements, ellipse MS, ellipse data file
not_used   ::  particle data file

"""

n = 3

lamVec = np.logspace(-10,1,n)
DmVec = np.logspace(-8,2,n)
LDVec = np.logspace(-1,2,n)
omegaVec = np.logspace(-6,-1,n)

r =  np.empty((n,n,n,n,400))

fk = 1.0E-3
fSs = 1.0E-5

infn = 'fractures_timeseries-pydrive.in'
exe = ['./ltaem',infn]
outfn = 'fractures_timeseries-pydrive.dat'
qfn = infn.replace('.in','.Q')

for i,j,k,l in product(range(n),repeat=4):
    
    lam = lamVec[i]
    Dm = DmVec[j]
    LD = LDVec[k]
    omega = omegaVec[l]
    
    print i,j,k,l,lam,Dm,LD,omega

    # compute matrix Ss from omega and fracture Ss
    matrixSs = fSs*(1/omega - 1)

    # variables for substituting into input file
    v = (outfn,fk,fSs,matrixSs,lam,Dm,LD)

    fh = open(infn,'w')
    fh.write(input_string % v)
    fh.close()

    p = subprocess.Popen(exe,stdout=PIPE,stderr=PIPE)
    stdout,stderr = p.communicate()

    r[i,j,k,l,:] = np.loadtxt(qfn,usecols=(1,))
    
print 'saving ...'
np.savez('r.npz',r)
