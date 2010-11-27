import numpy as np
import subprocess
import matplotlib.pyplot as plt
##import matplotlib.colors as clr

PIPE = subprocess.PIPE

nk = 12
nrows = 40
ncols = 40

kvec = np.logspace(-8,0,nk)[::-1] # permeability (large to small)
theta = np.pi/4.0

# want one corner of line to stay in same place as line grows
x0 = 0.21
y0 = 0.21

fvec = np.zeros_like(kvec)

xc = fvec*np.cos(theta) + x0 # x-component of center of line 
yc = fvec*np.sin(theta) + y0 # y-component of center of line 

x1 = 2*fvec*np.cos(theta) + x0
y1 = 2*fvec*np.sin(theta) + y0

pltx = np.zeros((nrows,ncols))
plty = np.zeros((nrows,ncols))
out = np.zeros((nrows,ncols,nk,3),np.float64)

# how to vary N and ms with kappa???
N = 10
cutoff = 1.0E-6

for j,(x,y,k) in enumerate(zip(xc,yc,kvec)):
   passed = False

   # start small, increase matrix size until stop failing for a given cutoff
   for ms in range(16,73,4):
      
      if not passed:
         fel = open('single_line_ellipses.in','w')
         fel.write("""%i :: N
22 :: M
%i  :: ms
2  :: ibnd
false :: calcin?
false :: storin?
0.0D-1  :: r
%.5e  :: x0
%.5e  :: y0
5.0D-1    :: f
%.22e  :: theta
1.0D-0  :: K
1.0D-0  :: Ss
1.5D-0  :: porosity
0.0D-1  :: area source strength
5.0D-1  :: boundary source strength
001  0.0D0  0.0D0 :: area time flag, tpar1, tpar2
001  0.0D0  0.0D0 :: boundary time flag, tpar1, tpar2
0  :: leaky type flag
1.0D0 :: K2
1.0D-4 :: Ss2
1.0D0 :: b2
False :: unconfined?
1.5D-1  :: Sy
2.0D-2  :: Kz
1.0D0  :: b
1.0D0 :: dimensionless skin
""" % (N,ms,x,y,theta))
         fel.close()

         fin = open('single_line.in','w')
         fin.write("""True  False  True  1  single_line_contour.xyz  dump.out  single_well.elem  single_well.geom    ::  re-compute coefficients?, particle?, contour?, output flag, outputFN, coeffFN (output), elementHierarchy (input), geometryFN (output)
2.0D-1  %.22e  1.0D-4   0  1.0D0  1.0D-4  1.0D0  %i  %.22e ::  por, k, Ss, leakFlag, K2, Ss2, b2, ellipse MS, ellipse cutoff
1.5D-1  2.0D0  False  1.0  :: Sy, Kz, unconfined?, b
40  40  1  :: nx, ny, nt
0.000000000000000000e+00 2.820512820512820901e-02 5.641025641025641801e-02 8.461538461538462008e-02 1.128205128205128360e-01 1.410256410256410520e-01 1.692307692307692402e-01 1.974358974358974561e-01 2.256410256410256721e-01 2.538461538461538880e-01 2.820512820512821039e-01 3.102564102564103199e-01 3.384615384615384803e-01 3.666666666666666963e-01 3.948717948717949122e-01 4.230769230769231282e-01 4.512820512820513441e-01 4.794871794871795601e-01 5.076923076923077760e-01 5.358974358974359919e-01 5.641025641025642079e-01 5.923076923076924238e-01 6.205128205128206398e-01 6.487179487179488557e-01 6.769230769230769607e-01 7.051282051282051766e-01 7.333333333333333925e-01 7.615384615384616085e-01 7.897435897435898244e-01 8.179487179487180404e-01 8.461538461538462563e-01 8.743589743589744723e-01 9.025641025641026882e-01 9.307692307692309042e-01 9.589743589743591201e-01 9.871794871794873361e-01 1.015384615384615552e+00 1.043589743589743657e+00 1.071794871794871984e+00 1.100000000000000089e+00 :: x
0.000000000000000000e+00 2.820512820512820901e-02 5.641025641025641801e-02 8.461538461538462008e-02 1.128205128205128360e-01 1.410256410256410520e-01 1.692307692307692402e-01 1.974358974358974561e-01 2.256410256410256721e-01 2.538461538461538880e-01 2.820512820512821039e-01 3.102564102564103199e-01 3.384615384615384803e-01 3.666666666666666963e-01 3.948717948717949122e-01 4.230769230769231282e-01 4.512820512820513441e-01 4.794871794871795601e-01 5.076923076923077760e-01 5.358974358974359919e-01 5.641025641025642079e-01 5.923076923076924238e-01 6.205128205128206398e-01 6.487179487179488557e-01 6.769230769230769607e-01 7.051282051282051766e-01 7.333333333333333925e-01 7.615384615384616085e-01 7.897435897435898244e-01 8.179487179487180404e-01 8.461538461538462563e-01 8.743589743589744723e-01 9.025641025641026882e-01 9.307692307692309042e-01 9.589743589743591201e-01 9.871794871794873361e-01 1.015384615384615552e+00 1.043589743589743657e+00 1.071794871794871984e+00 1.100000000000000089e+00 :: y
0.125  :: t
1.0D-6  1.0D-8  5  :: alpha, tolerance, M
0  file_not_used      :: number of circular elements, circular elements input file
1  single_line_ellipses.in      :: number of elliptical elements, elliptical elements input file
file_not_used     :: particle inputs file
""" % (k,ms,cutoff))
         fin.close()

         (stdout,stderr) = subprocess.Popen(['./ltaem','single_line.in'],
                                         stdout=PIPE,stderr=PIPE).communicate()

         fout = open('screen_j%3.3i_ms%3.3i_k%.2e.dbg' % (j,ms,k),'w')
         fout.write(stdout)
         fout.close()

         if ('not large enough to meet tolerance' in stdout):
            val = stdout.split('\n')[-1].split('(')[1].rstrip(')').split(',')
            q = np.complex(float(val[0]),float(val[1]))
         else:
            nq = int(stdout.split('\n')[0].split()[1])
            qlines = stdout.split('\n')[3:3+nq]
            qstrs = [line.lstrip('(').rstrip(')').split(',') for line in qlines]
            q = [np.complex(float(row[0]),float(row[1])) for row in qstrs]

         if ('CHECK_CUTOFF' not in stderr) and ('not large enough to meet tolerance' not in stdout):
            print  '%i,%.7e,%.7e,%i,%7e,passed' % (j,k,np.abs(q).max(),ms,cutoff)
   
            if 'gnuplot contour map style output' in stdout:
               if j == 0:
                  pltx = np.loadtxt('single_line_contour.xyz',skiprows=7)[:,0].reshape((nrows,ncols))
                  plty = np.loadtxt('single_line_contour.xyz',skiprows=7)[:,1].reshape((nrows,ncols))
   
               for kk in range(3):
                  out[:,:,j,kk] = np.loadtxt('single_line_contour.xyz',skiprows=7)[:,2+kk].reshape((nrows,ncols))
            else:
               print 'STDOUT::\n',stdout
               print 'STDERR::\n',stderr

            passed = True
         else:
            print '%i,%.7e,%.7e,%i,%7e,failed' % (j,k,np.abs(q).max(),ms,cutoff)
            passed = False

#nrm = clr.Normalize(vmin=0.0,vmax=out[:,:,0,0].max()) # use same colorscale across all plots
   
plt.figure(1,figsize=(28,21))
for r in range(3):
    for c in range(4):
        idx = 4*r + c

        plt.subplot(3,4,idx+1)
        plt.contourf(pltx,plty,out[:,:,idx,0],20)
        plt.plot([x0,x1[idx]],[y0,y1[idx]],'k-')
        plt.axis('image')
        plt.axis([0,1.1,0,1.1])
        plt.colorbar(shrink=0.5)
        plt.title('f=%.2e' % (fvec[idx]))

plt.savefig('compare-lines-head.png')
plt.close(1)

plt.figure(1,figsize=(26,21))
for r in range(3):
    for c in range(4):
        idx = 4*r + c

        plt.subplot(3,4,idx+1)
        plt.quiver(pltx,plty,out[:,:,idx,1],out[:,:,idx,2])
        ax = plt.axis()
        plt.plot([x0,x1[idx]],[y0,y1[idx]],'k-')
        plt.axis('image')
        plt.axis([0,1.1,0,1.1])
        plt.title('f=%.2e' % (fvec[idx]))

plt.savefig('compare-lines-vec.png')

