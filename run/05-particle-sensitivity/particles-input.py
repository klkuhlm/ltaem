import numpy as np
import subprocess 

nx,ny = (15,15)

delta = 1.0E-5

xv = np.linspace(-0.15,1.15,nx)
yv = np.linspace(-0.15,1.15,ny)

np.save('x.npy',xv)
np.save('y.npy',yv)

r = np.empty((ny,nx,2,2))  # deformation-gradient tensor (2x2 at each location)
C = np.empty((ny,nx,2,2))  # cauchy-green strain tensor (" ")
vec = np.empty((ny,nx,2,2))  # eigenvector field (strainlines + ?)
val = np.empty((ny,nx,2))

for i in range(ny):
    for j in range(nx):

        for k,(dx,dy) in enumerate(zip([0.0,delta,0.0],[0.0,0.0,delta])):
            
            fh = open('particles.in','w')
            fh.write("""1  1     :: #particles, streakSkip
T        :: forward tracking?
1        :: integration scheme (1=RKM,2=RK,3=analytic,4=fwd Euler)
1.0D-4   :: RKM tolerance
0.025    :: RKM max step length
5.0D-7   :: RKM min dt step size
5.0D-5   :: initial dt step size
%.7f     :: starting x
%.7f     :: starting y
1.15      :: starting time
1.25     :: ending time
F        :: beginning inside an element?
""" % (xv[j]+dx,yv[i]+dy))

            fh.close()

            subprocess.call(['rm','particles.dat'])

            fh = open('screen.out','w')
            subprocess.call(['./ltaem','particle_track.in'],stdout=fh)
            fh.close()

            p = np.loadtxt('particles.dat')
            # columns are : 0=time, 1=x_1, 2=y_1, 3=v_{x,1}, 4=v_{x,1}

            if k == 0:
                # no shift
                base = p[-1,1:3]  # get final x/y positions
                print i,j,xv[j],yv[i]
            elif k == 1:
                # x shift 
                r[i,j,0,0] = (p[-1,1] - base[0])/delta
                r[i,j,1,0] = (p[-1,2] - base[1])/delta
            else:
                # y shift 
                r[i,j,0,1] = (p[-1,1] - base[0])/delta
                r[i,j,1,1] = (p[-1,2] - base[1])/delta
                
                # compute Cauchy-Green strain tensor
                C[i,j,:,:] = np.dot(np.transpose(r[i,j,:,:]),r[i,j,:,:])
                val[i,j,:],vec[i,j,:,:] = np.linalg.eig(C[i,j,:,:])

np.save('def-grad-tensor.npy',r)
np.save('cauchy-green-tensor.npy',C)
np.save('eigenvector-field.npy',vec)
np.save('eigenvalue-field.npy',val)

