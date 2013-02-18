import numpy as np
import subprocess 

nx,ny = (10,10)

delta = 1.0E-4

xv = np.linspace(-0.5,1.5,nx)
yv = np.linspace(-0.5,1.5,ny)

r = np.empty((ny,nx,2,4))

for i in range(ny):
    for j in range(nx):

        for k,(dx,dy) in enumerate(zip([0.0,delta,0.0],[0.0,0.0,delta])):
            
            fh = open('particles.in','w')
            fh.write("""1  1     :: #particles, streakSkip
T        :: forward tracking?
1        :: integration scheme (1=RKM,2=RK,3=analytic,4=fwd Euler)
1.0D-4   :: RKM tolerance
0.025    :: RKM max step length
1.0D-6  :: RKM min dt step size
1.0D-4   :: initial dt step size
%.7f      :: starting x
%.7f      :: starting y
1.1      :: starting time
1.75      :: ending time
F        :: beginning inside an element?
""" % (xv[j]+dx,yv[j]+dy))

            fh.close()

            fh = open('screen.out','w')

            subprocess.call(['rm','particles.dat'])
            subprocess.call(['./ltaem','particle_track.in'],stdout=fh)
        
            fh.close()

            p = np.loadtxt('particles.dat')

            if k == 0:
                # no shift
                base = p[-1,1:]
                print i,j,xv[j],yv[i]
            elif k == 1:
                # x shift
                r[i,j,0,:] = (base - p[-1,1:])/delta
                #print 'x',base,':',p[-1,1:],':',r[i,j,0,:]
            else:
                # y shift
                r[i,j,1,:] = (base - p[-1,1:])/delta
                #print 'y',base,':',p[-1,1:],':',r[i,j,1,:]

np.save('r.npy',r)

        

