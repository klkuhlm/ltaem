import numpy as np
import subprocess 

nx,ny = (30,30)

xv = np.linspace(-0.5,1.5,nx)
yv = np.linspace(-0.5,1.5,ny)

r = np.empty((ny,nx))

for i in range(ny):
    for j in range(nx):

        fh = open('particles.in','w')
        fh.write("""1  1     :: #particles, streakSkip
T        :: forward tracking?
1        :: integration scheme
1.0D-4   :: RKM tolerance
0.025    :: RKM max step length
1.0D-6  :: RKM min dt step size
1.0D-4   :: initial dt step size
%.7f      :: starting x
%.7f      :: starting y
1.1      :: starting time
2.0      :: ending time
F        :: beginning inside an element?""" % (xv[j],yv[j]))

        fh.close()

        rtn = subprocess.call(['./ltaem','particle_track.in'])

        
        

