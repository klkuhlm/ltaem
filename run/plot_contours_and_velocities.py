import numpy as np
import matplotlib.pyplot as plt
from sys import argv

if argv > 1:
    basefn = argv[1]
else:
    print 'pass base filename as argument'
    exit

particle = False

t = np.atleast_1d(np.loadtxt(basefn+'_t.dat'))
nt = t.shape[0]
X = np.loadtxt(basefn+'_x.dat')
Y = np.loadtxt(basefn+'_y.dat')
ny = X.shape[0]
nx = X.shape[1]

print nx,ny,nt

h = np.zeros((ny,nx,nt))
v = np.zeros((ny,nx,2,nt))

for val in range(nt):
    h[:,:,val] = np.loadtxt('%s_head_%4.4i.dat' % (basefn,val+1))
    v[:,:,0,val] = np.loadtxt('%s_velx_%4.4i.dat' % (basefn,val+1))
    v[:,:,1,val] = np.loadtxt('%s_vely_%4.4i.dat' % (basefn,val+1))

# read in data defining boundaries of elements
fh = open(basefn.replace('_contour.xyz','.geom'),'r')
junk = fh.readline() # throw away first line
elements = []
for row in fh:
    if row[0] == '#':
        # comment character, beginning of new element
        elements.append([])
    else:
        elements[-1].append([float(x) for x in row.split()])

fh.close()
del elements[-1] # EOF comment at end

if particle:
    fh = open(basefn.replace('_contour.xyz','_particles.dat'),'r')
    particles = []
    junk = fh.readline()
    for row in fh:
        if row[0] == '#':
            particles.append([])
        else:
            particles[-1].append([float(x) for x in row.split()])
    fh.close()
    del particles[-1]

for val in range(nt):
    plt.figure(1)
    plt.subplot(121)
    plt.contour(X,Y,h[:,:,val],10)
    plt.axis('image')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.plot([0.0,-0.2],[-0.1,0.1],'k*') 
    for el in elements:
        e = np.array(el[:][:-2])
        plt.plot(e[:,0],e[:,1],'r-',lw=1)
    if particle:
        print 'len',len(particles)
        for part in particles:
            p = np.array(part[:][:-2])
            plt.plot(p[:,1],p[:,2],'k-',lw=0.5)
    
    plt.subplot(122)
    plt.contourf(X,Y,np.log10(np.abs(v[:,:,0,val]+v[:,:,1,val]*1j)),20)
##    plt.contour(X,Y,h[:,:,val],levels=[0])
    plt.quiver(X[::4,::4],Y[::4,::4],v[::4,::4,0,val],v[::4,::4,1,val],scale=1000)

    
    plt.axis('image')
    plt.xlabel('X')
    for el in elements:
        e = np.array(el[:][:-2])
        plt.plot(e[:,0],e[:,1],'r-',lw=1)
    
    plt.savefig('%s_%4.4i_.eps' %(basefn,val))
    plt.close(1)
