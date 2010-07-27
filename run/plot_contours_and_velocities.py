import numpy as np
import matplotlib.pyplot as plt

basefn = 'wedge_contour.xyz'
t = np.loadtxt(basefn+'_t.dat')
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

del elements[-1] # EOF comment at end

for val in range(nt):
    plt.figure(1)
#    plt.subplot(121)
    plt.contour(X,Y,h[:,:,val],50)
    plt.grid()
    plt.axis('image')
    plt.xlabel('X')
    plt.ylabel('Y')
#    plt.title('head')
    for el in elements:
        e = np.array(el[:][:-2])
        plt.plot(e[:,0],e[:,1],'r-',lw=3)
    
#    plt.subplot(122)
#    plt.quiver(X,Y,v[:,:,0,val],v[:,:,1,val])
#    plt.axis('image')
#    plt.xlabel('X')
#    plt.ylabel('Y')
#    plt.title('velocity')
#    for el in elements:
#        e = np.array(el[:][:-2])
#        plt.plot(e[:,0],e[:,1],'r-',lw=2)
    
    plt.savefig('%s_%4.4i_.eps' %(basefn,val))
    plt.close(1)
