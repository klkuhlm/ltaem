import matplotlib.pyplot as plt
import numpy as np

numcol = 6

# read contour map output
fn = 'pumping_well_contour.dat'

fh = open(fn,'r')
lines = fh.readlines()
fh.close()
header = lines[:4]
nt = int(header[1].split(':')[1])
nx = int(header[2].split(':')[1])
ny = int(header[3].split(':')[1])

###########################
# whole-file header: 5 lines
############################

# LT-AEM contour map output     -*-auto-revert-*-
# t: 4
# x: 40
# y: 40
# locations:1600

##########################
# each time header: 2 lines
##########################
 # t=   1.00000E-04
#      X           Y               head                velx                  vely                d(head)/d(log(t))

t = []
r = np.empty((nx*ny,numcol))

h = np.empty((nx,ny,nt))
v = np.empty((nx,ny,2,nt))
dh = np.empty((nx,ny,nt))

# initial header
idx = 5

for i in range(nt):
    # header at beginning of each time
    idx += 2
    f = []
    for j in xrange(nx*ny):
        f.append([float(x) for x in lines[idx].split()])
        idx += 1
    idx += 2 # empty rows between times
    r[:,:] = np.array(f)
    
    # array is flattened in Fortran style (left-most column changing fastest)
    h[:,:,i] =   r[:,2].reshape(nx,ny,order='F')
    v[:,:,0,i] = r[:,3].reshape(nx,ny,order='F')
    v[:,:,1,i] = r[:,4].reshape(nx,ny,order='F')
    dh[:,:,i] =  r[:,5].reshape(nx,ny,order='F')
    
X = r[:,0].reshape(nx,ny,order='F')
Y = r[:,1].reshape(nx,ny,order='F')

ts = 2

fig = plt.figure(figsize=(18,6))
plt.subplot(131)
plt.contourf(X,Y,h[:,:,ts],30)
plt.title('head')
plt.axis('equal')
plt.grid()
plt.subplot(132)
plt.quiver(X,Y,v[:,:,0,ts],v[:,:,1,ts])
plt.title('velocity vectors')
plt.axis('equal')
plt.grid()
plt.subplot(133)
plt.contourf(X,Y,dh[:,:,ts],30)
plt.title('head ln(t) derivative')
plt.axis('equal')
plt.grid()
plt.savefig('test-contours.png')

