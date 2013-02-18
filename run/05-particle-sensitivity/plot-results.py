import numpy as np
import matplotlib.pyplot as plt

x = np.load('x.npy')
y = np.load('y.npy')

nx = x.size
ny = y.size

X,Y = np.meshgrid(x,y)

#r = np.load('def-grad-tensor.npy')
#C = np.load('cauchy-green-tensor.npy')
val = np.load('eigenvalue-field.npy')
vec = np.load('eigenvector-field.npy')

sl = np.empty_like(val)

for i in range(ny):
    for j in range(nx):

        # select eigenvector associated with smallest eigenvalue
        # these are tangent to strainlines
        if val[i,j,0] < val[i,j,1]:
            sl[i,j,:] = vec[i,j,:,0]
        else:
            sl[i,j,:] = vec[i,j,:,1]

plt.figure(1,figsize=(14,14))
plt.quiver(X,Y,sl[:,:,0],sl[:,:,1])
plt.xlabel('X')
plt.ylabel('Y')
plt.savefig('strainlines.eps')
