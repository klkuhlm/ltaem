import numpy as np
import matplotlib.pyplot as plt

fn = 'multiple_lines_contour.xyz'
fh = open(fn,'r')

num = fh.readlines()[2:4] # get number x and number y
fh.seek(0)

nx = int(num[0].split(':')[-1])
ny = int(num[1].split(':')[-1])

# X,Y,head,velx,vely
# x column changes fastest
res = np.loadtxt(fh,skiprows=7)
fh.close()

# read in data defining boundaries of elements
fh = open(fn.replace('_contour.xyz','.geom'),'r')
junk = fh.readline() # throw away first line
elements = []
for row in fh:
    if row[0] == '#':
        # comment character, beginning of new element
        elements.append([])
    else:
        elements[-1].append([float(x) for x in row.split()])

del elements[-1] # EOF comment at end

X = res[:,0].reshape((ny,nx))
Y = res[:,1].reshape((ny,nx))
h = res[:,2].reshape((ny,nx))
vx = res[:,3].reshape((ny,nx))
vy = res[:,4].reshape((ny,nx))

plt.figure()
plt.subplot(121)
plt.contourf(X,Y,h,50)
plt.axis('image')
plt.xlabel('X')
plt.ylabel('Y')
plt.title('head')
for el in elements:
    e = np.array(el[:][:-2])
    plt.plot(e[:,0],e[:,1],'k-')

plt.subplot(122)
plt.quiver(X,Y,vx,vy)
plt.axis('image')
plt.xlabel('X')
plt.ylabel('Y')
plt.title('velocity')
for el in elements:
    e = np.array(el[:][:-2])
    plt.plot(e[:,0],e[:,1],'k-')

plt.savefig('test.png')
