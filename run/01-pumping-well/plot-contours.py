import matplotlib.pyplot as plt
import numpy as np

numcol = 6

# read contour map output
fn = "pumping_well_contour.dat"

with open(fn, "r", encoding="ascii") as fh:
    lines = fh.readlines()

header = lines[:4]
nt = int(header[1].split(":")[1])
nx = int(header[2].split(":")[1])
ny = int(header[3].split(":")[1])

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
#      X  Y  head  velx  vely  d(head)/d(log(t))

t = []
r = np.empty((nx * ny, numcol))

h = np.empty((nx, ny, nt))
v = np.empty((nx, ny, 2, nt))
dh = np.empty((nx, ny, nt))

# initial header
idx = 5

for i in range(nt):
    # header at beginning of each time
    t.append(float(lines[idx].strip().split("=")[1]))
    idx += 2
    f = []
    for j in range(nx * ny):
        f.append([float(x) for x in lines[idx].split()])
        idx += 1
    idx += 2  # empty rows between times
    r[:, :] = np.array(f)

    # array is flattened in Fortran style (left-most column changing fastest)
    h[:, :, i] = r[:, 2].reshape(nx, ny, order="F")
    v[:, :, 0, i] = r[:, 3].reshape(nx, ny, order="F")
    v[:, :, 1, i] = r[:, 4].reshape(nx, ny, order="F")
    dh[:, :, i] = r[:, 5].reshape(nx, ny, order="F")

X = r[:, 0].reshape(nx, ny, order="F")
Y = r[:, 1].reshape(nx, ny, order="F")

for j in range(nt):
    fig, ax = plt.subplots(1,3,num=1,figsize=(18, 6),constrained_layout=True)
    ax[0].contourf(X, Y, h[:, :, j], 30)
    ax[0].set_title(f"head (t={t[j]:.3E})")
    ax[0].set_aspect("equal", "box")
    ax[0].grid(True)
    ax[1].quiver(X, Y, v[:, :, 0, j], v[:, :, 1, j])
    ax[1].set_title("velocity vectors")
    ax[1].set_aspect("equal", "box")
    ax[1].grid(True)
    ax[2].contourf(X, Y, dh[:, :, j], 30)
    ax[2].set_title("head ln(t) derivative")
    ax[2].set_aspect("equal", "box")
    ax[2].grid(True)
    plt.savefig(f"test_contours_{j:03}.png")
    plt.close(1)
