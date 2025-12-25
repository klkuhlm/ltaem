import sys
from glob import glob
import numpy as np
import matplotlib.pyplot as plt

if len(sys.argv) > 1:
    basefn = sys.argv[1]
else:
    print("pass base filename as argument")
    sys.exit(1)

particle = False

t = np.atleast_1d(np.loadtxt(f"{basefn}_t.dat"))
nt = t.shape[0]
X = np.loadtxt(f"{basefn}_x.dat")
Y = np.loadtxt(f"{basefn}_y.dat")
ny = X.shape[0]
nx = X.shape[1]

print(nx, ny, nt)

h = np.zeros((ny, nx, nt))
v = np.zeros((ny, nx, 2, nt))

for val, headfn in enumerate(glob(f"{basefn}_head_*.dat")):
    h[:, :, val] = np.loadtxt(headfn)
    v[:, :, 0, val] = np.loadtxt(headfn.replace("head", "velx"))
    v[:, :, 1, val] = np.loadtxt(headfn.replace("head", "velx"))

# read in data defining boundaries of elements
with open(f"{basefn}.geom", "r", encoding="ascii") as fh:
    junk = fh.readline()  # throw away first line
    elements = []
    for row in fh:
        if row[0] == "#":
            # comment character, beginning of new element
            elements.append([])
        else:
            elements[-1].append([float(x) for x in row.split()])

del elements[-1]  # EOF comment at end

if particle:
    with open(f"{basefn}_particles.dat", "r", encoding="ascii") as fh:
        particles = []
        junk = fh.readline()
        for row in fh:
            if row[0] == "#":
                particles.append([])
            else:
                particles[-1].append([float(x) for x in row.split()])
    del particles[-1]

sc = [1000, 1000, 1000]

for val in range(nt):
    plt.figure(1)
    plt.subplot(121)
    plt.contour(X, Y, h[:, :, val], 20)
    plt.axis("image")
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.plot([0.0, -0.2], [-0.1, 0.1], "k*")
    for el in elements:
        e = np.array(el[:][:-2])
        plt.plot(e[:, 0], e[:, 1], "r-", lw=1)
    if particle:
        print(f"len # particles: {len(particles)}")
        for part in particles:
            p = np.array(part[:][:-2])
            mask = p[:, 0] < t[val]
            plt.plot(p[mask, 1], p[mask, 2], "k-", lw=0.5)

    plt.subplot(122)
    plt.contourf(X, Y, np.log10(np.abs(v[:, :, 0, val] + v[:, :, 1, val] * 1j)), 20)
    plt.quiver(
        X[::4, ::4],
        Y[::4, ::4],
        v[::4, ::4, 0, val],
        v[::4, ::4, 1, val],
        scale=sc[val],
    )

    plt.axis("image")
    plt.xlabel("X")
    for el in elements:
        e = np.array(el[:][:-2])
        plt.plot(e[:, 0], e[:, 1], "r-", lw=1)

    plt.savefig(f"plot_{basefn}_{val:4.4i}_.png", dpi=200)
    plt.close(1)
