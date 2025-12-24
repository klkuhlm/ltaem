import numpy as np
import matplotlib.pyplot as plt

fn = "pumping_well_contour.dat"
tvec = []
with open(fn, "r", encoding="ascii") as fh:
    lines = fh.readlines()
    nt = int(lines[1].strip().split(":")[1])
    nx = int(lines[2].strip().split(":")[1])
    ny = int(lines[3].strip().split(":")[1])

    for line in lines:
        if "t= " in line:
            tvec.append(float(line.strip().split("=")[1]))

# this skips blanks and comments
# reads in as one giant matrix
d = np.loadtxt(fn)

print(nt, nx, ny, tvec, d.shape)

t = np.array(tvec)

X = np.reshape(d[: nx * ny, 0], (nx, ny), order="F").transpose()
Y = np.reshape(d[: nx * ny, 1], (nx, ny), order="F").transpose()

h  = np.reshape(d[:, 2], (nx, ny, nt), order="F")
vx = np.reshape(d[:, 3], (nx, ny, nt), order="F")
vy = np.reshape(d[:, 4], (nx, ny, nt), order="F")
dh = np.reshape(d[:, 5], (nx, ny, nt), order="F")

data = [h,dh,vx,vy]

c = []
e = []
with open(fn.replace(".dat", ".geom"), "r", encoding="ascii") as fh:
    lines = fh.readlines()[1:]
    for j, line in enumerate(lines):
        if "#" in line:
            if "EOF" in line:
                break
            f = line.strip().split()
            t = f[1]
            num = int(f[2])
            m = int(f[4])
            dd = []
            for i in range(m):
                # x, y, local angle
                dd.append(
                    [float(x) for x in lines[j + i + 1].lstrip("#").strip().split()]
                )
            if t == "circle":
                c.append({"id": num, "m": m, "xt": np.array(dd)})
            elif t == "ellipse":
                e.append({"id": num, "m": m, "xt": np.array(dd)})
            else:
                print("invalid type", t)

titles = ["head","dh/d(lnt)","vx","vy"]
                
for j in range(nt):
    print(j)
    fig,axes = plt.subplots(2,2,num=1,figsize=(8,8),constrained_layout=True)
    for i,ax in enumerate(axes.ravel()):
        mmax = np.nanmax(data[i][:,:,j])
        mmin = np.nanmin(data[i][:,:,j])
        absmax = max(abs(mmax),abs(mmin))
        #print(j,i,mmax,mmin,absmax)
        CS = ax.contourf(X,Y,data[i][:,:,j].transpose(),levels=20,vmax=absmax, vmin=-absmax, cmap="RdBu_r")
        cbar = fig.colorbar(CS)
        ax.clabel(CS, fontsize=6)
        #ax.streamplot(X,Y,data[2][:,:,j].transpose(),data[3][:,:,j].transpose())
        #ax.quiver(X,Y,data[2][:,:,j].transpose(),data[3][:,:,j].transpose())
        ax.set_title(titles[i])
        ax.set_aspect("equal")

        for cc in c:
            xx = cc["xt"]
            ax.plot(xx[:,0],xx[:,1],"k-")
        for ee in e:
            xx = ee["xt"]
            ax.plot(xx[:,0],xx[:,1],"k-",lw=0.25)            
            
    fig.savefig(f"waves_contour_{j:03}.png",dpi=200)
    plt.close(1)
    
    
