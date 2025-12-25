import sys
import matplotlib.pyplot as plt
import numpy as np

# read timeseriesoutput
if not len(sys.argv) == 2:
    print("pass filename of .dat output to plot at command line")
    sys.exit(1)
else:
    fn = sys.argv[1]
    print(f"opening {fn} for reading hydrograph output")

with open(fn, "r", encoding="ascii") as fh:
    lines = fh.readlines()

header = lines[1]
outputtype = int(lines[0].split()[1])  # which output_flag was used
if outputtype < 10 or outputtype > 11:
    print("ERROR: only output_flag = {10,11} works with this script")
    sys.exit(1)
    
nloc = int(header.lstrip("#").split("locations")[0].strip())
nt = int(header.split("locations")[1].strip().rstrip("times").strip())


###########################
# whole-file header: 2 line
############################
# 10 LT-AEM time series output   -*-auto-revert-*-
# 8 locations 100 times

##########################
# each time header: 2 lines
##########################
# location: x= -1.2500E+00 y=  1.2500E+00   NW LOC
#     time              head                  velx                vely                 deriv

t = []

# containers for numpy arrays of output
mdv = []
hv = []
dhv = []
if outputtype == 10:
    vv = []


# initial header
idx = 2

for loc in range(nloc):
    # header at beginning of each time
    header = lines[idx]
    xloc = float(header.split("x=")[1].split("y=")[0].strip())
    yloc = float(header.split("y=")[1].split()[0].strip())
    name = " ".join(header.split()[6:])

    mdv.append({"xloc": xloc, "yloc": yloc, "name": name})
    idx += 2

    f = []
    for i in range(nt):
        f.append([float(x) for x in lines[idx].split()])
        idx += 1  # empty rows between locations
    idx += 2

    r = np.array(f)
    if loc == 0:
        t = r[:, 0]
    hv.append(r[:, 1])
    if outputtype == 10:
        vv.append(r[:, 2:4])
        dhv.append(r[:, 4])
    else:
        dhv.append(r[:, 2])

if outputtype == 10:
    nplotcol = 3
else:
    nplotcol = 2

fig, ax = plt.subplots(nplotcol, 1, figsize=(8, 8), constrained_layout=True)

colors = ["red", "green", "blue", "pink", "magenta", "cyan", "orange", "purple"]
ncol = len(colors)
lines = ["-", "--", "-.", ":"]

plt.suptitle(f"data from '{fn}'")

# loc = 2 # which of the locations to plot?
for loc in range(nloc):

    # two/three side-by-side timeseries plots
    ax[0].plot(
        t,
        (hv[loc]),
        linestyle=lines[int(loc / ncol)],
        color=colors[loc % ncol],
        label=(
            f"{mdv[loc]['name']} (x={mdv[loc]['xloc']:.5g}, y={mdv[loc]['yloc']:.5g})"
        ),
    )
    ax[0].set_xscale("log")
    ax[0].set_yscale("linear")
    # ax[0].set_ylim([0.001,3])
    #ax[0].set_xlabel("time")
    ax[0].set_ylabel("head")
    # ax[0].axis('equal')
    ax[0].grid(True)

    ax[1].plot(t, dhv[loc])
    ax[1].set_xscale("log")
    ax[1].set_yscale("linear")
    #ax[1].set_xlabel("time")
    ax[1].set_ylabel("$\\partial$ head/$\\ln t$")
    # ax[1].set_ylim([-20, 20])
    ax[1].grid(True)
    # ax[1].axis('equal')

    if outputtype == 10:
        if loc == 0:
            labelx = "$v_x$"
            labely = "$v_y$"
            labela = "$|v|$"
        else:
            labelx, labely, labela = (None, None, None)
        ax[2].plot(t, vv[loc][:, 0], "r:", label=labelx)
        ax[2].plot(t, vv[loc][:, 1], "g--", label=labely)
        ax[2].plot(t, np.abs(vv[loc][:, 0] + vv[loc][:, 1] * 1j), "k-", label=labela)
        ax[2].set_xscale("log")
        ax[2].set_yscale("linear")
        ax[2].set_xlabel("time")
        ax[2].set_ylabel("velocity")
        ax[2].grid(True)

ax[0].legend(loc=0, fontsize="x-small", ncol=2)
#ax[0].set_title("head")
# ax[1].legend(loc=0)
#ax[1].set_title("derivative")
if outputtype == 10:
    ax[2].legend(loc=0, fontsize="x-small")
    #ax[2].set_title("velocity")

plt.savefig(f"plot_{fn.replace('.dat','.png')}")
