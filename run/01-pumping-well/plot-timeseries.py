import matplotlib.pyplot as plt
import numpy as np

# read timeseriesoutput
fn = 'pumping_well_linearQ_hydrographs.dat'

fh = open(fn,'r')
lines = fh.readlines()
fh.close()
header = lines[1]
nloc = int(header.lstrip('#').split('locations')[0].strip())
nt =  int(header.split('locations')[1].strip().rstrip('times').strip())

###########################
# whole-file header: 2 line
############################
# LT-AEM time series output   -*-auto-revert-*-
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
vv = []
dhv = []

# initial header
idx = 2

for loc in range(nloc):
    # header at beginning of each time
    header = lines[idx]
    xloc = float(header.split('x=')[1].split('y=')[0].strip())
    yloc = float(header.split('y=')[1].split()[0].strip())
    name = ' '.join(header.split()[6:])
        
    mdv.append({'xloc':xloc,'yloc':yloc,'name':name})
    idx += 2

    f = []
    for i in range(nt):
        f.append([float(x) for x in lines[idx].split()])
        idx += 1 # empty rows between locations
    idx += 2

    r = np.array(f)
    if loc == 0:
        t = r[:,0]
    hv.append(r[:,1])
    vv.append(r[:,2:4])
    dhv.append(r[:,4])
    

loc = 2 # which of the locations to plot?

h = hv[loc]
v = vv[loc]
dh = dhv[loc]
md = mdv[loc]

# three side-by-side timeseries plots
fig = plt.figure(figsize=(12,6))
ax = fig.add_subplot(131)
ax.loglog(t,np.abs(h),label='head')
ax.set_xscale('linear')
ax.set_yscale('linear')
ax.legend(loc=0)
ax.set_xlabel('time')
ax.set_ylabel('head')
ax.axis('equal')
ax.grid()

ax = fig.add_subplot(132)
ax.loglog(t,dh,label='derivative')
ax.set_xscale('linear')
ax.set_yscale('linear')
ax.legend(loc=0)
ax.set_xlabel('time')
ax.set_ylabel('$\partial$ head/$\\ln t$')
ax.set_ylim([-20,20])
ax.grid()
#ax.axis('equal')

ax = fig.add_subplot(133)
ax.semilogx(t,v[:,0],'r:',label='$v_x$')
ax.semilogx(t,v[:,1],'g--',label='$v_y$')
ax.semilogx(t,np.abs(v[:,0]+v[:,1]*1j),'k-',label='$|v|$')
ax.set_xscale('linear')
ax.set_yscale('linear')
ax.legend(loc=0)
ax.set_xlabel('time')
ax.set_ylabel('velocity')
ax.grid()

plt.suptitle('%s (x=%.5g, y=%.5g)' % (md['name'],md['xloc'],md['yloc']))
plt.tight_layout()
plt.savefig('test-timeseries.png')

