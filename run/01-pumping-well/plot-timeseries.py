import matplotlib.pyplot as plt
import numpy as np

# read timeseriesoutput
fn = 'boise_timeseries.dat'

fh = open(fn,'r')
lines = fh.readlines()
fh.close()
header = lines[1]
outputtype = int(lines[0].split()[1])  # which output_flag was used
nloc = int(header.lstrip('#').split('locations')[0].strip())
nt =  int(header.split('locations')[1].strip().rstrip('times').strip())

if outputtype < 10 or outputtype > 11:
    print 'ERROR: only output_flag = {10,11} works with this script'

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
    if outputtype == 10:
        vv.append(r[:,2:4])
        dhv.append(r[:,4])
    else:
        dhv.append(r[:,2])

fig = plt.figure(figsize=(12,6))
first = True

colors = ['red','green','blue','pink','magenta','cyan','orange','purple']
ncol = len(colors)
lines = ['-','--','-.',':']

##plt.suptitle('%s (x=%.5g, y=%.5g)' % (md['name'],md['xloc'],md['yloc']))

#loc = 2 # which of the locations to plot?
for loc in range(nloc):
    
    if outputtype == 10:
        nplotcol = 3
    else:
        nplotcol = 2
    
    # two/three side-by-side timeseries plots
    if first:
        ax1 = fig.add_subplot(1,nplotcol,1)
    ax1.plot(t,(hv[loc]),linestyle=lines[int(loc/ncol)],color=colors[loc%ncol],
             label='%s (x=%.5g, y=%.5g)' % (mdv[loc]['name'],mdv[loc]['xloc'],mdv[loc]['yloc']))
    ax1.set_xscale('log')
    ax1.set_yscale('linear')
    #ax1.set_ylim([0.001,3])
    ax1.set_xlabel('time')
    ax1.set_ylabel('head')
    #ax1.axis('equal')
    ax1.grid()
    
    if first:
        ax2 = fig.add_subplot(1,nplotcol,2)
    ax2.plot(t,dhv[loc])
    ax2.set_xscale('log')
    ax2.set_yscale('linear')
    ax2.set_xlabel('time')
    ax2.set_ylabel('$\partial$ head/$\\ln t$')
    ax2.set_ylim([-20,20])
    ax2.grid()
    #ax2.axis('equal')
    
    if outputtype == 10:
        if first:
            ax3 = fig.add_subplot(133)
        ax3.plot(t,vv[loc][:,0],'r:',label='$v_x$')
        ax3.plot(t,vv[loc][:,1],'g--',label='$v_y$')
        ax3.plot(t,np.abs(vv[loc][:,0]+v[loc][:,1]*1j),'k-',label='$|v|$')
        ax3.set_xscale('linear')
        ax3.set_yscale('linear')
        ax3.set_xlabel('time')
        ax3.set_ylabel('velocity')
        ax3.grid()


#ax1.legend(loc=0)
ax1.set_title('head')
#ax2.legend(loc=0)
ax2.set_title('derivative')
if outputtype == 10:
#    ax3.legend(loc=0)
    ax3.set_title('velocity')

plt.tight_layout()
plt.savefig('test-timeseries.png')

