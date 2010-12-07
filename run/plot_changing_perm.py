import numpy as np
import matplotlib.pyplot as plt

# output is saved to csv with the following columns
# idx, K, q, MS, cutoff, pass/fail, N

fdata = open('changing-perm-06.csv','r')
lines = fdata.readlines()
fdata.close()

data = []

for line in lines:
    v = line.split(',')
    if len(v) > 1:
        if v[5].strip() == 'passed':
            data.append([int(v[0]), float(v[1]), float(v[2]), int(v[3]), float(v[4]), int(v[6])])

# set of unique cutoff values
cutoffs = list(set(zip(*data)[4]))
cutoffs.sort()

# set of unique max orders
orders = list(set(zip(*data)[5]))
orders.sort()
nn = len(orders)

p = np.array(data)
colors = ['r','g','b','c','m','y','k','r','g','b']

lt = ['-',':','--']

qq = np.logspace(np.log10(p[:,2].min()),np.log10(p[:,2].max()),100)

plt.figure(1)
for ii,n in enumerate(orders):
    C = (8.46 + 0.44*n)/(1.0 + 0.085*n)
    D = (0.240 + 0.0214*n)/(1.0 + 0.059*n)

    plt.loglog(qq, np.ceil(n + 3.0 + C*qq**D), 'k'+lt[ii], label='Shirts ($10^{-12}$) N=%i' % (n,))

    for jj,c in enumerate(cutoffs):
        if c == 1.0E-3 or c == 1.0E-9:
            mask = np.logical_and(p[:,4] == c, p[:,5] == n)
            plt.loglog(p[mask,2], p[mask,3], colors[jj]+lt[ii], label='$10^{%i}$, N=%i' % (np.log10(c),n))

plt.xlim([1.0E-1,1.3E+4])
plt.ylim([10,300])
plt.legend(loc=2)
plt.xlabel('max value of |q| needed')
plt.ylabel('required matrix size')  
plt.savefig('matrix-size-vs-q.eps')
plt.close(1)

