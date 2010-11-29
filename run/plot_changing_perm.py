import numpy as np
import matplotlib.pyplot as plt

# output is saved to csv with the following columns
# idx, K, q, MS, cutoff, pass/fail

fdata = open('changing-perm-04.csv','r')
lines = fdata.readlines()
fdata.close()

data = []

for line in lines:
    v = line.split(',')
    if len(v) > 1:
        if v[5].strip() == 'passed':
            data.append([int(v[0]), float(v[1]), float(v[2]), int(v[3]), float(v[4])])

# set of unique cutoff values
cutoffs = list(set(zip(*data)[4]))
cutoffs.sort()

p = np.array(data)
colors = ['r','g','b','c','m','y','k','r','g','b']

plt.figure(1)
for j,c in enumerate(cutoffs):
    mask = p[:,4] == c
    plt.semilogx(p[mask,2], p[mask,3], colors[j]+'.-', label='%.0e' % c)

plt.legend(loc=2)
plt.xlabel('max value of |q| required')
plt.ylabel('required matrix size for N=10')
plt.savefig('matrix-size-vs-q.eps')
plt.close(1)

