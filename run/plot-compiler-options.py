import matplotlib.pyplot as plt
import numpy as np

fhin = open('run-data.csv','r')

data = []

for line in fhin:
    data.append(line.split(','))

fhin.close()

n = zip(*data)[0]
t = zip(*data)[1]
m = zip(*data)[2]

name = [val.strip(' test.in').strip('./ltaem-').replace('native','n').replace('wholefile','w').replace('both','b') for val in n]
mem = [float(val) for val in m]
secs = [float(val.split(':')[0])*60 + float(val.split(':')[1]) for val in t]

del n,t,m

n = len(name)

plt.figure(1)
plt.subplot(211)
plt.bar(np.arange(n),100*(np.array(mem) - mem[0])/mem[0])
plt.ylabel('relative mem (%)')
plt.xticks(np.arange(n),[''*n])
plt.subplot(212)
plt.bar(np.arange(n),100*(np.array(secs) - secs[0])/secs[0])
plt.ylabel('relative time (%)')
plt.xticks(np.arange(n),name,rotation=35)
plt.savefig('compiler-options.png')

