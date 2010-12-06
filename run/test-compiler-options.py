import subprocess 
from glob import glob

PIPE = subprocess.PIPE


exes = glob('ltaem-*O[0,1,2,3]')
exes.sort()

fhout = open('run-data.csv','w')

for exe in exes:
    
    cmd = 'time -f "%C,%E,%M" ./'+exe+' test.in'
    print cmd

    p = subprocess.Popen(cmd,stdout=PIPE,stderr=PIPE,shell=True)
    (stdout,stderr) = p.communicate()
    
    # time output goes to stderr
    fhout.write(stderr.strip() + '\n')

fhout.close()
    
