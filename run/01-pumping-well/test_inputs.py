from subprocess import Popen,PIPE,check_output

names = ['porosity','K','Ss','leaky flag','K2','Ss2','b2','ellipse MS','ellipse cutoff','Sy','Kz','b','matrixSs','lambda','x','y','t','DH alpha','DH tol','DH M']
values = [1.0,1.0,0.1,0,1.0,1.0,1.0,20,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.00001,1.0E-6,1.0E-8,10]
nval = len(values)
subs = [-2,0,'junk']

check_output(["cp","pumping_well_circles.in","circles.in"])

for j in range(nval):
    for s in subs: 
        print names[j],'substitute',s
        v = values[:]
        v[j] = s
        P = Popen(["./ltaem","input.in"],stdout=PIPE,stderr=PIPE)
        fh = open('input.in','w')
        fh.write("""True  False  True  True  1  pumping_well_contour  dump.out  pumping_well.elem  pumping_well.geom    
%s  %s  %s  %s  %s  %s  %s  %s  %s  ::  por, k, Ss, leakFlag, K2, Ss2, b2, ellipse MS, ellipse cutoff
%s  %s  False  %s  :: Sy, Kz, unconfined?, b
%s  %s  False  :: matrixSs, matrix/fracture lambda, dualPorosityFlag?
4  4  4  :: nx, ny, nt
NAME 1|NAME 2|NAME 3  
%s -2.37179e+00 -2.24359e+00 -2.11538e+00  :: x 
%s -2.37179e+00 -2.24359e+00 -2.11538e+00  :: y
%s 0.001 0.01 0.1 :: t
%s  %s  %s  :: alpha, tolerance, M
1  circles.in     :: number of circular elements, circle data file
0  not_used    :: number of elliptical elements, ellipse data file 
not_used   ::  particle data file
""" % tuple(v))
        fh.close()

        stdout,stderr = P.communicate()
        print 'error >>',stderr.rstrip()
        print 'output >>'
        print stdout

P = Popen(["./ltaem","input.in"],stdout=PIPE,stderr=PIPE)
fh = open('input.in','w')
fh.write("""True  False  True  True  1  pumping_well_contour  dump.out  pumping_well.elem  pumping_well.geom    
%s  %s  %s  %s  %s  %s  %s  %s  %s  ::  por, k, Ss, leakFlag, K2, Ss2, b2, ellipse MS, ellipse cutoff
%s  %s  False  %s  :: Sy, Kz, unconfined?, b
%s  %s  False  :: matrixSs, matrix/fracture lambda, dualPorosityFlag?
4  4  4  :: nx, ny, nt
NAME 1|NAME 2|NAME 3  
%s -2.37179e+00 -2.24359e+00 -2.11538e+00  :: x 
%s -2.37179e+00 -2.24359e+00 -2.11538e+00  :: y
%s 0.001 0.01 0.1 :: t
%s  %s  %s  :: alpha, tolerance, M
1  circles.in     :: number of circular elements, circle data file
0  not_used    :: number of elliptical elements, ellipse data file 
not_used   ::  particle data file
""" % tuple(values))
fh.close()
stdout,stderr = P.communicate()

names = ['N','M','ibnd','r','x','y','K','Ss','porosity','leak flag','K2','Ss2','b2','Sy','Kz','b','mSs','lambda','skin','area Q','bdry Q','area flag','area P1','area P2','bdry flag','bdry P1','bdry P2']
values = [1,1,2, 0.1,2.0,3.0,1.0,1.0,1.0,1,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,2,0.0,0.0,2,0.0,0.0]
nval = len(values)
subs = [-2,0,'junk']
for j in range(nval):
    for s in subs: 
        print 'circle',names[j],'substitute',s
        v = values[:]
        v[j] = s
        P = Popen(["./ltaem","input.in"],stdout=PIPE,stderr=PIPE)

        fh = open('circles.in','w')
        fh.write("""%s  :: N
%s  :: M
%s  :: ibnd
F  :: calcin?
F  :: storin?
%s   :: r
%s    :: x0
%s    :: y0
%s    :: K
%s    :: Ss
%s    :: porosity
%s    ::  leaky type flag
%s  :: K2
%s  :: Ss2
%s  :: b2
F        :: unconfined?
%s   :: Sy
%s   :: Kz
%s   :: b
F      :: dual porosity flag
%s   :: matrix Ss
%s   :: fracture/matrix lambda
%s   :: dimensionless skin
%s    :: area source strength
%s    :: boundary source strength
%s  %s  %s  :: area time flag, tpar1, tpar2
%s  %s  %s  :: boundary time flag, tpar1, tpar2
"""% tuple(v))
        fh.close()

        stdout,stderr = P.communicate()
        print 'error >>',stderr.rstrip()
        print 'output >>'
        print stdout