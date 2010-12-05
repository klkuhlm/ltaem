EXTERNAL=cbessel.o
KKLOW=utility.o bessel_functions.o mathieu_functions2.o

##Work MacPro##############################
#LAPACKLIB=/Users/klkuhlm/src/lapack-3.2.1/lapack_MAC.a
#BLASLIB=/Users/klkuhlm/src/lapack-3.2.1/blas_MAC.a

##home linux desktop ##########
LAPACKLIB=/home/kris/src/lapack-3.3.0/lapack_LINUX.a
BLASLIB=/home/kris/src/lapack-3.3.0/blas_LINUX.a

##work linux desktop ##########
#LAPACKLIB=/home/klkuhlm/src/lapack-3.3.0/lapack_LINUX.a
#BLASLIB=/home/klkuhlm/src/lapack-3.3.0/blas_LINUX.a
#LAPACKLIB=/opt/intel/mkl/10.0.3.020/lib/em64t/libmkl_lapack.a

### for intel compiler
##################################
##BLASLIB=/home/klkuhlm/src/GotoBLAS2/libgoto2.a -pthread
##LAPACKLIB=/home/klkuhlm/src/lapack-3.2.1/lapack_LINUX.a

# main-program related
####################################
KKHIGH=kappa.o time.o ltaem_io.o invlap.o ellipse_mathieu_init.o circular_elements.o elliptical_elements.o calc_routines.o element_geometry.o solution.o particle_integration.o

MAIN=ltaem_main.o

TARG=ltaem
OBJS=$(EXTERNAL) constants.o $(KKLOW) type_definitions.o $(KKHIGH)

############## ifort 11.1 ##########
##F90=ifort
##OMP=-openmp -DOMP
##CPP=-cpp
##DFLAGS=-O0 -g -traceback -warn all -check bounds -check pointers -check uninit -check format -check output_conversion
###DFLAGS=-O0 -g -traceback -warn all -check all
##DFLAGS+=$(CPP) #-DDEBUG
##PFLAGS=$(CPP)  -xHOST -O3 -ipo # $(OMP)
##LINKFLAG=$(DFLAGS) $(OMP)
#######################################

#### gfortran #################
F90=gfortran-4.5
OMP=-fopenmp -DOMP
CPP=-cpp
DFLAGS=-O0 -Wall -Wextra -fbacktrace -g
DFLAGS+=-frange-check -fcheck=all -finit-real=snan -finit-integer=-99999
DFLAGS+=$(CPP) #-DDEBUG 
##DFLAGS+=-fdump-rtl-expand
PFLAGS=-O3 -march=native $(CPP) -fwhole-file  $(OMP)
LINKFLAG=$(PFLAGS) $(OMP)
###################################

# switch between performance or debugging flags
FLAGS=$(DFLAGS)

# the link step
main: $(OBJS) $(MAIN)
	$(F90) $(LINKFLAG) -o $(TARG) $(OBJS) $(MAIN) $(LAPACKLIB) $(BLASLIB)
	install $(TARG) ../run/

# this external library tends to produce spurious errors about implicit type conversions
cbessel.o : ./external/cbessel.f90
	$(F90) -c $(FLAGS) -o $@ $<

%.o: %.f90
	$(F90) -c $(FLAGS) -o $@ $<

%.o: %.f
	$(F90) -c $(FLAGS) -o $@ $<
clean:
	rm -f *.o *.mod $(TARG)

# dependency listing 
# "provides" : "depends on"
cbessel.o complex_bessel.mod : ./external/cbessel.f90
constants.o constants.mod : constants.f90
type_definitions.o type_definitions.mod : constants.mod mathieu_functions.mod type_definitions.f90
utility.o utility.mod : constants.mod
bessel_functions.o bessel_functions.mod : complex_bessel.mod constants.mod bessel_functions.f90
mathieu_functions2.o mathieu_functions.mod : constants.mod mathieu_functions2.f90
kappa.o kappa_mod.mod : constants.mod type_definitions.mod kappa.f90
time.o  time_mod.mod  : constants.mod type_definitions.mod utility.mod time.f90
ltaem_io.o file_ops.mod : type_definitions.mod constants.mod ltaem_io.f90
invlap.o inverse_laplace_transform.mod : constants.mod type_definitions.mod invlap.f90
solution.o solution_mod.mod : constants.mod type_definitions.mod circular_elements.mod elliptical_elements.mod solution.f90
ellipse_mathieu_init.o ellipse_mathieu_init.mod : constants.mod mathieu_functions.mod type_definitions.mod kappa_mod.mod ellipse_mathieu_init.f90
circular_elements.o circular_elements.mod : constants.mod kappa_mod.mod time_mod.mod utility.mod type_definitions.mod bessel_functions.mod circular_elements.f90
elliptical_elements.o elliptical_elements.mod : constants.mod kappa_mod.mod time_mod.mod utility.mod type_definitions.mod mathieu_functions.mod elliptical_elements.f90
calc_routines.o calc_routines.mod : type_definitions.mod constants.mod circular_elements.mod elliptical_elements.mod kappa_mod.mod time_mod.mod utility.mod calc_routines.f90
element_geometry.o geometry.mod : constants.mod type_definitions.mod file_ops.mod utility.mod element_geometry.f90
ltaem_main.o : constants.mod type_definitions.mod file_ops.mod inverse_laplace_transform.mod solution_mod.mod calc_routines.mod geometry.mod ellipse_mathieu_init.mod ltaem_main.f90
particle_integration.o particle_integrate.mod : constants.mod inverse_laplace_transform.mod type_definitions.mod calc_routines.mod

