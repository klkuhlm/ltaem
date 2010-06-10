# $Id: Makefile,v 1.1 2007/03/15 00:16:09 kris Exp kris $
EXTERNAL=cbessel.o
KKLOW=utility.o bessel_functions.o mathieu_functions2.o

# 
################################
LAPACKLIB=/Users/klkuhlm/src/lapack-3.2.1/lapack_MAC10.5.a
BLASLIB=/Users/klkuhlm/src/GotoBLAS2/libgoto2.a -lpthread

# for trunk gfortran version
################################
#BLASLIB=/home/kris/src/lapack-3.1.1/BLAS/SRC/blas_GFORTRAN4.4.a
#LAPACKLIB=/home/kris/src/lapack-3.1.1/lapack_GFORTRAN4.4.a


# main-program related
####################################
## matching.o not used
KKHIGH=kappa.o time.o ltaem_io.o invlap.o solultion.o ellipse_mathieu_init.o circular_elements.o elliptical_elements.o calc_routines.o particle_integration.o circular_geometry.o 
MAIN=ltaem_main.o
TARG=ltaem

OBJS=constants.o type_definitions.o $(EXTERNAL) $(KKLOW) $(KKHIGH)

########## ifort (work) related ##########
#F90=ifort
#DFLAGS=-O0 -g -traceback -cpp -warn all -check bounds -check pointers -check uninit
#PFLAGS=-O3 -cpp
#LINKFLAG=-O3
###################################

#### gfortran 4.4 (system version) ################# 
F90=gfortran
DFLAGS=-O0 -fbounds-check -frange-check -Wall -Wextra 
PFLAGS=-O3 -march=native
###################################

#### gfortran 4.6 (trunk) #################
#F90=gfortran-trunk
#DFLAGS=-W -g -O0 -fbounds-check -frange-check -static
#PFLAGS=-O3 -march=nocona -static
#LINKFLAG=-static
###################################

# switch between performance or debugging flags
FLAGS=$(DFLAGS)

# the link step
main: $(OBJS) $(MAIN)
	$(F90) $(LINKFLAG) -o $(TARG) $(OBJS) $(MAIN) $(LIBPATH) $(LAPACKLIB) $(BLASLIB)
	install $(TARG) ../run/

# warnings 140 and 141: implicit type conversion
cbessel.o : ./external/cbessel.f90
	$(F90) -c $(PFLAGS) $(W1FLAG) -o $@ $<

%.o: %.f
	$(F90) -c $(PFLAGS) -o $@ $<

%.o: %.f90
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
mathieu_functions2.o mathieu_functions.mod : constants.mod mathieu_functions.f90
kappa.o kappa_mod.mod : constants.mod type_definitions.mod kappa.f90
time.o  time_mod.mod  : constants.mod type_definitions.mod time.f90
ltaem_io.o file_ops.mod : type_definitions.mod constants.mod ltaem_io.f90
invlap.o inverse_laplace_transform.mod : constants.mod type_definitions.mod invlap.f90
solultion.o solution.mod : constants.mod type_definitions.mod circular_elements.mod elliptical_elements.mod solution.f90
ellipse_mathieu_init.o ellipse_mathieu_init.mod : constants.mod mathieu_functions.mod type_definitions.mod kappa_mod.mod ellipse_mathieu_init.f90
circular_elements.o circular_elements.mod : constants.mod kappa_mod.mod time_mod.mod utility.mod type_definitions.mod bessel_functions.mod circular_elements.f90
elliptical_elements.o elliptical_elements.mod : constants.mod kappa_mod.mod time_mod.mod utility.mod type_definitions.mod mathieu_functions.mod elliptical_elements.f90
calc_routines.o 
element_geometry.o geometry.mod : constants.mod type_definitions.mod file_ops.mod utility.mod element_geometry.f90
ltaem_main.o : constants.mod type_definitions.mod file_ops.mod inverse_laplace_transform.mod solution.mod calc_routines.mod element_geometry.mod ellipse_mathieu_init.mod ltaem_main.f90



