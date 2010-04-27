# $Id: Makefile,v 1.1 2007/03/15 00:16:09 kris Exp kris $
EXTERNAL=cbessel.o
KKLOW=special_functions.o bessel_functions.o invlap.o file_ops.o calcshareddata.o shared_matching_data.o

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
KKHIGH=wells.o circ_passive_elements.o matching.o calc_routines.o particle_integration.o circular_geometry.o 
MAIN=ltaem_main.o
TARG=ltaem

OBJS=constants.o element_specs.o $(EXTERNAL) $(KKLOW) $(KKHIGH)

# give warnings for underflow from Bessel functions
DEBUG=-Dbessuf 

# for the dependency-generation
F90OBJS=constants.o element_specs.o $(KKLOW) $(KKHIGH)
F90SRC=./external/cbessel.f90 $(patsubst %.o,%.f90,$(F90OBJS) $(MAIN)) 

########### g95 ##########
#F90=g95
#DFLAGS=-O0 -Wall -Wextra -fbounds-check -freal=nan -finteger=-999999 -flogical=false -ftrace=full
#WUNFLAG=-Wunused-module-vars -Wunused-module-procs -Wunused-parameter
#W1FLAG=-Wno=140,141
#W2FLAG=-Wno=159
#W3FLAG=-Wno=136
#W5FLAG=-Wno=165
#PFLAGS=-O3 -march=nocona -funroll-loops -fomit-frame-pointer
################################

########## ifort (work) related ##########
#F90=ifort
#DFLAGS=-O0 -g -traceback -cpp -warn all -check bounds -check pointers -check uninit
#PFLAGS=-O3 -cpp
#LINKFLAG=-O3
###################################

#### gfortran 4.2.X (system version) ################# 
F90=gfortran
DFLAGS=-W -O0 -fbounds-check -frange-check -Wall
PFLAGS=-O3 -march=nocona
###################################

#### gfortran 4.4 (trunk) #################
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

# warning 136: unused module variable
file_ops.o :  file_ops.f90
	$(F90) -c $(FLAGS) $(W3FLAG) -o $@ $<

## include output generated using "make dep" below
cbessel.o complex_bessel.mod: ./external/cbessel.f90
constants.o constants.mod: constants.f90
element_specs.o element_specs.mod: constants.mod element_specs.f90
special_functions.o matrix_inverse.mod bessel_functions.mod error_handler.mod: constants.mod complex_bessel.mod special_functions.f90
invlap.o inverse_laplace_transform.mod: constants.mod invlap.f90
file_ops.o file_ops.mod: constants.mod error_handler.mod element_specs.mod file_ops.f90
calcshareddata.o calc_shared_data.mod: constants.mod calcshareddata.f90
shared_matching_data.o shared_matching_data.mod: constants.mod shared_matching_data.f90
wells.o wells.mod: bessel_functions.mod element_specs.mod constants.mod wells.f90
circ_passive_elements.o circ_passive_elements.mod: element_specs.mod wells.mod bessel_functions.mod constants.mod shared_matching_data.mod circ_passive_elements.f90
matching.o matching_old.mod: circ_passive_elements.mod wells.mod matrix_inverse.mod shared_matching_data.mod bessel_functions.mod element_specs.mod constants.mod matching.f90
#match_matrix.o matching_new.mod: wells.mod matrix_inverse.mod shared_matching_data.mod bessel_functions.mod element_specs.mod constants.mod match_matrix.f90
calc_routines.o calc_routines.mod: wells.mod bessel_functions.mod element_specs.mod constants.mod calc_routines.f90
particle_integration.o particle_integrate.mod: error_handler.mod calc_routines.mod calc_shared_data.mod element_specs.mod inverse_laplace_transform.mod constants.mod particle_integration.f90
circular_geometry.o circular_geometry.mod: file_ops.mod shared_matching_data.mod element_specs.mod constants.mod circular_geometry.f90
ltaem_main.o: particle_integrate.mod circular_geometry.mod calc_routines.mod matching_old.mod inverse_laplace_transform.mod error_handler.mod file_ops.mod element_specs.mod calc_shared_data.mod constants.mod ltaem_main.f90
bessel_functions.o bessel_functions_deriv.mod: constants.mod complex_bessel.mod bessel_functions.f90

%.o: %.f
	$(F90) -c $(PFLAGS) -o $@ $<

%.o: %.f90
	$(F90) -c $(FLAGS) -o $@ $<



clean:
	rm -f *.o *.mod $(TARG)

# create list of dependencies using compiler
dep:
	g95 -M -cpp $(F90SRC) > dep.out && echo -e "\n\t***\ndependencies output to dep.out\n\t***\n"

