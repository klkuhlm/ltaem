# $Id: Makefile,v 1.1 2007/03/15 00:16:09 kris Exp kris $
EXTERNAL=cbessel.o
KKLOW=utility.o bessel_functions.o mathieu_functions2.o

##Work MacPro##############################
LAPACKLIB=/Users/klkuhlm/src/lapack-3.2.1/lapack_MAC.a
BLASLIB=/Users/klkuhlm/src/lapack-3.2.1/blas_MAC.a

##Work linux laptop / home linux desktop ##########
#LAPACKLIB=/home/kris/src/lapack-3.2.1/lapack_LINUX.a
#BLASLIB=/home/kris/src/lapack-3.2.1/blas_LINUX.a

### for intel compiler
##################################
##BLASLIB=/home/klkuhlm/src/GotoBLAS2/libgoto2.a -pthread
##LAPACKLIB=/home/klkuhlm/src/lapack-3.2.1/lapack_LINUX.a

# main-program related
####################################
KKHIGH=kappa.o time.o ltaem_io.o invlap.o ellipse_mathieu_init.o circular_elements.o elliptical_elements.o calc_routines.o element_geometry.o solution.o particle_integration.o

KKHIGHUNSAT=unsat_io.o unsat_circular_elements.o unsat_elliptical_elements.o unsat_calc_routines.o unsat_element_geometry.o unsat_solution.o

MAIN=ltaem_main.o
MAINUNSAT=unsat_main.o

TARG=ltaem
TARGUNSAT=unsat
OBJS=$(EXTERNAL) constants.o $(KKLOW) type_definitions.o $(KKHIGH)
OBJSUNSAT=$(EXTERNAL) constants.o $(KKLOW) unsat_type_definitions.o $(KKHIGHUNSAT)

############ ifort 11.1 ##########
##F90=ifort
###DFLAGS=-O0 -g -traceback -cpp -DDEBUG -warn all -check bounds -check pointers -check uninit -check format -check output_conversion
##DFLAGS=-O0 -g -traceback -cpp -DDEBUG -warn all -check all
##PFLAGS=-O0 -cpp
##LINKFLAG=-O0
#####################################

#### gfortran #################
F90=gfortran-4.6
DFLAGS=-O0 -Wall -Wextra -fbacktrace -g
DFLAGS+=-frange-check -fcheck=all -finit-real=snan -finit-integer=-99999 -fbacktrace #-fwhole-file
DFLAGS+=-cpp -DDEBUG 
#DFLAGS+=-fdump-rtl-expand
PFLAGS=-O3 -march=native -cpp -fwhole-file  -fopenmp
LINKFLAG=-$(PFLAGS)
###################################

# switch between performance or debugging flags
FLAGS=$(PFLAGS)

# the link step
main: $(OBJS) $(MAIN)
	$(F90) $(LINKFLAG) -o $(TARG) $(OBJS) $(MAIN) $(LAPACKLIB) $(BLASLIB)
	install $(TARG) ../run/

unsat: $(OBJSUNSAT) $(MAINUNSAT)
	$(F90) $(LINKFLAG) -o $(TARGUNSAT) $(OBJSUNSAT) $(MAINUNSAT) $(LAPACKLIB) $(BLASLIB)
	install $(TARGUNSAT) ../run/

# this external library tends to produce spurious errors about implicit type conversions
cbessel.o : ./external/cbessel.f90
	$(F90) -c $(FLAGS) -o $@ $<

%.o: %.f
	$(F90) -c $(FLAGS) -o $@ $<

%.o: %.f90
	$(F90) -c $(FLAGS) -o $@ $<

clean:
	rm -f *.o *.mod $(TARG) $(TARGUNSAT)

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

unsat_type_definitions.o unsat_type_definitions.mod : constants.mod mathieu_functions.mod unsat_type_definitions.f90
unsat_io.o unsat_file_ops.mod : unsat_type_definitions.mod constants.mod unsat_io.f90
unsat_solution.o unsat_solution_mod.mod : constants.mod unsat_type_definitions.mod unsat_circular_elements.mod unsat_elliptical_elements.mod unsat_solution.f90
unsat_circular_elements.o unsat_circular_elements.mod : constants.mod utility.mod unsat_type_definitions.mod bessel_functions.mod unsat_circular_elements.f90
unsat_elliptical_elements.o unsat_elliptical_elements.mod : constants.mod utility.mod unsat_type_definitions.mod mathieu_functions.mod unsat_elliptical_elements.f90
unsat_calc_routines.o unsat_calc_routines.mod : unsat_type_definitions.mod constants.mod unsat_circular_elements.mod unsat_elliptical_elements.mod utility.mod unsat_calc_routines.f90
unsat_element_geometry.o unsat_geometry.mod : constants.mod unsat_type_definitions.mod unsat_file_ops.mod utility.mod unsat_element_geometry.f90
unsat_main.o : constants.mod unsat_type_definitions.mod unsat_file_ops.mod unsat_solution_mod.mod unsat_calc_routines.mod unsat_geometry.mod unsat_main.f90



