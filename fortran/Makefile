#
# Copyright (c) 2011-2019 Kristopher L. Kuhlman (klkuhlm at sandia dot gov)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#

EXTERNAL=cbessel.o
KKLOW=type_definitions.o utility.o bessel_functions.o mathieu_functions2.o

ARBLOW=bessel_arb_wrapper.o cbessel_arb.o

####### use system libraries 
LAPACKLIB=-llapack
BLASLIB=-lblas

## if no system libraries put statically linked libraries in this directory
#FBLASLAPACKDIR=/home/klkuhlm/src/petsc/tgcc-c-debug/externalpackages/fblaslapack-3.4.2
#LAPACKLIB=${FBLASLAPACKDIR}/libflapack.a
#BLASLIB=${FBLASLAPACKDIR}/libfblas.a


####### use Intel math kernel library
##MKLROOT=/opt/intel/mkl/lib/intel64
##LAPACKLIB= -Wl,--start-group  \
##	$(MKLROOT)/libmkl_intel_lp64.a \
##	$(MKLROOT)/libmkl_sequential.a \
##	$(MKLROOT)/libmkl_core.a \
##	-Wl,--end-group -lpthread -lm

# main-program related
####################################
KKHIGH=kappa.o time.o ellipse_mathieu_init.o ltaem_io.o invlap.o  \
	circular_elements.o elliptical_elements.o calc_routines.o \
	element_geometry.o solution.o particle_integration.o

MAIN=ltaem_main.o

TARG=ltaem
ifdef ARB
  OBJS=constants.o $(ARBLOW) $(KKLOW) geometryUtil.o $(KKHIGH)
else
  OBJS=$(EXTERNAL) constants.o $(KKLOW) geometryUtil.o $(KKHIGH)
endif

SRCS=$(OBJS:.o=.f90) $(MAIN:.o=.f90)

############################### ifort >= 12.0 ##########
##F90=ifort
##OMP=-openmp
##CPP=-cpp
#######DFLAGS=-O0 -g -traceback -warn all -check bounds -check pointers -check uninit -check format -check output_conversion
######DFLAGS=-O0 -g -traceback -check all -warn all
##PFLAGS=-xHOST -O3 -fast  
########################################################

#### gfortran >= 4.6.1 #################
CC=gcc
F90=gfortran
#CC=/usr/local/stow/gcc-trunk/bin/gcc
#F90=/usr/local/stow/gcc-trunk/bin/gfortran
##OMP=-fopenmp
##CPP=-cpp
###DFLAGS=-O2 -march=native -fprofile-generate -flto
DFLAGS=-Og -Wall -Wextra  -g -fcheck=all #-ffpe-trap=overflow,underflow 
###DFLAGS+=-fdump-rtl-expand 
##LINKOPT=-fbacktrace -B/usr/lib/x86_64-linux-gnu  # needed to make dev version of gcc to work?
##PFLAGS=-O2 -march=native -pipe -funroll-loops -ftree-vectorize
PFLAGS=-O3 -march=native -flto=auto
#PFLAGS=-Og -g -fcheck=all
#################################

#PFLAGS+=$(OMP)

# switch between performance or debugging flags
# based on whethere the environment variable PERF is defined
ifdef PERF
  FLAGS=$(PFLAGS)
else
  FLAGS=$(DFLAGS)
endif

LINKFLAG=$(FLAGS) $(LINKOPT) -flto=auto

# the link step
ifdef ARB
  main: $(OBJS) $(MAIN)
	$(F90) $(LINKFLAG) -o $(TARG) $(OBJS) $(MAIN) $(LAPACKLIB) $(BLASLIB) -larb -lflint
	install $(TARG) ../run/
else
  main: $(OBJS) $(MAIN)
	$(F90) $(LINKFLAG) -o $(TARG) $(OBJS) $(MAIN) $(LAPACKLIB) $(BLASLIB)
	install $(TARG) ../run/
endif

# dependency listing
ifdef ARB
  cbessel_arb.o complex_bessel.mod: cbessel_arb.f90
	$(F90) -c $(FLAGS) -o $@ $<
  bessel_arb_wrapper.o: bessel_arb_wrapper.c
	$(CC) -c $(FLAGS) -o $@ $<
else
  cbessel.o complex_bessel.mod: cbessel.f90
	$(F90) -c $(FLAGS) -o $@ $<
endif

%.o: %.f90
	$(F90) -c $(FLAGS) -o $@ $<

constants.o constants.mod: constants.f90
type_definitions.o type_definitions.mod: type_definitions.f90 \
 constants.mod
utility.o utility.mod: utility.f90 constants.mod 
bessel_functions.o bessel_functions.mod: bessel_functions.f90 \
 constants.mod complex_bessel.mod complex_bessel.mod 
mathieu_functions2.o mathieu_functions.mod: mathieu_functions2.f90 \
 constants.mod type_definitions.mod utility.mod complex_bessel.mod  
geometryUtil.o geomconv.mod: geometryUtil.f90 constants.mod constants.mod \
 type_definitions.mod constants.mod type_definitions.mod utility.mod \
 constants.mod type_definitions.mod constants.mod constants.mod \
 type_definitions.mod constants.mod type_definitions.mod utility.mod \
 constants.mod type_definitions.mod
kappa.o kappa_mod.mod: kappa.f90 constants.mod type_definitions.mod \
 constants.mod type_definitions.mod
time.o time_mod.mod: time.f90 constants.mod type_definitions.mod \
 utility.mod constants.mod type_definitions.mod
ellipse_mathieu_init.o ellipse_mathieu_init.mod: ellipse_mathieu_init.f90 \
 constants.mod mathieu_functions.mod type_definitions.mod kappa_mod.mod \
 constants.mod
ltaem_io.o file_ops.mod: ltaem_io.f90 constants.mod type_definitions.mod \
 type_definitions.mod constants.mod ellipse_mathieu_init.mod
invlap.o inverse_laplace_transform.mod: invlap.f90 constants.mod \
 type_definitions.mod
circular_elements.o circular_elements.mod: circular_elements.f90 \
 constants.mod kappa_mod.mod time_mod.mod utility.mod \
 type_definitions.mod bessel_functions.mod constants.mod kappa_mod.mod \
 time_mod.mod utility.mod type_definitions.mod bessel_functions.mod \
 constants.mod kappa_mod.mod time_mod.mod type_definitions.mod \
 bessel_functions.mod constants.mod kappa_mod.mod type_definitions.mod \
 bessel_functions.mod constants.mod kappa_mod.mod time_mod.mod \
 type_definitions.mod bessel_functions.mod constants.mod kappa_mod.mod \
 time_mod.mod type_definitions.mod bessel_functions.mod
elliptical_elements.o elliptical_elements.mod: elliptical_elements.f90 \
 constants.mod kappa_mod.mod time_mod.mod type_definitions.mod \
 mathieu_functions.mod utility.mod constants.mod type_definitions.mod \
 time_mod.mod mathieu_functions.mod constants.mod type_definitions.mod \
 time_mod.mod mathieu_functions.mod constants.mod time_mod.mod \
 type_definitions.mod mathieu_functions.mod
calc_routines.o calc_routines.mod: calc_routines.f90 type_definitions.mod \
 constants.mod circular_elements.mod elliptical_elements.mod \
 kappa_mod.mod time_mod.mod type_definitions.mod constants.mod \
 circular_elements.mod elliptical_elements.mod kappa_mod.mod time_mod.mod \
 utility.mod type_definitions.mod constants.mod type_definitions.mod \
 constants.mod constants.mod type_definitions.mod geomconv.mod \
 constants.mod
element_geometry.o geometry.mod: element_geometry.f90 constants.mod \
 type_definitions.mod file_ops.mod geomconv.mod constants.mod \
 type_definitions.mod utility.mod geomconv.mod
solution.o solution_mod.mod: solution.f90 constants.mod \
 type_definitions.mod circular_elements.mod elliptical_elements.mod
particle_integration.o particle_integrate.mod: particle_integration.f90 \
 constants.mod inverse_laplace_transform.mod type_definitions.mod \
 calc_routines.mod utility.mod constants.mod \
 inverse_laplace_transform.mod type_definitions.mod calc_routines.mod \
 constants.mod inverse_laplace_transform.mod type_definitions.mod \
 calc_routines.mod constants.mod utility.mod \
 inverse_laplace_transform.mod type_definitions.mod calc_routines.mod \
 constants.mod constants.mod constants.mod type_definitions.mod \
 constants.mod utility.mod type_definitions.mod
ltaem_main.o: ltaem_main.f90 constants.mod type_definitions.mod \
 file_ops.mod inverse_laplace_transform.mod particle_integrate.mod \
 solution_mod.mod calc_routines.mod geometry.mod ellipse_mathieu_init.mod

## other targets
clean:
	rm -f *.o *.mod $(TARG)

# to make awful dependency listing above (doesn't work for ARB)
deps:
	gfortran -cpp -MM $(SRCS) > ltaem.deps

# build python module of stripped-down version of mathieu libraries using f2py
pymodule fmathieu.so : cbessel.f90 mathieu.f90
	f2py --fcompiler=gnu95 $(LAPACKLIB) $(BLASLIB) -c fmathieu.pyf ./external/cbessel.f90 mathieu.f90


