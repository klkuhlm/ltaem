#
# Copyright (c) 2011 Kristopher L. Kuhlman (klkuhlm at sandia dot gov)
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

# provide a link to the static lapack/blas libraries in the source directory
LAPACKLIB=lapack.a
BLASLIB=blas.a

# main-program related
####################################
KKHIGH=kappa.o time.o ltaem_io.o invlap.o ellipse_mathieu_init.o \
	circular_elements.o elliptical_elements.o calc_routines.o \
	element_geometry.o solution.o particle_integration.o

MAIN=ltaem_main.o

TARG=ltaem
OBJS=$(EXTERNAL) constants.o $(KKLOW) geometryUtil.o $(KKHIGH)

################## ifort >= 12.0 ##########
##F90=ifort
##OMP=-openmp
##CPP=-cpp
###DFLAGS=-O0 -g -traceback -warn all -check bounds -check pointers -check uninit -check format -check output_conversion
##DFLAGS=-O0 -g -traceback -warn all -check all
##PFLAGS=-xHOST -O3 $(OMP)
##LINKFLAG=$(DFLAGS) $(OMP)
###########################################

## gfortran #################
F90=gfortran-4.6 # must be >= 4.6.1
#F90=/opt/local/bin/dragonegg-3.0-gfortran 
OMP=-fopenmp
CPP=-cpp
DFLAGS=-O0 -Wall -Wextra -g -std=f2008
DFLAGS+=-frange-check -fbounds-check -fcheck=all
##DFLAGS+=-fdump-rtl-expand
PFLAGS=-Ofast -march=native -flto
#################################

BFLAGS=$(PFLAGS)
PFLAGS+=$(PFLAGEXTRA) $(OMP)

# switch between performance or debugging flags
FLAGS=$(PFLAGS)

LINKFLAG=$(FLAGS) -fbacktrace $(OMP)

# the link step
main: $(OBJS) $(MAIN)
	$(F90) $(LINKFLAG) -o $(TARG) $(OBJS) $(MAIN) $(LAPACKLIB) $(BLASLIB)
	install $(TARG) ../run/

# this external library tends to produce spurious errors about implicit type conversions
cbessel.o : ./external/cbessel.f90
	$(F90) -c $(BFLAGS) -o $@ $<

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
utility.o utility.mod : constants.mod utility.f90
geometryUtil.o geomConv.mod : constants.mod type_definitions.mod utility.mod geometryUtil.f90
bessel_functions.o bessel_functions.mod : complex_bessel.mod constants.mod bessel_functions.f90
mathieu_functions2.o mathieu_functions.mod : constants.mod mathieu_functions2.f90
kappa.o kappa_mod.mod : constants.mod type_definitions.mod kappa.f90
time.o  time_mod.mod  : constants.mod type_definitions.mod utility.mod time.f90
ltaem_io.o file_ops.mod : type_definitions.mod constants.mod ltaem_io.f90
invlap.o inverse_laplace_transform.mod : constants.mod type_definitions.mod invlap.f90
solution.o solution_mod.mod : constants.mod type_definitions.mod circular_elements.mod \
	elliptical_elements.mod solution.f90
ellipse_mathieu_init.o ellipse_mathieu_init.mod : constants.mod mathieu_functions.mod \
	type_definitions.mod kappa_mod.mod ellipse_mathieu_init.f90
circular_elements.o circular_elements.mod : constants.mod kappa_mod.mod time_mod.mod \
	utility.mod type_definitions.mod bessel_functions.mod circular_elements.f90
elliptical_elements.o elliptical_elements.mod : constants.mod kappa_mod.mod time_mod.mod \
	utility.mod type_definitions.mod mathieu_functions.mod elliptical_elements.f90
calc_routines.o calc_routines.mod : type_definitions.mod constants.mod circular_elements.mod \
	elliptical_elements.mod kappa_mod.mod time_mod.mod utility.mod calc_routines.f90
element_geometry.o geometry.mod : constants.mod type_definitions.mod file_ops.mod \
	utility.mod element_geometry.f90
ltaem_main.o : constants.mod type_definitions.mod file_ops.mod inverse_laplace_transform.mod \
	solution_mod.mod calc_routines.mod geometry.mod ellipse_mathieu_init.mod ltaem_main.f90
particle_integration.o particle_integrate.mod : constants.mod inverse_laplace_transform.mod \
	type_definitions.mod calc_routines.mod

