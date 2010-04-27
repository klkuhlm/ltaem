# -*-makefile-*-
# $Id: Makefile.arc,v 1.1 2007/07/18 02:01:54 kris Exp kris $
EXTERNAL=cbessel.o invlap.o
KK=error_handler.o file_ops.o arc_test.o

LIBPATH=-L/usr/lib64
LAPACKLIB=-llapack
BLASLIB=-lpthread

# main-program related
####################################
MAIN=arc_source_test.o
TARG=arc_source

OBJS=constants.o element_specs.o $(EXTERNAL) $(KK)

# for the dependency-generation
F90SRC=$(patsubst %.o,%.f90,$(OBJS) $(MAIN)) ./external/cbessel.f90

########### g95 (work) related ##########
F90=g95
DFLAGS=-Wall -Wextra -g -O0 -fbounds-check -ftrace=full -freal=NaN -finteger=-999999 -flogical=false
#DFLAGS+=-Wunused-module-vars -Wunused-module-procs -Wunused-parameter
W1FLAG=-Wno=140,141
W2FLAG=-Wno=159
W3FLAG=-Wno=136
W5FLAG=-Wno=165
PFLAGS=-O3 -march=nocona -fshort-circuit
##################################

######  gfortran 4.1.2 ##########
#F90=gfortran
#DFLAGS=-W -O0 -g -fbounds-check
#BLASLIB=~/src/GotoBLAS/libgoto.a -lpthread
#PFLAGS=-O3 -march=nocona
##################################

# switch between performance or debugging flags
FLAGS=$(DFLAGS)

# the link & install step for testing the line source
main: $(OBJS) $(MAIN)
	$(F90) $(LINKFLAG) -o $(TARG) $(OBJS) $(MAIN) \
  $(LIBPATH) $(LAPACKLIB) $(BLASLIB) $(BLASFLAG)
	install $(TARG) ../run/

# warnings 140 and 141: implicit type conversion
cbessel.o : ./external/cbessel.f90
	$(F90) -c $(PFLAGS) $(W1FLAG) -o $@ ./external/cbessel.f90

# warning 136: unused module variable
file_ops.o : file_ops.f90
	$(F90) -c $(FLAGS) $(W3FLAG) -o $@ $<


%.o: %.f90
	$(F90) -c $(FLAGS) -o $@ $<


clean:
	rm -f *.o *.mod $(TARG)

cleankk:
	rm -f $(MAIN) element_specs.o $(KKLOW) $(KKHIGH) $(TARG)

# create list of dependencies using compiler
dep:
	g95 -M $(F90SRC) > dep.out && echo -e "\n\t***\ndependencies output to dep.out\n\t***\n"

