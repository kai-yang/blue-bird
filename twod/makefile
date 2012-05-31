#FF=ifort
#F90=ifort -traceback -CB

FF=gfortran
F90=$(FF)
CC=g++

#F77FLAGS=-O3 -xT
#F90FLAGS=-O3 -xT

F77FLAGS=-O3
F90FLAGS=-O3
CCFLAGS=-O3 -I../../ext/boost/Linux/include  -I../../ext/python/Linux/include

LIB_PATH= -Wl,-rpath,-L/usr/local/lib
#LIB_FFT= -ldrfftw_mpi -ldfftw_mpi -ldrfftw -ldfftw
LIB_LAPACK=-llapack -lblas
LIB_MPI =  
#LIB_PIM = -lpim_double
 
CLIBS = ../../ext/boost/Linux/lib/libboost_python.a ../../ext/python/Linux/lib/libpython.so
LIBS = $(LIB_PATH) $(LIB_FFT) $(LIB_LAPACK) $(LIB_PIM) $(LIB_MPI)

F77SR= dgetrf.f dgetri.f

F90SR= global_com.f90 global_fast.f90 global_geom.f90 global_dim.f90\
quadratures.f90 misc_dbl.f90 parallel.f90\
var_arrays.f90 intg_rtns.f90 layers.f90 fielde_dbl.f90 dir_field.f90\
mat_vec_mult.f90 geoinp.f90 tfqmr_arif.f90 interfaces.f90 

CCSR = python_bindings.cc 
CCMAIN = export_python.cc
MAIN = frequency.f90

MOM=$(F90SR:.f90=.o)
MOM1=$(F77SR:.f=.o)
CCO=$(CCSR:.cc=.o)
CCMAINO=$(CCMAIN:.cc=.o)
MAINO=$(MAIN:.f90=.o)


.SUFFIXES:      .o .f .c .C .h .f90

%.o: %.f90
	$(F90) $(F90FLAGS) -c $< -o $@
%.o: %.f
	$(FF) $(F77FLAGS) -c $< -o $@

%.o: %.cc
	$(CC) $(CCFLAGS) -c $< -o $@

all:  UT-MOM1 utfs

libutfs_basic.a: $(MOM) $(MOM1)
	ar r $@ $(MOM) $(MOM1) 

libutfs.a: $(MOM) $(MOM1) $(CCO)
	ar r $@ $(MOM) $(MOM1) $(CCO)

UT-MOM1: libutfs_basic.a $(MAINO) 
	$(F90) $(MAINO) libutfs_basic.a $(LIBS) -o $@

utfs: libutfs.a $(CCMAINO) 
	$(CC) $(CCFLAGS) $(CCMAINO) $(LIBS) -lgfortran $(CLIBS) libutfs.a -o $@

clean:
	rm -f *.o
	rm -f *.out
	rm -f *.time
	rm -f *.mod
	rm -f core
	rm -f UT-MOM1
