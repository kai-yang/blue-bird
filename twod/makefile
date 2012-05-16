#FF=ifort
#F90=ifort -traceback -CB

FF=gfortran
F90=$(FF)

#F77FLAGS=-O3 -xT
#F90FLAGS=-O3 -xT

F77FLAGS=-O3
F90FLAGS=-O3

LIB_PATH= -Wl,-rpath,-L$TACC_MKL_LIB -L$(HOME)/local/lib -L/usr/local/lib
#LIB_FFT= -ldrfftw_mpi -ldfftw_mpi -ldrfftw -ldfftw
LIB_LAPACK=-llapack -lblas
LIB_MPI =  
#LIB_PIM = -lpim_double
 
LIBS = $(LIB_PATH) $(LIB_FFT) $(LIB_LAPACK) $(LIB_PIM) $(LIB_MPI)

F77SR= dgetrf.f dgetri.f

F90SR= global_com.f90 global_fast.f90 global_geom.f90 global_dim.f90\
quadratures.f90 misc_dbl.f90 parallel.f90\
var_arrays.f90 intg_rtns.f90 layers.f90 fielde_dbl.f90 dir_field.f90\
mat_vec_mult.f90 geoinp.f90 tfqmr_arif.f90 interfaces.f90 frequency.f90

MOM=$(F90SR:.f90=.o)
MOM1=$(F77SR:.f=.o)

.SUFFIXES:      .o .f .c .C .h .f90

.f90.o: 
	$(F90) $(F90FLAGS) -c $<
.f.o:
	$(FF) $(F77FLAGS) -c $<

all:  UT-MOM1
UT-MOM1: $(MOM) $(MOM1)
	$(F90) $(MOM) $(MOM1) $(LIBS) -o $@

clean:
	rm -f *.o
	rm -f *.out
	rm -f *.time
	rm -f *.mod
	rm -f core
	rm -f UT-MOM1
