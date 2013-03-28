PLATFORM := $(shell uname)
BUILD := $(BUILD)

FF = /usr/local/bin/gfortran
F90 = $(FF)
AR = /usr/bin/ar

LIB_PATH = -Wl,-rpath,-L/usr/local/lib
LIB_FFT =
LIB_LAPACK = -llapack -lblas
LIB_PIM =
LIB_MPI =
LIBS = $(LIB_PATH) $(LIB_FFT) $(LIB_LAPACK) $(LIB_PIM) $(LIB_MPI)

#COMMON_FLAGS = -fbounds-check -cpp
COMMON_FLAGS = -cpp
DEBUG_FLAGS = -g -ffpe-trap=invalid,zero,overflow
#DEBUG_FLAGS = -g
RELEASE_FLAGS = -march=native -ffast-math -funroll-loops -O3 -fno-protect-parens -flto -fstack-arrays -finline-limit=600 -fwhole-program -fexternal-blas -mfpmath=sse -msse2
#MAIN_PROG_FLAGS = -DENABLED_FOR_MAIN_PROGRAM


F77SRCS = dgetrf.f \
	  dgetri.f

F77OBJS = $(F77SRCS:.f=.o) 

F90SRCS = global_com.f90 \
	  global_fast.f90 \
	  global_geom.f90 \
	  global_dim.f90 \
	  quadratures.f90 \
	  misc_dbl.f90 \
	  var_arrays.f90 \
	  intg_rtns.f90 \
	  fielde_dbl.f90 \
	  fielde_diel_dbl.f90 \
	  mat_vec_mult.f90 \
	  geoinp.f90

F90OBJS = $(F90SRCS:.f90=.o)

API_SRCS = interfaces.f90 \
	   tfqmr_arif.f90 \
	   dir_field.f90

API_OBJS = $(API_SRCS:.f90=.o)

MAIN_SRCS = frequency.f90
MAIN_OBJS = $(MAIN_SRCS:.f90=.o)


all:	UT-MOM1


%.o: %.f90
ifeq ($(BUILD),Release)
	$(F90) $(COMMON_FLAGS) $(MAIN_PROG_FLAGS) $(RELEASE_FLAGS) -c $< -o $@
else
	$(F90) $(COMMON_FLAGS) $(MAIN_PROG_FLAGS) $(DEBUG_FLAGS) -c $< -o $@
endif

%.o: %.f
ifeq ($(BUILD),Release)
	$(FF) $(COMMON_FLAGS) $(MAIN_PROG_FLAGS) $(RELEASE_FLAGS) -c $< -o $@
else
	$(FF) $(COMMON_FLAGS) $(MAIN_PROG_FLAGS) $(DEBUG_FLAGS) -c $< -o $@
endif


#UT-MOM1: MAIN_PROG_FLAGS = -DENABLED_FOR_MAIN_PROGRAM
UT-MOM1: libutfs_core.a $(MAIN_OBJS) 
ifeq ($(BUILD),Release)
	$(F90) $(COMMON_FLAGS) $(RELEASE_FLAGS) $(MAIN_OBJS) libutfs_core.a $(LIBS) -o $@
#	strip -s $@
else
	$(F90) $(COMMON_FLAGS) $(DEBUG_FLAGS) $(F90OBJS) $(F77OBJS) $(API_OBJS) $(MAIN_OBJS) $(LIBS) -o $@
endif


libutfs_core.a: clean $(F90OBJS) $(F77OBJS) $(API_OBJS)
	ar r $@ $(F90OBJS) $(F77OBJS) $(API_OBJS)
	-mkdir -p $(PLATFORM)/lib
	cp -f $@ $(PLATFORM)/lib/$@


clean:
	rm -f *.o *.out *.time *.mod core UT-MOM1
