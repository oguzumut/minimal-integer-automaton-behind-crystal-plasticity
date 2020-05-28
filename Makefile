PROGS= ~/bin/peng_fourier.exe
OBJS=   mkl_dfti.o real_to_cmplx.o dyn_array.o  main_amorphous_v2.o
INCLUDES=

FFLAGS1= -O3  -fno-range-check
#-openmp
FFLAGS1 = -O3

F90 = ifort


  FFTWLIBS= 
 STDLIBS  =    
 INCLUDES =    
#MKLROOT = /opt/intel/compilers_and_libraries_2017.0.102/mac/mkl
STDLIBS  =   ${MKLROOT}/lib/libmkl_blas95_lp64.a ${MKLROOT}/lib/libmkl_lapack95_lp64.a -L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl

INCLUDES =      -I${MKLROOT}/include/intel64/lp64 -I${MKLROOT}/include   

STDLIBS=   -L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl

INCLUDES =  -m64 -I${MKLROOT}/include

SDLIBS =   ${MKLROOT}/lib/intel64/libmkl_blas95_lp64.a ${MKLROOT}/lib/intel64/libmkl_lapack95_lp64.a -liomp5 -lpthread -lm -ldl

INCLUDES =   -I${MKLROOT}/include/intel64/lp64 -mkl=parallel 

all: $(PROGS)



$(PROGS):$(OBJS)
	$(F90) $(FFLAGS1) $(INCLUDES) -o $(PROGS)  $(OBJS) $(FFTWLIBS) $(STDLIBS)



.SUFFIXES:      .f90 .o

.f90.o:
	$(F90) $(FFLAGS1) $(INCLUDES)  $(FOPTS) $(USER_DEFS) -c $<

clean:
	\rm -f *.bak *.o 

clean_all:
	\rm -f *.bak *.o $(PROGS)

