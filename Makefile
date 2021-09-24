F90C	= ifort
FFLAGS = -O2 -cm -w -vec_report0 -sox -openmp -openmp_report0
LAPACKLIB = -L/usr/local/lapack-3.4.1 -llapack -lblas -ltmglib

F90FLAGS = -DMATRIX_SINGLE $(FFLAGS)

DISTFILES =  utils.o inifile.o Matrix_utils.o settings.o GetDist.o

OBJFILES= modelfit_in.o fftlog.o cdgamma.o drfft.o utils.o inifile.o Matrix_utils.o settings.o \
	paramdef.o propose.o calclike.o MCMC.o driver.o

default: cosmomc

all : cosmomc getdist

.f.o:
	$(F90C) $(F90FLAGS) -c $<

%.o: %.f90
	$(F90C) $(F90FLAGS) -c $*.f90

%.o: %.F90
	$(F90C) $(F90FLAGS) -c $*.F90


cosmomc: $(OBJFILES) 	
	$(F90C) -o cosmomc $(OBJFILES) $(LAPACKLIB) $(F90FLAGS)

run:
	mpirun -np 1 ./cosmomc params.ini

clean:
	rm -f *.o *.mod *.d *.pc *.obj core

getdist: $(DISTFILES)
	$(F90C) -o getdist $(DISTFILES) $(LAPACKLIB) $(F90FLAGS) 
