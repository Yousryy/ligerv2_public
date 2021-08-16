EXEC   = BC

OBJS   = main.o src/Input_read.o src/IO.o src/create_galaxy_density.o
#src/power_spec_funcs.o

INCL   = src/vardef.h


#OPT   += -DMEASURETIME

#PAR   += -DNOOPENMP


#OPTIONS =  $(OPT) $(PAR)

#FFTW_INCL = -I/users/aludlow/InstallPrograms/fftw-2.1.5/include/
#FFTW_LIBS = -L/users/aludlow/InstallPrograms/fftw-2.1.5/lib/

CC       = gcc #mpicc

OPTIMIZE =   -O3 -Wall #-pg    # optimization and warning flags (default)


FFTW_LIB =  -lfftw3_threads  -lfftw3
#ifeq (SINGLEFFTW,$(findstring SINGLEFFTW,$(OPT)))
#  FFTW_LIB = $(FFTW_LIBS) -lsrfftw_mpi -lsfftw_mpi -lsrfftw -lsfftw
#else
#  FFTW_LIB = $(FFTW_LIBS) -ldrfftw_mpi -ldfftw_mpi -ldrfftw -ldfftw
#endif

LIBS   =   -lm  $(FFTW_LIB)

CFLAGS =   $(OPTIONS)  $(OPTIMIZE)


GSL_LIB =  -lgsl -lgslcblas
LIBS+= -lgsl
LIBS+= -lgslcblas
$(EXEC): $(OBJS)
	$(CC) $(LPATHS) $(CPATHS) $(OPTIMIZE) $(OBJS) $(LIBS)   -o  $(EXEC)
$(OBJS): $(INCL)
clean:
	rm -f $(OBJS) $(EXEC)
