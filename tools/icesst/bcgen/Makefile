# Set up special characters
null  :=

# Check for netcdf locations

# FFLAGS and LDFLAGS will have at least these values in all situations
# netCDF/4.9.0-iimpi-2022a
FFLAGS  = -c 
FFLAGS += $(shell nc-config --fflags)
LDFLAGS = $(shell nc-config --flibs)

# Determine platform 
UNAMES := $(shell uname -s)

DEBUG=TRUE

#------------------------------------------------------------------------
# Linux
#------------------------------------------------------------------------

ifeq ($(UNAMES),Linux)
  ifeq ($(USER_FC),pgf90)
    FC = pgf90
    ifeq ($(DEBUG),TRUE)
      FFLAGS += -g -Mbounds
    else
      FFLAGS += -fast
    endif
  else
    FC = ifort
    ifeq ($(DEBUG),TRUE)
      FFLAGS += -gdwarf-4
    else
      FFLAGS += -O3
    endif
  endif
endif

#------------------------------------------------------------------------
# SGI
#------------------------------------------------------------------------

ifeq ($(UNAMES),IRIX64)
  FC = f90
  FFLAGS  += -64 
  LDFLAGS += -64
  ifeq ($(DEBUG),TRUE)
    FFLAGS  += -C -g
    LDFLAGS += -g
  else
    FFLAGS += -O2
  endif
endif

#------------------------------------------------------------------------
# IBM
#------------------------------------------------------------------------

ifeq ($(UNAMES),AIX)
  FC = xlf90
  FFLAGS += -qsuffix=f=f90
  ifeq ($(DEBUG),TRUE)
    FFLAGS  += -C -g
    LDFLAGS += -g
  else
    FFLAGS += -O2
  endif
endif

OBJS = calcclim.o calcfull.o driver.o bcgen.o output_dateinfo.o prec.o \
        setup_outfile.o solver.o wrap_nf.o types.o

bcgen: $(OBJS)
	$(FC) -o $@ $(OBJS) $(LDFLAGS)

clean:
	$(RM) bcgen *.o *.mod *.stb *.MOD

.SUFFIXES:
.SUFFIXES: .f90 .o

.f90.o:
	$(FC) $(FFLAGS) $*.f90

calcclim.o: prec.o solver.o
calcfull.o: prec.o solver.o
driver.o: prec.o
bcgen.o: prec.o types.o
output_dateinfo.o: prec.o
prec.o:
setup_outfile.o: prec.o types.o
solver.o: prec.o
wrap_nf.o: prec.o
types.o: prec.o
