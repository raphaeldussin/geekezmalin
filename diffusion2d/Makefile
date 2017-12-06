### Makefile

### OSX/Linux with gfortran
F90=gfortran 
NCDF=-I$(shell nf-config --includedir) $(shell nf-config --flibs)
FFLAGS= -Wall -O $(NCDF) -ffree-line-length-none -fconvert=big-endian

EXEC = diffusion2d

all: $(EXEC)

diffusion2d: diffusion2d.f90
	$(F90) diffusion2d.f90 -o diffusion2d $(FFLAGS)

clean:
	\rm diffusion2d

run: diffusion2d
	./diffusion2d namelist
