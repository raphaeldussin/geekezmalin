# Makefile for large analysis tools ROMS
# (Raf Dussin)

F90=gfortran
NC_INC=$(shell nc-config --includedir)
NC_LIB=$(shell nc-config --flibs)
FFLAGS= -Wall -O -I$(NC_INC) $(NC_LIB) 

EXEC = compute_mean_Tsquared
       
all: $(EXEC)

# rules
compute_mean_Tsquared: compute_mean_Tsquared.f90
	$(F90) compute_mean_Tsquared.f90 -o compute_mean_Tsquared $(FFLAGS)

clean:
	\rm compute_mean_Tsquared
