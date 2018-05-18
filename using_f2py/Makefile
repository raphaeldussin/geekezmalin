
all: lib_watermass.so creeping_sea.so

lib_watermass.so: lib_watermass.f90
	f2py -c -m lib_watermass lib_watermass.f90

creeping_sea.so: creeping_sea.f90
	f2py -c -m creeping_sea creeping_sea.f90

clean:
	\rm -Rf *.so *.dSYM

getdata:
	wget https://data.nodc.noaa.gov/woa/WOA13/DATA/temperature/netcdf/5564/1.00/woa13_5564_t00_01.nc
	wget https://data.nodc.noaa.gov/woa/WOA13/DATA/salinity/netcdf/5564/1.00/woa13_5564_s00_01.nc

run: lib_watermass.so
	./example_watermass.py
