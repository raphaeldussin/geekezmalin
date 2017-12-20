#!/usr/bin/env python

import netCDF4 as nc
import lib_watermass
import numpy as np

#---------------------------------------------------------------
# sample data set (use make getdata to download WOA annual data)
fnameT = 'woa13_5564_t00_01.nc'
fnameS = 'woa13_5564_s00_01.nc'

fid        = nc.Dataset(fnameT,'r')
lon        = fid.variables['lon'][:].squeeze()
lat        = fid.variables['lat'][:].squeeze()
depth      = fid.variables['depth'][:].squeeze()
lon_bnds   = fid.variables['lon_bnds'][:].squeeze()
lat_bnds   = fid.variables['lat_bnds'][:].squeeze()
depth_bnds = fid.variables['depth_bnds'][:].squeeze()
temp       = fid.variables['t_an'][:].squeeze()
spval      = fid.variables['t_an']._FillValue
fid.close()

fid        = nc.Dataset(fnameS,'r')
salt       = fid.variables['s_an'][:].squeeze()
fid.close()

#---------------------------------------------------------------
# compute metrics (dx,dy,dz) of ocean grid cells and mask
Rearth = 6378e+3
nx = lon.shape[0]
ny = lat.shape[0]
nz = depth.shape[0]

dz = np.empty((nz,ny,nx))

dxflat = (2*np.pi*Rearth/360) * (lon_bnds[:,1] - lon_bnds[:,0])
dyflat = (2*np.pi*Rearth/360) * (lat_bnds[:,1] - lat_bnds[:,0])
dzcolumn = depth_bnds[:,1] - depth_bnds[:,0]

dxflat2d, dyflat2d = np.meshgrid(dxflat,dyflat)
lon2d, lat2d = np.meshgrid(lon,lat)

dx = dxflat2d * np.cos(2*np.pi*lat2d/360)
dy = dyflat2d

for jj in np.arange(ny):
	for ji in np.arange(nx):
		dz[:,jj,ji] = dzcolumn

mask = np.ones(temp.shape)
mask[np.where(temp.mask == True)] = 0

#---------------------------------------------------------------
# sanity checks
# compute water volume between 50C and 55C (never observed)
wmass_testimpossible = lib_watermass.volume_watermass_from_ts(dx.T,dy.T,dz.T,temp.T,salt.T,50.,55.,0.,40.)

# compute all possible volume
wmass_testall = lib_watermass.volume_watermass_from_ts(dx.T,dy.T,dz.T,temp.T,salt.T,-10.,100.,0.,50.)

# and compare to volume from metrics
wmass_from_metrics = (dz * dx * dy * mask).sum()

print '---Control checks---'
print 'Impossible water mass, volume =', wmass_testimpossible
print 'All possible water, volume =', wmass_testall
print 'volume from scale factors = ', wmass_from_metrics

#---------------------------------------------------------------
# now let's learn some stuff
# compute volume of water 20C < T < 40C
wmass_20to40C = lib_watermass.volume_watermass_from_ts(dx.T,dy.T,dz.T,temp.T,salt.T,20.,40.,0.,40.)
# compute volume of water 10C < T < 20C
wmass_10to20C = lib_watermass.volume_watermass_from_ts(dx.T,dy.T,dz.T,temp.T,salt.T,10.,20.,0.,40.)
# compute volume of water 0C < T < 10C
wmass_0to10C = lib_watermass.volume_watermass_from_ts(dx.T,dy.T,dz.T,temp.T,salt.T,0.,10.,0.,40.)

print '---Volume of ocean in temperature ranges---'
print 'The percentage of ocean waters 20C < T < 40C is ', 100 * wmass_20to40C / wmass_from_metrics, '%'
print 'The percentage of ocean waters 10C < T < 20C is ', 100 * wmass_10to20C / wmass_from_metrics, '%'
print 'The percentage of ocean waters  0C < T < 10C is ', 100 * wmass_0to10C / wmass_from_metrics, '%'

