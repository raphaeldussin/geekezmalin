#!/usr/bin/env python

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
import netCDF4 as nc
import matplotlib.cm as cm
import matplotlib.colors as cl

def readnc(filein,varin):
        ''' read data from netcdf file '''
        fid = nc.Dataset(filein,'r')
        out = fid.variables[varin][:].squeeze()
        fid.close()
        return out

def setup_map_ccs(plt_topo=True,topofile='etopo_ccs.nc'):
        ''' set the map for california current system, with etopo optionally '''

	plt.figure(figsize=[8.,8.])
	ax = plt.axes(projection=ccrs.PlateCarree())
	ax.set_extent([-142, -108, 18, 51])
	ax.coastlines(resolution='50m')
	ax.gridlines(draw_labels=False,xlocs=np.arange(-150,-90,10),ylocs=np.arange(10,70,10))
	ax.set_xticks(np.arange(-140,-100,10))
	ax.set_yticks(np.arange(20,60,10))

        if plt_topo:
                paltopo = cm.binary
                lon_topo = readnc(topofile,'topo_lon') + 360.
                lat_topo = readnc(topofile,'topo_lat')
                lon_topo2, lat_topo2 = np.meshgrid(lon_topo,lat_topo)
                topo = readnc(topofile,'topo')
                topomin=0 ; topomax=4000
                topo = np.ma.masked_less_equal(topo,0.)
                normtopo = cl.Normalize(vmin=topomin, vmax=topomax)
                T = ax.contourf(lon_topo2,lat_topo2,topo,100,cmap=paltopo,norm=normtopo)
	else:
		ax.stock_img()
        return ax


#--- read data from WOA annual temperature ---
lon = readnc('./woa13_5564_t00_01.nc','lon')
lat = readnc('./woa13_5564_t00_01.nc','lat')
temp = readnc('./woa13_5564_t00_01.nc','t_an')[0,:,:]

#--- set parameters for contour plot ---
pal = cm.gist_ncar
cbarfmt = "%01g"
norm = cl.Normalize(vmin=6, vmax=26)
contours=np.arange(6,26+0.2,0.2)
ticks=np.arange(6.,26+3,3.)

#--- try a first plot using cartopy default background ---
ax = setup_map_ccs(plt_topo=False)
C = ax.contourf(lon,lat,temp,contours,cmap=pal,norm=norm)
plt.colorbar(C,format=cbarfmt,shrink=0.8,ticks=ticks)
plt.show()
plt.close()

# UGLY ! right ?
# the only background image is at downsampled too much and
# we don't want ocean background where WOA is not present

#--- now let's plot the topography from the ETOPO 5min file ---
ax = setup_map_ccs()
C = ax.contourf(lon,lat,temp,contours,cmap=pal,norm=norm)
plt.colorbar(C,format=cbarfmt,shrink=0.8,ticks=ticks)
plt.show()
plt.close()

# feeling better ? so do I :)
