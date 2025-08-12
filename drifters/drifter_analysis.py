#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 28 15:28:59 2025

@author: cassandra

For examples see: https://opendrift.github.io/gallery/index.html
"""

from opendrift.models.oceandrift import OceanDrift
from opendrift.readers.reader_netCDF_CF_generic import Reader
from datetime import datetime, timedelta
import numpy as np
import geojson
from shapely.geometry import Polygon, Point
import matplotlib.pyplot as plt
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeature


# Create OpenDrift model instance
o = OceanDrift(loglevel=20)

# Load forecasts as readers
ocean = Reader('./data/currents.nc')
wind = Reader('./data/wind2.nc')
waves = Reader('./data/waves.nc')
o.add_reader([ocean, wind, waves])

# Simulation details
release_time = datetime(2025, 8, 12, 0, 0)
duration = timedelta(days=10)

# Configure some things about the drifters
# No idea what these should really be! 
# Someone needs to read the docs and figure it out...
o.set_config('seed:terminal_velocity', 0.0)       # neutral buoyancy, no sinking/rising
o.set_config('seed:wind_drift_factor', 0.5)       # wind effect -- not sure what this should be
o.set_config('seed:current_drift_factor', 1.0)    # full current effect
o.set_config('drift:stokes_drift', True)          # enable Stokes drift from waves
o.set_config('seed:z', 0)                          # release at surface

# Some stuff to represent uncertainty... 
o.set_config('drift:horizontal_diffusivity', 10)  # typical units: m^2/s; try 5–50 to start
o.set_config('drift:current_uncertainty', 0.1)  # fraction of current speed, e.g., 0.1 = 10%
o.set_config('drift:wind_uncertainty', 0.1)  # 20% variation

# Seed elemetns in a grid
lat_min = 31.2
lon_min = 30.0    
lat_max = 33.3
lon_max = 33.4

for lat in np.linspace(lat_min, lat_max, 50):
    for lon in np.linspace(lon_min, lon_max, 50):
        o.seed_elements(lon=lon, lat=lat, number=10, time=release_time)
        
# Run & save simulation
o.run(duration=duration, time_step=3600, outfile='./simulation.nc')


"""
Read results
"""
# Plot results -- Takes a long time with many drifters! 
o.plot(fast=True, filename='./drifter_waves_wind.png')

# Open simulation
nc = xr.open_dataset('./simulation.nc')

# SIMPLE Analysis
start_lons = nc.lon[:,0].values
start_lats = nc.lat[:,0].values
end_lats = nc.lat[:,-1].values
end_lons = nc.lon[:,-1].values

# Code to check if something reached the Gaza coast
gaza_coast = geojson.load(open('/home/cassandra/Arish/data/gaza_coast.geojson', 'r'))
gaza_coast_polygon_coords = gaza_coast['features'][0]['geometry']['coordinates'][0]
gaza_coast_polygon = Polygon(gaza_coast_polygon_coords)
gx,gy = gaza_coast_polygon.exterior.xy
def reached_gaza(lat, lon):
    return gaza_coast_polygon.contains(Point(lon, lat))

# Evalute if they reached gaza
success = np.array([reached_gaza(la, lo) for la, lo in zip(end_lats, end_lons)])

# Convert to a probability at each point
lats = np.linspace(lat_min, lat_max, 50)
lons = np.linspace(lon_min, lon_max, 50)
p_success = []
for lat in lats:
    lat_close = np.around(np.abs(start_lats - lat), 5) == 0
    p_success.append([])
    for lon in lons:
        lon_close = np.around(np.abs(start_lons - lon), 5) == 0
        drifter_index = np.logical_and(lat_close, lon_close)
        p_success[-1].append(np.sum(success[drifter_index])/len(success[drifter_index]))
p_success = np.array(p_success)    


"""
Plot Result
"""
# Plot some representation of paths and success...
fig = plt.figure(figsize=(6,4))
fig.suptitle("Drifter Fate By Release Location")

# Plot some map lines
sp = ccrs.Stereographic(central_longitude=lon_min, central_latitude=lat_min)
ax = fig.add_subplot(1, 1, 1, projection=sp)
corners = [29.5, 36, 31, 34]
f = cfeature.GSHHSFeature(scale='intermediate', levels=[1])
ax.add_geometries(f.geometries(), ccrs.PlateCarree(),
                  facecolor=cfeature.COLORS['land'],
                  edgecolor='black')
ax.set_extent(corners, crs=ccrs.PlateCarree())
gl = ax.gridlines(ccrs.PlateCarree(), draw_labels=True)

# Plot the probability map
#lolo, lala = np.meshgrid(lons, lats)
m = ax.pcolormesh(lons, lats, 100*p_success, shading='nearest',
                  transform=ccrs.PlateCarree(), cmap='viridis_r')
cb = plt.colorbar(m)
cb.ax.set_ylabel("Probability to reach\nGaza after 10 days")

# Plot Gaza location
plt.plot(gx,gy, label='Gaza Coast', color='magenta',
         transform=ccrs.PlateCarree())
plt.legend(framealpha=1)

plt.savefig('./probability.png')
plt.show()
