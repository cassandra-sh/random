#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 28 15:39:40 2025

@author: cassandra
"""
import copernicusmarine
from datetime import datetime, timedelta
import pytz
import xarray as xr
import numpy as np

# Note to self: don't store passwords in plain text!
# This generates a file. Need only run once, I think... 
copernicusmarine.login(username='', 
                       password='')

# Get today's date and make a datetime
# For some reason datetime.today() doesn't work without
# this step
d = datetime.today().date()
day, month, year = d.day, d.month, d.year
today = datetime(day=day, year=year, month=month)

# Today + 10 days for forecast
date_start = today
date_end = today + timedelta(days=10)

# a 20 day subset, 2 days ago, for hindcast
wind_start = today - timedelta(days=22)
wind_end   = today - timedelta(days=2)

# Bounds for eastern Mediterranean
lat_min = 29.88
lon_min = 27.07
lat_max = 35.38
lon_max = 36.78

# surface currents (2 components)
ds_cur = copernicusmarine.subset(
    dataset_id="cmems_mod_med_phy-cur_anfc_4.2km_PT15M-i",
    variables=["uo", "vo"],
    minimum_longitude=lon_min, maximum_longitude=lon_max,
    minimum_latitude=lat_min, maximum_latitude=lat_max,
    start_datetime=date_start, end_datetime=date_end,
    output_directory="./data", output_filename="currents.nc")

# significant wave height & stokes drift
ds_wav = copernicusmarine.subset(
    dataset_id="cmems_mod_med_wav_anfc_4.2km_PT1H-i",
    variables=["VHM0", "VSDX", "VSDY"],  
    minimum_longitude=lon_min, maximum_longitude=lon_max,
    minimum_latitude=lat_min, maximum_latitude=lat_max,
    start_datetime=date_start, end_datetime=date_end,
    output_directory="./data", output_filename="waves.nc")

# 10m wind components eastward and northward
# note -- wind forecast not available
ds_wind = copernicusmarine.subset(
    dataset_id="cmems_obs-wind_glo_phy_nrt_l4_0.125deg_PT1H",
    variables=["eastward_wind", "northward_wind"],  
    minimum_longitude=lon_min, maximum_longitude=lon_max,
    minimum_latitude=lat_min, maximum_latitude=lat_max,
    start_datetime=wind_start, end_datetime=wind_end,
    output_directory="./data", output_filename="wind.nc")

"""
Wind Forecast
-- just project the average day, into the future... 
"""
# Generate a forecast based on the hindcast. 
# Its going to be extremely simple, just taking the average
# day (at each point) and repeating it for the forecast extent

# Load the hindcast
wind = xr.open_dataset('./data/wind.nc')

# Generate a 1 day timedelta timestamp, with the same dt as the hindcast
t = wind.time.values
dt = np.nanmean(np.diff(t)).astype('timedelta64[m]')
t_day = np.arange(0, np.timedelta64(24, 'h'), dt)
day = np.timedelta64(1, 'D').astype(dt.dtype)

# Convert the timestamp to an index for each hour of the day
# Need some fancy indexing math...
time_elapsed = t-t[0]
time_elapsed_per_day = time_elapsed.astype('timedelta64[m]').astype(int) % day.astype(int)
per_day_index = np.around(time_elapsed_per_day / dt.astype(int)).astype(int)

# Now for each step throughout the average day, evaluate the average
east_avg  = []
north_avg = []
for i in np.unique(per_day_index):
    east_avg.append( np.nanmean(wind.eastward_wind[ np.where(per_day_index == i)[0]], axis=0))
    north_avg.append(np.nanmean(wind.northward_wind[np.where(per_day_index == i)[0]], axis=0))

# Now project the average day into the future (lets say, 15 days)
# Start with a timestamp
t_future = np.arange(t[-1], t[-1] + 15 * np.timedelta64(24, 'h'), dt)
# Repeat the average into the future
east_future = np.concat([east_avg for d in range(15)])
north_future = np.concat([north_avg for d in range(15)])

# Create new DataArray with new time and values
d = {"coords": {"time": t_future,
                 "longitude": wind.longitude.values,
                 "latitude": wind.latitude.values},
     "attrs": wind.attrs,
     "data_vars": {"eastward_wind": (["time", "latitude", "longitude"], east_future),
                    "northward_wind": (["time", "latitude", "longitude"], north_future)}}
new_wind = xr.Dataset(data_vars=d["data_vars"],
                      coords=d["coords"],
                      attrs=d["attrs"])

# Save
new_wind.to_netcdf("./data/wind2.nc")







