import sys
sys.path.append('/home/jwalker/dynamics/python/atmos-tools')
sys.path.append('/home/jwalker/dynamics/python/atmos-read')

import os
import numpy as np
import matplotlib.pyplot as plt
import xray
import merra
import atmos as atm

datadir = '/home/jennifer/datastore/merra/daily/'

filestr = 'merra_uv200_40E-120E_60S-60N_'
year = 1980
lon1, lon2 = 60, 100
lat1, lat2 = 10, 30

filename = '%s%s%d.nc' % (datadir, filestr, year)
ds = atm.ncload(filename)

iplot = {'U' : 1, 'V' : 2, 'rel_vort' : 3, 'Ro' : 4}

plt.figure(figsize=(12,9))
for var in ds.data_vars:
    plt.subplot(2, 2, iplot[var])
    atm.pcolor_latlon(ds[var].mean(axis=0))
    plt.title(var)

dsbar = xray.Dataset()
for var in ds.data_vars:
    print(var)
    dsbar[var] = atm.mean_over_geobox(ds[var], lat1, lat2, lon1, lon2)

plt.figure(figsize=(12, 9))
for var in dsbar.data_vars:
    plt.subplot(2, 2, iplot[var])
    plt.plot(dsbar.Day, dsbar[var])
    plt.title(var)
