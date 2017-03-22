import sys
sys.path.append('/home/jwalker/dynamics/python/atmos-tools')
sys.path.append('/home/jwalker/dynamics/python/atmos-read')

import os
import numpy as np
import matplotlib.pyplot as plt
import merra
import atmos as atm
from merra import load_daily_season

datadir = '/home/jwalker/datastore/merra/daily/'
savedir = datadir
#plev = 200
plev = 850
varnm = 'T'
years = np.arange(1979, 2015)
season = 'ann'
lon1, lon2 = 40, 120
lat1, lat2 = -60, 60
nperday = 8

def pathstr(var, plev):
    return datadir + 'merra_' + var + str(plev) + '_'

def outfile(year, var, plev):
    lats = atm.latlon_labels([lat1, lat2], 'lat', '%.0f', deg_symbol=False)
    lons = atm.latlon_labels([lon1, lon2], 'lon', '%.0f', deg_symbol=False)
    subset = '%s-%s_%s-%s' % (lons[0], lons[1], lats[0], lats[1])
    return datadir + 'merra_%s%d_%s_%d.nc' % (var, plev, subset, year)

for year in years:
    print('Loading data')
    data = load_daily_season(pathstr(varnm, plev), year, season, varnm,
                             lat1, lat2, lon1, lon2)

    print('Calculating daily means from 3-hourly data')
    days = np.arange(1, data.shape[0]/nperday + 1)
    data = atm.daily_from_subdaily(data, nperday, dayname='Day', dayvals=days)

    savefile = outfile(year, varnm, plev)
    print('Saving to ' + savefile)
    atm.save_nc(savefile, data)
