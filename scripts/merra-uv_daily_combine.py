import sys
sys.path.append('/home/jwalker/dynamics/python/atmos-tools')
sys.path.append('/home/jwalker/dynamics/python/atmos-read')

import os
import numpy as np
import matplotlib.pyplot as plt
import merra
import atmos as atm
from merra import load_daily_season

datadir = '/home/jwalker/eady/datastore/merra/daily/'
savedir = datadir
#plev = 200
plev = 850
years = np.arange(1979, 2015)
season = 'ann'
lon1, lon2 = 40, 120
lat1, lat2 = -60, 60
nperday = 8

def pathstr(var, plev):
    return datadir + 'merra_' + var + str(plev) + '_'

def outfile(year, plev):
    lats = atm.latlon_labels([lat1, lat2], 'lat', '%.0f', deg_symbol=False)
    lons = atm.latlon_labels([lon1, lon2], 'lon', '%.0f', deg_symbol=False)
    subset = '%s-%s_%s-%s' % (lons[0], lons[1], lats[0], lats[1])
    return datadir + 'merra_uv%d_%s_%d.nc' % (plev, subset, year)

for year in years:
    print('Loading U')
    u = load_daily_season(pathstr('u', plev), year, season, 'U',
                          lat1, lat2, lon1, lon2)
    print('Loading V')
    v = load_daily_season(pathstr('v', plev), year, season, 'V',
                          lat1, lat2, lon1, lon2)
    print('Calculating vorticity and Rossby number')
    rel_vort, _ , _ = atm.vorticity(u, v)
    Ro = atm.rossby_num(u, v)

    print('Calculating daily means from 3-hourly data')
    days = np.arange(1, u.shape[0]/nperday + 1)
    u = atm.daily_from_subdaily(u, nperday, dayname='Day', dayvals=days)
    v = atm.daily_from_subdaily(v, nperday, dayname='Day', dayvals=days)
    rel_vort = atm.daily_from_subdaily(rel_vort, nperday, dayname='Day',
                                       dayvals=days)
    Ro = atm.daily_from_subdaily(Ro, nperday, dayname='Day', dayvals=days)

    print('Saving to ' + outfile(year, plev))
    atm.save_nc(outfile(year, plev), u, v, rel_vort, Ro)
