import sys
sys.path.append('/home/jwalker/dynamics/python/atmos-tools')
sys.path.append('/home/jwalker/dynamics/python/atmos-read')

import os
import numpy as np
import matplotlib.pyplot as plt
import merra
import atmos as atm
from merra import load_daily_season

datadir = '/home/jennifer/datastore/merra/daily/'
#year = 1995
year = 1979

ustr = datadir + 'merra_u200_'
vstr = datadir + 'merra_v200_'

season = 'jan'
u = load_daily_season(ustr, year, season, 'U')
plt.figure()
atm.pcolor_latlon(u.mean(dim='TIME'))

season = 'jja'
lon1, lon2 = 0, 100
lat1, lat2 = -20, 50
u = load_daily_season(ustr, year, season, 'U', lat1, lat2, lon1, lon2)
plt.figure()
atm.pcolor_latlon(u.mean(dim='TIME'))

ds = load_daily_season(ustr, year, season, None, lat1, lat2, lon1, lon2)
u2 = ds['U']
print((u == u2).any())

season = 'ann'
lon1, lon2 = 20, 120
lat1, lat2 = -60, 60
n = 8   # Number of time points per day
days = np.arange(1, 366)
u = load_daily_season(ustr, year, season, 'U', lat1, lat2, lon1, lon2)
udaily = atm.daily_from_subdaily(u, n, dayvals=days)
plt.figure()
atm.pcolor_latlon(udaily.mean(axis=0), axlims=(lat1,lat2,lon1,lon2))

lat1, lat2 = 10, 30
lon1, lon2 = 60, 100
ubar = atm.mean_over_geobox(u, lat1, lat2, lon1, lon2)
ubar_daily = atm.mean_over_geobox(udaily, lat1, lat2, lon1, lon2)
plt.figure(figsize=(7,8))
plt.subplot(211)
plt.plot(ubar)
plt.subplot(212)
plt.plot(ubar_daily)
