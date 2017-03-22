import sys
sys.path.append('/home/jwalker/dynamics/python/atmos-tools')
sys.path.append('/home/jwalker/dynamics/python/atmos-read')

import os
import xray
import numpy as np
import matplotlib.pyplot as plt
import merra
import atmos as atm

# ----------------------------------------------------------------------

datadir = atm.homedir() + 'datastore/merra/daily/'
years = np.arange(1979, 2015)
months = [4, 5, 6, 7, 8, 9]
lon1, lon2 = 40, 120
lat1, lat2 = -60, 60
varnm = 'T'
filestr = 'merra_T200_'
filestr2 = '40E-120E_60S-60N_apr-sep_%d.nc'
datafiles = [datadir + filestr + filestr2 % year for year in years]
nperday = 8

for y, year in enumerate(years):
    for m, mon in enumerate(months):
        filn = datadir + filestr + '%d%02d.nc' % (year, mon)
        dayvals = atm.season_days(atm.month_str(mon), atm.isleap(year))
        print('Loading ' + filn)
        with xray.open_dataset(filn) as ds:
            var_in = atm.subset(ds[varnm], 'lat', lat1, lat2, 'lon', lon1, lon2)
            var_in = atm.squeeze(var_in)
        # Compute daily means of 3-hourly data
        var_in = atm.daily_from_subdaily(var_in, nperday, dayname='Day',
                                         dayvals=dayvals)
        if m == 0:
            var = var_in
        else:
            var = xray.concat((var, var_in), dim='Day')

    savefile = datafiles[y]
    print('Saving to ' + savefile)
    atm.save_nc(savefile, var)
