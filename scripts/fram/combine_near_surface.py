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
#months = [4, 5, 6, 7, 8, 9]
months = range(1, 13)
lon1, lon2 = 40, 120
lat1, lat2 = -60, 60
varnms = ['T', 'H', 'QV']
filestr = 'merra_%s850_40E-120E_60S-60N_'
#varnms = ['T', 'H', 'QV', 'V']
#filestr = 'merra_%s950_40E-120E_60S-60N_'
#filestr2 = 'apr-sep_%d.nc'
filestr2 = '%d.nc'
nperday = 8

subset_dict = {'lat' : (lat1, lat2), 'lon' : (lon1, lon2)}

for varnm in varnms:
    datafiles = [datadir + filestr % varnm + filestr2 % year for year in years]
    for y, year in enumerate(years):
        for m, mon in enumerate(months):
            filn = datadir + filestr % varnm + '%d%02d.nc' % (year, mon)
            dayvals = atm.season_days(atm.month_str(mon), atm.isleap(year))
            print('Loading ' + filn)
            with xray.open_dataset(filn) as ds:
                var_in = atm.subset(ds[varnm], subset_dict)
                var_in = atm.squeeze(var_in)
            if m == 0:
                var = var_in
            else:
                var = xray.concat((var, var_in), dim='Day')

        savefile = datafiles[y]
        print('Saving to ' + savefile)
        atm.save_nc(savefile, var)
