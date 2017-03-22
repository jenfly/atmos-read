import sys
sys.path.append('/home/jwalker/dynamics/python/atmos-tools')
sys.path.append('/home/jwalker/dynamics/python/atmos-read')

import os
import numpy as np
import merra
import atmos as atm

# ----------------------------------------------------------------------
# Save single pressure-level daily data for each month
savedir = atm.homedir() + 'datastore/merra/daily/'
#plev = 950
plev = 850
lon1, lon2 = 40, 120
lat1, lat2 = -60, 60
filestr = savedir + 'merra_%s%d_40E-120E_60S-60N_%d%02d.nc'
#varnms = ['T', 'H', 'QV', 'V']
varnms = ['T', 'H', 'QV']
#years = np.arange(1979, 2015)
#months = np.arange(4, 10)
years = np.arange(1979, 2000)
months = np.arange(1, 13)
nperday = 8

def filename(varname, plev, year, month, filestr):
    filen = filestr % (varname, plev, year, month)
    print('Saving to ' + filen)
    return filen

subset_dict = {'plev' : (plev, plev), 'lon' : (lon1, lon2), 'lat' : (lat1, lat2)}

for varnm in varnms:
    for year in years:
        for m, month in enumerate(months):
            datestr = '%d%02d' % (year, month)        
            print varnm
            var = merra.read_daily(varnm, year, month, subset_dict=subset_dict)

            # Compute daily means of 3-hourly data
            dayvals = atm.season_days(atm.month_str(month), atm.isleap(year))
            var = atm.daily_from_subdaily(var, nperday, dayname='Day',
                                          dayvals=dayvals)
            # Save to file
            atm.save_nc(filename(varnm, plev, year, month, filestr), var)
