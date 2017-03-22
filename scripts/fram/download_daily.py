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
plev = 200
varnm = 'T'
years = np.arange(1979, 2015)
#months = np.arange(1, 13)
months = [1, 2, 3, 10, 11, 12]

def filename(varname, plev, datestr, savedir):
    filen = savedir + 'merra_%s%d_%s.nc' % (varname, plev, datestr)
    print('Saving to ' + filen)
    return filen

for year in years:
    for month in months:
        datestr = '%d%02d' % (year, month)

        var = merra.read_daily(varnm, year, month, subset1=('plev', plev, plev))
        atm.save_nc(filename(varnm, plev, datestr, savedir), var)
