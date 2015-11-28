import sys
sys.path.append('/home/jwalker/dynamics/python/atmos-tools')
sys.path.append('/home/jwalker/dynamics/python/atmos-read')

import os
import numpy as np
import merra
import atmos as atm

# ----------------------------------------------------------------------
# Save single pressure-level u, v daily data for each month
savedir = atm.homedir() + 'datastore/merra/daily/'
plev = 850
years = np.arange(1979, 2015)
#months = np.arange(1, 13)
months = np.arange(4, 10)

def filename(varname, plev, datestr, savedir):
    filen = savedir + 'merra_%s%d_%s.nc' % (varname, plev, datestr)
    print('Saving to ' + filen)
    return filen

for year in years:
    for month in months:
        datestr = '%d%02d' % (year, month)

        u = merra.read_daily('u', year, month, subset1=('plev', plev, plev))
        atm.save_nc(filename('u', plev, datestr, savedir), u)

        v = merra.read_daily('v', year, month, subset1=('plev', plev, plev))
        atm.save_nc(filename('v', plev, datestr, savedir), v)
