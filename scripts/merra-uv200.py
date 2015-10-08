import sys
sys.path.append('/home/jwalker/dynamics/python/atmos-tools')
sys.path.append('/home/jwalker/dynamics/python/atmos-read')

import os
import merra
import atmos as atm

# ----------------------------------------------------------------------
# Save 200mb u, v daily data for each month
savedir = '/home/jennifer/datastore/merra/daily/'

def filename(varname, datestr, savedir):
    filen = savedir + 'merra_' + varname + '_' + datestr + '.nc'
    print('Saving to ' + filen)
    return filen

for year in range(1988, 1996):
    for month in range(1, 13):

        datestr = '%d%02d' % (year, month)

        v = merra.load_daily(year, month, 'v', subset1=('plev', 200, 200))
        atm.save_nc(filename('v200', datestr, savedir), v)
