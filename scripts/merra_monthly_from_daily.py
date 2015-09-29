import sys
sys.path.append('/home/jwalker/dynamics/python/atmos-tools')
sys.path.append('/home/jwalker/dynamics/python/atmos-read')

import os
import merra
import atmos as atm

def filename(varname, datestr):
    nm = '/home/jwalker/datastore/merra/monthly/' + varname + datestr
    print('Saving to ' + nm)
    return nm

for year in [1979]:
    for month in range(1, 13):
        datestr = '_%d%02d.nc' % (year, month)
        u, v, uv = merra.monthly_from_daily(year, month, 'u', 'v')
        atm.save_nc(filename('u', datestr), u)
        atm.save_nc(filename('u', datestr), v)
        atm.save_nc(filename('uv', datestr), uv)
