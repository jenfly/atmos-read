import sys
sys.path.append('/home/jwalker/dynamics/python/atmos-tools')
sys.path.append('/home/jwalker/dynamics/python/atmos-read')

import os
import merra
import atmos as atm

def filename(varname, datestr):
    return '/home/jwalker/eady/merra/scratch/' + varname + datestr

year = 1979
month = 1

datestr = '_%d%02d.nc' % (year, month)

# One variable, surface
precip = merra.monthly_from_daily(year, month, 'precip')
atm.save_nc(filename('precip', datestr), precip)

# One variable, pressure-level
u = merra.monthly_from_daily(year, month, 'u')
atm.save_nc(filename('u', datestr), u)

# Two variables, surface
precip, evap, precip_evap = merra.monthly_from_daily(year, month, 'precip',
                                                     'evap')
atm.save_nc(filename('precip_evap', datestr), precip, evap, precip_evap)

# Two variables, pressure-level
u, v, uv = merra.monthly_from_daily(year, month, 'u', 'v')
atm.save_nc(filename('uv', datestr), u, v, uv)
