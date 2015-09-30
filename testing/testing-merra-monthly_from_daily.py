import sys
sys.path.append('/home/jwalker/dynamics/python/atmos-tools')
sys.path.append('/home/jwalker/dynamics/python/atmos-read')

import os
import merra
import atmos as atm

def filename(varname, datestr):
    savedir = '/net/eady/data1/jwalker/datastore/merra/monthly/'
    filen = savedir + varname + datestr
    print('Saving to ' + filen)
    return filen

year = 1979
month = 1

datestr = '_%d%02d.nc' % (year, month)

# One variable, surface
precip = merra.monthly_from_daily(year, month, 'precip')
atm.save_nc(filename('precip', datestr), precip)

# One variable, pressure-level
T = merra.monthly_from_daily(year, month, 'T')
atm.save_nc(filename('T', datestr), T)

# Variable and fluxes, pressure-level
q, uq, vq = merra.monthly_from_daily(year, month, 'q', fluxes=True)
atm.save_nc(filename('q_flx', datestr), q, uq, vq)
