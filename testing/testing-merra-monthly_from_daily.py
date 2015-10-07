import sys
sys.path.append('/home/jwalker/dynamics/python/atmos-tools')
sys.path.append('/home/jwalker/dynamics/python/atmos-read')

import merra
import atmos as atm

scratchdir = '/net/eady/data1/jwalker/datastore/scratch/'

def filename(varname, datestr):
    savedir = '/net/eady/data1/jwalker/datastore/merra/monthly/'
    filen = savedir + varname + datestr
    print('Saving to ' + filen)
    return filen

year = 1979
month = 2

datestr = '_%d%02d.nc' % (year, month)

# One variable, surface
precip = merra.monthly_from_daily(year, month, 'precip', fluxes=False,
                                  scratchdir=scratchdir)
atm.save_nc(filename('precip', datestr), precip)

# One variable, pressure-level
T = merra.monthly_from_daily(year, month, 'T', fluxes=False,
                             scratchdir=scratchdir)
atm.save_nc(filename('T', datestr), T)

# Variable and fluxes, pressure-level
ds = merra.monthly_from_daily(year, month, 'q', fluxes=True,
                              scratchdir=scratchdir)
ds.to_netcdf(filename('q_flx', datestr))
