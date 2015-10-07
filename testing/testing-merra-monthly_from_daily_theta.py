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

theta = merra.monthly_from_daily(year, month, 'theta', fluxes=False,
                                 scratchdir=scratchdir)
atm.save_nc(filename('theta', datestr), theta)

theta_e = merra.monthly_from_daily(year, month, 'theta_e', fluxes=False,
                                   scratchdir=scratchdir)
atm.save_nc(filename('theta_e', datestr), theta_e)

# Fluxes
ds = merra.monthly_from_daily(year, month, 'theta', fluxes=True,
                              scratchdir=scratchdir)
ds.to_netcdf(filename('theta_flx', datestr))

ds = merra.monthly_from_daily(year, month, 'theta_e', fluxes=True,
                              scratchdir=scratchdir)
ds.to_netcdf(filename('theta_e_flx', datestr))
