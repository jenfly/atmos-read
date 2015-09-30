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

theta = merra.monthly_from_daily(year, month, 'theta')
atm.save_nc(filename('theta', datestr), theta)

theta_e = merra.monthly_from_daily(year, month, 'theta_e')
atm.save_nc(filename('theta_e', datestr), theta_e)

# Fluxes
theta, utheta, vtheta = merra.monthly_from_daily(year, month, 'theta',
                                                 fluxes=True)
atm.save_nc(filename('theta_flx', datestr), theta, utheta, vtheta)

theta_e, utheta_e, vtheta_e = merra.monthly_from_daily(year, month, 'theta_e',
                                                       fluxes=True)
atm.save_nc(filename('theta_e_flx', datestr), theta_e, utheta_e, vtheta_e)
