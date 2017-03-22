import sys
sys.path.append('/home/jwalker/dynamics/python/atmos-tools')
sys.path.append('/home/jwalker/dynamics/python/atmos-read')

import os
import merra
import atmos as atm

def filename(varname, datestr):
    nm = '~/datastore/merra/monthly/' + varname + datestr
    print('Saving to ' + nm)
    return nm

flux_vars = ['u', 'q', 'T', 'hgt']
noflux_vars = ['v', 'omega', 'ps']

for year in [1979]:
    for month in range(1, 13):
        datestr = '_%d%02d.nc' % (year, month)

        for var in flux_vars:
            print('*******' + var + '*********')
            var, u_var, v_var = merra.monthly_from_daily(year, month, var,
                                                         fluxes=True)
            atm.save_nc(filename(var + '_flx', datestr), var, u_var, v_var)

        for var in noflux_vars:
            print('*******' + var + '*********')
            var  = merra.monthly_from_daily(year, month, var, fluxes=False)
            atm.save_nc(filename(var, datestr), var)
