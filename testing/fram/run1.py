import sys
sys.path.append('/home/jwalker/dynamics/python/atmos-tools')
sys.path.append('/home/jwalker/dynamics/python/atmos-read')

import atmos as atm
import merra
from merra import calc_fluxes

scratchdir = '/net/eady/data1/jwalker/datastore/scratch/'

def filename(varname, datestr):
    savedir = '/net/eady/data1/jwalker/datastore/merra/monthly/'
    filen = savedir + varname + datestr
    print('Saving to ' + filen)
    return filen

year = 1979
month = 3

datestr = '_%d%02d.nc' % (year, month)

ds = calc_fluxes(year, month, scratchdir=scratchdir)
ds.to_netcdf(filename('fluxes', datestr))
