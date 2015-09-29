import sys
sys.path.append('/home/jwalker/dynamics/python/atmos-tools')
sys.path.append('/home/jwalker/dynamics/python/atmos-read')

import os
import merra
import atmos as atm

# ----------------------------------------------------------------------
# prevdir = os.getcwd()
# os.chdir('..')
# ----------------------------------------------------------------------

# Save 200mb u, v daily data for each month

savedir = '/home/jennifer/datastore/merra/daily/'

def filename(varname, datestr, savedir):
    return savedir + 'merra_' + varname + '_' + datestr + '.nc'

for year in range(1979, 1980):
    for month in range(7, 13):
#
# for year in range(1979, 2015):
#     for month in range(1, 13):

        datestr = '%d%02d' % (year, month)

        u = merra.load_daily(year, month, 'u', subset1=('plev', 200, 200))
        atm.save_nc(filename('u200', datestr, savedir), u)

        v = merra.load_daily(year, month, 'v', subset1=('plev', 200, 200))
        atm.save_nc(filename('v200', datestr, savedir), v)


# ----------------------------------------------------------------------
# os.chdir(prevdir)
# ----------------------------------------------------------------------
