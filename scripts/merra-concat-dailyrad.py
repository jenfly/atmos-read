import sys
sys.path.append('/home/jwalker/dynamics/python/atmos-tools')
sys.path.append('/home/jwalker/dynamics/python/atmos-read')
sys.path.append('/home/jwalker/dynamics/python/monsoon-onset')

import xarray as xray
import numpy as np
import collections
import pandas as pd
import scipy

import atmos as atm

years = range(1980, 2016)
datadir = atm.homedir() + 'datastore/merra2/dailyrad/'
savefiles = {yr : datadir + 'merra2_RAD_%d.nc4' % yr for yr in years}
# ----------------------------------------------------------------------

def get_filenames(datadir, year):

    prod_dict = {yr : '100' for yr in range(1980, 1992)}
    for yr in range(1992, 2001):
        prod_dict[yr] = '200'
    for yr in range(2001, 2011):
        prod_dict[yr] = '300'
    for yr in range(2011, 2016):
        prod_dict[yr] = '400'
    prod = prod_dict[year]

    filestr = (datadir + '%d/' % year + 'MERRA2_' + prod
               + '.tavg1_2d_rad_Nx.%d%02d%02d.SUB.nc4')
    files = []
    for mon in range(1, 13):
        days = range(1, atm.days_this_month(year, mon) + 1)
        for day in days:
            files.append(filestr % (year, mon, day))

    return files

for year in years:
    files = get_filenames(datadir, year)
    savefile = savefiles[year]
    ds = atm.load_concat(files, concat_dim='time')
    ds = ds.rename({'time' : 'day'})
    ds = ds.drop(['time_bnds', 'bnds'])
    ds['day'] = range(1, len(ds['day']) + 1)
    print('Saving to ' + savefile)
    ds.to_netcdf(savefile)
