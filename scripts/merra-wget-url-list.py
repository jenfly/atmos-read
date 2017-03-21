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
savedir = atm.homedir() + 'dynamics/python/atmos-read/scripts/merra_urls/'
savefiles = {yr : savedir + 'merra2_rad_%d.txt' % yr for yr in years}
# ----------------------------------------------------------------------

def get_url(year, mon, day, lat1='-65', lat2='65', lon1='40', lon2='120'):
    yrstr = '%d' % year
    mstr = '%02d' % mon
    dstr = '%02d' % day

    prod_dict = {yr : '100' for yr in range(1980, 1992)}
    for yr in range(1992, 2001):
        prod_dict[yr] = '200'
    for yr in range(2001, 2011):
        prod_dict[yr] = '300'
    for yr in range(2011, 2016):
        prod_dict[yr] = '400'
    prod = prod_dict[year]

    url = ('http://goldsmr4.gesdisc.eosdis.nasa.gov/daac-bin/OTF/' +
           'HTTP_services.cgi?FILENAME=%2Fdata%2FMERRA2%2FM2T1NXRAD.5.12.4%2F' +
            yrstr + '%2F' + mstr + '%2FMERRA2_' + prod + '.tavg1_2d_rad_Nx.' +
            yrstr + mstr + dstr + '.nc4&FORMAT=bmM0Lw&BBOX=' + lat1 +
            '%2C' + lon1 + '%2C' + lat2 + '%2C' + lon2 + '&LABEL=MERRA2_' +
            prod + '.tavg1_2d_rad_Nx.' + yrstr + mstr + dstr +
            '.SUB.nc4&FLAGS=1&SHORTNAME=M2T1NXRAD&SERVICE=SUBSET_MERRA2' +
            '&LAYERS=&VERSION=1.02&VARIABLES=lwgnt%2Clwtup%2Cswgnt%2Cswtnt')
    return url

for year in years:
    print(year)
    with open(savefiles[year], 'w') as f:
        for mon in range(1, 13):
            print('  %d' % mon)
            days = range(1, atm.days_this_month(year, mon) + 1)
            for day in days:
                print('    %d' % day)
                url = get_url(year, mon, day)
                f.write(url + '\n')
