import sys
sys.path.append('/home/jwalker/dynamics/python/atmos-tools')
sys.path.append('/home/jwalker/dynamics/python/atmos-read')

import xray
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
import pandas as pd
import atmos as atm
import precipdat
import merra


def filename(year, varnm, subset_str):
    savedir = atm.homedir() + 'datastore/merra/daily/'
    filn = savedir + 'merra_%s_%s_%d.nc' % (varnm, subset_str, year)
    print('Saving to ' + filn)
    return filn

# MERRA daily surface fluxes
varnms = ['EFLUX', 'HFLUX', 'EVAP']
years = range(1979, 2015)
months = range(1, 13)
lon1, lon2 = 40, 120
subset_dict = {'lon' : (lon1, lon2)}
subset_str = '%dE-%dE_90S-90N' % (lon1, lon2)
nperday = 24

for varnm in varnms:
    for yr in years:
        for mon in months:
            var1 = merra.read_daily(varnm, yr, mon, subset_dict=subset_dict)
            dayvals = atm.season_days(atm.month_str(mon), atm.isleap(yr))
            var1 = atm.daily_from_subdaily(var1, nperday, dayvals=dayvals)
            if mon == 1:
                var = var1
            else:
                var = xray.concat([var, var1], dim='day')

        # Save daily mean data for this year
        atm.save_nc(filename(yr, varnm, subset_str), var)
