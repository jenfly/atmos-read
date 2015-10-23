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


def filename(year):
    savedir = atm.homedir() + 'datastore/merra/daily/'
    filn = savedir + 'merra_precip_%d.nc' % year
    print('Saving to ' + filn)
    return filn

# MERRA daily precip
years = range(1979, 2015)
months = range(1, 13)
nperday = 24

for yr in years:
    for mon in months:
        prec1 = merra.read_daily('precip', yr, mon)
        dayvals = atm.season_days(atm.month_str(mon), atm.isleap(yr))
        prec1 = atm.daily_from_subdaily(prec1, nperday, dayvals=dayvals)
        if mon == 1:
            precip = prec1
        else:
            precip = xray.concat([precip, prec1], dim='day')

    # Save daily mean data for this year
    atm.save_nc(filename(yr), precip)
