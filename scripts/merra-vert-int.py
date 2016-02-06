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

# ----------------------------------------------------------------------
# Download precipitable water

datadir = atm.homedir() + 'datastore/merra/daily/'
years = np.arange(1979, 2015)
months = np.arange(1, 13)

varnms = ['TQV']
subset_dict = {'lon' : (40, 120)}
subset = '_40E-120E_90S-90N'
nperday = 24

def get_url(year, month, day):
    url = ('http://goldsmr2.sci.gsfc.nasa.gov/opendap/MERRA/MAI1NXINT.5.2.0/'
           '%d/%02d/MERRA100.prod.assim.inst1_2d_int_Nx.%d%02d%02d.hdf')
    url = url % (year, month, year, month, day)
    return url

def savefile(datadir, varnm, year, month, subset):
    return '%smerra_%s%s_%d%02d.nc' % (datadir, varnm, subset, year, month)


for varnm in varnms:
    for year in years:
        for month in months:
            days = range(1, atm.days_this_month(year, month) + 1)
            jdays = atm.season_days(atm.month_str(month), atm.isleap(year))
            urls = [get_url(year, month, day) for day in days]
            var = atm.load_concat(urls, varnm, concat_dim='TIME',
                                  subset_dict=subset_dict)
            var = atm.daily_from_subdaily(var, nperday, dayname='day',
                                          dayvals=jdays)
            filenm = savefile(datadir, varnm, year, month, subset)
            print('Saving to ' + filenm)
            atm.save_nc(filenm, var)
