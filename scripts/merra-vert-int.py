import sys
sys.path.append('/home/jwalker/dynamics/python/atmos-tools')
sys.path.append('/home/jwalker/dynamics/python/atmos-read')

import os
import xray
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
import pandas as pd
import atmos as atm
import precipdat
import merra

# ----------------------------------------------------------------------
# Download vertically integrated fluxes

version, years = 'merra2', np.arange(1980, 2016)
# version, years = 'merra', np.arange(1979, 2016)

varnms = ['UFLXQV', 'VFLXQV', 'DQVDT_ANA']
plev = None

datadir = atm.homedir() + 'datastore/' + version + '/daily/'
months = np.arange(1, 13)
lon1, lon2 = 40, 120
lat1, lat2 = -90, 90

subset_dict = {'lon' : (lon1, lon2), 'lat' : (lat1, lat2)}
if plev is not None:
    subset_dict['plev'] = plev
lonstr = atm.latlon_labels([lon1, lon2], 'lon', deg_symbol=False, join_str='-')
latstr = atm.latlon_labels([lat1, lat2], 'lat', deg_symbol=False, join_str='-')
subset = '_%s_%s' % (lonstr, latstr)
if plev is not None:
    subset = '%d' % plev + subset
time_dim = {'merra' : 'TIME', 'merra2' : 'time'}[version]


def monthlyfile(datadir, varnm, year, month, subset):
    return '%smerra_%s%s_%d%02d.nc' % (datadir, varnm, subset, year, month)

def yrlyfile(datadir, varnm, year, subset):
    return '%smerra_%s%s_%d.nc' % (datadir, varnm, subset, year)

for varnm in varnms:
    for year in years:
        for month in months:
            url_dict = merra.get_urls(year, month, version, varnm)
            days = range(1, atm.days_this_month(year, month) + 1)
            jdays = atm.season_days(atm.month_str(month), atm.isleap(year))
            urls = [url_dict['%d%02d%02d' % (year, month, day)] for day in days]
            var = atm.load_concat(urls, varnm, concat_dim=time_dim,
                                  subset_dict=subset_dict)
            nperday = len(var[time_dim]) / len(days)
            var = atm.daily_from_subdaily(var, nperday, dayname='day',
                                          dayvals=jdays)
            filenm = monthlyfile(datadir, varnm, year, month, subset)
            print('Saving to ' + filenm)
            var.to_dataset().to_netcdf(filenm)

        # Consolidate monthly files into yearly files
        files = [monthlyfile(datadir, varnm, year, m, subset) for m in months]
        var = atm.load_concat(files, varnm, concat_dim='day')
        filenm = yrlyfile(datadir, varnm, year, subset)
        print('Saving to ' + filenm)
        var.to_dataset().to_netcdf(filenm)
        print('Deleting monthly files:')
        for filenm in files:
            print(filenm)
            os.remove(filenm)
